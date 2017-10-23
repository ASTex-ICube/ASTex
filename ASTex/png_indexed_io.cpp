#include "png_indexed_io.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>

#ifdef WIN32
#include <itk_png.h>
#else
#include <png.h>
#endif // WIN32



namespace ASTex
{

int ASTEX_API read_png_palette(const char *filename, int& width, int& height, int& nb_colors, unsigned char*& palette, unsigned char*& data)
{

	png_byte color_type;
	png_byte bit_depth;
	png_bytep *row_pointers;
	png_bytep buffer_img;


	FILE *fp = fopen(filename, "rb");
	if(!fp)
	{
		std::cerr << "can not open file "<<filename<<std::endl;
		return 1;// can not open file
	}

	png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if(!png)
	{
		std::cerr << "out of memory while reading file "<<filename<<std::endl;
		return 4;	// out of memory
	}

	png_infop info = png_create_info_struct(png);
	if(!info)
	{
		std::cerr << "out of memory while reading png file "<<filename<<std::endl;
		free(png);
		return 4;	// out of memory
	}

	if(setjmp(png_jmpbuf(png)))
	{
		std::cerr << "problem while reading png file "<<filename<<std::endl;
		png_destroy_read_struct(&png, &info, NULL);
		return 2	;
	}

	png_init_io(png, fp);
	png_read_info(png, info);

	width      = png_get_image_width(png, info);
	height     = png_get_image_height(png, info);
	color_type = png_get_color_type(png, info);
	bit_depth  = png_get_bit_depth(png, info);

	if( color_type != PNG_COLOR_TYPE_PALETTE )
	{
		std::cerr << "no color palette in png file "<<filename<<std::endl;
		return 2;
	}

	png_read_update_info(png, info);

	png_uint_32 rowbytes = png_uint_32(png_get_rowbytes(png, info));

	// alloc image buffer
	buffer_img = (png_byte*)malloc(rowbytes*height);
	if (buffer_img == NULL)
	{
		png_destroy_read_struct(&png, &info, NULL);
		return 4; // out of memory
	}

	// alloc pointers table of lines
	row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
	// compute
	for(int y = 0; y < height; y++)
		row_pointers[y] = buffer_img + y*rowbytes;

	// read image
	png_read_image(png, row_pointers);


	// translate in simple unsigned char
	data = new unsigned char[width*height];
	unsigned char* ptr = data;

	switch(bit_depth)
	{
		case 1:
			for(int y = 0; y < height; y++)
			{
				for(int x = 0; x < width/8; x++)
				{
					png_byte p = row_pointers[y][x];
					*ptr++ = p>>7;
					*ptr++ = (p>>6)&0x1;
					*ptr++ = (p>>5)&0x1;
					*ptr++ = (p>>4)&0x1;
					*ptr++ = (p>>3)&0x1;
					*ptr++ = (p>>2)&0x1;
					*ptr++ = (p>>1)&0x1;
					*ptr++ = p&0x1;
				}
			}
			break;
		case 2:
			for(int y = 0; y < height; y++)
			{
				for(int x = 0; x < width/4; x++)
				{
					png_byte p = row_pointers[y][x];
					*ptr++ = p>>6;
					*ptr++ = (p>>4)&0x3;
					*ptr++ = (p>>2)&0x3;
					*ptr++ = p&0x3;
				}
			}
			break;
		case 4:
			for(int y = 0; y < height; y++)
			{
				for(int x = 0; x < width/2; x++)
				{
					png_byte p = row_pointers[y][x];
					*ptr++ = p>>4;
					*ptr++ = p&0xf;
				}
			}
			break;
		case 8:
			for(int y = 0; y < height; y++)
			{
				memcpy(ptr,row_pointers[y],width);
				ptr += width;
			}
			break;
	}

	png_colorp locPal;
	png_uint_32 res = png_get_PLTE(png, info, &locPal, &nb_colors);

	if (res)
	{
		palette = new unsigned char[3*nb_colors];
		ptr = palette;

		memcpy(ptr,locPal,3*nb_colors);

		for(int i=0; i<nb_colors; ++i)
		{
			memcpy(ptr,&locPal[i],sizeof(png_color));
			ptr+= sizeof(png_color);
		}
	}

	// clean & free png alloc & close file
	png_read_end(png,info);
	png_destroy_read_struct(&png, &info, NULL);
	free(buffer_img);
	free(row_pointers);
	fclose(fp);


	if( res != PNG_INFO_PLTE )
	{
		std::cerr << "problem loading palette in png file "<<filename<<std::endl;
		return 5;
	}

	return 0; // OK
}




int ASTEX_API write_png_palette(const char *filename, int width, int height, int nb_colors, unsigned char* palette, unsigned char* data)
{
	png_byte bit_depth;
	png_bytep *row_pointers;
	png_bytep buffer_img;

	FILE *fp = fopen(filename, "wb");
	if(!fp) abort();

	png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png)
		return 4;


	png_infop info = png_create_info_struct(png);
	if(!info)
	{
		free(png);
		return 4;	// out of memory
	}

	if(setjmp(png_jmpbuf(png)))
	{
		png_destroy_write_struct(&png, &info);
		return 2	;
	}

	png_init_io(png, fp);

//	determnie bit_depth function of palette size
	bit_depth = 8;
	if (nb_colors<=2)
		bit_depth = 1;
	else if (nb_colors<=4)
		bit_depth = 2;
	else if (nb_colors<=16)
		bit_depth = 4;

	png_set_IHDR(
				png,
				info,
				width, height,
				bit_depth,
				PNG_COLOR_TYPE_PALETTE,
				PNG_INTERLACE_NONE,
				PNG_COMPRESSION_TYPE_DEFAULT,
				PNG_FILTER_TYPE_DEFAULT
				);

	// fill palette
	// first alloc
	png_colorp locPal;
	locPal = (png_colorp)malloc(nb_colors*sizeof(png_color));
	// traverse and memcopy
	unsigned char* ptr = palette;
	for(int i=0; i<nb_colors; ++i)
	{
		memcpy(&locPal[i],ptr,sizeof(png_color));
		ptr+= sizeof(png_color);
	}

	// set in png struct
	png_set_PLTE(png,info,locPal,nb_colors);
	png_write_info(png, info);

	// free local allocated palette
	free(locPal);

	png_uint_32 rowbytes = png_uint_32(png_get_rowbytes(png, info));
	// alloc image buffer
	buffer_img = (png_byte*)malloc(rowbytes*height);
	if (buffer_img == NULL)
	{
		png_destroy_read_struct(&png, &info, NULL);
		return 4; // out of memory
	}

	// alloc pointers table of lines
	row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
	// compute
	for(int y = 0; y < height; y++)
		row_pointers[y] = buffer_img + y*rowbytes;

	// reuse ptr for data traversal
	ptr = data;

	switch(bit_depth)
	{
		case 1:
			for(int y = 0; y < height; y++)
			{
				for(int x = 0; x < width/4; x++)
				{
					png_byte p = (*ptr++)<<7;
					p |= (*ptr++)<<6;
					p |= (*ptr++)<<5;
					p |= (*ptr++)<<4;
					p |= (*ptr++)<<3;
					p |= (*ptr++)<<2;
					p |= (*ptr++)<<1;
					p |= (*ptr++);
					row_pointers[y][x] = p;
				}
			}
			break;
		case 2:
			for(int y = 0; y < height; y++)
			{
				for(int x = 0; x < width/4; x++)
				{
					png_byte p = (*ptr++)<<6;
					p |= (*ptr++)<<4;
					p |= (*ptr++)<<2;
					p |= (*ptr++);
					row_pointers[y][x] = p;
				}
			}
			break;
		case 4:
			for(int y = 0; y < height; y++)
			{
				for(int x = 0; x < width/2; x++)
				{
					png_byte p = (*ptr++)<<4;
					p |= (*ptr++);
					row_pointers[y][x] = p;
				}
			}
			break;
		case 8:
			unsigned char* ptr = data;
			for(int y = 0; y < height; y++)
			{
				memcpy(row_pointers[y],ptr,width);
				ptr += width;
			}
			break;
	}


	// write png and finish
	png_write_image(png, row_pointers);
	png_write_end(png, info);
	png_destroy_write_struct(&png, &info);

	// free allocated buffers
	free(buffer_img);
	free(row_pointers);

	// and close file
	fclose(fp);

	return 0;
}


}


