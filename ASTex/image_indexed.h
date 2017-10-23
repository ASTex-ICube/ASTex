/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/



#ifndef __ASTEX_IMAGE_PALETTE__
#define __ASTEX_IMAGE_PALETTE__

#include <ASTex/image_gray.h>
#include <ASTex/png_indexed_io.h>
#include <ASTex/image_rgb.h>


namespace ASTex
{


// some small function to easy float/double -> uint8_t
// and uint8_t -> double/float

template <typename CHANNEL_TYPE>
uint8_t toU8(CHANNEL_TYPE v)
{
	 return (uint8_t)(v*255);
}

template <>
inline uint8_t toU8(uint8_t v)
{
	 return v;
}

template <typename CHANNEL_TYPE>
CHANNEL_TYPE fromU8(uint8_t v)
{
	 return CHANNEL_TYPE(v)/255;
}

template <>
inline uint8_t fromU8<uint8_t>(uint8_t v)
{
	 return v;
}


/**
 * @brief The Image indexed class
 * @tparam COLOR_CHANNEL_TYPE (uchar/float/double)
 */
template<typename COLOR_CHANNEL_TYPE>
class ImageIndexed: public ImageGrayu8
{
public:
	typedef typename itk::RGBPixel<COLOR_CHANNEL_TYPE> PaletteColorType;


protected:
	typename std::vector< PaletteColorType > palette_;

public:

	/**
	 * @brief load an indexed PNG file
	 * @param filename
	 * @return
	 */
	bool loadIndexedPNG(const std::string& filename)
	{
		int width;
		int height;
		int nb_colors;
		uint8_t* palette = nullptr;
		uint8_t* data = nullptr;
		if (read_png_palette(filename.c_str(), width, height, nb_colors, palette, data) == 0)
		{

			initItk(width,height);

			this->itk_img_->GetPixelContainer()->SetImportPointer(data,width*height,true);

			palette_.resize(nb_colors);
			for (int i=0;i<nb_colors;++i)
			{
				palette_[i][0] = toU8<COLOR_CHANNEL_TYPE>(palette[3*i]);
				palette_[i][1] = toU8<COLOR_CHANNEL_TYPE>(palette[3*i+1]);
				palette_[i][2] = toU8<COLOR_CHANNEL_TYPE>(palette[3*i+2]);
			}
			delete[] palette;
			return true;

		}

		std::cerr << "Error loading "<< filename << std::endl;
		return false;
	}

	/**
	 * @brief save in an indexed PNG file
	 * @param filename
	 */
	void saveIndexedPNG(const std::string& filename) /*const*/
	{

		uint8_t* localPalette = new uint8_t[3*palette_.size()];
		for (size_t i=0; i<palette_.size(); ++i)
		{
			localPalette[3*i  ] = fromU8<COLOR_CHANNEL_TYPE>(palette_[i][0]);
			localPalette[3*i+1] = fromU8<COLOR_CHANNEL_TYPE>(palette_[i][1]);
			localPalette[3*i+2] = fromU8<COLOR_CHANNEL_TYPE>(palette_[i][2]);
		}

		write_png_palette(filename.c_str(), this->width(), this->height(), palette_.size(), localPalette, this->itk_img_->GetPixelContainer()->GetBufferPointer());

		delete[] localPalette;
	}


	/**
	 * @brief get ref the pixel color
	 * @param ind index in palette
	 * @return
	 */
	PaletteColorType& color(uint8_t ind)
	{
		return palette_[ind];
	}


	/**
	 * @brief get the palette
	 * @return ref on vector of itk::RGBPixel<COLOR_CHANNEL_TYPE>
	 */
	typename std::vector< PaletteColorType >& palette()
	{
		return palette_;
	}

	/**
	 * @brief get the palette (const version)
	 * @return const ref vector of itk::RGBPixel<COLOR_CHANNEL_TYPE>
	 */
	const typename std::vector< PaletteColorType >& palette() const
	{
		return palette_;
	}

	/**
	 * @brief create a RGB image from indexed image
	 * @return pointer on the created RGB image
	 */
	ImageRGB<COLOR_CHANNEL_TYPE> createRGB()
	{
		ImageRGB<COLOR_CHANNEL_TYPE> dst;
		dst.initItk();

		dst.itk()->SetRegions(this->itk_img_->GetLargestPossibleRegion());
		dst.itk()->Allocate();


		IteratorIndexed inputIt = this->beginIteratorIndexed();
		typename ImageRGB<COLOR_CHANNEL_TYPE>::IteratorIndexed outputIt = dst.beginIteratorIndexed();

		while( !inputIt.IsAtEnd() )
		{
			outputIt.Set(palette_[inputIt.Get()]);
			++inputIt;
			++outputIt;
		}

		return  dst;
	}



};



typedef ImageIndexed< uint8_t >	ImageIndexedu8;

typedef ImageIndexed< float >			ImageIndexedf;

typedef ImageIndexed< double >			ImageIndexedd;


}

#endif



