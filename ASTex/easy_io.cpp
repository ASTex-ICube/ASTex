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




#include "easy_io.h"
#include "colorspace_filters.h"

namespace ASTex
{

namespace IO
{

void ASTEX_API save01_in_u8(const ImageGrayd& img, const std::string& filename)
{
	auto filter = ColorSpace::FilterGray01To255< ImageGrayd::ItkImg,ImageGrayu8::ItkImg>::New();
	filter->SetInput(img.itk());
	ImageGrayu8 out(filter->GetOutput());
	out.save(filename);
}

void ASTEX_API save01_in_u8(const ImageGrayf& img, const std::string& filename)
{
	auto filter = ColorSpace::FilterGray01To255< ImageGrayf::ItkImg,ImageGrayu8::ItkImg>::New();
	filter->SetInput(img.itk());
	ImageGrayu8 out(filter->GetOutput());
	out.save(filename);
}


void ASTEX_API save01_in_u8(const ImageRGBd& img, const std::string& filename)
{
	auto filter = ColorSpace::FilterRGB01To255< ImageRGBd::ItkImg,ImageRGBu8::ItkImg>::New();
	filter->SetInput(img.itk());
	ImageRGBu8 out(filter->GetOutput());
	out.save(filename);
}

void ASTEX_API save01_in_u8(const ImageRGBf& img, const std::string& filename)
{
	auto filter = ColorSpace::FilterRGB01To255< ImageRGBf::ItkImg,ImageRGBu8::ItkImg>::New();
	filter->SetInput(img.itk());
	ImageRGBu8 out(filter->GetOutput());
	out.save(filename);
}

void ASTEX_API save01_in_u8(const ImageRGBAd& img, const std::string& filename)
{
	auto filter = ColorSpace::FilterRGB01To255< ImageRGBAd::ItkImg,ImageRGBAu8::ItkImg>::New();
	filter->SetInput(img.itk());
	ImageRGBAu8 out(filter->GetOutput());
	out.save(filename);
}

void ASTEX_API save01_in_u8(const ImageRGBAf& img, const std::string& filename)
{
	auto filter = ColorSpace::FilterRGB01To255< ImageRGBAf::ItkImg,ImageRGBAu8::ItkImg>::New();
	filter->SetInput(img.itk());
	ImageRGBAu8 out(filter->GetOutput());
	out.save(filename);
}



bool ASTEX_API loadu8_in_01(ImageGrayd& img, const std::string& filename)
{
	ImageGrayu8 in;
	bool ok = in.load(filename);
	if (!ok)
		return false;

	auto filter = ColorSpace::FilterGray255To01<ImageGrayu8::ItkImg,ImageGrayd::ItkImg>::New();
	filter->SetInput(in.itk());
	filter->Update();
	img = ImageGrayd(filter->GetOutput());
	return true;
}


bool ASTEX_API loadu8_in_01(ImageGrayf& img, const std::string& filename)
{
	ImageGrayu8 in;
	bool ok = in.load(filename);
	if (!ok)
		return false;

	auto filter = ColorSpace::FilterGray255To01<ImageGrayu8::ItkImg,ImageGrayf::ItkImg>::New();
	filter->SetInput(in.itk());
	filter->Update();
	img = ImageGrayf(filter->GetOutput());
	return true;
}


bool ASTEX_API loadu8_in_01(ImageRGBd& img, const std::string& filename)
{
	ImageRGBu8 in;
	bool ok = in.load(filename);
	if (!ok)
		return false;

	auto filter = ColorSpace::FilterRGB255To01<ImageRGBu8::ItkImg,ImageRGBd::ItkImg>::New();
	filter->SetInput(in.itk());
	filter->Update();
	img = ImageRGBd(filter->GetOutput());
	return true;
}



bool ASTEX_API loadu8_in_01(ImageRGBf& img, const std::string& filename)
{
	ImageRGBu8 in;
	bool ok = in.load(filename);
	if (!ok)
		return false;

	auto filter = ColorSpace::FilterRGB255To01<ImageRGBu8::ItkImg,ImageRGBf::ItkImg>::New();
	filter->SetInput(in.itk());
	filter->Update();
	img = ImageRGBf(filter->GetOutput());
	return true;
}


bool ASTEX_API loadu8_in_01(ImageRGBAd& img, const std::string& filename)
{
	ImageRGBAu8 in;
	bool ok = in.load(filename);
	if (!ok)
		return false;

	auto filter = ColorSpace::FilterRGB255To01<ImageRGBAu8::ItkImg,ImageRGBAd::ItkImg>::New();
	filter->SetInput(in.itk());
	filter->Update();
	img = ImageRGBAd(filter->GetOutput());
	return true;
}



bool ASTEX_API loadu8_in_01(ImageRGBAf& img, const std::string& filename)
{
	ImageRGBAu8 in;
	bool ok = in.load(filename);
	if (!ok)
		return false;

	auto filter = ColorSpace::FilterRGB255To01<ImageRGBAu8::ItkImg,ImageRGBAf::ItkImg>::New();
	filter->SetInput(in.itk());
	filter->Update();
	img = ImageRGBAf(filter->GetOutput());
	return true;
}



void ASTEX_API save01_in_u16(const ImageGrayd& img, const std::string& filename)
{
    auto filter = ColorSpace::FilterGray01To65535< ImageGrayd::ItkImg,ImageGrayu16::ItkImg>::New();
    filter->SetInput(img.itk());
    ImageGrayu16 out(filter->GetOutput());
    out.save(filename);
}

void ASTEX_API save01_in_u16(const ImageGrayf& img, const std::string& filename)
{
    auto filter = ColorSpace::FilterGray01To65535< ImageGrayf::ItkImg,ImageGrayu16::ItkImg>::New();
    filter->SetInput(img.itk());
    ImageGrayu16 out(filter->GetOutput());
    out.save(filename);
}


void ASTEX_API save01_in_u16(const ImageRGBd& img, const std::string& filename)
{
    auto filter = ColorSpace::FilterRGB01To65535< ImageRGBd::ItkImg,ImageRGBu16::ItkImg>::New();
    filter->SetInput(img.itk());
    ImageRGBu16 out(filter->GetOutput());
    out.save(filename);
}

void ASTEX_API save01_in_u16(const ImageRGBf& img, const std::string& filename)
{
    auto filter = ColorSpace::FilterRGB01To65535< ImageRGBf::ItkImg,ImageRGBu16::ItkImg>::New();
    filter->SetInput(img.itk());
    ImageRGBu16 out(filter->GetOutput());
    out.save(filename);
}

void ASTEX_API save01_in_u16(const ImageRGBAd& img, const std::string& filename)
{
    auto filter = ColorSpace::FilterRGB01To65535< ImageRGBAd::ItkImg,ImageRGBAu16::ItkImg>::New();
    filter->SetInput(img.itk());
    ImageRGBAu16 out(filter->GetOutput());
    out.save(filename);
}

void ASTEX_API save01_in_u16(const ImageRGBAf& img, const std::string& filename)
{
    auto filter = ColorSpace::FilterRGB01To65535< ImageRGBAf::ItkImg,ImageRGBAu16::ItkImg>::New();
    filter->SetInput(img.itk());
    ImageRGBAu16 out(filter->GetOutput());
    out.save(filename);
}

bool ASTEX_API loadu16_in_01(ImageGrayd& img, const std::string& filename)
{
    ImageGrayu16 in;
    bool ok = in.load(filename);
    if (!ok)
        return false;

    auto filter = ColorSpace::FilterGray65535To01<ImageGrayu16::ItkImg,ImageGrayd::ItkImg>::New();
    filter->SetInput(in.itk());
    filter->Update();
    img = ImageGrayd(filter->GetOutput());
    return true;
}



bool ASTEX_API loadu16_in_01(ImageGrayf& img, const std::string& filename)
{
    ImageGrayu16 in;
    bool ok = in.load(filename);
    if (!ok)
        return false;

    auto filter = ColorSpace::FilterGray65535To01<ImageGrayu16::ItkImg,ImageGrayf::ItkImg>::New();
    filter->SetInput(in.itk());
    filter->Update();
    img = ImageGrayf(filter->GetOutput());
    return true;
}


bool ASTEX_API loadu16_in_01(ImageRGBd& img, const std::string& filename)
{
    ImageRGBu16 in;
    bool ok = in.load(filename);
    if (!ok)
        return false;

    auto filter = ColorSpace::FilterRGB65535To01<ImageRGBu16::ItkImg,ImageRGBd::ItkImg>::New();
    filter->SetInput(in.itk());
    filter->Update();
    img = ImageRGBd(filter->GetOutput());
    return true;
}



bool ASTEX_API loadu16_in_01(ImageRGBf& img, const std::string& filename)
{
    ImageRGBu16 in;
    bool ok = in.load(filename);
    if (!ok)
        return false;

    auto filter = ColorSpace::FilterRGB65535To01<ImageRGBu16::ItkImg,ImageRGBf::ItkImg>::New();
    filter->SetInput(in.itk());
    filter->Update();
    img = ImageRGBf(filter->GetOutput());
    return true;
}


bool ASTEX_API loadu16_in_01(ImageRGBAd& img, const std::string& filename)
{
    ImageRGBAu16 in;
    bool ok = in.load(filename);
    if (!ok)
        return false;

    auto filter = ColorSpace::FilterRGB65535To01<ImageRGBAu16::ItkImg,ImageRGBAd::ItkImg>::New();
    filter->SetInput(in.itk());
    filter->Update();
    img = ImageRGBAd(filter->GetOutput());
    return true;
}



bool ASTEX_API loadu16_in_01(ImageRGBAf& img, const std::string& filename)
{
    ImageRGBAu16 in;
    bool ok = in.load(filename);
    if (!ok)
        return false;

    auto filter = ColorSpace::FilterRGB65535To01<ImageRGBAu16::ItkImg,ImageRGBAf::ItkImg>::New();
    filter->SetInput(in.itk());
    filter->Update();
    img = ImageRGBAf(filter->GetOutput());
    return true;
}



bool ASTEX_API loadRGB2gray(const std::string& filename, ImageGrayd& output)
{

	ImageRGBd input_color;

	input_color.load(filename);
	output.initItk(input_color.width(), input_color.height(),true);

	typedef ImageRGBd::ItkImg IMG_DBL;

	// first transform [0,255] double -> [0,1] double
	ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::Pointer filter0 =
			ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::New();
	filter0->SetInput(input_color.itk());

	filter0->Update();

	input_color.itk() = filter0->GetOutput() ;

	output.initItk(input_color.width(),input_color.height(),true);

	for (int j = 0; j < input_color.height() ; ++j)
		for (int i = 0; i < input_color.width(); ++i)
		{
			if (input_color.pixelAbsolute(i,j)[0]  == input_color.pixelAbsolute(i,j)[0])
			{
				output.pixelAbsolute(i,j)= input_color.pixelAbsolute(i,j)[0] ;
			}
			else
			{
				output.pixelAbsolute(i,j)=0;
			}
		}

	return true;
}



bool ASTEX_API loadRGB2luminance(const std::string& filename, ImageGrayd& output)
{
	typedef ImageRGBd::ItkImg IMG_DBL;

	ImageRGBd input_color;
	input_color.load(filename);

	const int im_w = input_color.width();
	const int im_h = input_color.height();
	output.initItk(im_w, im_h,true);


	// first transform [0,255] double -> [0,1] double
	ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::Pointer filter0 =
			ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::New();
	filter0->SetInput(input_color.itk());

		// RGB double -> XYZ float
	ColorSpace::FilterRGBtoXYZ<IMG_DBL,IMG_DBL>::Pointer filter1 =
			ColorSpace::FilterRGBtoXYZ<IMG_DBL,IMG_DBL>::New();
	filter1->SetInput(filter0->GetOutput());
	filter1->Update();

	ImageRGBd input_XYZ;
	input_XYZ.itk()= filter1->GetOutput() ;

	for (int j = 0; j < im_h ; ++j)
		for (int i = 0; i < im_w; ++i)
		{
			if(input_XYZ.pixelAbsolute(i,j)[1]  == input_XYZ.pixelAbsolute(i,j)[1])
			{
				output.pixelAbsolute(i,j)= input_XYZ.pixelAbsolute(i,j)[1] ;
			}
			else
			{
				output.pixelAbsolute(i,j)=0;
			}
		}

	return true;
}


bool ASTEX_API load2lightness(const std::string& filename, ImageGrayd& output)
{

	typedef ImageRGBd::ItkImg IMG_DBL;

	ImageRGBd input_color;
	input_color.load(filename);

	const int im_w = input_color.width();
	const int im_h = input_color.height();
	output.initItk(im_w, im_h,true);


	// first transform [0,255] double -> [0,1] double
	ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::Pointer filter0 =
			ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::New();
	filter0->SetInput(input_color.itk());

		// RGB double -> XYZ float
	ColorSpace::FilterRGBtoXYZ<IMG_DBL,IMG_DBL>::Pointer filter1 =
			ColorSpace::FilterRGBtoXYZ<IMG_DBL,IMG_DBL>::New();
	filter1->SetInput(filter0->GetOutput());

	// XYZ float -> LUV float
	ColorSpace::FilterXYZtoLUV<IMG_DBL,IMG_DBL>::Pointer filter2 =
			ColorSpace::FilterXYZtoLUV<IMG_DBL,IMG_DBL>::New();
	filter2->SetInput(filter1->GetOutput());

	filter2->Update();

	ImageRGBd input_LUV;
	input_LUV.itk()= filter2->GetOutput() ;

	for (int j = 0; j < im_h ; ++j)
		for (int i = 0; i < im_w; ++i)
		{
			if(input_LUV.pixelAbsolute(i,j)[0]  == input_LUV.pixelAbsolute(i,j)[0])
			{
				output.pixelAbsolute(i,j)= input_LUV.pixelAbsolute(i,j)[0] / 100.0 ;
			}
			else
			{
				output.pixelAbsolute(i,j)=0;
			}

		}

	return true;
}





}

}
