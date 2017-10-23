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



#include <ASTex/mask.h>

#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

namespace ASTex
{


void Mask::export_binary_image (ImageGrayu8& im) 
{
	im.for_all_pixels([&] (ImageGrayu8::PixelType& p, int i, int j)
	{
		p = 255*this->operator()(i,j);
	});
}


MaskBool::MaskBool(int width, int height):
	MaskImage<uint8_t>(width,height)
{
	img_mask_.for_all_pixels([] (MaskPixelType& p)
	{
		p = 255; // ???
	});
}

MaskBool::MaskBool(int width, int height, bool val):
	MaskImage<uint8_t>(width,height)
{
	if (val)
		img_mask_.for_all_pixels([] (MaskPixelType& p)
		{
			p = 255;
		});
	else
		img_mask_.for_all_pixels([] (MaskPixelType& p)
		{
			p = 0;
		});

}

void MaskBool::clone(const MaskBool& mb)
{
	if(!img_mask_.is_initialized_as(mb.img_mask_))
	{
		img_mask_.initItk(mb.img_mask_.width(),mb.img_mask_.height());
	}

	int nb = img_mask_.width()*img_mask_.height();
	const uint8_t* ptr_src = mb.img_mask_.getPixelsPtr();
	uint8_t* ptr_dst = img_mask_.getPixelsPtr();
	for(int i=0;i<nb;++i)
		*ptr_dst++ = *ptr_src++;
}


void MaskBool::random(double proportion)
{
	srand(time(NULL));
	img_mask_.for_all_pixels([&] (MaskPixelType& p)
	{
		p = (double(rand()) / double(RAND_MAX)) <= proportion;
	});
}

void MaskBool::clear(bool val)
{
	MaskPixelType v = 255*MaskPixelType(val);
	img_mask_.for_all_pixels([v](MaskPixelType& P)
	{
		P = v;
	});
}


void MaskBool::neg()
{
	int nb = img_mask_.width()*img_mask_.height();
	uint8_t* ptr = img_mask_.getPixelsPtr();
	for(int i=0;i<nb;++i)
	{
		*ptr = ~(*ptr);
		++ptr;
	}
}


MaskBool MaskBool::operator~() const
{
	MaskBool mb(img_mask_.width(),img_mask_.height());

	int nb = img_mask_.width()*img_mask_.height();
	const uint8_t* ptr_src = img_mask_.getPixelsPtr();
	uint8_t* ptr_dst = mb.img_mask_.getPixelsPtr();
	for(int i=0;i<nb;++i)
		*ptr_dst++ = ~(*ptr_src++);
	return mb;
}


MaskBool MaskBool::erosion(int radius) const
{
	typedef itk::BinaryBallStructuringElement<MaskPixelType, 2>  StructuringElementType;
	StructuringElementType structuringElement;
	structuringElement.SetRadius(radius);
	structuringElement.CreateStructuringElement();

	typedef itk::BinaryErodeImageFilter <MaskImageType::ItkImg, MaskImageType::ItkImg, StructuringElementType> 	BinaryErodeImageFilterType;
	BinaryErodeImageFilterType::Pointer erodeFilter	= BinaryErodeImageFilterType::New();

	erodeFilter->SetErodeValue(MaskPixelType(255));
	erodeFilter->SetInput(this->img_mask_.itk());
	erodeFilter->SetKernel(structuringElement);
	erodeFilter->Update();
	return MaskBool(erodeFilter->GetOutput());
}

int MaskBool::count_true() const
{
	int co=0;
	int nb = img_mask_.width()*img_mask_.height();
	const uint8_t* ptr = img_mask_.getPixelsPtr();
	for(int i=0;i<nb;++i)
		if (*ptr++) co++;
	return co;
}

bool MaskBool::full_false() const
{
	int nb = img_mask_.width()*img_mask_.height();
	const uint8_t* ptr = img_mask_.getPixelsPtr();
	for(int i=0;i<nb;++i)
		if (*ptr++ != 0) return false;
	return true;
}



MaskDist::MaskDist(int32_t width, int32_t height, double _dist, std::vector<std::vector<int32_t>> _coord_seeds ):
	MaskBool(width,height),
	dist(_dist),coord_seeds(_coord_seeds)
{
	update();
}

double MaskDist::get_ratio_true()
{
	int32_t W = this->img_mask_.width();
	int32_t H = this->img_mask_.height();
	double count_pix=W*H;
	double count_true=0;

	for (int32_t y = 0; y < H; ++y)
		for (int32_t x = 0; x < W; ++x)
			if (img_mask_.pixelAbsolute(x,y)==1)
				count_true++;

	return count_true/count_pix;
}

void MaskDist::update_dist(double _dist)
{
	dist=_dist;
	update();
}

void MaskDist::update_seeds(std::vector<std::vector<int32_t>> new_seeds)
{
	coord_seeds.clear();
	coord_seeds.swap(new_seeds);
	update();
}


void MaskDist::update()
{
	img_mask_.for_all_pixels([&] (MaskPixelType& p)
	{
		p = 1;
	});

	int32_t W = this->img_mask_.width();
	int32_t H = this->img_mask_.height();

	for (int32_t y = 0; y < H; ++y)
	{
		for (int32_t x = 0; x < W; ++x)
		{
			for (std::size_t k = 0; k < coord_seeds.size(); ++k)
			{
				double current_dist = std::sqrt((x-coord_seeds[k][0])*(x-coord_seeds[k][0])+(y-coord_seeds[k][1])*(y-coord_seeds[k][1]));
				if (current_dist<dist)
				{
					img_mask_.pixelAbsolute(x,y)=0;
					break;
				}
			}
		}
	}
}



} // namespace ASTex
