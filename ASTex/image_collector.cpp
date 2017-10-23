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



#include <ASTex/image_collector.h>
#include <iostream>

namespace ASTex
{

void ImageCollector::add (const ImageGrayd& im, double mean_shift)
{
    m_images.push_back(im);
    m_shifts.push_back(mean_shift);
}

ImageGrayd ImageCollector::collect ()
{
	int W = 0, H = 0;

    for (uint32_t i=0; i<m_images.size(); ++i)
    {
		W += m_images[i].width();
		H = std::max(H,m_images[i].height());
    }

	ImageGrayd image(W,H,true);

    uint32_t dx=0;
    for (uint32_t i=0; i<m_images.size(); ++i)
    {
		uint32_t Wi = m_images[i].width();
		uint32_t Hi = m_images[i].height();

		uint32_t dxi = m_images[i].itk()->GetOrigin()[0];
		uint32_t dyi = m_images[i].itk()->GetOrigin()[1];

        double shift = m_shifts[i];

        for (uint32_t x=0 ; x < Wi ; ++x)
        {
            for (uint32_t y=0 ; y < Hi ; ++y)
            {
				image.pixelAbsolute(x+dx,y) = m_images[i].pixelAbsolute(x+dxi,y+dyi) + shift;
            }
        }
        dx += Wi;
    }
    return image;
}

} // namespace
