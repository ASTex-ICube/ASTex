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



#ifndef __ASTEX_IMG__COLLECTOR__
#define __ASTEX_IMG__COLLECTOR__

#include <ASTex/image_gray.h>

namespace ASTex
{

/**
 * \class image_collector io.h
 * \brief collects images and builds a patchwork
 */
class ASTEX_API ImageCollector
{
public :
	/**
	 * \brief add an image in the collection
     * \param [in] im scalar image
	 * \param [in] mean_shift (default = 0.0) value for the mean value shifting
	 */
	void add (const ImageGrayd& im, double mean_shift = 0.0);
	/**
	 * \brief computes a patchwork image from the collection
	 * \return patchwork image
	 */
	ImageGrayd collect ();

private :
	std::vector<ImageGrayd> m_images;
    std::vector<double> m_shifts;
};


}
#endif
