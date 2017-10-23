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




#include <ASTex/easy_io.h>


namespace ASTex
{

int import_texton(const std::string& file_path, ImageRGBd& texton );

/**
 * @brief decomposition periodic + smooth
 * @param input input grayd image
 * @return a pair of image (periodic,smooth)
 */
std::pair<ImageGrayd,ImageGrayd> decompo(const ImageGrayd& input);

/**
 * @brief corrige les images ?
 * @param in input
 * @param per periodic (clamped to [0,1])
 * @param smo smooth (shifted +0.5 & clamped to [0,1])
 */
void corrige(const ImageGrayd& in, ImageGrayd& per, ImageGrayd& smo);

/**
 * @brief analyse_texton
 * @param input_color
 * @param images_alpha_visu
 * @param images_alpha_visu_bounded
 */
void analyse_texton( const ASTex::ImageRGBd& input_color,
					 std::vector<ASTex::ImageRGBd>& images_alpha_visu,
					 std::vector<ASTex::ImageRGBd>& images_alpha_visu_bounded);


}
