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




#include <fstream>
#include <cmath>
#include <algorithm>

#include <ASTex/local_spectrum.h>
#include <ASTex/dll.h>


namespace ASTex
{


void ASTEX_API bilateral_filter(const ImageGrayd& input, ImageGrayd& output, int filtre_size, double id, double cd);

void ASTEX_API bilateral_filter(const ImageRGBd& input, ImageRGBd& output, int filtre_size, double id, double cd);



namespace Fourier
{

/**
 * @brief Filtre bilateral avec descripteur de frequence
 */
void ASTEX_API frequency_bilateral_filter(const ImageGrayd& input, ImageGrayd& output, int filtre_size, LocalSpectrum& lsp , double id, double cd);

void ASTEX_API frequency_bilateral_filter(const ImageRGBd& input_color, /*const ImageGrayd& input,*/ ImageRGBd& output, int filtre_size, LocalSpectrum& lsp , double id, double cd);

void ASTEX_API frequency_joint_bilateral_filter(const ImageGrayd& input, ImageGrayd& output,const std::vector<ImageGrayd>& guidance_maps, int filtre_size, double id, double cd);


} //namespace Fourier


} //namespace ASTex
