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



#ifndef __ASTEX_PNG_INDEXED_IO__
#define __ASTEX_PNG_INDEXED_IO__

#include "dll.h"

namespace ASTex
{

/**
 * @brief read indexed png image with palette
 * @param filename
 * @param width
 * @param height
 * @param nb_colors
 * @param palette
 * @param data
 * @return error code 0 = ok
 */
int ASTEX_API read_png_palette(const char *filename, int& width, int& height, int& nb_colors, unsigned char*& palette, unsigned char*& data);

/**
 * @brief write_png_palette
 * @param filename
 * @param width
 * @param height
 * @param nb_colors
 * @param palette
 * @param data
 * @return error code 0 = ok
 */
int ASTEX_API write_png_palette(const char *filename, int width, int height, int nb_colors, unsigned char* palette, unsigned char* data);

}

#endif
