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



#include <string>

/**
 * @brief biscalenoisepatchexg
 * @param filename_source
 * @param base_dir
 * @param NCLUSTERS
 * @param CONTENTS
 * @param NSCALES
 * @param NROTANGLES
 * @param TSCALE
 * @param SZ_MULT
 */
void biscalenoisepatchexg(const std::string& filename_source, const std::string& base_dir, int NCLUSTERS, int CONTENTS, int NSCALES, int NROTANGLES, float TSCALE, int SZ_MULT);
