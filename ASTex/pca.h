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



#ifndef _ASTEX_PCA_H_
#define _ASTEX_PCA_H_

#include <ASTex/mask.h>
#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

namespace ASTex
{

class ASTEX_API PCA
{
protected:
	ImageRGBd m_im;
	// eigenvectors
	double m_v1 [3];
	double m_v2 [3];
	double m_v3 [3];
	// mean values
	double m_meancolor [3];
	// covariance matrix
	double m_covar [3][3];


	double dot(double v[3], double w[3]);

	double normalize(double v[3]);

	void cross(double v[3], double w[3], double res[3]);


public :
	PCA (const ImageRGBd& input);

	void computeMeanColor(const Mask& mask);

	void computeCovariance(const Mask& mask);

	void mat_mult(double vin[3], double vout[3]);

	void computeEigenVectors();

	void computePCA (const Mask& m);

	void exportToTXT(std::string savefile);

	void project(ImageGrayd& res1, ImageGrayd& res2, ImageGrayd& res3) const;

	void back_project(const ImageGrayd& coord1, const ImageGrayd& coord2, const ImageGrayd& coord3, ImageRGBd& res) const;
};

} //ASTex

#endif
