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

#include <Eigen/Eigen>
#include <cmath>

#include <ASTex/mask.h>
#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

namespace ASTex
{

template<typename T>
class ASTEX_API PCA
{
protected:

	using ImageType = ImageRGB<T>;
	using PixelType = typename ImageType::PixelType;
	using DataType	= typename ImageType::DataType;

	ImageType m_im;
	// eigenvectors
	double m_v1 [3];
	double m_v2 [3];
	double m_v3 [3];
	// mean values
	double m_meancolor [3];
	// covariance matrix
	double m_covar [3][3];
	double m_A_inverse[3][3]; // inverse transpose
	double m_det;
	double m_A[3][3];


	double dot(double v[3], double w[3]);

	double normalize(double v[3]);

	void cross(double v[3], double w[3], double res[3]);


public :
	PCA (const ImageType& input);

	void computeMeanColor(const Mask& mask);

	void computeCovariance(const Mask& mask);

	void mat_mult(double vin[3], double vout[3]);

	void computeEigenVectors();

	void computePCA (const Mask& m);

	void exportToTXT(std::string savefile);

	void project(ImageType &res);

	void back_project(const ImageType &coord, ImageType& res) const;

	Eigen::Vector3d eigenVector(unsigned int i);
};

template<typename T>
PCA<T>::PCA (const ImageType& input)
{
	m_im.share(input);
}

template<typename T>
void PCA<T>::computeMeanColor(const Mask& mask)
{
	m_meancolor[0] = 0;
	m_meancolor[1] = 0;
	m_meancolor[2] = 0;
	int nbpix = 0;
	for (int x = 0; x < m_im.width(); ++x)
	{
		for (int y = 0; y < m_im.height(); ++y)
		{
			if (mask(x,y))
			{
				++nbpix;
				m_meancolor[0] += m_im.pixelAbsolute(x,y)[0];
				m_meancolor[1] += m_im.pixelAbsolute(x,y)[1];
				m_meancolor[2] += m_im.pixelAbsolute(x,y)[2];
			}
		}
	}
	m_meancolor[0] /= double(nbpix);
	m_meancolor[1] /= double(nbpix);
	m_meancolor[2] /= double(nbpix);
}

template<typename T>
void PCA<T>::computeCovariance(const Mask& mask)
{
	for (int z = 0; z < 3; ++z)
	{
		for (int t = 0; t < 3; ++t)
		{
			m_covar[z][t] = 0;
		}

	}

	int nbpix = 0;
	for (int x = 0; x < m_im.width(); ++x)
	{
		for (int y = 0; y < m_im.height(); ++y)
		{
			if (mask(x,y))
			{
				++nbpix;
				for (int z = 0; z < 3; ++z)
				{
					for (int t = 0; t < 3; ++t)
					{
						m_covar[z][t] += ( m_im.pixelAbsolute(x,y)[z] - m_meancolor[z] ) * ( m_im.pixelAbsolute(x,y)[t] - m_meancolor[t] );
					}
				}
			}
		}
	}

	for (int z = 0; z < 3; ++z)
	{
		for (int t = 0; t < 3; ++t)
		{
			m_covar[z][t] /= double(nbpix);
		}
	}
}

template<typename T>
void PCA<T>::mat_mult(double vin[3], double vout[3])
{
	for (int z = 0; z < 3; ++z)
	{
		vout[z] = m_covar[z][0]*vin[0] + m_covar[z][1]*vin[1] + m_covar[z][2]*vin[2];
	}
}

template<typename T>
double PCA<T>::dot(double v[3], double w[3])
{
	return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

template<typename T>
double PCA<T>::normalize(double v[3])
{
	double norm = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
	v[0] /= norm;
	v[1] /= norm;
	v[2] /= norm;
	return norm;
}

template<typename T>
void PCA<T>::cross(double v[3], double w[3], double res[3])
{
	res[0] = v[1]*w[2] - v[2]*w[1];
	res[1] = v[2]*w[0] - v[0]*w[2];
	res[2] = v[0]*w[1] - v[1]*w[0];
}

template<typename T>
void PCA<T>::computeEigenVectors()
{
	/*** computing first EigenVector v1 = (a1 = theta, a2= phi in polar) ***/
	int i;
	double v [3]{};
	double w [3];
	double a1=0.0, a2=0.0;
	for(i=2;i<20;i++)
	{
		m_v1[0] = std::cos(a1);
		m_v1[1] = std::sin(a1)*std::cos(a2);
		m_v1[2] = std::sin(a1)*std::sin(a2); // cartesian coordinates
		mat_mult(m_v1,v); // matrix-vector mult by covar matrix
		w[0]=-std::sin(a1);
		w[1]=std::cos(a1)*std::cos(a2);
		w[2]=std::cos(a1)*std::sin(a2);
		double d1 = dot(v,w); // dot product
		w[0]=0.0;
		w[1]=-std::sin(a1)*std::sin(a2);
		w[2]=std::sin(a1)*std::cos(a2);
		double d2= dot(v,w); // dot product
		a1+=M_PI*(d1<0.0 ? (-1.0) : (1.0))/pow(2.0,(double)i);                /* adjust ai */
		a2+=M_PI*(d2<0.0 ? (-1.0) : (1.0))/pow(2.0,(double)i);
	}
	m_v1[0]=std::cos(a1); m_v1[1]=std::sin(a1)*std::cos(a2); m_v1[2]=std::sin(a1)*std::sin(a2); // first eigen vector

	/*** computing 2nd & 3rd EigenVectors (a1 = theta in polar in plane (v2,v3) ) ***/
	m_v2[0]=-m_v1[1]; m_v2[1]=m_v1[0]; m_v2[2]=0.0;
	normalize(m_v2);
	cross(m_v1,m_v2,m_v3);
	mat_mult(m_v2,v);
	double k1=-dot(v,m_v2);
	double k2=2.0*dot(v,m_v3);
	mat_mult(m_v3,v);
	k1 += dot(v,m_v3); /*   k1= v3.C.v3-v2.C.v2 ; k2= 2 * v3.C.v2   */
	a1=0.0; // ons single angle around v1
	for(i=2;i<20;i++)
	{
		double d1=std::sin(2.0*a1)*k1+std::cos(2.0*a1)*k2;
		a1 += M_PI*(d1<0.0 ? (-1.0) : (1.0))/pow(2.0,(double)i);
	}
	v[0]=m_v2[0]*std::cos(a1);
	v[1]=m_v2[1]*std::cos(a1);
	v[2]=m_v2[2]*std::cos(a1);
	w[0]=m_v3[0]*std::sin(a1);
	w[1]=m_v3[1]*std::sin(a1);
	w[2]=m_v3[2]*std::sin(a1);
	m_v2[0]=v[0]+w[0];
	m_v2[1]=v[1]+w[1];
	m_v2[2]=v[2]+w[2];
	normalize(m_v2);
	cross(m_v1,m_v2,m_v3);
	//    hvMat3<T> res(v1,v2,v3);
	//    res.transpose(); // TO DO : check with JMD if vectors are rows or columns
}

template<typename T>
void PCA<T>::computePCA (const Mask& m)
{
	computeMeanColor(m);
	computeCovariance(m);
	computeEigenVectors();
}

template<typename T>
void PCA<T>::exportToTXT(std::string savefile)
{
	std::ofstream fichier(savefile, std::ios::out | std::ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert
	fichier << "mean color"<< std::endl;
	fichier << m_meancolor[0] << "\t"<< m_meancolor[1] << "\t"<< m_meancolor[2] << std::endl;
	fichier << "eigen vectors"<< std::endl;
	fichier << m_v1[0] << "\t"<< m_v1[1] << "\t"<< m_v1[2] << std::endl;
	fichier << m_v2[0] << "\t"<< m_v2[1] << "\t"<< m_v2[2] << std::endl;
	fichier << m_v3[0] << "\t"<< m_v3[1] << "\t"<< m_v3[2] << std::endl;
	fichier.close();
}

template<typename T>
void PCA<T>::project(ImageType &res)
{
	res.initItk(m_im.width(), m_im.height());

	for (int x = 0; x < m_im.width(); ++x)
	{
		for (int y = 0; y < m_im.height(); ++y)
		{
			res.pixelAbsolute(x,y)[0] = (m_im.pixelAbsolute(x,y)[0]-m_meancolor[0]) * m_v1[0] + (m_im.pixelAbsolute(x,y)[1]-m_meancolor[1]) * m_v1[1] + (m_im.pixelAbsolute(x,y)[2]-m_meancolor[2]) * m_v1[2];
			res.pixelAbsolute(x,y)[1] = (m_im.pixelAbsolute(x,y)[0]-m_meancolor[0]) * m_v2[0] + (m_im.pixelAbsolute(x,y)[1]-m_meancolor[1]) * m_v2[1] + (m_im.pixelAbsolute(x,y)[2]-m_meancolor[2]) * m_v2[2];
			res.pixelAbsolute(x,y)[2] = (m_im.pixelAbsolute(x,y)[0]-m_meancolor[0]) * m_v3[0] + (m_im.pixelAbsolute(x,y)[1]-m_meancolor[1]) * m_v3[1] + (m_im.pixelAbsolute(x,y)[2]-m_meancolor[2]) * m_v3[2];
		}
	}

	m_A[0][0] = m_v1[0];
	m_A[0][1] = m_v2[0];
	m_A[0][2] = m_v3[0];
	m_A[1][0] = m_v1[1];
	m_A[1][1] = m_v2[1];
	m_A[1][2] = m_v3[1];
	m_A[2][0] = m_v1[2];
	m_A[2][1] = m_v2[2];
	m_A[2][2] = m_v3[2];

	m_det = m_A[0][0]*m_A[1][1]*m_A[2][2] + m_A[0][1]*m_A[1][2]*m_A[2][0] + m_A[0][2]*m_A[1][0]*m_A[2][1] - m_A[0][0]*m_A[1][2]*m_A[2][1] - m_A[0][1]*m_A[1][0]*m_A[2][2] - m_A[0][2]*m_A[1][1]*m_A[2][0];
	m_A_inverse[0][0] = (m_A[1][1] * m_A[2][2] - m_A[2][1] * m_A[1][2]) / m_det;

	m_A_inverse[0][1] = - (m_A[1][0] * m_A[2][2] - m_A[2][0] * m_A[1][2]) / m_det;
	m_A_inverse[1][0] = - (m_A[0][1] * m_A[2][2] - m_A[2][1] * m_A[0][2]) / m_det;

	m_A_inverse[0][2] = (m_A[1][0] * m_A[2][1] - m_A[2][0] * m_A[1][1]) / m_det;
	m_A_inverse[2][0] = (m_A[0][1] * m_A[1][2] - m_A[1][1] * m_A[0][2]) / m_det;

	m_A_inverse[1][1] = (m_A[0][0] * m_A[2][2] - m_A[2][0] * m_A[0][2]) / m_det;

	m_A_inverse[1][2] = - (m_A[0][0] * m_A[2][1] - m_A[2][0] * m_A[0][1]) / m_det;
	m_A_inverse[2][1] = - (m_A[0][0] * m_A[1][2] - m_A[1][0] * m_A[0][2]) / m_det;

	m_A_inverse[2][2] = (m_A[0][0] * m_A[1][1] - m_A[1][0] * m_A[0][1]) / m_det;
}

template<typename T>
void PCA<T>::back_project(const ImageType& coord, ImageType& res) const
{
	res.initItk(coord.width(), coord.height());
	assert(res.width() == coord.width() && res.height() == coord.height());

	for (int x = 0; x < res.width(); ++x)
	{
		for (int y = 0; y < res.height(); ++y)
		{
			res.pixelAbsolute(x,y)[0] = m_meancolor[0] + m_A_inverse[0][0]*coord.pixelAbsolute(x,y)[0] + m_A_inverse[0][1]*coord.pixelAbsolute(x,y)[1] + m_A_inverse[0][2]*coord.pixelAbsolute(x,y)[2];
			res.pixelAbsolute(x,y)[1] = m_meancolor[1] + m_A_inverse[1][0]*coord.pixelAbsolute(x,y)[0] + m_A_inverse[1][1]*coord.pixelAbsolute(x,y)[1] + m_A_inverse[1][2]*coord.pixelAbsolute(x,y)[2];
			res.pixelAbsolute(x,y)[2] = m_meancolor[2] + m_A_inverse[2][0]*coord.pixelAbsolute(x,y)[0] + m_A_inverse[2][1]*coord.pixelAbsolute(x,y)[1] + m_A_inverse[2][2]*coord.pixelAbsolute(x,y)[2];
		}
	}
}

template<typename T>
Eigen::Vector3d PCA<T>::eigenVector(unsigned int i)
{
	Eigen::Vector3d v;
	switch(i)
	{
		case 0:
			std::memcpy(&v, m_v1, sizeof(Eigen::Vector3d));
			break;
		case 1:
			std::memcpy(&v, m_v2, sizeof(Eigen::Vector3d));
			break;
		case 2:
			std::memcpy(&v, m_v3, sizeof(Eigen::Vector3d));
			break;
		default:
			break;
	}
	return v;
}

} //ASTex

#endif
