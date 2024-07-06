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



#ifndef __CONTENTEXCHG__POISSONDISKSAMPLER__
#define __CONTENTEXCHG__POISSONDISKSAMPLER__




#include <ASTex/image_gray.h>
#include "Eigen/Core"


namespace ASTex
{

namespace ContentExchg
{


using Point = Eigen::Vector2d;


class NeighborhoodSearchGrid
{
private:
	Eigen::Vector2d                 cellSize_;
	ASTex::ImageGray32              cells_;
	int                             width_;
	int                             height_;
	double                          minDistSq_;
	std::vector<Eigen::Vector2d>    points_;

public:
	NeighborhoodSearchGrid( int width, int height, double minDist );

	void                     put( const Eigen::Vector2d &coords );
	void                     put( itk::Index<2> cell, const Eigen::Vector2d &coords );

	const Eigen::Vector2d*   get( const Eigen::Vector2d &coords ) const;
	const Eigen::Vector2d*   get( itk::Index<2> cell ) const;

	void                     clear();

	bool                     isNeighborhoodFree( const Eigen::Vector2d &coords ) const;

	itk::Index<2>            gridCoord( const Eigen::Vector2d &coords ) const;
};





class PoissonDiskSampler
{
	static Point generateRandomPointAround( Eigen::Vector2d &point, double minDist );

public:
	/** Generates N seeds over the domain [0,width[x[0,height[ with a Poisson disk sampling.
	 *
	 *  \param width        Domain size along the X axis.
	 *  \param height       Domain size along the Y axis.
	 *  \param samples      Array where to store generated samples.
	 *  \param distMin      Minimum distance allowed between samples.
	 *  \param nNewPoints   Number of new points to generate and check around valid samples for the generation of new ones.
	 */
	static void generateSamples( int width, int height, std::vector< Eigen::Vector2d > &samples, double minDist, int nNewPoints = 100 );
};




inline void NeighborhoodSearchGrid::clear()
{
	points_.clear();
	cells_.for_all_pixels( [] (ASTex::ImageGray32::PixelType &p){ p = -1; });
}

inline itk::Index<2> NeighborhoodSearchGrid::gridCoord( const Eigen::Vector2d &coords ) const
{
	itk::Index<2> cellCoord;
	cellCoord[0] = (int)(coords[0] / cellSize_[0]);
	cellCoord[1] = (int)(coords[1] / cellSize_[1]);
	return cellCoord;
}

inline void NeighborhoodSearchGrid::put( const Eigen::Vector2d &coords )
{
	put( gridCoord( coords ), coords );
}


inline void NeighborhoodSearchGrid::put( itk::Index<2> cell, const Eigen::Vector2d &coords )
{
	cell[0] = (cell[0] + cells_.width ()) % cells_.width ();
	cell[1] = (cell[1] + cells_.height()) % cells_.height();

	cells_.pixelAbsolute( cell ) = (int32_t) points_.size();
	points_.push_back( coords );
}


inline const Eigen::Vector2d* NeighborhoodSearchGrid::get( const Eigen::Vector2d &coords ) const
{
	return get( gridCoord( coords ) );
}


inline const Eigen::Vector2d* NeighborhoodSearchGrid::get( itk::Index<2> cell ) const
{
	cell[0] = (cell[0] + cells_.width ()) % cells_.width ();
	cell[1] = (cell[1] + cells_.height()) % cells_.height();

	int32_t pointId = cells_.pixelAbsolute( cell );
	return (pointId >= 0)?  &points_[pointId]  :  NULL;
}




} // namespace ContentExchg


} // ASTex

#endif // __CONTENTEXCHG__POISSONDISKSAMPLER__
