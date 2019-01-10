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



#include "old_PoissonDiskSampler.h"


namespace ASTex
{

namespace ContentExchg
{


NeighborhoodSearchGrid::NeighborhoodSearchGrid( int width, int height, double minDist ) :
	width_( width ),
	height_( height ),
	minDistSq_( minDist*minDist )
{
	cellSize_.setConstant( minDist / std::sqrt(2.0) );

	cells_ = ASTex::ImageGray32( (int) std::ceil(width  / cellSize_[0]),
								 (int) std::ceil(height / cellSize_[1]) );

	cellSize_[0] = double(width ) / cells_.width ();
	cellSize_[1] = double(height) / cells_.height();

	clear();
}



bool NeighborhoodSearchGrid::isNeighborhoodFree( const Eigen::Vector2d &coords ) const
{
	itk::Index<2> cell = gridCoord( coords );

	for( int y=-2; y<=2; y++ )
		for( int x=-2; x<=2; x++ )
		{
			itk::Index<2> neighb;
			neighb[0] = cell[0] + x;
			neighb[1] = cell[1] + y;

			auto point = get( neighb );
			if( !point )
				continue;

			Eigen::Vector2d delta( 0.0, 0.0 );

			if( neighb[0] < 0 )
				delta[0] = -width_;
			else if( neighb[0] >= cells_.width() )
				delta[0] =  width_;

			if( neighb[1] < 0 )
				delta[1] = -height_;
			else if( neighb[1] >= cells_.height() )
				delta[1] =  height_;

			double sqDist = (coords - *point - delta).squaredNorm();
			if( sqDist < minDistSq_ )
				return false;
		}

	return true;
}






void PoissonDiskSampler::generateSamples( int width, int height, std::vector< Eigen::Vector2d > &samples, double minDist, int nNewPoints )
{
	// Creation of the grid used for neighborhood query.

	NeighborhoodSearchGrid grid( width, height, minDist );

	std::vector<Point> processingList;


	// Creation of the first point.

	samples.clear();
	samples.push_back( Point( double(width ) * rand() / RAND_MAX,
							  double(height) * rand() / RAND_MAX ) );

	grid.put( samples.back() );

	processingList.push_back( samples.back() );


	while( !processingList.empty() )
	{
		// Extract a random candidate from the processing list.

		int candidateId = rand() % processingList.size();
		Point candidate = processingList[candidateId];

		processingList[candidateId] = processingList.back();
		processingList.pop_back();

		//

		for( int i=0; i<nNewPoints; ++i )
		{
			Point newPoint = generateRandomPointAround( candidate, minDist );

			if( newPoint[0] < 0.0 )
				newPoint[0] += width;
			else if( newPoint[0] >= width )
				newPoint[0] -= width;

			if( newPoint[1] < 0.0 )
				newPoint[1] += height;
			else if( newPoint[1] >= height )
				newPoint[1] -= height;

			if( grid.isNeighborhoodFree( newPoint ) )
			{
				grid.put( newPoint );
				processingList.push_back( newPoint );
				samples.push_back( newPoint );
			}
		}
	}
}


Point PoissonDiskSampler::generateRandomPointAround( Eigen::Vector2d &point, double minDist )
{
	double rnd1 = double(rand()) / RAND_MAX;
	double rnd2 = double(rand()) / RAND_MAX;

	double radius = minDist * (rnd1 + 1.0);

	double angle = 2.0 * M_PI * rnd2;

	return point + radius * Point( std::cos(angle), std::sin(angle) );
}

}
}
