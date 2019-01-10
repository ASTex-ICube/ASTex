#include "sampler.h"

namespace ASTex
{

namespace Stamping
{

std::vector<Eigen::Vector2f> SamplerRegular::generate()
{
	std::vector<Eigen::Vector2f> SamplePoints;
	SamplePoints.reserve(m_nbPoints);
	float step = float(1.0/float(m_nbPoints));
	for (unsigned i = 0; i < m_nbPoints; ++i)
	{
		for (unsigned j = 0; j < m_nbPoints; ++j)
		{
			Eigen::Vector2f tmp_point(float(i)*step,float(j)*step);
			SamplePoints.push_back(tmp_point);
		}
	}

	return SamplePoints;
}

std::vector<Eigen::Vector2f> SamplerOrigin::generate()
{
	std::vector<Eigen::Vector2f> samplePoints;
	samplePoints.reserve(m_nbPoints);
	for(unsigned i=0; i<m_nbPoints; ++i)
		samplePoints.push_back(Eigen::Vector2f(0, 0));
	return samplePoints;
}

std::vector<Eigen::Vector2f> SamplerUniform::generate()
{
	std::vector<Eigen::Vector2f> SamplePoints;
	SamplePoints.reserve(m_nbPoints);
	srand(time(NULL));
	for (unsigned i = 0; i < m_nbPoints; ++i)
	{
		Eigen::Vector2f tmp_point(double(rand())/RAND_MAX, double(rand())/RAND_MAX);
		SamplePoints.push_back(tmp_point);
	}
	return SamplePoints;
}

std::vector<Eigen::Vector2f> SamplerPoisson::generate()
{
	float minDistance=m_minDistance;
	if ( minDistance < 0.0f )
	{
		minDistance = sqrt( float(m_nbPoints) ) / m_nbPoints;
	}

	std::vector<Eigen::Vector2f> SamplePoints;
	std::vector<Eigen::Vector2f> ProcessList;

	// create the grid
	float CellSize = minDistance / sqrt( 2.0f );

	int GridW = ( int )ceil( 1.0f / CellSize );
	int GridH = ( int )ceil( 1.0f / CellSize );

	sGrid Grid( GridW, GridH, CellSize );

	while(SamplePoints.size()<m_nbPoints)
	{
		Eigen::Vector2f FirstPoint;
		do {
			FirstPoint = Eigen::Vector2f( m_generator.RandomFloat(), m_generator.RandomFloat() );
		} while (!(m_generateInCircle ? IsInCircle(FirstPoint) : IsInRectangle(FirstPoint)));

		// update containers
		ProcessList.push_back( FirstPoint );
		SamplePoints.push_back( FirstPoint );
		Grid.Insert( FirstPoint );

		// generate new points for each point in the queue
		while ( !ProcessList.empty() && SamplePoints.size() < m_nbPoints )
		{
			Eigen::Vector2f Point = PopRandom( ProcessList, m_generator );

			for ( int i = 0; i < m_newPointsCount; i++ )
			{
				Eigen::Vector2f NewPoint = GenerateRandomPointAround( Point, minDistance, m_generator );

				bool Fits = m_generateInCircle ? IsInCircle(NewPoint) : IsInRectangle(NewPoint);

				if ( Fits && !Grid.IsInNeighbourhood( NewPoint, minDistance, CellSize ) )
				{
					ProcessList.push_back( NewPoint );
					SamplePoints.push_back( NewPoint );
					Grid.Insert( NewPoint );
					continue;
				}
			}
		}
	}

	return SamplePoints;
}

} //namsepace Stamping

} //namespace ASTex
