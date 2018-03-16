#include "sampler.h"

namespace ASTex
{

namespace Stamping
{

std::vector<Eigen::Vector2f> RegularSampler::generate()
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

std::vector<Eigen::Vector2f> UniformSampler::generate()
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

std::vector<Eigen::Vector2f> PoissonSampler::generate()
{
    if ( m_minDistance < 0.0f )
    {
        m_minDistance = sqrt( float(m_nbPoints) ) / float(m_nbPoints);
    }

    std::vector<Eigen::Vector2f> SamplePoints;
    std::vector<Eigen::Vector2f> ProcessList;

    // create the grid
    float CellSize = m_minDistance / sqrt( 2.0f );

    int GridW = ( int )ceil( 1.0f / CellSize );
    int GridH = ( int )ceil( 1.0f / CellSize );

    sGrid Grid( GridW, GridH, CellSize );

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
#if POISSON_PROGRESS_INDICATOR
        // a progress indicator, kind of
        if ( SamplePoints.size() % 100 == 0 ) std::cout << ".";
#endif // POISSON_PROGRESS_INDICATOR

        Eigen::Vector2f Point = PopRandom( ProcessList, m_generator );

        for ( int i = 0; i < m_newPointsCount; i++ )
        {
            Eigen::Vector2f NewPoint = GenerateRandomPointAround( Point, m_minDistance, m_generator );

            bool Fits = m_generateInCircle ? IsInCircle(NewPoint) : IsInRectangle(NewPoint);

            if ( Fits && !Grid.IsInNeighbourhood( NewPoint, m_minDistance, CellSize ) )
            {
                ProcessList.push_back( NewPoint );
                SamplePoints.push_back( NewPoint );
                Grid.Insert( NewPoint );
                continue;
            }
        }
    }

#if POISSON_PROGRESS_INDICATOR
    std::cout << std::endl << std::endl;
#endif // POISSON_PROGRESS_INDICATOR

    return SamplePoints;
}

} //namsepace Stamping

} //namespace ASTex
