#ifndef __SAMPLER__H__
#define __SAMPLER__H__

#include <vector>
#include <iostream>
#include <fstream>
#include <memory.h>

#include <cmath>
#include <Eigen/Core>

namespace ASTex
{

namespace Stamping
{

class SamplerBase
{
public:

    SamplerBase() {}

    virtual std::vector<Eigen::Vector2f> generate()=0;
};

class RegularSampler : public SamplerBase
{
public:
    RegularSampler(int nbPoints=0) :
        SamplerBase(),
        m_nbPoints(nbPoints)
    {}

    void setNbPoints(unsigned nbPoints) {m_nbPoints = nbPoints;}

    std::vector<Eigen::Vector2f> generate();

private:
    unsigned int m_nbPoints;
};

class UniformSampler : public SamplerBase
{
public:
    UniformSampler(int nbPoints=0) :
        SamplerBase(),
        m_nbPoints(nbPoints)
    {}

    void setNbPoints(unsigned nbPoints) {m_nbPoints = nbPoints;}

    std::vector<Eigen::Vector2f> generate();

private:
    unsigned int m_nbPoints;
};

class PoissonSampler: public SamplerBase
{
public:
    PoissonSampler(int nbPoints=0) :
        SamplerBase(),
        m_nbPoints(nbPoints),
        m_generator(),
        m_newPointsCount(30),
        m_generateInCircle(false),
        m_minDistance(-1.0f)
    {}

    /**
        Return a vector of generated points

        NewPointsCount - refer to bridson-siggraph07-poissondisk.pdf for details (the value 'k')
        Circle  - 'true' to fill a circle, 'false' to fill a rectangle
        MinDist - minimal distance estimator, use negative value for default
    **/

    void setNbPoints(unsigned nbPoints) {m_nbPoints = nbPoints;}
    void setNewPointsCount(int count) {m_newPointsCount = std::max(0, count);}
    void setGenerateInCircle(bool b) {m_generateInCircle = b;}
    void setMinDistance(float d) {m_minDistance = d;}

    std::vector<Eigen::Vector2f> generate()
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

private:
    class DefaultPRNG
    {
    public:
        DefaultPRNG()
        : m_Gen( std::random_device()() )
        , m_Dis( 0.0f, 1.0f )
        {
            // prepare PRNG
            m_Gen.seed( time( nullptr ) );
        }

        explicit DefaultPRNG( uint32_t seed )
        : m_Gen( seed )
        , m_Dis( 0.0f, 1.0f )
        {
        }

        float RandomFloat()
        {
            return static_cast<float>( m_Dis( m_Gen ) );
        }

        int RandomInt( int Max )
        {
            std::uniform_int_distribution<> DisInt( 0, Max );
            return DisInt( m_Gen );
        }

    private:
        std::mt19937 m_Gen;
        std::uniform_real_distribution<float> m_Dis;
    };

    struct sGridPoint
    {
        sGridPoint( int X, int Y )
            : x( X )
            , y( Y )
        {}
        int x;
        int y;
    };


    struct sGrid
    {
        sGrid( int W, int H, float CellSize )
            : m_W( W )
            , m_H( H )
            , m_CellSize( CellSize )
        {
            m_Grid.resize( m_H );

            for ( auto i = m_Grid.begin(); i != m_Grid.end(); i++ ) { i->resize( m_W ); }
        }
        sGridPoint ImageToGrid( const Eigen::Vector2f& P, float CellSize )
        {
            return sGridPoint( ( int )( P[0] / CellSize ), ( int )( P[1] / CellSize ) );
        }
        void Insert( const Eigen::Vector2f& P )
        {
            sGridPoint G = ImageToGrid( P, m_CellSize );
            m_Grid[ G.x ][ G.y ] = P;
        }
        bool IsInNeighbourhood( Eigen::Vector2f Point, float MinDist, float CellSize )
        {
            sGridPoint G = ImageToGrid( Point, CellSize );

            // number of adjucent cells to look for neighbour points
            const int D = 5;

            // scan the neighbourhood of the point in the grid
            for ( int i = G.x - D; i < G.x + D; i++ )
            {
                for ( int j = G.y - D; j < G.y + D; j++ )
                {
                    if ( i >= 0 && i < m_W && j >= 0 && j < m_H )
                    {
                        Eigen::Vector2f P = m_Grid[ i ][ j ];

                        Eigen::Vector2f np = P-Point;
                        float X = np[0];
                        float Y = np[1];
                        float d = std::sqrt( X*X + Y*Y );

                        if ( d < MinDist ) { return true; }
                    }
                }
            }


            return false;
        }

    private:
        int m_W;
        int m_H;
        float m_CellSize;

        std::vector< std::vector<Eigen::Vector2f> > m_Grid;
    };

    bool IsInRectangle(Eigen::Vector2f p) const
    {
        return p[0] >= 0 && p[1] >= 0 && p[0] <= 1 && p[1] <= 1;
    }
    //
    bool IsInCircle(Eigen::Vector2f p) const
    {
        float fx = p[0] - 0.5f;
        float fy = p[1] - 0.5f;
        return ( fx*fx + fy*fy ) <= 0.25f;
    }

    Eigen::Vector2f PopRandom( std::vector<Eigen::Vector2f>& Points, DefaultPRNG& Generator )
    {
        const int Idx = Generator.RandomInt( Points.size()-1 );
        const Eigen::Vector2f P = Points[ Idx ];
        Points.erase( Points.begin() + Idx );
        return P;
    }

    Eigen::Vector2f GenerateRandomPointAround( const Eigen::Vector2f& P, float MinDist, DefaultPRNG& Generator )
    {
        // start with non-uniform distribution
        float R1 = Generator.RandomFloat();
        float R2 = Generator.RandomFloat();

        // radius should be between MinDist and 2 * MinDist
        float Radius = MinDist * ( R1 + 1.0f );

        // random angle
        float Angle = 2 * 3.141592653589f * R2;

        // the new point is generated around the point (x, y)
        float X = P[0] + Radius * cos( Angle );
        float Y = P[1] + Radius * sin( Angle );

        return Eigen::Vector2f( X, Y );
    }

    int m_nbPoints;             ///< Number of points to be generated at the end of the process. No default, must be set.
    DefaultPRNG m_generator;    ///< The random number generator. The user can't change that, as we aren't preoccupied with low-level considerations.
    int m_newPointsCount;       ///< Defines the number of points the process tries to pick in a single loop... 30 by default. TODO: refer to bridson-siggraph07-poissondisk.pdf for details (the value 'k')
    int m_generateInCircle;     ///< Defines whether you want to generate the points inside circles of diameter 1 or a square of dimension 1x1. Off by default. TODO: check with circle on.
    float m_minDistance;        ///< Defines the minimum distance between generated points. If lower than 0, the process determines an optimal distance automatically. -1.0f by default.
};

} //namespace Stamping

} //namespace ASTex


#endif //__SAMPLER__H__
