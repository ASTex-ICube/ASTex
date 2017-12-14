#ifndef _TEXTON_CPU_H_
#define _TEXTON_CPU_H_

#include "vector"

#include <vector>
#include <iostream>
#include <fstream>
#include <memory.h>

#include <cmath>
#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>
#include <ASTex/fourier.h>
#include <ASTex/local_spectrum.h>
#include <ASTex/utils.h>
#include <ASTex/distances_maps.h>
#include <Eigen/Core>

namespace ASTex
{

typedef Eigen::Vector2f vec2;

// Geoffrey sampling class 

class PointSet
{
public:
    std::vector<vec2> Generate(size_t);
};

class RegularSampling: public PointSet
{
    public:
        RegularSampling() {;}

        /**
            Return a vector of generated points

            NewPointsCount - refer to bridson-siggraph07-poissondisk.pdf for details (the value 'k')
            Circle  - 'true' to fill a circle, 'false' to fill a rectangle
            MinDist - minimal distance estimator, use negative value for default
        **/
        std::vector<vec2> Generate(size_t NumPoints)
        {
            std::vector<vec2> SamplePoints;
            SamplePoints.reserve(NumPoints);
            float step = float(1.0/float(NumPoints));
            for (int i = 0; i < NumPoints; ++i) 
            {
                for (int j = 0; j < NumPoints; ++j) 
                {
                    vec2 tmp_point(float(i)*step,float(j)*step);
                    SamplePoints.push_back(tmp_point);
                }
            }

            return SamplePoints;
        }
};

class RandomSampling: public PointSet
{
    public:
        RandomSampling() {;}

        /**
            Return a vector of generated points

            NewPointsCount - refer to bridson-siggraph07-poissondisk.pdf for details (the value 'k')
            Circle  - 'true' to fill a circle, 'false' to fill a rectangle
            MinDist - minimal distance estimator, use negative value for default
        **/
        std::vector<vec2> Generate(size_t NumPoints)
        {
            std::vector<vec2> SamplePoints;
            SamplePoints.reserve(NumPoints);
            srand(time(NULL));
            for (int i = 0; i < NumPoints; ++i) 
            {
                    vec2 tmp_point((rand()%100)/100.f ,(rand()%100)/100.f );
                    SamplePoints.push_back(tmp_point);
            }
            return SamplePoints;
        }
};

class PoissonSampling: public PointSet
{
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

    DefaultPRNG Generator;

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
        sGridPoint ImageToGrid( const vec2& P, float CellSize )
        {
            return sGridPoint( ( int )( P[0] / CellSize ), ( int )( P[1] / CellSize ) );
        }
        void Insert( const vec2& P )
        {
            sGridPoint G = ImageToGrid( P, m_CellSize );
            m_Grid[ G.x ][ G.y ] = P;
        }
        bool IsInNeighbourhood( vec2 Point, float MinDist, float CellSize )
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
                        vec2 P = m_Grid[ i ][ j ];

                        vec2 np = P-Point;
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

        std::vector< std::vector<vec2> > m_Grid;
    };

    bool IsInRectangle(vec2 p) const
    {
        return p[0] >= 0 && p[1] >= 0 && p[0] <= 1 && p[1] <= 1;
    }
    //
    bool IsInCircle(vec2 p) const
    {
        float fx = p[0] - 0.5f;
        float fy = p[1] - 0.5f;
        return ( fx*fx + fy*fy ) <= 0.25f;
    }

public:
    PoissonSampling() {
        Generator = DefaultPRNG();
    }

    template <typename PRNG>
    vec2 PopRandom( std::vector<vec2>& Points, PRNG& Generator )
    {
        const int Idx = Generator.RandomInt( Points.size()-1 );
        const vec2 P = Points[ Idx ];
        Points.erase( Points.begin() + Idx );
        return P;
    }

    template <typename PRNG>
    vec2 GenerateRandomPointAround( const vec2& P, float MinDist, PRNG& Generator )
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

        return vec2( X, Y );
    }

    /**
        Return a vector of generated points

        NewPointsCount - refer to bridson-siggraph07-poissondisk.pdf for details (the value 'k')
        Circle  - 'true' to fill a circle, 'false' to fill a rectangle
        MinDist - minimal distance estimator, use negative value for default
    **/
    std::vector<vec2> Generate(
        size_t NumPoints,
        int NewPointsCount = 30,
        bool Circle = false,
        float MinDist = -1.0f
    )
    {
        if ( MinDist < 0.0f )
        {
            MinDist = sqrt( float(NumPoints) ) / float(NumPoints);
        }

        std::vector<vec2> SamplePoints;
        std::vector<vec2> ProcessList;

        // create the grid
        float CellSize = MinDist / sqrt( 2.0f );

        int GridW = ( int )ceil( 1.0f / CellSize );
        int GridH = ( int )ceil( 1.0f / CellSize );

        sGrid Grid( GridW, GridH, CellSize );

        vec2 FirstPoint;
        do {
            FirstPoint = vec2( Generator.RandomFloat(), Generator.RandomFloat() );
        } while (!(Circle ? IsInCircle(FirstPoint) : IsInRectangle(FirstPoint)));

        // update containers
        ProcessList.push_back( FirstPoint );
        SamplePoints.push_back( FirstPoint );
        Grid.Insert( FirstPoint );

        // generate new points for each point in the queue
        while ( !ProcessList.empty() && SamplePoints.size() < NumPoints )
        {
    #if POISSON_PROGRESS_INDICATOR
            // a progress indicator, kind of
            if ( SamplePoints.size() % 100 == 0 ) std::cout << ".";
    #endif // POISSON_PROGRESS_INDICATOR

            vec2 Point = PopRandom<DefaultPRNG>( ProcessList, Generator );

            for ( int i = 0; i < NewPointsCount; i++ )
            {
                vec2 NewPoint = GenerateRandomPointAround( Point, MinDist, Generator );

                bool Fits = Circle ? IsInCircle(NewPoint) : IsInRectangle(NewPoint);

                if ( Fits && !Grid.IsInNeighbourhood( NewPoint, MinDist, CellSize ) )
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
};

//////////////STAMPER///////////////


class StampBase
{
public:
    StampBase();

    virtual ImageRGBd::PixelType pixelAbsolute(int x, int y)=0;

protected:
    size_t m_width;
    size_t m_height;
};

class StampDiscrete : public StampBase
{
public:
    StampDiscrete();
    StampDiscrete(const ImageRGBd &stamp);

    ImageRGBd::PixelType pixelAbsolute(int x, int y);

private:

    ImageRGBd m_stamp;
};

//class StampFunctionalInt : public StampBase
//{
//public:

//    StampFunctionalInt();
//    StampFunctionalInt(const function<ImageRGBd::PixelType (int x, int y)> &f);

//    ImageRGBd::PixelType pixelAbsolute(int x, int y);

//    void setWidth(size_t width) {m_width=width;}
//    void setHeight(size_t height) {m_height=height;}

//private:
//    function<ImageRGBd::PixelType (int x, int y)> m_f;
//};

class Stamper
{
    protected:
        std::vector<vec2> m_pointArray;

        ImageRGBd m_stamp;

    public:

        Stamper(const std::vector<vec2> &pointArray, const ImageRGBd &tampon);

        virtual ImageRGBd generate(int imageWidth, int imageHeight) = 0;
};

class BombingStamper : public Stamper
{
public:

    BombingStamper(const std::vector<vec2> &pointArray, const ImageRGBd &tampon);

    ImageRGBd generate(int imageWidth, int imageHeight);
};

class TextonStamper : public Stamper
{
public:

    TextonStamper(const std::vector<vec2> &pointArray, const ImageRGBd &tampon);

    ImageRGBd generate(int imageWidth, int imageHeight);

    //get

    double ratioX() {return m_ratioX;}
    double ratioY() {return m_ratioY;}
    bool periodicity() {return m_periodicity;}
    bool bilinearInterpolation() {return m_bilinearInterpolation;}

    //set

    void setRatioX(double ratioX) {m_ratioX=ratioX; assert(ratioX>0 && "TextonStamper::setRatioX: ratioX must be > 0");}
    void setRatioY(double ratioY) {m_ratioY=ratioY; assert(ratioY>0 && "TextonStamper::setRatioY: ratioY must be > 0");}
    void setPeriodicity(bool periodicity) {m_periodicity = periodicity;}
    void setBilinearInterpolation(bool bi) {m_bilinearInterpolation = bi;}

private:

    double m_ratioX; //< these define the resolution of the texton relative to the resolution of the texture (in X and Y)
    double m_ratioY;

    bool m_periodicity; //< if the stamping process is periodic or not.

    bool m_bilinearInterpolation; //< if the stamping process is allowed to stamp in between pixels or not.
};

} //namespace ASTex


#endif //_TEXTON_CPU_H_
