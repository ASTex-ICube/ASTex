#ifndef __ASTEX__REGION_GROWING_H__
#define __ASTEX__REGION_GROWING_H__




#include "neighborhood.h"
#include <ASTex/image_gray.h>
#include <array>




namespace ASTex
{


class GrowingStrategyRegular
{
public:
	template <typename TPixel>
    inline TPixel operator()( std::vector< TPixel > &processingQueue )
    {
		auto current = processingQueue.back();
        processingQueue.pop_back();
        return current;
    }
};


class GrowingStrategyRandom
{
public:
	template <typename TPixel>
    inline TPixel operator()( std::vector< TPixel > &processingQueue )
    {
		int candidate = rand() % processingQueue.size();
		auto current = processingQueue[candidate];

        processingQueue[candidate] = processingQueue.back();
        processingQueue.pop_back();

        return current;
    }
};


/** Operator for the extraction of regions by a growing algorithm from seeds.
 */
template < typename TNeighborhood >
class RegionGrowingOperator
{
    using TPixel = typename TNeighborhood::PixelType;

    ASTex::ImageGray< bool >        isVisited_;
    std::vector< TPixel >           processingQueue_;
    std::vector< itk::Index<2> >    currentRegionBoundary_;

public:
    /** Creates a region growing operator on the basis of the provided image size.
     *  \param img          Image from which neighborhood queries have to be done.
     */
    template <typename TImage>
    inline RegionGrowingOperator( const TImage &img ) :
        isVisited_( img.width(), img.height() )
    {
        clear();
    }

    /** Creates a region growing operator on the basis of the provided image size.
     *  \param width        Width of the image from which region growing has to be done.
     *  \param height       Height of the image from which region growing has to be done.
     */
    inline RegionGrowingOperator( int width, int height ) :
        isVisited_( width, height )
    {
        clear();
    }

    /** Checks if the operator already marked the given pixel as assigned to a region
     *  during some previous call to the newRegionFromSeed() function.
     *  \param  p           Pixel to check.
     *  \return             True if the pixel has already been assigned, false otherwise. 
     */
    inline bool isAssignedToARegion( const itk::Index<2> &p )
    {
        return isVisited_.pixelAbsolute( p );
    }

    /** Reinitializes the region growing operator by marking all pixels as not assigned to any region.
     */
    inline void clear()
    {
        isVisited_.parallel_for_all_pixels( [](bool &p)
        {
            p = false;
        });
    }

    /** Checks if the operator already marked the given pixel as assigned to a region
     *  during some previous call to the newRegionFromSeed() function.
     *  \tparam TGrowingStrategy        Strategy to use for next pixel selection among all candidate during the growing algorithm.
     *  \param  seed                    Coordinates of the seed pixel from which region growing must start.
     *  \param  growingConditionFunc    Function defining the condition for a pixel to be accepted in the region.
     *  \param  addToRegionFunc         Function defining what to do once a pixel is accepted into the region.
     */
	template <typename TGrowingStrategy,
              typename TGrowingConditionFunc,
              typename TAddToRegionFunc>
    void newRegionFromSeed( const itk::Index<2>         &seed,
                            const TGrowingConditionFunc &growingConditionFunc,
                            const TAddToRegionFunc      &addToRegionFunc )
    {
        // Add the seed to the processing queue.

        processingQueue_.push_back( TPixel(seed) );
        isVisited_.pixelAbsolute( seed ) = true;

		TGrowingStrategy getNextPixel;
        TNeighborhood neighborhood( isVisited_ );
		typename TNeighborhood::NeighborsList neighbors( 4 );

        // Region growing loop.

        do
	    {
            // Take randomly one of the candidate pixels of the processing queue.
		    auto current = getNextPixel( processingQueue_ );

            // Check if this pixel fulfils the given region growing condition.
            if( !growingConditionFunc(current) )
            {
                currentRegionBoundary_.push_back( current );
                continue;
            }

            // If yes, it is added to the region.
            neighborhood.get4( current, neighbors );
            addToRegionFunc( current, neighbors );

		    // Its four neighbors that have not been visited yet are then added to the processing queue.
            for( auto &nn : neighbors )
                if( !isVisited_.pixelAbsolute(nn) )
                {
                    isVisited_.pixelAbsolute(nn) = true;
                    processingQueue_.push_back(nn);
                }
	    } while( !processingQueue_.empty() );

        // Clear the "visited" flag only for pixels that have been visited but not added to the region.

        for( auto &p : currentRegionBoundary_ )
            isVisited_.pixelAbsolute( p ) = false;

        currentRegionBoundary_.clear();
    }
};


} // namespace ASTex




#endif // __ASTEX__REGION_GROWING_H__
