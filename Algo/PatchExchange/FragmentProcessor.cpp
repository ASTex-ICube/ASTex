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



#include "FragmentProcessor.h"


namespace ASTex
{

namespace ContentExchg
{


FragmentProcessor::FragmentProcessor( const ASTex::ImageRGBu8 &image ) :
	isFragmented_( false ),
	sourceImage_( image ),
	idMap_( image.width(), image.height() )
{
}


float FragmentProcessor::colorSquareDist( const ColorF &c1, const ColorF &c2 )
{
    ColorF diffC = c1 - c2;
    return diffC.GetRed  ()*diffC.GetRed  () +
           diffC.GetGreen()*diffC.GetGreen() +
           diffC.GetBlue ()*diffC.GetBlue ();
}


int FragmentProcessor::fragmentCount() const
{
	return (int) allFragments_.size();
}


int FragmentProcessor::fragmentIdAt( int x, int y ) const
{
    assert( isFragmented_ );
	return idMap_.pixelAbsolute(x,y);
}


const Fragment& FragmentProcessor::fragmentById( int fid ) const
{
    assert( fid < (int) allFragments_.size() );
	return allFragments_[fid];
}


const Fragment& FragmentProcessor::fragmentAt( int x, int y ) const
{
    return fragmentById( fragmentIdAt(x,y) );
}


const std::vector<Fragment>& FragmentProcessor::fragments() const
{
    return allFragments_;
}


const ASTex::ImageGrayu32& FragmentProcessor::idMap() const
{
    return idMap_;
}


unsigned int FragmentProcessor::createFragments( unsigned int fragmentMaxSize, unsigned int threshold, ASTex::ImageGrayu8 *labelMap )
{
    allFragments_.clear();
    isFragmented_ = false;
    const float tolerance = threshold * threshold;

	// Fragmentation loop.

    ASTex::RegionGrowingOperator< ASTex::NeighborhoodPeriodic > regionOp( sourceImage_ );
    auto seed = idMap_.beginConstIterator();

	do
	{
        // Create a new fragment.

        allFragments_.push_back( Fragment() );
	    Fragment &currentFragment = allFragments_.back();

        currentFragment.id = (uint32_t) allFragments_.size() - 1;

        unsigned char currentFragmentLabel = 0;
	    if( labelMap )
            currentFragmentLabel = labelMap->pixelAbsolute( seed.GetIndex() );

        ColorU32 fragmentColorSum( 0U );

        // Build the fragment by region growing.

        regionOp.newRegionFromSeed< ASTex::GrowingStrategyRandom >(
            seed.GetIndex(),

            // Condition function for accepting a pixel to the current region.
            [&]( const ASTex::NeighborhoodPeriodic::PixelType &pixel )
            {
                // Check fragment size.
                if( currentFragment.pixels.size() == fragmentMaxSize )
                    return false;

                // Check pixel label.
	    		if( labelMap  &&  labelMap->pixelAbsolute(pixel) != currentFragmentLabel )
                    return false;

                // Check color divergence.
                ColorU8  col = sourceImage_.pixelAbsolute( pixel );
                ColorU32 updatedColorSum = fragmentColorSum + col;
                ColorF   meanColor = ColorF(updatedColorSum) * (1.0f / (currentFragment.pixels.size()+1));

				float sqDistToMean = FragmentProcessor::colorSquareDist( ColorF(col), meanColor );
                return sqDistToMean < tolerance;
            },

            // Processing function applied to accepted pixels.
			[&]( const ASTex::NeighborhoodPeriodic::PixelType &pixel, const ASTex::NeighborhoodPeriodic::NeighborsList& /*neighbors*/ )
            {
                idMap_.pixelAbsolute( pixel ) = currentFragment.id;
			    currentFragment.pixels.push_back( pixel );
                fragmentColorSum += sourceImage_.pixelAbsolute( pixel );
            }
        );

        // Pick the next unassigned pixel in the image and use it as seed for the new fragment.

        while( !seed.IsAtEnd()  &&  regionOp.isAssignedToARegion( seed.GetIndex() ) )
            ++ seed;
    } while( !seed.IsAtEnd() );

    isFragmented_ = true;

	return (int) allFragments_.size();
}


unsigned int FragmentProcessor::cleanupSmallestFragments( unsigned int fragmentMinSize )
{
	if( !isFragmented_ )
		return 0;

	// Create a list ordered by increasing size of all fragments which are smaller than the given threshold.

    struct FragmentInfo
    {
        uint32_t    id;
        size_t      size;
        inline bool operator<( const FragmentInfo &fi ) const
        {
            return size < fi.size  ||  (size == fi.size  &&  id < fi.id);
        }
    };

    std::set< FragmentInfo > orderedFragments;

    for( auto &frag : allFragments_ )
        if( frag.pixels.size() < fragmentMinSize )
	    {
		    FragmentInfo fragmentInfo;
            fragmentInfo.id   = frag.id;
            fragmentInfo.size = frag.pixels.size();
    	    orderedFragments.insert( fragmentInfo );
	    }

    // Process this list by starting with the smallest fragment.

    ASTex::NeighborhoodPeriodic::NeighborsList neighbors( 4 );
    ASTex::NeighborhoodPeriodic neighborhood( sourceImage_ );

    while( !orderedFragments.empty() )
	{
        Fragment &fragment = allFragments_[ orderedFragments.begin()->id ];
		orderedFragments.erase( orderedFragments.begin() );

		// Create a queue containing all its pixels.

		std::list<PixelPos> pixels;
		for( auto &p : fragment.pixels )
			pixels.push_back( p );

		// Process the queue until it's empty.

		while( !pixels.empty() )
		{
            // Get the first pixel in the queue.

            PixelPos pixel = pixels.front();
			pixels.pop_front();

			// Recover its 4-neighbors that belong to different fragments.

            neighborhood.get4( pixel, neighbors );
            int nNeighbors = 4;

            for( int i=0; i<nNeighbors; )
                if( idMap_.pixelAbsolute(neighbors[i]) != idMap_.pixelAbsolute(pixel) )
                    ++ i;
                else
                    neighbors[i] = neighbors[--nNeighbors];

			// If no such neighbor exists, the pixel is put at the end of the queue for futur processing.

			if( nNeighbors == 0 )
            {
				pixels.push_back( pixel );
                continue;
            }

            // Otherwise, the neighbor which is the closest to the the current pixel (in terms of Euclidean distance in RGB space) is looked for.

			ColorF pixelColor = sourceImage_.pixelAbsolute( pixel );

			float minColorDist = std::numeric_limits<float>::max();
			PixelPos bestNeighbor;

            for( int i=0; i<nNeighbors; ++i )
			{
                ColorF neighborColor = sourceImage_.pixelAbsolute( neighbors[i] );
                float colorSqDist = colorSquareDist( neighborColor, pixelColor );

				if( colorSqDist < minColorDist )
				{
					minColorDist = colorSqDist;
					bestNeighbor = neighbors[i];
				}
			}

            // The fragment this neighbor belongs to is the one to which the current pixel must be added.

            uint32_t bestNeighborId = idMap_.pixelAbsolute( bestNeighbor );
            Fragment &bestNeighborFragment = allFragments_[ bestNeighborId ];

			bestNeighborFragment.pixels.push_back( pixel );
			idMap_.pixelAbsolute( pixel ) = bestNeighborId;

            // Update the increasing size ordered list of fragments, if needed.

            FragmentInfo neighborFragmentInfo;
            neighborFragmentInfo.id = bestNeighborId;
			neighborFragmentInfo.size = (unsigned int) bestNeighborFragment.pixels.size() - 1;

            if( neighborFragmentInfo.size < fragmentMinSize )
            {
                orderedFragments.erase( neighborFragmentInfo );

                neighborFragmentInfo.size ++;
                if( neighborFragmentInfo.size < fragmentMinSize )
                    orderedFragments.insert( neighborFragmentInfo );
            }
		}

		// The fragment which has been emptied is eventually cleared and removed from the ordered list.

        fragment.pixels.clear();
    }

	// ID map is updated according to the new list of fragments.

    int beg = 0;
    int end = allFragments_.size() - 1;

    do
    {
        if( !allFragments_[beg].pixels.empty() )
            ++ beg;
        else if( allFragments_[end].pixels.empty() )
            -- end;
        else
        {
            allFragments_[beg] = allFragments_[end];
            allFragments_[beg].id = beg;

            for( auto &p : allFragments_[beg].pixels )
                idMap_.pixelAbsolute(p) = beg;

            ++ beg;
            -- end;
        }
    } while( beg < end );

    allFragments_.resize( beg );

	return (unsigned int) allFragments_.size();
}


void FragmentProcessor::updateFragmentsAttributes()
{
    std::vector<bool> isNeighborOfCurrentFragment;
    isNeighborOfCurrentFragment.assign( allFragments_.size(), false );

    // Computes the remaining attributes of fragments, namely their lists of neighboring fragments and their centroids.

    ASTex::RegionGrowingOperator< ASTex::NeighborhoodPeriodicRelative > regionOp( sourceImage_ );

    for( auto &currentFragment : allFragments_ )
    {
        // Initialize the centroid and the set of neighboring fragments.

        currentFragment.centroid.setConstant( 0.0 );
        currentFragment.neighbors.clear();

        // Perform the seed fill.

        regionOp.newRegionFromSeed< ASTex::GrowingStrategyRegular >(
            currentFragment.pixels.front(),

            // Condition function for accepting a pixel to the current region.
            [&]( const ASTex::NeighborhoodPeriodicRelative::PixelType &pixel )
            {
                return idMap_.pixelAbsolute( pixel.absolute ) == currentFragment.id;
            },

            // Processing function applied to accepted pixels.
            [&]( const ASTex::NeighborhoodPeriodicRelative::PixelType &pixel, const ASTex::NeighborhoodPeriodicRelative::NeighborsList &neighbors )
            {
                currentFragment.centroid += Eigen::Vector2d( pixel.relative[0], pixel.relative[1] );

                for( auto &nn : neighbors )
                {
                    auto &neighbFragmentId = idMap_.pixelAbsolute( nn );
                    if( neighbFragmentId != currentFragment.id  &&  !isNeighborOfCurrentFragment[neighbFragmentId] )
                    {
                        isNeighborOfCurrentFragment[neighbFragmentId] = true;
                        currentFragment.neighbors.push_back( neighbFragmentId );
                    }
                }
            }
        );

        // Finalize the centroid computation.

        currentFragment.centroid /= (double) currentFragment.pixels.size();

        for( int i=0; i<2; ++i )
            if( currentFragment.centroid[i] < 0.0 )
                currentFragment.centroid[i] += sourceImage_.size()[i];
            else if( currentFragment.centroid[i] >= sourceImage_.size()[i] )
                currentFragment.centroid[i] -= sourceImage_.size()[i];

        // Clear the "isNeighbor" flag.

        for( auto &nn : currentFragment.neighbors )
            isNeighborOfCurrentFragment[nn] = false;
    }
}


}
}
