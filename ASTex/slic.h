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



#ifndef __ASTEX__SLIC_H__
#define __ASTEX__SLIC_H__

#include <Eigen/Core>

#include <ASTex/region_growing/connected_components.h>
#include <itkImageToImageFilter.h>
#include <ASTex/image_rgb.h>
#include <ASTex/colorspace_filters.h>
#include <mutex>




namespace ASTex
{


struct SuperPixel
{
    uint32_t                        label;
    Eigen::Vector2d                 coord;
    Eigen::Vector3d                 color;
    std::vector< PixelAbsolute >    pixels;
    std::vector< uint32_t >         neighbors;
};

struct SuperPixelSearchRegion
{
    Eigen::Vector2d centerOffset;
    itk::Index<2>   min;
    itk::Index<2>   max;
};

struct SuperPixelCenter
{
    SuperPixel                              stats;
    unsigned int                            pixelCount;
    std::vector< SuperPixelSearchRegion >   regions;
};


/**
 *  Filter performing the labelling of an image according to the SLIC superpixel algorithm.
 */
template < typename TNeighborhood >
class SLICSuperPixelFilter : public itk::ImageToImageFilter< ASTex::ImageRGBu8::ItkImg, ASTex::ImageGrayu32::ItkImg >
{
public:
	// Filter is not template so define input & output types for more simple writing
	using TInputImage  = ASTex::ImageRGBu8::ItkImg;
	using TOutputImage = ASTex::ImageGrayu32::ItkImg;
	using TInputPixel  = TInputImage::PixelType;
	using TOutputPixel = TOutputImage::PixelType;

	// Standard class typedefs & macros for ikt
	typedef SLICSuperPixelFilter                                    Self;
	typedef itk::ImageToImageFilter< TInputImage, TOutputImage >    Superclass;
	typedef itk::SmartPointer< Self >                               Pointer;
	typedef itk::SmartPointer< const Self >                         ConstPointer;
	itkNewMacro( Self );
	itkTypeMacro( SLICSuperPixelFilter, ImageToImageFilter );


    struct BoundaryCondition
    {
        inline static void CheckRegions( SuperPixelCenter &superPixelCenter,
                                         std::vector<SuperPixelSearchRegion> &searchRegions,
                                         const ASTex::Size &imgSize )
        {
            const bool TEMPLATE_NEIGHBORHOOD_TYPE_IS_INVALID = false;
            assert( TEMPLATE_NEIGHBORHOOD_TYPE_IS_INVALID );
			std::ignore = superPixelCenter;
			std::ignore = searchRegions;
			std::ignore = imgSize;
			std::ignore = TEMPLATE_NEIGHBORHOOD_TYPE_IS_INVALID;
        }

        inline static void AdjustCentroids( SuperPixelCenter &superPixelCenter,
                                            const ASTex::ImageGrayu32 &labelMap )
        {
            const bool TEMPLATE_NEIGHBORHOOD_TYPE_IS_INVALID = false;
            assert( TEMPLATE_NEIGHBORHOOD_TYPE_IS_INVALID );
			std::ignore = superPixelCenter;
			std::ignore = labelMap;
			std::ignore = TEMPLATE_NEIGHBORHOOD_TYPE_IS_INVALID;
		}
    };


	/**
	 * @brief set center for filter computation
	 * @param c
	 */
	void setCenter(const itk::Offset<2>& c)
	{
		center_ = c;
	}

	/**
	 * @brief set radius for filter computation
	 * @param r
	 */
	void setRadius(long r)
	{
		radius_ = r;
	}

    void setSuperPixelSize( int s )
    {
        superPixelSize_ = s;
    }

    void setIterationNumber( unsigned int nIter )
    {
        nIterations_ = nIter;
    }

    unsigned int getSuperPixelCount() const
    {
        return superPixelCenters_.size();
    }

    const SuperPixel& getSuperPixel( unsigned int i ) const
    {
        return superPixelCenters_[i].stats;
    }

    /** Range from 1 to 40 (default = 10).
     *  When m is large, spatial proximity is more important and the resulting superpixels are more compact (i.e. they have
     *  a lower area to perimeter ratio). When m is small, the resulting superpixels adhere more tightly to image boundaries,
     *  but have less regular size and shape.
     */
    void setCompactness( double compactness )
    {
        compactness_ = compactness;
    }

protected:
    using CC = ASTex::ConnectedComponent< ASTex::ImageGrayu32::PixelType >;

    struct CCAggregationGraphEdge
    {
        CC          *cc;
        uint32_t    label;
        double      colorDist;

        inline bool operator<( const CCAggregationGraphEdge &scc ) const
        {
            return colorDist < scc.colorDist  ||
                   (colorDist == scc.colorDist  &&  (label < scc.label  ||
                   (label == scc.label  &&  cc < scc.cc)));
        }
    };

    /** Functor converting a pixel from RGB colorspace (ranging [0-255]) to LAB colorspace.
     */
    template< class TInput, class TOutput>
	class functorRGBtoLAB : ColorSpace::fonctorRGB255To01<TInput,TOutput>
    {
    public:
    	inline TOutput operator()( const TInput & p ) const
    	{
		    TOutput q = ASTex::ColorSpace::fonctorRGB255To01<TInput,TOutput>::operator()( p );
            q = ASTex::ColorSpace::colorRGBtoXYZ<TOutput,TOutput>( q );
            q = ASTex::ColorSpace::colorXYZtoLAB<TOutput,TOutput>( q );
            return q;
	    }
    };

    /** Speedup structure specifically designed for the case of SLIC algorithm.
     *
     *  A 2D grid is built to store superpixels. Its cells have size bigger than the spacing S between two superpixels,
     *  so as to ensure that the search regions of two superpixels located at at least 2-cells away from each other do
     *  not overlap. In order to avoid concurrent accesses, parallel processing can be done in four steps, by inpsecting
     *  cells in the following order:
     *    1 - cells of even rows and even columns,
     *    2 -    ''    even    ''    odd    ''   ,
     *    3 -    ''    odd     ''    even   ''   ,
     *    4 -    ''    odd     ''    odd    ''   .
     *  Morever, to avoid problems at image boundaries, the size of the grid is forced to be even in each dimension.
     */
    class SLICParallelProcGrid
    {
    public:
        using Cell = std::list<SuperPixelCenter*>;

    private:
        enum ProcessingState
        {
            EVEN_ROWS_EVEN_COLS = 0,
            EVEN_ROWS_ODD_COLS  ,
            ODD_ROWS_EVEN_COLS  ,
            ODD_ROWS_ODD_COLS   ,
            FINISHED            ,
        };

        int                 imgWidth_;
        int                 imgHeight_;

        Eigen::Vector2i     gridSize_;
        Eigen::Vector2d     cellSize_;
        std::vector<Cell>   cells_;
        bool                isPopulated_;

        /** Recover the grid cell corresponding to the given pixel coordinates.
         *  \param      coord       Pixel coordinates.
         *  \return                 Reference to the corresponding grid cell.
         */
        inline Cell& getCell( const Eigen::Vector2d &coord )
        {
            unsigned int x = std::min( (int)(coord[0] / cellSize_[0]), gridSize_[0]-1 );
            unsigned int y = std::min( (int)(coord[1] / cellSize_[1]), gridSize_[1]-1 );
            return cells_[ y*gridSize_[0] + x ];
        }

        /** Recover the bounds of the image pixel region corresponding to the given grid cell.
         *  \param[in]  cell        Grid cell.
         *  \param[out] xMin        Region lower bound along the X axis (first region column).
         *  \param[out] xMax        Region upper bound along the X axis (column after the last region column).
         *  \param[out] yMin        Region lower bound along the Y axis (first region row).
         *  \param[out] yMax        Region upper bound along the Y axis (row after the last region row).
         */
        inline void getCellImageRegion( const Eigen::Vector2i &cell, int &xMin, int &xMax, int &yMin, int &yMax ) const
        {
            xMin = int( cell[0] * cellSize_[0] );
            yMin = int( cell[1] * cellSize_[1] );

            if( cell[0] == gridSize_[0]-1 )
                xMax = imgWidth_;
            else
                xMax = int( (cell[0]+1) * cellSize_[0] );

            if( cell[1] == gridSize_[1]-1 )
                yMax = imgHeight_;
            else
                yMax = int( (cell[1]+1) * cellSize_[1] );
        }

        /** Considering the current processing state, return the row id where next step must be started (odd or even rows).
         *  \param  state           Current processing state.
         *  \return                 Starting row id.
         */
        inline int getStartingRowId( ProcessingState state ) const
        {
            return ((state >> 1) & 1);
        }

        /** Considering the current processing state, return the column id where the next row starts (odd or even columns).
         *  \param  state           Current processing state.
         *  \return                 Starting column id.
         */
        inline int getStartingColumnId( ProcessingState state ) const
        {
            return (state & 1);
        }

        /** Considering the current processing state, move the given cell to the next one to be processed.
         *  \param[in,out]  cell            Current cell location.
         *  \param[in,out]  state           Current processing state.
         */
        inline bool goToNextCell( Eigen::Vector2i &cell, ProcessingState &state ) const
        {
            bool hasStateChanged = false;

            cell[0] += 2;
            if( cell[0] >= gridSize_[0] )       // go to next row.
            {
                cell[1] += 2;
                if( cell[1] >= gridSize_[1] )   // go to next step.
                {
                    state = static_cast<ProcessingState>( int(state) + 1 );
                    hasStateChanged = true;
                    cell[1] = getStartingRowId( state );
                }
                cell[0] = getStartingColumnId( state );
            }

            return hasStateChanged;
        }

    public:
        /** Constructor.
         *  \param  imgWidth            Width of the image to process.
         *  \param  imgHeight           Height of the image to process.
         *  \param  superPixelSpacing   Value indicating the radius of the search region around each superpixel center.
         */
        inline SLICParallelProcGrid( int imgWidth, int imgHeight, double superPixelSpacing ) :
            imgWidth_ ( imgWidth  ),
            imgHeight_( imgHeight ),
            isPopulated_( false )
        {
            gridSize_[0] = imgWidth  / (4*superPixelSpacing);
            gridSize_[1] = imgHeight / (4*superPixelSpacing);

            gridSize_[0] = ((gridSize_[0] >> 1) << 1);  // Ensure that grid size is even.
            gridSize_[1] = ((gridSize_[1] >> 1) << 1);

            cellSize_[0] = (double) imgWidth  / gridSize_[0];
            cellSize_[1] = (double) imgHeight / gridSize_[1];

            cells_.resize( gridSize_[0] * gridSize_[1] );
        }

        /** Initialize the speedup grid content with the given list of superpixel centers.
         *  \param superPixelCenters    List of superpixel centers.
         */
        inline void populate( std::vector<SuperPixelCenter> &superPixelCenters, bool firstIteration )
        {
            if( isPopulated_ )
                for( auto &c : cells_ )
                    c.clear();

            if( firstIteration )
            {
                for( auto &center : superPixelCenters )
                    getCell( center.stats.coord ).push_back( &center );
            }
            else
            {
                for( auto &center : superPixelCenters )
                    if( center.pixelCount )
                        getCell( center.stats.coord ).push_back( &center );
            }

            isPopulated_ = true;
        }

        /** Parallel processing of the grid content.
         *  \param parallelProcFunc     Processing function to apply to each grid cell.
         */
        template <typename TParallelProcFunc>
        inline void parallelProcessing( const TParallelProcFunc &parallelProcFunc )
        {
            std::mutex m;
            std::condition_variable condition;

            unsigned int nThreads = std::thread::hardware_concurrency();
            std::vector< std::thread* > threads( nThreads );

            unsigned int nThreadsFinished = 0;

            Eigen::Vector2i nextCell( 0, 0 );
            ProcessingState processingState = EVEN_ROWS_EVEN_COLS;
            bool hasStateChanged = false;

            for( unsigned int i=0; i<nThreads; ++i )
                threads[i] = new std::thread( [&]()
                {
                    do
                    {
                        std::unique_lock<std::mutex> lock( m );

                        if( hasStateChanged )
                        {
                            if( processingState == FINISHED )
                                break;

                            if( (++nThreadsFinished) < nThreads )
                                condition.wait( lock );
                            else
                            {
                                hasStateChanged = false;
                                nThreadsFinished = 0;
                                condition.notify_all();
                            }

                            continue;
                        }

                        auto &gridCell = cells_[ nextCell[1]*gridSize_[0] + nextCell[0] ];

                        int xMin, xMax, yMin, yMax;
                        getCellImageRegion( nextCell, xMin, xMax, yMin, yMax );

                        hasStateChanged = goToNextCell( nextCell, processingState );

                        lock.unlock();

                        parallelProcFunc( gridCell, xMin, xMax, yMin, yMax );
                    } while( processingState != FINISHED );
                });

            for( unsigned int i=0; i<nThreads; ++i )
            {
                threads[i]->join();
                delete threads[i];
            }
        }
    };


    /// filter data
	itk::Offset<2>                  center_;
	long                            radius_;

    int                             superPixelSize_;
	unsigned int                    nIterations_;
    double                          compactness_;

    ASTex::Size                     inputImageSize_;
    double                          superPixelSpacing_;
    std::vector< SuperPixelCenter > superPixelCenters_;
    double                          colorDistStrength_;


	/// protected constructor (forbid usage of new and variable declaration)
	SLICSuperPixelFilter() :
		center_({{0,0}}),
		radius_(1),
        superPixelSize_(0),
        nIterations_(10),
        compactness_(10)
	{}

	virtual ~SLICSuperPixelFilter() {}

    /** Initialization SLIC superpixel centers by generating seeds on a regular grid,
     *  ensuring that none of them is located on a sharp feature.
     *
     *  \param[in]      inputImageLAB       Input image converted to LAB colorspace.
     */
    inline void generateSuperPixelSeeds( const ASTex::ImageRGBf &inputImageLAB )
    {
        // Generate equally spaced center coordinates by following a regular grid pattern.

        superPixelSpacing_ = std::sqrt( (double) superPixelSize_ );

        double sx = inputImageSize_[0] / std::floor(inputImageSize_[0]/superPixelSpacing_);
        double sy = inputImageSize_[1] / std::floor(inputImageSize_[1]/superPixelSpacing_);

        superPixelCenters_.clear();
        superPixelCenters_.reserve( std::ceil(inputImageSize_[0]/sx) * std::ceil(inputImageSize_[1]/sy) );

        for( double y=0.5*sy; y<inputImageSize_[1]; y+=sy )
            for( double x=0.5*sx; x<inputImageSize_[0]; x+=sx )
            {
                // Create a new superpixel center.

                superPixelCenters_.push_back( SuperPixelCenter() );
                auto &c = superPixelCenters_.back();

                c.stats.label = superPixelCenters_.size() - 1;

                c.stats.color.setZero();

                c.stats.coord[0] = std::floor( x );
                c.stats.coord[1] = std::floor( y );

                // Ensure that its coordinates are at least 2 pixels far away from image boundaries (required by the next step).

                for( int dim=0; dim<2; ++dim )
                    if( c.stats.coord[dim] < 2 )
                        c.stats.coord[dim] = 2;
                    else if( c.stats.coord[dim] > inputImageSize_[dim]-3 )
                        c.stats.coord[dim] = inputImageSize_[dim] - 3;

                // Move the center to the pixel of its 3x3 neighborhood for which the gradient value is the smallest,
                // so as to ensure that it is not located over a salient edge.

                double smallestGradient = std::numeric_limits<double>::max();
				int bestX=0;
				int bestY=0;

                for( int ny = int(c.stats.coord[1])-1; ny <= int(c.stats.coord[1])+1; ++ny )
                    for( int nx = int(c.stats.coord[0])-1; nx <= int(c.stats.coord[0])+1; ++nx )
                    {
                        auto gx = inputImageLAB.pixelAbsolute( nx+1, ny ) -
                                  inputImageLAB.pixelAbsolute( nx-1, ny );

                        auto gy = inputImageLAB.pixelAbsolute( nx, ny+1 ) -
                                  inputImageLAB.pixelAbsolute( nx, ny-1 );

                        double g = gx[0]*gx[0] + gx[1]*gx[1] + gx[2]*gx[2] +
                                   gy[0]*gy[0] + gy[1]*gy[1] + gy[2]*gy[2];

                        if( g < smallestGradient )
                        {
                            smallestGradient = g;
                            bestX = nx;
                            bestY = ny;
                        }
                    }

                c.stats.coord[0] = bestX;
                c.stats.coord[1] = bestY;
            }
    }



//    template <template <typename T> class NT>
//    inline void searchRegionBoundaryCheck( SuperPixelCenter &superPixelCenter,
//                                           std::vector<SuperPixelSearchRegion> &searchRegions )
//    {
//        const bool TEMPLATE_NEIGHBORHOOD_TYPE_IS_INVALID = false;
//        assert( TEMPLATE_NEIGHBORHOOD_TYPE_IS_INVALID );
//    }
	//    // In case of a non-periodic image, this initial region is simply clamped wrt. image boundaries.
	//    template <>
	//	inline void searchRegionBoundaryCheck<NeighborhoodClamped>( SuperPixelCenter &superPixelCenter,
	//                                                                       std::vector<SuperPixelSearchRegion> &searchRegions )
	//    {
	//        for( int i=0; i<2; ++i )
	//        {
	//            if( searchRegions[0].min[i] < 0 )
	//                searchRegions[0].min[i] = 0;
	//            else if( searchRegions[0].max[i] > (int) inputImageSize_[i]-1 )
	//                searchRegions[0].max[i] = inputImageSize_[i]-1;
	//        }
	//    }

	//    // In case of a periodic image, this region is subdivided wrt. to every possibly crossed image boundary,
	//    // resulting to a maximum of 4 different regions, if the initial one is located onto an image corner.
	//    template <>
	//	inline void searchRegionBoundaryCheck<NeighborhoodPeriodic>( SuperPixelCenter &superPixelCenter,
	//                                                                        std::vector<SuperPixelSearchRegion> &searchRegions )
	//    {
	//        searchRegions[0].min[0] = (searchRegions[0].min[0] + inputImageSize_[0]) % inputImageSize_[0];
	//        searchRegions[0].max[0] = (searchRegions[0].max[0] + inputImageSize_[0]) % inputImageSize_[0];
	//        searchRegions[0].min[1] = (searchRegions[0].min[1] + inputImageSize_[1]) % inputImageSize_[1];
	//        searchRegions[0].max[1] = (searchRegions[0].max[1] + inputImageSize_[1]) % inputImageSize_[1];

	//        for( int dim=0; dim<2; ++dim )
	//        {
	//            size_t currentRegionCount = searchRegions.size();
	//            for( size_t n=0; n<currentRegionCount; ++n )
	//                if( searchRegions[n].max[dim] < searchRegions[n].min[dim] )
	//                {
	//                    searchRegions.push_back( searchRegions[n] );

	//                    searchRegions[n].min[dim] = 0;
	//                    searchRegions.back().max[dim] = inputImageSize_[dim]-1;

	//                    if( superPixelCenter.stats.coord[dim] > 0.5*inputImageSize_[dim] )
	//                        searchRegions[n].centerOffset[dim] = -(double) inputImageSize_[dim];
	//                    else
	//                        searchRegions.back().centerOffset[dim] = (double) inputImageSize_[dim];
	//                }
	//        }

	//        if( searchRegions.size() > 1 )
	//            superPixelCenter.regions = searchRegions;
	//        else
	//            superPixelCenter.regions.clear();
	//    }


/*
	// In case of a non-periodic image, this initial region is simply clamped wrt. image boundaries.
	template <typename N>
	auto searchRegionBoundaryCheck( SuperPixelCenter &superPixelCenter, std::vector<SuperPixelSearchRegion> &searchRegions )
						-> typename std::enable_if<std::is_same<N,NeighborhoodClamped<PixelAbsolute>>::value ||
							std::is_same<N,NeighborhoodClamped<PixelRelative>>::value,void>::type
	{
		for( int i=0; i<2; ++i )
		{
			if( searchRegions[0].min[i] < 0 )
				searchRegions[0].min[i] = 0;
			else if( searchRegions[0].max[i] > (int) inputImageSize_[i]-1 )
				searchRegions[0].max[i] = inputImageSize_[i]-1;
		}
	}

	// In case of a periodic image, this region is subdivided wrt. to every possibly crossed image boundary,
	// resulting to a maximum of 4 different regions, if the initial one is located onto an image corner.
	template <typename N>
	auto searchRegionBoundaryCheck( SuperPixelCenter &superPixelCenter,
				  std::vector<SuperPixelSearchRegion> &searchRegions )
	-> typename std::enable_if<std::is_same<N,NeighborhoodPeriodic<PixelAbsolute>>::value ||
		std::is_same<N,NeighborhoodPeriodic<PixelRelative>>::value,void>::type
	{
		searchRegions[0].min[0] = (searchRegions[0].min[0] + inputImageSize_[0]) % inputImageSize_[0];
		searchRegions[0].max[0] = (searchRegions[0].max[0] + inputImageSize_[0]) % inputImageSize_[0];
		searchRegions[0].min[1] = (searchRegions[0].min[1] + inputImageSize_[1]) % inputImageSize_[1];
		searchRegions[0].max[1] = (searchRegions[0].max[1] + inputImageSize_[1]) % inputImageSize_[1];

		for( int dim=0; dim<2; ++dim )
		{
			size_t currentRegionCount = searchRegions.size();
			for( size_t n=0; n<currentRegionCount; ++n )
				if( searchRegions[n].max[dim] < searchRegions[n].min[dim] )
				{
					searchRegions.push_back( searchRegions[n] );

					searchRegions[n].min[dim] = 0;
					searchRegions.back().max[dim] = inputImageSize_[dim]-1;

					if( superPixelCenter.stats.coord[dim] > 0.5*inputImageSize_[dim] )
						searchRegions[n].centerOffset[dim] = -(double) inputImageSize_[dim];
					else
						searchRegions.back().centerOffset[dim] = (double) inputImageSize_[dim];
				}
		}

		if( searchRegions.size() > 1 )
			superPixelCenter.regions = searchRegions;
		else
			superPixelCenter.regions.clear();
	}


*/

    /** Recover the rectangular regions of pixels for which distance to the given superpixel center must be computed.
     *
     *  In case a periodic image is processed, this function may produce multiple regions, if image boundaries are
     *  crossed. Otherwise, only one clamped region is created.
     *
     *  \param[in,out]  superPixelCenter        The superpixel center for which search regions must be recovered.
     *  \param[in,out]  searchRegions           The resulting rectangles of pixels for which distance to the center must be computed.
     */
    inline void getSearchRegions( SuperPixelCenter &superPixelCenter,
                                  std::vector<SuperPixelSearchRegion> &searchRegions )
    {
        // An initial rectangular region is created, that does not account for image boundaries.

        searchRegions.resize( 1 );
        searchRegions[0].centerOffset.setZero();
        searchRegions[0].min[0] = (int) std::floor( superPixelCenter.stats.coord[0] - superPixelSpacing_ );
        searchRegions[0].max[0] = (int) std::ceil ( superPixelCenter.stats.coord[0] + superPixelSpacing_ );
        searchRegions[0].min[1] = (int) std::floor( superPixelCenter.stats.coord[1] - superPixelSpacing_ );
        searchRegions[0].max[1] = (int) std::ceil ( superPixelCenter.stats.coord[1] + superPixelSpacing_ );

        BoundaryCondition::CheckRegions( superPixelCenter, searchRegions, inputImageSize_ );
    }

    /** Determine the closest superpixel center for each image pixel, based on a metric combining spatial and chromatic distances.
     *
     *  \param[in,out]  parallelProcGrid    The speedup structure for efficient parallel processing.
     *  \param[in]      inputImageLAB       Input image converted to LAB colorspace.
     *  \param[out]     distanceMap         Map storing the per-pixel distance to the closest superpixel center.
     *  \param[out]     labelMap            Map storing the per-pixel label of the closest superpixel center (its number in the superPixelCenters_ array).
     */
    inline void assignCentersToPixels( SLICParallelProcGrid &parallelProcGrid,
									   const ImageRGBf &inputImageLAB,
									   ImageGrayd &distanceMap,
									   ImageGrayu32 &labelMap )
    {
        // Reset distance map to inifinite (or at least very high) value.

		distanceMap.parallel_for_all_pixels( []( ImageGrayd::PixelType &p )
        {
            p = std::numeric_limits<double>::max();
        });

        // For each superpixel center, recover the impacted regions of pixels and compute pixel-to-center distances.
        // Each pixel then keeps track of the superpixel center for which this distance is the smallest.

		parallelProcGrid.parallelProcessing( [&](typename SLICParallelProcGrid::Cell &cell, int /*xMin*/, int /*xMax*/, int /*yMin*/, int /*yMax*/ )
        {
            std::vector<SuperPixelSearchRegion> searchRegions;
            searchRegions.reserve( 4 );

            for( auto center : cell )
            {
                getSearchRegions( *center, searchRegions );

                for( auto &region : searchRegions )
                    for( int y=region.min[1]; y<=region.max[1]; ++y )
                        for( int x=region.min[0]; x<=region.max[0]; ++x )
                        {
                            double &pixelDist = distanceMap.pixelAbsolute( x, y );

                            // First, check the spatial distance only. If it is greater than the best distance recorded so far, stop immediately.

                            double sqDistSpatial = ( center->stats.coord + region.centerOffset - Eigen::Vector2d(x,y) ).squaredNorm();
                            if( sqDistSpatial >= pixelDist )
                                continue;

                            // Otherwise, compute the final SLIC metric, combining spatial and color distances.

                            auto pixelColor = inputImageLAB.pixelAbsolute(x,y);
                            double sqDistColor = ( center->stats.color - Eigen::Vector3d(pixelColor[0],pixelColor[1],pixelColor[2]) ).squaredNorm();

                            double d = sqDistColor * colorDistStrength_ + sqDistSpatial;

                            // If it is better than the recorded one, store this superpixel center as the closest one.

                            if( d < pixelDist  ||  (d == pixelDist && center->stats.label < labelMap.pixelAbsolute(x,y)) )
                            {
                                pixelDist = d;
                                labelMap.pixelAbsolute( x, y ) = center->stats.label;
                            }
                        }
            }
        });
    }


//    template <template <typename T> class NT>
//    inline void applyCentroidCorrection( SuperPixelCenter &superPixelCenter,
//										 const ImageGrayu32 &labelMap )
//    {
//        const bool TEMPLATE_NEIGHBORHOOD_TYPE_IS_INVALID = false;
//        assert( TEMPLATE_NEIGHBORHOOD_TYPE_IS_INVALID );
//    }
/*
	template <typename N>
	auto applyCentroidCorrection( SuperPixelCenter &superPixelCenter, const ImageGrayu32 &labelMap )
	-> typename std::enable_if<std::is_same<N,NeighborhoodClamped<PixelAbsolute>>::value ||
		std::is_same<N,NeighborhoodClamped<PixelRelative>>::value,void>::type
	{
	}

	template <typename N>
	auto applyCentroidCorrection( SuperPixelCenter &superPixelCenter, const ImageGrayu32 &labelMap )
	-> typename std::enable_if<std::is_same<N,NeighborhoodPeriodic<PixelAbsolute>>::value ||
		std::is_same<N,NeighborhoodPeriodic<PixelRelative>>::value,void>::type

	{
		for( auto &region : superPixelCenter.regions )
			for( int y=region.min[1]; y<=region.max[1]; ++y )
				for( int x=region.min[0]; x<=region.max[0]; ++x )
					if( labelMap.pixelAbsolute(x,y) == superPixelCenter.stats.label )
						superPixelCenter.stats.coord -= region.centerOffset;
	}
*/

    enum
    {
        STORE_PIXELS        ,
        DO_NOT_STORE_PIXELS ,
    };

//    template <int T>
//    inline void managePixelStoring( SuperPixelCenter &center, int x, int y )
//    {
//        const bool PIXEL_STORING_POLICY_IS_INVALID = false;
//        assert( PIXEL_STORING_POLICY_IS_INVALID );
//    }

	template <int SP>
	inline void managePixelStoring( SuperPixelCenter &center, int x, int y )
    {
		if (SP==STORE_PIXELS)
			center.stats.pixels.push_back( {x,y} );
		else
			if (SP==DO_NOT_STORE_PIXELS)
			{}
			else
			{
				const bool PIXEL_STORING_POLICY_IS_INVALID = false;
				assert( PIXEL_STORING_POLICY_IS_INVALID );
				std::ignore = PIXEL_STORING_POLICY_IS_INVALID;
			}
    }


    /** Update superpixels' mean colors and centroids based on the labelling given in the provided label map.
     *
     *  \param[in,out]  parallelProcGrid    The speedup structure for efficient parallel processing.
     *  \param[in]      inputImageLAB       Input image converted to LAB colorspace.
     *  \param[in]      labelMap            Map storing the per-pixel label of the closest superpixel center.
     */
    template <int TPixelStoringPolicy>
    inline void updateSuperPixelCenters( SLICParallelProcGrid &parallelProcGrid,
										 const ImageRGBf &inputImageLAB,
										 const ImageGrayu32 &labelMap )
    {
        // Initialize mean color, centroid and pixel count of each superpixel center to zero.

        for( auto &center : superPixelCenters_ )
        {
            center.stats.color.setZero();
            center.stats.coord.setZero();
            center.pixelCount = 0;
        }

        // Count the number of pixels belonging to each superpixel, and sum their color and position values.

		parallelProcGrid.parallelProcessing( [&](typename SLICParallelProcGrid::Cell& /*cell*/, int xMin, int xMax, int yMin, int yMax )
        {
            for( int y=yMin; y<yMax; ++y )
                for( int x=xMin; x<xMax; ++x )
                {
                    uint32_t label = labelMap.pixelAbsolute( x, y );

                    auto &center = superPixelCenters_[label];
                    auto &color = inputImageLAB.pixelAbsolute( x, y );

                    center.stats.coord += Eigen::Vector2d( x, y );
                    center.stats.color += Eigen::Vector3d( color.GetRed(), color.GetGreen(), color.GetBlue() );
                    center.pixelCount ++;

                    managePixelStoring<TPixelStoringPolicy>( center, x, y );
                }
        });

        // Update mean color and centroid, and ensure that the latter is located inside image bounds.

        for( auto &center : superPixelCenters_ )
            if( center.pixelCount )
            {
                // In the case of a periodic image, ensure that summing for centroid computation is correct even for superpixels crossing image boundaries.

                BoundaryCondition::AdjustCentroids( center, labelMap );

                // Compute mean color and centroid.

                center.stats.color /= (double) center.pixelCount;
                center.stats.coord /= (double) center.pixelCount;

                // Ensure that centroid coordinates are inside image boundaries.

                for( int i=0; i<2; ++i )
                    if( center.stats.coord[i] < 0.0 )
                        center.stats.coord[i] += inputImageSize_[i];
                    else if( center.stats.coord[i] >= inputImageSize_[i] )
                        center.stats.coord[i] -= inputImageSize_[i];
            }
    }

    /** Cleaning step that build a label map with a single connected component per superpixel.
     *  This is done by a relabelling strategy, based on components neighborhood.
     *
     *  \param[in,out]  labelMap        Superpixel map to process for relabelling.
     */
    inline void cleanIsolatedComponents( std::vector<CC> &connectedComponents,
										 ImageGrayu32 &labelMap ) const
    {
        // Only the biggest connected component associated to each superpixel is retained...

        std::vector< CC* > superPixelComponents;
        superPixelComponents.assign( superPixelCenters_.size(), NULL );

        for( auto &cc : connectedComponents )
            if( superPixelComponents[cc.label] == NULL  ||  superPixelComponents[cc.label]->pixels.size() < cc.pixels.size() )
                superPixelComponents[cc.label] = &cc;

        // ... and used as a seed to initialize the aggregation algorithm.

        std::set< CCAggregationGraphEdge > processingQueue;

        for( auto cc : superPixelComponents )
            if( cc )
            {
                CCAggregationGraphEdge candidate;
                candidate.label     = cc->label;
                candidate.cc        = cc;
                candidate.colorDist = -std::numeric_limits<double>::max();  // force seeds to be processed first.
                processingQueue.insert( candidate );
            }

        // Aggregation: the processing queue contains couples of neighbor components, one already aggregated (the destination)
        // and another one which is the candidate for aggregation. All these couples are ordered by increasing distance between
        // the two components' mean LAB colors. Couples with the smallest distances are aggregated first.

        std::vector< bool > isComponentLabelled;
        isComponentLabelled.assign( connectedComponents.size(), false );

        while( !processingQueue.empty() )
        {
            // Get next candidate component for aggregation.
            auto candidate = *processingQueue.begin();
            processingQueue.erase( processingQueue.begin() );

            auto cc = candidate.cc;

            // If already aggregated, ignore it.
            if( isComponentLabelled[cc->id] )
                continue;

            isComponentLabelled[cc->id] = true;

            // Perform relabelling, if required.
            if( cc->label != candidate.label )
            {
                cc->label = candidate.label;
                for( auto &p : cc->pixels )
                    labelMap.pixelAbsolute( p ) = cc->label;
            }

            // Add new component couples into the processing queue, one for each of its unlabelled neighbors.
            CCAggregationGraphEdge graphEdge;
            graphEdge.label = cc->label;

            for( auto nn : cc->neighbors )
                if( !isComponentLabelled[nn->id] )
                {
                    graphEdge.cc = nn;
                    auto colorDiff = superPixelCenters_[cc->label].stats.color - superPixelCenters_[nn->label].stats.color;
                    graphEdge.colorDist = colorDiff[0]*colorDiff[0] + colorDiff[1]*colorDiff[1] + colorDiff[2]*colorDiff[2];
                    processingQueue.insert( graphEdge );
                }
        }
    }

	//
	// overriden method that generate the output
	//
	void GenerateData() ITK_OVERRIDE
	{
		//  we are responsible of allocating the images buffers
		// in this simple case (output of same size than inputs) just call:
		this->AllocateOutputs();

        // Convert input picture from RGB to LAB color space.

		typedef itk::UnaryFunctorImageFilter< ImageRGBu8::ItkImg, ImageRGBf::ItkImg,
			functorRGBtoLAB<ImageRGBu8::PixelType, ImageRGBf::PixelType> > FilterRGBtoLABType;

		typename FilterRGBtoLABType::Pointer toLABFilter = FilterRGBtoLABType::New();
        toLABFilter->SetInput( this->GetInput() );
        toLABFilter->Update();
		ImageRGBf inputImageLAB( toLABFilter->GetOutput() );

        inputImageSize_ = GetInput()->GetLargestPossibleRegion().GetSize();

        // Generate superpixel seeds.

        generateSuperPixelSeeds( inputImageLAB );

        // Apply iteration, by alternating a phase where each pixel is assigned to the closest superpixel
        // and a phase where superpixel centroids and mean colors are updated accordingly.

		ImageGrayu32 labelMap( this->GetOutput() );
		ImageGrayd distanceMap( inputImageSize_[0], inputImageSize_[1] );

        colorDistStrength_ = superPixelSize_ / (compactness_*compactness_);

        SLICParallelProcGrid parallelProcGrid( inputImageLAB.width(), inputImageLAB.height(), superPixelSpacing_ );

        for( unsigned int i=0; i<nIterations_; ++i )
        {
            parallelProcGrid.populate( superPixelCenters_, i==0 );
            assignCentersToPixels( parallelProcGrid, inputImageLAB, distanceMap, labelMap );
            updateSuperPixelCenters<DO_NOT_STORE_PIXELS>( parallelProcGrid, inputImageLAB, labelMap );
        }

        // Recover the list of connected components.

        std::vector< CC > connectedComponents;
        connectedComponents.reserve( superPixelCenters_.size() );
        ASTex::getConnectedComponentsWithNeighborhood< TNeighborhood >( labelMap, connectedComponents );

        // Cleaning step, where isolated connected components are aggregated to valid superpixels.

        cleanIsolatedComponents( connectedComponents, labelMap );

        // On the basis of the clean label map, update superpixel centroids, mean colors and lists of pixels.

        updateSuperPixelCenters<STORE_PIXELS>( parallelProcGrid, inputImageLAB, labelMap );

        // Recover neighborhoods of superpixels from the neighborhood of connected components extracted before relabelling.

        std::vector< std::set<uint32_t> > labelNeighborhood( superPixelCenters_.size() );

        for( auto &cc : connectedComponents )
            for( auto nn : cc.neighbors )
                if( nn->label != cc.label )
                    labelNeighborhood[cc.label].insert( nn->label );

        for( size_t i=0; i<superPixelCenters_.size(); ++i )
        {
            superPixelCenters_[i].stats.neighbors.reserve( labelNeighborhood[i].size() );
            for( auto nn : labelNeighborhood[i] )
                superPixelCenters_[i].stats.neighbors.push_back( nn );
        }
    }

private:
	// to avoid filter object copy
	SLICSuperPixelFilter( const Self & );
	void operator=( const Self & );
};



// In case of a non-periodic image, this initial region is simply clamped wrt. image boundaries.
template <>
inline void SLICSuperPixelFilter<ASTex::NeighborhoodClamped>::BoundaryCondition::CheckRegions( SuperPixelCenter& /*superPixelCenter*/,
                                                                                               std::vector<SuperPixelSearchRegion> &searchRegions,
                                                                                               const ASTex::Size &imgSize )
{
    for( int i=0; i<2; ++i )
    {
        if( searchRegions[0].min[i] < 0 )
            searchRegions[0].min[i] = 0;
        else if( searchRegions[0].max[i] > (int) imgSize[i]-1 )
            searchRegions[0].max[i] = imgSize[i]-1;
    }
}

// In case of a periodic image, this region is subdivided wrt. to every possibly crossed image boundary,
// resulting to a maximum of 4 different regions, if the initial one is located onto an image corner.
template <>
inline void SLICSuperPixelFilter<ASTex::NeighborhoodPeriodic>::BoundaryCondition::CheckRegions( SuperPixelCenter &superPixelCenter,
                                                                                                std::vector<SuperPixelSearchRegion> &searchRegions,
                                                                                                const ASTex::Size &imgSize )
{
    searchRegions[0].min[0] = (searchRegions[0].min[0] + imgSize[0]) % imgSize[0];
    searchRegions[0].max[0] = (searchRegions[0].max[0] + imgSize[0]) % imgSize[0];
    searchRegions[0].min[1] = (searchRegions[0].min[1] + imgSize[1]) % imgSize[1];
    searchRegions[0].max[1] = (searchRegions[0].max[1] + imgSize[1]) % imgSize[1];

    for( int dim=0; dim<2; ++dim )
    {
        size_t currentRegionCount = searchRegions.size();
        for( size_t n=0; n<currentRegionCount; ++n )
            if( searchRegions[n].max[dim] < searchRegions[n].min[dim] )
            {
                searchRegions.push_back( searchRegions[n] );

                searchRegions[n].min[dim] = 0;
                searchRegions.back().max[dim] = imgSize[dim]-1;

                if( superPixelCenter.stats.coord[dim] > 0.5*imgSize[dim] )
                    searchRegions[n].centerOffset[dim] = -(double) imgSize[dim];
                else
                    searchRegions.back().centerOffset[dim] = (double) imgSize[dim];
            }
    }

    if( searchRegions.size() > 1 )
        superPixelCenter.regions = searchRegions;
    else
        superPixelCenter.regions.clear();
}

template <>
inline void SLICSuperPixelFilter<ASTex::NeighborhoodPeriodicRelative>::BoundaryCondition::CheckRegions( SuperPixelCenter &superPixelCenter,
                                                                                                        std::vector<SuperPixelSearchRegion> &searchRegions,
                                                                                                        const ASTex::Size &imgSize )
{
    SLICSuperPixelFilter<ASTex::NeighborhoodPeriodic>::BoundaryCondition::CheckRegions( superPixelCenter, searchRegions, imgSize );
}




template <>
inline void SLICSuperPixelFilter<ASTex::NeighborhoodClamped>::BoundaryCondition::AdjustCentroids( SuperPixelCenter& /*superPixelCenter*/,
																								  const ASTex::ImageGrayu32& /*labelMap*/ )
{
}

template <>
inline void SLICSuperPixelFilter<ASTex::NeighborhoodPeriodic>::BoundaryCondition::AdjustCentroids( SuperPixelCenter &superPixelCenter,
                                                                                                   const ASTex::ImageGrayu32 &labelMap )
{
    for( auto &region : superPixelCenter.regions )
        for( int y=region.min[1]; y<=region.max[1]; ++y )
            for( int x=region.min[0]; x<=region.max[0]; ++x )
                if( labelMap.pixelAbsolute(x,y) == superPixelCenter.stats.label )
                    superPixelCenter.stats.coord -= region.centerOffset;
}

template <>
inline void SLICSuperPixelFilter<ASTex::NeighborhoodPeriodicRelative>::BoundaryCondition::AdjustCentroids( SuperPixelCenter &superPixelCenter,
                                                                                                           const ASTex::ImageGrayu32 &labelMap )
{
    SLICSuperPixelFilter<ASTex::NeighborhoodPeriodic>::BoundaryCondition::AdjustCentroids( superPixelCenter, labelMap );
}




} // namespace ASTex




#endif // __ASTEX__SLIC_H__
