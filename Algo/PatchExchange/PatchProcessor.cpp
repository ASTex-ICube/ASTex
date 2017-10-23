#include "PatchProcessor.h"
#include "PoissonDiskSampler.h"
#include <itkResampleImageFilter.h>
#include <itkScaleTransform.h>


namespace ASTex
{

namespace ContentExchg
{


PatchProcessor::PatchProcessor( const FragmentProcessor &fragmentProc ) :
    fragmentProc_( fragmentProc )
{
}


int PatchProcessor::patchCount() const
{
    return (int) allPatches_.size();
}


Patch& PatchProcessor::patchById( int pid )
{
    assert( pid < (int) allPatches_.size() );
    return allPatches_[pid];
}


const Patch& PatchProcessor::patchById( int pid ) const
{
    assert( pid < (int) allPatches_.size() );
    return allPatches_[pid];
}


std::vector<Patch>& PatchProcessor::patches()
{
    return allPatches_;
}


const std::vector<Patch>& PatchProcessor::patches() const
{
    return allPatches_;
}


void PatchProcessor::createPatches( int requiredPatchNumber )
{
    struct FragmentAggregationInfo
    {
        double      dist;
        uint32_t    fragmentId;
        uint32_t    patchId;
        inline bool operator<( const FragmentAggregationInfo &f ) const { return dist < f.dist; }
    };

    std::multiset<FragmentAggregationInfo> aggregationQueue;

    auto &fragMap = fragmentProc_.idMap();
    idMap_ = ASTex::ImageGrayu32( fragMap.width(), fragMap.height() );
    std::vector<bool> isFragmentAssigned( fragmentProc_.fragmentCount() );

    // Determine a seed fragment for every patch.

    double distanceBetweenPatches = std::min( fragMap.width(), fragMap.height() ) / std::sqrt( (double) 2.0 * requiredPatchNumber );
    bool areAllPatchesInitialized = true;

    do
    {
        areAllPatchesInitialized = true;

		for( auto&& flag : isFragmentAssigned )
            flag = false;

        aggregationQueue.clear();

        // Patch centers are genereted by a Poisson disk sampling.

        std::vector<Point> patchSeeds;
        PoissonDiskSampler::generateSamples( fragMap.width(), fragMap.height(), patchSeeds, distanceBetweenPatches );

        // The number of centers defines the patch count.

        allPatches_.resize( patchSeeds.size() );

        // For each of them, try to find a seed fragment.

		for( std::size_t pid=0; pid<patchSeeds.size() && areAllPatchesInitialized; ++pid )
        {
            allPatches_[pid].id = pid;

            // First, the fragment intersecting the center is considered. If this fragment has already
            // been assigned to another patch, try with its neighbors.

            std::list<const Fragment*> candidateFragments;
            candidateFragments.push_back( &fragmentProc_.fragmentAt( patchSeeds[pid][0], patchSeeds[pid][1] ) );
            for( auto &neighbId : candidateFragments.back()->neighbors )
                candidateFragments.push_back( &fragmentProc_.fragmentById(neighbId) );

            bool found = false;

            for( auto f : candidateFragments )
                if( !isFragmentAssigned[f->id] )
                {
                    isFragmentAssigned[f->id] = true;

                    allPatches_[pid].centroid[0] = (itk::IndexValueType) patchSeeds[pid][0];
                    allPatches_[pid].centroid[1] = (itk::IndexValueType) patchSeeds[pid][1];

                    FragmentAggregationInfo agregInfo;
                    agregInfo.patchId = pid;
                    agregInfo.fragmentId = f->id;
                    agregInfo.dist = 0.0;
                    aggregationQueue.insert( agregInfo );

                    found = true;

                    break;
                }

            // If none of all these fragments is free, the process is started again with new centers.

            areAllPatchesInitialized &= found;
        }

        // To avoid infinite loop due to a distance between patches lower than the average size of fragment
        // (which may infinitely prevent some patches to find an unassigned fragment among the candidates),
        // the generation process will consider an increasing distance between patch centers.

        distanceBetweenPatches *= 1.2;

    } while( !areAllPatchesInitialized );

    for( auto &it : aggregationQueue )
        isFragmentAssigned[it.fragmentId] = false;

    // Agregate fragments 

    while( !aggregationQueue.empty() )
    {
        uint32_t patchId = aggregationQueue.begin()->patchId;
        uint32_t fragmentId = aggregationQueue.begin()->fragmentId;
        aggregationQueue.erase( aggregationQueue.begin() );

        if( isFragmentAssigned[fragmentId] )
            continue;

        isFragmentAssigned[fragmentId] = true;
        allPatches_[patchId].fragments.push_back( fragmentId );

        for( auto &pixel : fragmentProc_.fragmentById(fragmentId).pixels )
            idMap_.pixelAbsolute( pixel ) = patchId;

        for( auto neighb : fragmentProc_.fragmentById(fragmentId).neighbors )
            if( !isFragmentAssigned[neighb] )
            {
                FragmentAggregationInfo agregInfo;
                agregInfo.patchId = patchId;
                agregInfo.fragmentId = neighb;
                agregInfo.dist = std::numeric_limits<double>::max();

                for( int yy=-fragMap.height(); yy<=fragMap.height(); yy+=fragMap.height() )
                    for( int xx=-fragMap.width(); xx<=fragMap.width(); xx+=fragMap.width() )
                    {
                        Eigen::Vector2d centroid( allPatches_[patchId].centroid[0], allPatches_[patchId].centroid[1] );
                        double d = (fragmentProc_.fragmentById(neighb).centroid + Eigen::Vector2d(xx,yy) - centroid).squaredNorm();
                        if( d < agregInfo.dist )
                            agregInfo.dist = d;
                    }

                aggregationQueue.insert( agregInfo );
            }
    }
}


void PatchProcessor::computePatchBoundaries()
{
    static const uint32_t VISITED_FLAG = 0x80000000;

    // For each patch...

    ASTex::NeighborhoodPeriodic::NeighborsList neighbors( 8 );
    ASTex::NeighborhoodPeriodic neighborhood( idMap_ );

    for( auto &patch : allPatches_ )
    {
        // ... check every pixel belonging to it.

        for( auto &frag : patch.fragments )
            for( auto &pixel : fragmentProc_.fragmentById(frag).pixels )
            {
                // Recover its 8-neighbors and check each of them.

                neighborhood.get8( pixel, neighbors );

                bool isPixelOnBoundary = false;

                // If one neighbor belongs to another patch and is visited for the first time, it is added
                // to the patch boundary (this is how patch's "external boundary" pixels are recovered).

                for( auto &nn : neighbors )
                {
                    uint32_t &neighbPatchId = idMap_.pixelAbsolute( nn );
                    if( neighbPatchId != patch.id )
                    {
                        isPixelOnBoundary = true;
                        if( !(neighbPatchId & VISITED_FLAG) )
                        {
                            patch.boundary.push_back( nn );
                            neighbPatchId |= VISITED_FLAG;
                        }
                    }
                }

                // The current pixel is also added to the boundary if at least one of the 8-neighbors has been detected
                // as belonging to another patch (this is how patch's "internal boundary" pixels are recovered).

                if( isPixelOnBoundary )
                    patch.boundary.push_back( pixel );
            }

        // Clear the "visited" flag.

        for( auto &pixel : patch.boundary )
            idMap_.pixelAbsolute( pixel ) &= (~VISITED_FLAG);
    }
}


void PatchProcessor::getErosionMap( ASTex::ImageGrayu32 &erosionMap, uint32_t thickness )
{
    assert( !allPatches_.empty() );
    assert( !allPatches_[0].boundary.empty() );

    uint32_t currentDist = 1;

    std::vector<PixelPos> queues[2];
    auto processingQueue = &queues[0];
    auto nextQueue = &queues[1];


    erosionMap.parallel_for_all_pixels( [&] (ASTex::ImageGrayu32::PixelType& p)
	{
        p = 0;
	});

    for( auto &patch : allPatches_ )
        for( auto &pixel : patch.boundary )
            if( erosionMap.pixelAbsolute(pixel) == 0 )
            {
                erosionMap.pixelAbsolute( pixel ) = currentDist;
                processingQueue->push_back( pixel );
            }

    nextQueue->reserve( processingQueue->size() );

    ASTex::NeighborhoodPeriodic::NeighborsList neighbors( 8 );
    ASTex::NeighborhoodPeriodic neighborhood( idMap_ );

    while( !processingQueue->empty()  &&  currentDist < thickness )
    {
        ++ currentDist;

        while( !processingQueue->empty() )
        {
            PixelPos pixel = processingQueue->back();
            processingQueue->pop_back();

            neighborhood.get8( pixel, neighbors );

            for( auto &nn : neighbors )
                if( erosionMap.pixelAbsolute(nn) == 0 )
                {
                    erosionMap.pixelAbsolute( nn ) = currentDist;
                    nextQueue->push_back( nn );
                }
        }

        auto tmp = nextQueue;
        nextQueue = processingQueue;
        processingQueue = tmp;
    }

    // Fill all remaining internal pixels with the thickness value.

    if( currentDist == thickness )
    {
        erosionMap.parallel_for_all_pixels( [&] (ASTex::ImageGrayu32::PixelType& p)
	    {
            if( p == 0 )
                p = thickness;
        });
    }
}


void PatchProcessor::getPatchMap( ASTex::ImageRGBAf &patchMap, uint32_t thickness )
{
    ASTex::ImageGrayu32 erosionMap( idMap_.width(), idMap_.height() );
    getErosionMap( erosionMap, thickness );

    patchMap = ASTex::ImageRGBAf( idMap_.width(), idMap_.height() );

    patchMap.parallel_for_all_pixels( [&] (ASTex::ImageRGBAf::PixelType &pixel, int x, int y)
    {
        PixelPos currentPos = { x, y };

        auto &patch = patchById( idMap_.pixelAbsolute(x,y) );
        PixelPos shiftedCentroid = patch.centroid;

        for( int i=0; i<2; ++i )
        {
            int closestCentroidDist = std::numeric_limits<int>::max();
            for( int c=-1; c<=1; ++c )
            {
                int centroidCoord = patch.centroid[i] + c*patchMap.size()[i];
                int centroidDist = (centroidCoord - currentPos[i])*(centroidCoord - currentPos[i]);
                if( centroidDist < closestCentroidDist )
                {
                    closestCentroidDist = centroidDist;
                    shiftedCentroid[i] = centroidCoord;
                }
            }
        }

        pixel.SetRed  ( (int32_t) patch.id );
        pixel.SetGreen( (int32_t) x - shiftedCentroid[0] );
        pixel.SetBlue ( (int32_t) y - shiftedCentroid[1] );
        pixel.SetAlpha( erosionMap.pixelAbsolute(x,y) / float(thickness) );
    });
}


template <class T>
void BuildImagePyramid( const T &img, std::vector<T> &pyramid, itk::SizeValueType targetSize )
{
    ASTex::Size reducedSize = img.size();
    reducedSize[0] >>= 1;
    reducedSize[1] >>= 1;

    itk::SizeValueType *currentSize = (img.size()[0] < img.size()[1])?  &reducedSize[0] : &reducedSize[1];

    auto transform = itk::ScaleTransform< double, T::ItkImg::ImageDimension >::New();
    transform->SetScale( 2.0 );

	auto interpolator = itk::LinearInterpolateImageFunction<typename T::ItkImg, double >::New();

	typedef itk::ResampleImageFilter<typename T::ItkImg, typename T::ItkImg > ResamplerType;

    auto inputImage = img.itk();

    pyramid.push_back( img );

    while( *currentSize > targetSize )
    {
		typename ResamplerType::Pointer resampler = ResamplerType::New();
        resampler->SetTransform( transform );
        resampler->SetOutputOrigin( 0.5 );
        resampler->SetInterpolator( interpolator );
        resampler->SetInput( inputImage );
        resampler->SetSize( reducedSize );
        resampler->Update();

        inputImage = resampler->GetOutput();
        pyramid.push_back( T(inputImage) );

        reducedSize[0] >>= 1;
        reducedSize[1] >>= 1;
    }
}


double PatchProcessor::computeErrorAt( const ASTex::ImageRGBu8 &image,
                                                     const std::vector<PixelPos> &boundaryPixels,
                                                     const std::vector<Eigen::Vector2d> &transformedBoundaryPixels,
                                                     const Eigen::Vector2d &offset )
{
    double summedError = 0.0;

    for( unsigned int i=0; i<boundaryPixels.size(); ++i )
    {
        const PixelPos &pixel = boundaryPixels[i];
        Eigen::Vector2d transformedPixel = transformedBoundaryPixels[i] + offset + Eigen::Vector2d(0.5,0.5);

        int transformedX = ((int) transformedPixel[0] + 4*image.width ()) % image.width ();
        int transformedY = ((int) transformedPixel[1] + 4*image.height()) % image.height();

        itk::RGBPixel<double> pixelColor = image.pixelAbsolute( pixel );
        itk::RGBPixel<double> transformedColor = image.pixelAbsolute( transformedX, transformedY );

        auto colorDifference = transformedColor - pixelColor;
        double pixelError = colorDifference.GetRed  ()*colorDifference.GetRed  () +
                            colorDifference.GetGreen()*colorDifference.GetGreen() +
                            colorDifference.GetBlue ()*colorDifference.GetBlue ();

        summedError += pixelError;
    }

    return summedError;
}


void PatchProcessor::computeErrorMap( const ASTex::ImageRGBu8 &image,
                                                    const std::vector<PixelPos> &boundaryPixels,
                                                    const Eigen::Matrix2d &transform,
                                                    const Eigen::Vector2d &centroid,
                                                    ASTex::ImageGrayd &errorMap )
{
    std::vector< Eigen::Vector2d > transformedBoundaryPixels;
    transformedBoundaryPixels.reserve( boundaryPixels.size() );

    for( auto &pixel : boundaryPixels )
    {
        Eigen::Vector2d transformedPixel = (transform * (Eigen::Vector2d(pixel[0],pixel[1]) - centroid) );
        transformedBoundaryPixels.push_back( transformedPixel );
    }

    errorMap.parallel_for_all_pixels( [&] (ASTex::ImageGrayd::PixelType &summedError, int offsetX, int offsetY)
    {
        Eigen::Vector2d offset( offsetX, offsetY );
        summedError = computeErrorAt( image, boundaryPixels, transformedBoundaryPixels, offset );
    });
}


void PatchProcessor::extractLocalMinima( const ASTex::ImageGrayd &errorMap,
                                                       double angle,
                                                       double scale,
                                                       unsigned int count,
                                                       std::multiset<PatchContent,std::greater<PatchContent>> &orderedMinima )
{
    ASTex::NeighborhoodPeriodic::NeighborsList neighbors( 8 );
    ASTex::NeighborhoodPeriodic neighborhood( errorMap );

    errorMap.for_all_pixels( [&] (const ASTex::ImageGrayd::PixelType &errVal, int x, int y) {
        PixelPos pos = { x, y };

        neighborhood.get8( pos, neighbors );

        for( auto nn : neighbors )
            if( errorMap.pixelAbsolute(nn) < errVal )
                return;

        PatchContent content;
        content.error  = errVal;
        content.offset = pos;
        content.angle  = angle;
        content.scale  = scale;

        orderedMinima.insert( content );
        if( orderedMinima.size() > count )
            orderedMinima.erase( orderedMinima.begin() );
    });
}


void PatchProcessor::downSamplePatch( const std::vector<ASTex::ImageRGBu8> &imagePyramid,
													const Patch &patch,
                                                    unsigned int downsampleLevel,
                                                    Eigen::Vector2d &downsampledCentroid,
                                                    std::vector<PixelPos> &downsampledBoundary )
{
    if( downsampleLevel == 0 )
    {
        downsampledCentroid[0] = patch.centroid[0];
        downsampledCentroid[1] = patch.centroid[1];
        downsampledBoundary = patch.boundary;
    }
    else
    {
        double downsampleFactor( 1 << downsampleLevel );

        downsampledCentroid[0] = patch.centroid[0] / downsampleFactor;
        downsampledCentroid[1] = patch.centroid[1] / downsampleFactor;

        downsampledBoundary.clear();

        ASTex::ImageGrayu8 flagMap( imagePyramid[downsampleLevel].width(), imagePyramid[downsampleLevel].height() );
        flagMap.parallel_for_all_pixels( [] (ASTex::ImageGrayu8::PixelType &p) { p = 0; } );

        for( auto &pixel : patch.boundary )
        {
            PixelPos downsampledPixel = {
                pixel[0] >> downsampleLevel,
                pixel[1] >> downsampleLevel
            };

            if( flagMap.pixelAbsolute(downsampledPixel) == 0 )
            {
                flagMap.pixelAbsolute(downsampledPixel) = 1;
                downsampledBoundary.push_back( downsampledPixel );
            }
        }
    }
}


void PatchProcessor::findAlternativeContents( /*const ASTex::ImageRGBu8 &sourceImage,*/
                                                            unsigned int nContents,
                                                            unsigned int downSamplingMinSize )
{
    std::vector<double> rotationAngles;
    rotationAngles.push_back( 0.0 );

    std::vector<double> scalingFactors;
    scalingFactors.push_back( 1.0 );

	findAlternativeContents( /*sourceImage,*/ rotationAngles, scalingFactors, nContents, downSamplingMinSize );
}


void PatchProcessor::findAlternativeContents( /*const ASTex::ImageRGBu8 &sourceImage,*/
                                                            const std::vector<double> &rotationAngles,
                                                            const std::vector<double> &scalingFactors,
                                                            unsigned int nContents,
                                                            unsigned int downSamplingMinSize )
{
    std::vector< ASTex::ImageRGBu8 > imagePyramid;
    BuildImagePyramid( fragmentProc_.sourceImage(), imagePyramid, downSamplingMinSize );


    ASTex::ImageGrayd errorMap( imagePyramid.back().width(), imagePyramid.back().height() );
    NeighborhoodSearchGrid safeRadiusSearchGrid( errorMap.width(), errorMap.height(), 4.0f );

    Eigen::Vector2d downsampledCentroid;
    std::vector<PixelPos> downsampledBoundary;

    unsigned int downsampleLevel = (int) imagePyramid.size() - 1;
    unsigned int downsampleFactor = (1 << downsampleLevel);

    // The process starts from the coarsest level of the image pyramid. An exhaustive search is performed by computing
    // a SSD along the patch boundary for every possible location in the image, and every provided rotation and scale.
    // Only a given number of local minima are retained as candidates for alternative contents.

std::cout << std::endl << "  * Dense search at lowest resolution:   0%" << std::flush;
unsigned int ppp = 0;
    for( auto &patch : allPatches_ )
    {
        downSamplePatch( imagePyramid, patch, downsampleLevel, downsampledCentroid, downsampledBoundary );

        patch.contents.clear();
        std::multiset<PatchContent,std::greater<PatchContent>> orderedMinima;

        for( auto angle : rotationAngles )
        {
            const double cosAngle = std::cos( -angle );
            const double sinAngle = std::sin( -angle );

            for( auto scale : scalingFactors )
            {
                Eigen::Matrix2d transform;
                transform(0,0) =  cosAngle / scale;
                transform(1,0) =  sinAngle / scale;
                transform(0,1) = -transform(1,0);
                transform(1,1) =  transform(0,0);

                computeErrorMap( imagePyramid.back(), downsampledBoundary, transform, downsampledCentroid, errorMap );

                extractLocalMinima( errorMap, angle, scale, 4*scalingFactors.size()*rotationAngles.size()*nContents, orderedMinima );
printf( "\b\b\b\b%3i%%", 100 * (++ppp) / int(patchCount() * rotationAngles.size() * scalingFactors.size()) );
            }
        }

        safeRadiusSearchGrid.clear();

        for( auto m=orderedMinima.rbegin(); m!=orderedMinima.rend(); ++m )
        {
            Eigen::Vector2d minimumLocation( m->offset[0], m->offset[1] );
            if( safeRadiusSearchGrid.isNeighborhoodFree( minimumLocation ) )
            {
                safeRadiusSearchGrid.put( minimumLocation );
                patch.contents.push_back( *m );
            }
        }
    }
printf( "\n" );


    // These candidates are then refined by a coarse-to-fine approach: each of them is up-sampled and a better solution
    // is looked for by computing its error again in a small neighborhood in the lower pyramid level. Refinement is
    // repeated until the lowest level is reached.

std::cout << "  * Hierarchical refinement:   0%" << std::flush;
ppp = 0;
    while( downsampleLevel > 0 )
    {
        -- downsampleLevel;
        downsampleFactor >>= 1;

        ASTex::NeighborhoodPeriodic::NeighborsList neighbors( 8 );
        ASTex::NeighborhoodPeriodic neighborhood( imagePyramid[downsampleLevel] );

        for( auto &patch : allPatches_ )
        {
            downSamplePatch( imagePyramid, patch, downsampleLevel, downsampledCentroid, downsampledBoundary );

            std::vector< Eigen::Vector2d > transformedBoundaryPixels;
            transformedBoundaryPixels.reserve( downsampledBoundary.size() );

            for( auto &m : patch.contents )
            {
                Eigen::Matrix2d transform;
                transform(0,0) =  std::cos(m.angle) / m.scale;
                transform(1,0) =  std::sin(m.angle) / m.scale;
                transform(0,1) = -transform(1,0);
                transform(1,1) =  transform(0,0);

                transformedBoundaryPixels.clear();
                for( auto &pixel : downsampledBoundary )
                {
                    Eigen::Vector2d transformedPixel = (transform * (Eigen::Vector2d(pixel[0],pixel[1]) - downsampledCentroid) );
                    transformedBoundaryPixels.push_back( transformedPixel );
                }

                m.offset[0] <<= 1;
                m.offset[1] <<= 1;
                m.error = computeErrorAt( imagePyramid[downsampleLevel], downsampledBoundary, transformedBoundaryPixels, Eigen::Vector2d(m.offset[0],m.offset[1]) );

                neighborhood.get8( m.offset, neighbors );

                for( auto &nn : neighbors )
                {
                    double nnError = computeErrorAt( imagePyramid[downsampleLevel], downsampledBoundary, transformedBoundaryPixels, Eigen::Vector2d(nn[0],nn[1]) );
                    if( nnError < m.error )
                    {
                        m.error = nnError;
                        m.offset = nn;
                    }
                }
            }
printf( "\b\b\b\b%3i%%", 100 * (++ppp) / int(patchCount()*(imagePyramid.size()-1)) );
        }
    }
printf( "\n" );

    
    // Based on the updated error values, only the desired number of candidates is eventually retained.

    for( auto &patch : allPatches_ )
    {
        std::partial_sort( patch.contents.begin(), patch.contents.begin()+nContents, patch.contents.end() );
        patch.contents.resize( nContents );
    }
}

}
}
