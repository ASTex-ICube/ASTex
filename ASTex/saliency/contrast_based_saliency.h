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



#ifndef __ASTEX__CONTRAST_BASED_SALIENCY_H__
#define __ASTEX__CONTRAST_BASED_SALIENCY_H__




#include <ASTex/slic.h>
#include <ASTex/saliency/permutohedral_lattice.h>


namespace ASTex
{


/**
 *  Filter performing the labelling of an image according to the SLIC superpixel algorithm.
 */
template < typename TOutputImage >
class ContrastBasedSaliencyFilter : public itk::ImageToImageFilter< ASTex::ImageRGBu8::ItkImg, TOutputImage >
{
public:
	// Filter is not template so define input & output types for more simple writing
	using TInputImage  = ASTex::ImageRGBu8::ItkImg;
	using TInputPixel  = TInputImage::PixelType;
	using TOutputPixel = typename TOutputImage::PixelType;

	// Standard class typedefs & macros for ikt
	typedef ContrastBasedSaliencyFilter                             Self;
	typedef itk::ImageToImageFilter< TInputImage, TOutputImage >    Superclass;
	typedef itk::SmartPointer< Self >                               Pointer;
	typedef itk::SmartPointer< const Self >                         ConstPointer;
	itkNewMacro( Self );
	itkTypeMacro( ContrastBasedSaliencyFilter, ImageToImageFilter );

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

    void setUniquenessVariance( double variance )
    {
        uniquenessVariance_ = variance;
    }

    void setDistributionVariance( double variance )
    {
        distributionVariance_ = variance;
    }

    void setExponentialScaling( double k )
    {
        k_ = k;
    }

    void setAlphaBeta( double alpha, double beta )
    {
        alpha_ = alpha;
        beta_  = beta;
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
    template <typename PixelType>
    inline void CheckOutputType()
    {
		static_assert(std::is_floating_point<PixelType>::value,"PixelType must be float or double");
//        const bool OUTPUT_IMAGE_FORMAT_MUST_BE_GRAYSCALE_FLOAT_OR_DOUBLE = false;
//        assert( OUTPUT_IMAGE_FORMAT_MUST_BE_GRAYSCALE_FLOAT_OR_DOUBLE );
    }

//    template <> inline void CheckOutputType<float >() {}
//    template <> inline void CheckOutputType<double>() {}


    /** Functor converting a pixel from RGB colorspace (ranging [0-255]) to LAB colorspace.
     */
    template< class TInputOutput>
    class functorLABtoRGB : ASTex::ColorSpace::fonctorRGB01To255<TInputOutput,TInputOutput>
    {
    public:
    	inline TInputOutput operator()( const TInputOutput &p ) const
    	{
            TInputOutput q = ASTex::ColorSpace::colorLABtoXYZ<TInputOutput,TInputOutput>( p );
            q = ASTex::ColorSpace::colorXYZtoRGB<TInputOutput,TInputOutput>( q );
		    q = ASTex::ColorSpace::fonctorRGB01To255<TInputOutput,TInputOutput>::operator()( q );
            return q;
	    }
    };

    using SLICFilterType    = ASTex::SLICSuperPixelFilter< ASTex::NeighborhoodClamped >;
    template <int D>
    using BlurData          = Eigen::Matrix< double, D, 1 >;

    /// filter data
	itk::Offset<2>          center_;
	long                    radius_;

    int                     superPixelSize_;
    unsigned int            nIterations_;
    double                  compactness_;

    double                  uniquenessVariance_;
    double                  distributionVariance_;
    double                  k_;
    double                  alpha_;
    double                  beta_;

	/// protected constructor (forbid usage of new and variable declaration)
	ContrastBasedSaliencyFilter() :
		center_({{0,0}}),
		radius_(1),
        superPixelSize_(0),
        nIterations_(10),
        compactness_(10),
        uniquenessVariance_(0.25),
        distributionVariance_(20.0),
        k_(3.0),
        alpha_( 1.0 / 30.0 ),
        beta_( 1.0 / 30.0 )
	{}

	virtual ~ContrastBasedSaliencyFilter() {}


    inline SLICFilterType::Pointer getImageAbstraction()
    {
        SLICFilterType::Pointer slicFilter = SLICFilterType::New();
		slicFilter->SetInput( this->GetInput() );
        slicFilter->setSuperPixelSize( superPixelSize_ );
        slicFilter->setIterationNumber( nIterations_ );
        slicFilter->setCompactness( compactness_ );
        slicFilter->Update();
        return slicFilter;
    }

    inline void getUniquenessValues( const SLICFilterType::Pointer slicFilter,
                                     std::vector<double> &uniqueness )
    {
        // Allocate uniqueness filter value buffer.

        const unsigned int nSuperPixels = slicFilter->getSuperPixelCount();
        uniqueness.resize( nSuperPixels );

        // Compute the sample position scaling factor based on the gaussian kernel variance.

        itk::Size<2> imgSize = slicFilter->GetOutput()->GetLargestPossibleRegion().GetSize();
        const double uniquenessScaling  = 1.0 / (uniquenessVariance_ * std::max(imgSize[0],imgSize[1]));

        // Create the permutohedral lattice used for gaussian blurring speedup.

        PermutohedralLattice< double, 2, 5 > pLattice;
        pLattice.reserveRecordings( nSuperPixels );

        // 1st step: splatting data to lattice vertices.

        for( unsigned int i=0; i<nSuperPixels; ++i )
        {
            auto &sp = slicFilter->getSuperPixel( i );

            BlurData<5> value;
            value[0] = 1.0;
            value[1] = sp.color[0];
            value[2] = sp.color[1];
            value[3] = sp.color[2];
            value[4] = sp.color.dot( sp.color );

            pLattice.splatAndDeclare( sp.coord * uniquenessScaling, value );
        }

        // 2nd step: blurring of values at lattice vertices by applying a serie of separable 1D gaussian kernels.

        pLattice.blur();

        // 3rd step: slicing lattice vertex blurred values at data locations.

        for( unsigned int i=0; i<nSuperPixels; ++i )
        {
            BlurData<5> blurred = pLattice.slice( i );
            auto &c = slicFilter->getSuperPixel( i ).color;

            uniqueness[i] = blurred[0]*c.dot(c)
                          - 2.0*c.dot( Eigen::Vector3d(blurred[1],blurred[2],blurred[3]) )
                          + blurred[4];
        }

        // Normalize values.

        normalize( uniqueness );
    }

    inline void getDistributionValues( const SLICFilterType::Pointer slicFilter,
                                       std::vector<double> &distribution )
    {
        // Allocate distribution filter value buffer.

        const unsigned int nSuperPixels = slicFilter->getSuperPixelCount();
        distribution.resize( nSuperPixels );

        // Compute the color scaling factor based on the gaussian kernel variance.

        const double distributionScaling  = 1.0 / distributionVariance_;

        // Create the permutohedral lattice used for gaussian blurring speedup.

        PermutohedralLattice< double, 3, 4 > pLattice;
        pLattice.reserveRecordings( nSuperPixels );

        // 1st step: splatting data to lattice vertices.

        for( unsigned int i=0; i<nSuperPixels; ++i )
        {
            auto &sp = slicFilter->getSuperPixel( i );

            BlurData<4> value;
            value[0] = 1.0;
            value[1] = sp.coord[0];
            value[2] = sp.coord[1];
            value[3] = sp.coord.dot( sp.coord );

            pLattice.splatAndDeclare( sp.color * distributionScaling, value );
        }

        // 2nd step: blurring of values at lattice vertices by applying a serie of separable 1D gaussian kernels.

        pLattice.blur();

        // 3rd step: slicing lattice vertex blurred values at data locations.

        for( unsigned int i=0; i<nSuperPixels; ++i )
        {
            BlurData<4> blurred = pLattice.slice( i );
            //distribution[i] = blurred[3] / blurred[0]
            //                - (blurred[1]*blurred[1] + blurred[2]*blurred[2]) / (blurred[0]*blurred[0]);
            distribution[i] = blurred[3]*blurred[0] - blurred[1]*blurred[1] - blurred[2]*blurred[2];
        }

        // Normalize values.

        normalize( distribution );
    }

    inline void normalize( std::vector<double> &values )
    {
        double minValue = values[0], maxValue = values[0];

        for( auto v : values )
        {
            if( v < minValue )
                minValue = v;
            else if( v > maxValue )
                maxValue = v;
        }

        double scale = 1.0 / (maxValue - minValue);

        for( auto &v : values )
            v = (v - minValue) * scale;
    }

    inline void assignSaliency( SLICFilterType::Pointer slicFilter,
                                std::vector<double> &uniqueness,
                                std::vector<double> &distribution )
    {
		auto inputImg = this->GetInput();
        itk::Size<2> imgSize = inputImg->GetLargestPossibleRegion().GetSize();

        auto superPixelMap = slicFilter->GetOutput();

        // Convert the mean colors of superpixels from LAB to RGB colorspace.

        std::vector< Eigen::Vector3d > superPixelRGBColors( slicFilter->getSuperPixelCount() );
		functorLABtoRGB< itk::RGBPixel<double> > convertLABtoRGB;

        for( unsigned int i=0; i<slicFilter->getSuperPixelCount(); ++i )
        {
            auto &c = slicFilter->getSuperPixel(i).color;
            itk::RGBPixel<double> cRGB;
            cRGB.Set( c[0], c[1], c[2] );
            cRGB = convertLABtoRGB( cRGB );
            superPixelRGBColors[i] = Eigen::Vector3d( cRGB[0], cRGB[1], cRGB[2] ) * beta_;
        }

        // Compute saliency values for each superpixel.

        auto &superPixelSaliency = uniqueness;

        for( unsigned int i=0; i<superPixelSaliency.size(); ++i )
            superPixelSaliency[i] = uniqueness[i] * std::exp( -k_ * distribution[i] );

        // Splat superpixel saliency map pixel values.

        PermutohedralLattice< double, 5, 2 > pLattice;

        BlurData<5> point;
        BlurData<2> value;
        value[0] = 1.0;

		for( uint32_t y=0; y<imgSize[1]; ++y )
			for( uint32_t x=0; x<imgSize[0]; ++x )
            {
                auto spId = superPixelMap->GetPixel( {x,y} );
                auto &c = superPixelRGBColors[spId];

                point[0] = x * alpha_;
                point[1] = y * alpha_;
                point[2] = c[0];
                point[3] = c[1];
                point[4] = c[2];

                value[1] = superPixelSaliency[spId];

                pLattice.splatOnly( point, value );
            }

        // Declare the position where resulting values will be red in order to enforce the construction lattice vertices in their neighborhood.

        pLattice.reserveRecordings( imgSize[0]*imgSize[1] );

		for( uint32_t y=0; y<imgSize[1]; ++y )
			for( uint32_t x=0; x<imgSize[0]; ++x )
            {
                auto &c = inputImg->GetPixel( {x,y} );

                point[0] = x * alpha_;
                point[1] = y * alpha_;
                point[2] = c[0] * beta_;
                point[3] = c[1] * beta_;
                point[4] = c[2] * beta_;

                pLattice.declareOnly( point );
            }

        // Perform blurring.

        pLattice.blur();

        // Recover blurred values at each pixel of the input image.

		auto saliencyMap = this->GetOutput();

        auto saliencyMin =  std::numeric_limits<TOutputPixel>::max();
        auto saliencyMax = -std::numeric_limits<TOutputPixel>::max();

		for( uint32_t y=0, n=0; y<imgSize[1]; ++y )
			for( uint32_t x=0; x<imgSize[0]; ++x, ++n )
            {
                BlurData<2> blurred = pLattice.slice( n );
                TOutputPixel saliency = TOutputPixel( blurred[1] / (blurred[0] + 1e-10) );

                if( saliency < saliencyMin )
                    saliencyMin = saliency;
                if( saliency > saliencyMax )
                    saliencyMax = saliency;

                saliencyMap->SetPixel( {x,y}, saliency );
            }

		for( uint32_t y=0, n=0; y<imgSize[1]; ++y )
			for( uint32_t x=0; x<imgSize[0]; ++x, ++n )
            {
                TOutputPixel &saliency = saliencyMap->GetPixel( {x,y} );
                saliency = (saliency - saliencyMin) / (saliencyMax - saliencyMin);
            }
   }


	//
	// overriden method that generate the output
	//
	void GenerateData() ITK_OVERRIDE
	{
        CheckOutputType<TOutputPixel>();
        
        //  we are responsible of allocating the images buffers
		// in this simple case (output of same size than inputs) just call:
		this->AllocateOutputs();

        SLICFilterType::Pointer slicFilter = getImageAbstraction();

        std::vector<double> uniqueness;
        getUniquenessValues( slicFilter, uniqueness );

        std::vector<double> distribution;
        getDistributionValues( slicFilter, distribution );

        assignSaliency( slicFilter, uniqueness, distribution );
    }

private:
	// to avoid filter object copy
	ContrastBasedSaliencyFilter( const Self & );
	void operator=( const Self & );
};


} // namespace ASTex




#endif // __ASTEX__CONTRAST_BASED_SALIENCY_H__
