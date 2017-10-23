#ifndef __CONTENTEXCHG__FRAGMENTPROCESSOR__
#define __CONTENTEXCHG__FRAGMENTPROCESSOR__




#include <ASTex/region_growing/region_growing.h>
#include <ASTex/image_rgb.h>
#include <Eigen/Eigen>



namespace ASTex
{

namespace ContentExchg
{


using PixelPos = itk::Index<2>;

struct Fragment
{
    uint32_t                id;
    Eigen::Vector2d         centroid;
    std::vector<PixelPos>   pixels;
    std::vector<uint32_t>   neighbors;
};


/** \class FragmenterProcessor managing image fragmentation given three main parameters:
 *	min, max fragment size and threshold value defining when pixels should be separated into two fragments.
 */
class FragmentProcessor
{
private:
    using ColorF   = itk::RGBPixel<float>;
    using ColorU32 = itk::RGBPixel<uint32_t>;
    using ColorU8  = itk::RGBPixel<uint8_t>;

	bool                            isFragmented_;

    const ASTex::ImageRGBu8         &sourceImage_;
	std::vector<Fragment>           allFragments_;
	ASTex::ImageGrayu32             idMap_;

	static float                    colorSquareDist( const ColorF &c1, const ColorF &c2 );

public:
	/** \brief Constructor stores the given image and initializes internal parameters
	 *  \param image            The RGBu8 image that must be fragmented.
	 */
	FragmentProcessor( const ASTex::ImageRGBu8 &image );

    /** \brief Recover the source image to which the fragment processor has to work on.
     *  \return                 Source image of the fragment processor.
     */
    const ASTex::ImageRGBu8&        sourceImage() const { return sourceImage_; }
    /** \brief Recover the number of created fragments (zero if the createFragments() function has not been called yet).
     *  \return                 Number of fragments.
     */
	int                             fragmentCount() const;
    /** \brief Recover the ID of the fragment the given pixel belongs to.
     *  \param x                Image horizontal coordinate of the pixel.
     *  \param y                Image vertical coordinate of the pixel.
     *  \return                 ID of the correponding fragment.
     */
	int                             fragmentIdAt( int x, int y ) const;
    /** \brief Recover the fragment corresponding to the provided ID.
     *  \param fid              ID of the desired fragment.
     *  \return                 Fragment corresponding to the given ID.
     */
	const Fragment&                 fragmentById( int fid ) const;
    /** \brief Recover the fragment the given pixel belongs to.
     *  \param x                Image horizontal coordinate of the pixel.
     *  \param y                Image vertical coordinate of the pixel.
     *  \return                 Correponding fragment.
     */
	const Fragment&                 fragmentAt( int x, int y ) const;
    /** \brief Recover the map of fragment IDs.
     *  \return                 Image of fragment IDs.
     */
    const ASTex::ImageGrayu32&      idMap() const;
    /** \brief Recover the set of created fragments.
     *  \return                 Vector of all created fragments.
     */
    const std::vector<Fragment>&    fragments() const;

	/** \brief Classify the given image into a list of fragments by a region-growing algorithm.
	 *
	 *  \param fragmentMaxSize  The maximal fragment size.
	 *  \param threshold        The color distance over which neighboring pixels will be separated: between 30 and 100.
     *  \parma labelMap         A map defining per-pixel labels. If provided, fragments a created so that they cannot span over different labels.
	 *
	 *  \return                 The number of created fragments.
	 */
    unsigned int                    createFragments( unsigned int fragmentMaxSize, unsigned int threshold, ASTex::ImageGrayu8 *labelMap = NULL );

	/** \brief Removes fragments that are too small by dispatching their pixels into adjacent fragments.
     *  Warning : no check is done on previously defined max size used for fragmentation.
	 *
	 *  \param fragmentMinSize  The minimum size allowed for a fragment.
     *
	 *  \return                 The new number of fragments.
	 */
    unsigned int                    cleanupSmallestFragments( unsigned int fragmentMinSize );

    /** \brief Computes additional fragment attributes, such as centroid and neighborhood relationship. */
    void                            updateFragmentsAttributes();
};



} // namespace ContentExchg
}



#endif // __CONTENTEXCHG__FRAGMENTPROCESSOR__
