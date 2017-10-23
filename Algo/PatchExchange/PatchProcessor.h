#ifndef __CONTENTEXCHG__PATCHPROCESSOR__
#define __CONTENTEXCHG__PATCHPROCESSOR__




#include "FragmentProcessor.h"
#include "ASTex/image_rgba.h"

namespace ASTex
{

namespace ContentExchg
{


struct PatchContent
{
    double                      error;
    PixelPos                    offset;
    double                      angle;
    double                      scale;
    inline bool operator<( const PatchContent &m ) const  { return error < m.error; }
    inline bool operator>( const PatchContent &m ) const  { return error > m.error; }
};

struct Patch
{
    PixelPos                    centroid;
    uint32_t                    id;
    std::vector<uint32_t>       fragments;
    std::vector<PixelPos>       boundary;
    std::vector<PatchContent>   contents;
};

class PatchProcessor
{
    const FragmentProcessor&    fragmentProc_;
    std::vector<Patch>          allPatches_;
	ASTex::ImageGrayu32         idMap_;

    double                      computeErrorAt( const ASTex::ImageRGBu8 &image,
                                                const std::vector<PixelPos> &boundaryPixels,
                                                const std::vector<Eigen::Vector2d> &transformedBoundaryPixels,
                                                const Eigen::Vector2d &offset );

    void                        computeErrorMap( const ASTex::ImageRGBu8 &image,
                                                 const std::vector<PixelPos> &boundaryPixels,
                                                 const Eigen::Matrix2d &transform,
                                                 const Eigen::Vector2d &centroid,
                                                 ASTex::ImageGrayd &errorMap );

    void                        extractLocalMinima( const ASTex::ImageGrayd &errorMap,
                                                    double angle,
                                                    double scale,
                                                    unsigned int count,
                                                    std::multiset<PatchContent,std::greater<PatchContent>> &orderedMinima );

    void                        downSamplePatch( const std::vector<ASTex::ImageRGBu8> &imagePyramid,
                                                 const ContentExchg::Patch &patch,
                                                 unsigned int downsampleLevel,
                                                 Eigen::Vector2d &downsampledCentroid,
                                                 std::vector<PixelPos> &downsampledBoundary );

    void                        getErosionMap( ASTex::ImageGrayu32 &erosionMap, uint32_t thickness );

public:
    PatchProcessor( const FragmentProcessor &fragmentProc );

	int                         patchCount() const;
	Patch&                      patchById( int pid );
	const Patch&                patchById( int pid ) const;
    std::vector<Patch>&         patches();
    const std::vector<Patch>&   patches() const;

    void                        createPatches( int requiredPatchNumber );
    void                        computePatchBoundaries();
	void                        findAlternativeContents( unsigned int nContents,
                                                         unsigned int downSamplingMinSize );
	void                        findAlternativeContents( const std::vector<double> &rotationAngles,
                                                         const std::vector<double> &scalingFactors,
                                                         unsigned int nContents,
                                                         unsigned int downSamplingMinSize );

    void                        getPatchMap( ASTex::ImageRGBAf &patchMap, uint32_t thickness );
};

} // namespace ContentExchg
}



#endif //__CONTENTEXCHG__PATCHPROCESSOR__
