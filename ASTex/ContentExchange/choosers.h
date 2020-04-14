#ifndef __CTEXCH_CHOOSER_H__
#define __CTEXCH_CHOOSER_H__

#include "patchProcessor.h"

namespace ASTex
{

namespace ContentExchange
{

//
//Base classes
//

template<typename I>
class PatchProcessor<I>::ChooserBase
{
public:
    ChooserBase()
    {}

protected:

    virtual void generate() = 0;
};

template<typename I>
class PatchProcessor<I>::ChooserPatches : public PatchProcessor<I>::ChooserBase
{
    using Base = PatchProcessor<I>::ChooserBase;

public:
    ChooserPatches() {}

    //Parameters

    void setNbPatches(unsigned requestedNbPatches) {m_requestedNbPatches = requestedNbPatches;}

    //Generator

    virtual void generate() = 0;

protected:

    unsigned m_requestedNbPatches;
};

template<typename I>
class PatchProcessor<I>::ChooserContents : public PatchProcessor<I>::ChooserBase
{
public:
    using Base = PatchProcessor<I>::ChooserBase;

    ChooserContents() {}

    //Parameters

    void setNbContents(unsigned nbContents) {nbContents = nbContents;}

    //Generator

    virtual void generate() = 0;

protected:

    unsigned m_nbContents;
};

template<typename I>
class PatchProcessor<I>::ChooserPatchesAndContents : public PatchProcessor<I>::ChooserBase
{
public:
    using Base = PatchProcessor<I>::ChooserBase;

    ChooserPatchesAndContents() {}

    //Parameters

    void setNbPatches(unsigned requestedNbPatches) {m_requestedNbPatches = requestedNbPatches;}
    void setNbContents(unsigned nbContents) {nbContents = nbContents;}

    //Generator

    virtual void generate() = 0;

protected:

    unsigned m_requestedNbPatches;
    unsigned m_nbContents;
};

//
//Overrides of ChooserPatches
//

//ChooserPatches_FromImage

template<typename I>
template <typename IPatchmap>
class PatchProcessor<I>::ChooserPatches_FromImage : public PatchProcessor<I>::ChooserPatches
{
public:
    using Base = PatchProcessor<I>::ChooserPatches;

    ChooserPatches_FromImage() :
        PatchProcessor<I>::ChooserPatches(),
        m_patchmap()
    {}

    //Parameters

    void setPatchImage(const IPatchmap &patchmap);

    //Generator

    void generate() override;

private:
    IPatchmap m_patchmap;
};

template<typename I>
template<typename IPatchmap>
void PatchProcessor<I>::ChooserPatches_FromImage<IPatchmap>::setPatchImage(const IPatchmap &patchmap)
{
    m_patchmap = patchmap;
}

template<typename I>
template<typename IPatchmap>
void PatchProcessor<I>::ChooserPatches_FromImage<IPatchmap>::generate()
{
    assert(this &&
           "Base::generate: a patchProcessor must have been set (it can only be set by the PatchProcessor class)");
    assert(m_patchmap.size() == this->m_texture.size() &&
           "ChooserPatches_OverlappingCircles::generate: patchImage should have the same size as texture");
    Histogram<IPatchmap> histogramPI(m_patchmap);
    assert(histogramPI.binsNumber()<=CTEXCH_MAX_PATCHES &&
           "PatchProcessor::debug_setPatchFromImageRGBd: there should be up to 64 different colors in the given image");
    word64 w=0x1;
    ImageMask64 patchMask;
    patchMask.initItk(m_patchmap.width(), m_patchmap.height(), true);
    for(auto it=histogramPI.begin(); it!=histogramPI.end(); ++it)
    {
        m_patchmap.for_all_pixels([&] (const typename IPatchmap::PixelType &pix, int x, int y)
        {
            if((*it).first==pix)
            {
                word64 &wPixel=reinterpret_cast<word64&>(patchMask.pixelAbsolute(x, y));
                wPixel |= w;
            }
        });
        w*=2;
    }
    this->m_patchMapMipmap.setTexture(patchMask);
    this->m_patchMapMipmap.setMode(this->m_mipmapMode);
    this->m_patchMapMipmap.generate();
}

//ChooserPatches_RegularGrid

template<typename I>
class PatchProcessor<I>::ChooserPatches_RegularGrid : public PatchProcessor<I>::ChooserPatches
{
public:
    using Base = PatchProcessor<I>::ChooserPatches;

    ChooserPatches_RegularGrid() :
        PatchProcessor<I>::ChooserPatches() {}

    //Parameters

    //

    //Generator

    void generate() override;
};

template<typename I>
void PatchProcessor<I>::ChooserPatches_RegularGrid::generate()
{
    assert(this &&
           "Base::generate: a patchProcessor must have been set (it can only be set by the PatchProcessor class)");
    assert(this->m_texture.is_initialized() &&
           "PatchProcessor::initializePatches: texture uninitialized (use PatchProcessor::setTexture with an initialized texture)");
    assert(m_requestedNbPatches <= CTEXCH_MAX_PATCHES &&
           "PatchProcessor::initializePatches: synthesis cannot allow that many patches (max given by CTEXCH_MAX_PATCHES)");
    srand(this->m_seed);
    double nbPixelsInWidthPerPatch, nbPixelsInHeightPerPatch;
    int nbPatchesPerSide = int(std::sqrt(Base::m_requestedNbPatches));
    ImageMask64 patchMap;
    patchMap.initItk(this->m_texture.width(), this->m_texture.height());
    nbPixelsInWidthPerPatch = this->m_texture.width()/double(nbPatchesPerSide);
    nbPixelsInHeightPerPatch = this->m_texture.height()/double(nbPatchesPerSide);

    patchMap.for_all_pixels([&] (ImageMask64::PixelType &pix, int x, int y)
    {
        int xId = x/nbPixelsInWidthPerPatch;
        int yId = y/nbPixelsInHeightPerPatch;
        int id = yId * nbPatchesPerSide + xId;
        pix = uint64_t(std::pow(2, id));
    });

    this->m_patchMapMipmap.setTexture(patchMap);
    this->m_patchMapMipmap.setMode(this->m_mipmapMode);
    this->m_patchMapMipmap.generate();
}

//ChooserPatches_OverlappingCircles

template<typename I>
class PatchProcessor<I>::ChooserPatches_OverlappingCircles : public PatchProcessor<I>::ChooserPatches
{
public:
    using Base = PatchProcessor<I>::ChooserPatches;

    ChooserPatches_OverlappingCircles() :
        PatchProcessor<I>::ChooserPatches() {}

    //Parameters

    //

    //Generator

    void generate() override;
};

template<typename I>
void PatchProcessor<I>::ChooserPatches_OverlappingCircles::generate()
{
    assert(this &&
           "Base::generate: a patchProcessor must have been set (it can only be set by the PatchProcessor class)");
    assert(this->m_texture.is_initialized() &&
           "PatchProcessor: texture uninitialized (use PatchProcessor::setTexture with an initialized texture)");
    assert(m_requestedNbPatches <= CTEXCH_MAX_PATCHES &&
           "ChooserPatches_OverlappingCircles::generate: synthesis cannot allow that many patches (max given by CTEXCH_MAX_PATCHES)");
    srand(this->m_seed);
    Stamping::SamplerPoisson sampler(Base::m_requestedNbPatches);
    std::vector<Eigen::Vector2f> centroids = sampler.generate();
    ImageMask64 patchMap;
    patchMap.initItk(this->m_texture.width(), this->m_texture.height());
    float r = this->m_texture.width() * this->m_texture.height() / float(Base::m_requestedNbPatches*Base::m_requestedNbPatches*2);
    int i = 0;
    for(std::vector<Eigen::Vector2f>::const_iterator cit = centroids.begin(); cit != centroids.end(); ++cit, ++i)
    {
        word64 w=0x1;
        w <<= i;
        patchMap.for_all_pixels([&] (ImageMask64::PixelType &pix, int x, int y)
        {
            int x0, y0;
            x0 = (*cit)[0]*this->m_texture.width();
            y0 = (*cit)[1]*this->m_texture.height();
            if(std::sqrt((x0-x) * (x0-x) + (y0-y) * (y0-y)) < r)
                reinterpret_cast<word64&>(pix) |= w;
        });
    }

    this->m_patchMapMipmap.setTexture(patchMap);
    this->m_patchMapMipmap.setMode(this->m_mipmapMode);
    this->m_patchMapMipmap.generate();
}

//ChooserPatches_RandomPoisson

template<typename I>
class PatchProcessor<I>::ChooserPatches_RandomPoisson : public PatchProcessor<I>::ChooserPatches
{
public:
    using Base = PatchProcessor<I>::ChooserPatches;

    ChooserPatches_RandomPoisson() :
        ChooserPatches() {}

    //Parameters

    //

    //Generator

    void generate() override;
};

template<typename I>
void PatchProcessor<I>::ChooserPatches_RandomPoisson::generate()
{
    assert(this &&
           "Base::generate: a patchProcessor must have been set (it can only be set by the PatchProcessor class)");
    //patch mask initialization
    assert(m_requestedNbPatches <= CTEXCH_MAX_PATCHES);
    static ImageMask64::PixelType zero;
    ImageMask64 patchMap;
    patchMap.initItk(this->m_texture.width(), this->m_texture.height(), true);
    patchMap.for_all_pixels([&] (ImageMask64::PixelType &pix)
    {
        pix = zero;
    });

    //sampler poisson but with the correct number of samples
    std::vector< Eigen::Vector2d > samples;
    std::vector< Eigen::Vector2d > samplesValidated;
    double minDist = std::sqrt(this->m_texture.width()*this->m_texture.height())/(2.0*std::sqrt(Base::m_requestedNbPatches) + 1);
    ContentExchg::PoissonDiskSampler::generateSamples(this->m_texture.width(), this->m_texture.height(), samples, minDist, Base::m_requestedNbPatches);
    std::minstd_rand minstdRandGenerator;
    minstdRandGenerator.seed(this->m_seed+1);
    assert(samples.size()>=Base::m_requestedNbPatches);
    for(size_t i=0; i<Base::m_requestedNbPatches; ++i)
    {
        std::uniform_int_distribution<size_t> distribution(0, samples.size()-1);
        size_t k = distribution(minstdRandGenerator);
        samplesValidated.push_back(samples[k]);
        samples.erase(samples.begin()+int(k));
    }

    //region growing
    unsigned *grownTexelCountdownVector = new unsigned[Base::m_requestedNbPatches]();
    unsigned totalNumberOfMarkedPixels = 0;
    for(unsigned i=0; i<Base::m_requestedNbPatches; ++i)
    {
        unsigned maxTexelsForPatchi = ((0.5+(rand()/double(RAND_MAX)))*this->m_texture.width()*this->m_texture.height())/Base::m_requestedNbPatches;
        grownTexelCountdownVector[i] = maxTexelsForPatchi;
    }
    do
    {
        word64 w = 0x1;
        int i = 0;
        for(auto &seed : samplesValidated)
        {
            ASTex::RegionGrowingOperator< ASTex::NeighborhoodPeriodic > regionOp( patchMap );
            ImageMask64::PixelType seedTex;
            seedTex = w;
            itk::Index<2> seedInt;
            seedInt[0]=int(seed[0]);
            seedInt[1]=int(seed[1]);
            regionOp.newRegionFromSeed< ASTex::GrowingStrategyRandom >(
                        seedInt,

                        // Condition function for accepting a pixel to the current region.
                        [&]( const ASTex::NeighborhoodPeriodic::PixelType &texelId )
            {
                return (patchMap.pixelAbsolute(texelId) == 0x0 || patchMap.pixelAbsolute(texelId)==seedTex)
                        && grownTexelCountdownVector[i] > 0;
            },

            // Processing function applied to accepted pixels.
            [&]( ASTex::NeighborhoodPeriodic::PixelType &texelId, const ASTex::NeighborhoodPeriodic::NeighborsList& /*neighbors*/ )
            {
                ImageMask64::PixelType &pix = patchMap.pixelAbsolute(texelId);
                if(pix == 0x0)
                {
                    --grownTexelCountdownVector[i];
                    ++totalNumberOfMarkedPixels;
                    pix = seedTex;
                }
            });
            w*=2;
            grownTexelCountdownVector[i++]+=64; //whatever average number is fine, this is used to add more pixels in case the image isn't filled
        }
    }
    while(totalNumberOfMarkedPixels < this->m_texture.width()*this->m_texture.height());

    //Debug only
    ASTex::ImageGrayu8 patches;
    patches.initItk(this->m_texture.width(), this->m_texture.height());
    word64 w=0x1;
    for(unsigned i=0; i<Base::m_requestedNbPatches; ++i)
    {
        patchMap.for_all_pixels([&] (const word64 &pix, int x, int y)
        {
            if(pix == w)
                patches.pixelAbsolute(x, y)=uint8_t((i*255)/(Base::m_requestedNbPatches-1));
        });
        w *= 2;
    }
    ChooserPatches_FromImage<ASTex::ImageGrayu8> cpfi;
    cpfi.setPatchProcessor(this);
    cpfi.setPatchImage(patches);

    delete[] grownTexelCountdownVector;

    return cpfi.generate();
}

//ChooserPatches_HighContrast

template<typename I>
class PatchProcessor<I>::ChooserPatches_HighContrast : PatchProcessor<I>::ChooserPatches
{
public:
    using Base = PatchProcessor<I>::ChooserPatches;

    ChooserPatches_HighContrast() :
        Base(),
        m_fragmentMinSize(20),
        m_fragmentMaxSize(100),
        m_fragmentColorThreshold(40),
        m_downsamplingMinSize(128)
      {}

    //Parameters

    void setFragmentMinSize(unsigned fragmentMinSize) {m_fragmentMinSize = fragmentMinSize;}
    void setFragmentMaxSize(unsigned fragmentMaxSize) {m_fragmentMaxSize = fragmentMaxSize;}
    void setFragmentColorThreshold(unsigned fragmentColorThreshold) {m_fragmentColorThreshold = fragmentColorThreshold;}
    void setDownsamplingMinSize(unsigned downsamplingMinSize) {m_downsamplingMinSize = downsamplingMinSize;}

    //Generator

    void generate() override;

private:
    unsigned m_fragmentMinSize;
    unsigned m_fragmentMaxSize;
    unsigned m_fragmentColorThreshold;
    unsigned m_downsamplingMinSize;
};

template<typename I>
void PatchProcessor<I>::ChooserPatches_HighContrast::generate()
{
    assert(this &&
           "Base::generate: a patchProcessor must have been set (it can only be set by the PatchProcessor class)");
    srand(this->m_seed);

    ContentExchg::FragmentProcessor fProc( this->m_texture );
    fProc.createFragments( m_fragmentMaxSize, m_fragmentColorThreshold );
    fProc.cleanupSmallestFragments( m_fragmentMinSize );
    fProc.updateFragmentsAttributes();

    // Patches generation from old content exchange

    ContentExchg::PatchProcessor pProc( fProc );
    pProc.createPatches( int(Base::m_requestedNbPatches) );
    pProc.computePatchBoundaries();

    ASTex::ImageRGBu8::PixelType *idColors = new ASTex::ImageRGBu8::PixelType [ fProc.fragmentCount() ];
    for( int i=0; i<pProc.patchCount(); ++i )
    {
        idColors[i].SetRed  ( uint8_t((i * 255)/pProc.patchCount()-1) );
        idColors[i].SetGreen( uint8_t((i * 255)/pProc.patchCount()-1) );
        idColors[i].SetBlue ( uint8_t((i * 255)/pProc.patchCount()-1) );
    }

    ASTex::ImageRGBu8 *imgPatch = new ASTex::ImageRGBu8( this->m_texture.width(), this->m_texture.height() );
    for( auto &p : pProc.patches() )
        for( auto &frag : p.fragments )
            for( auto &pixel : fProc.fragmentById(int(frag)).pixels )
                imgPatch->pixelAbsolute(pixel) = idColors[p.id];

    //content generation from old content exchange

    std::vector<double> rotations;
    rotations.push_back( 0.0 );

    std::vector<double> scales;
    scales.push_back( 1.0000 );

    //translation to new content exchange
    PatchProcessor<I>::ChooserPatches_FromImage<ASTex::ImageRGBu8> cpfi(this);
    cpfi.setPatchImage(*imgPatch);
    return cpfi.generate();
}

//Overrides of ChooserContents

template<typename I>
class PatchProcessor<I>::ChooserContents_RandomPoisson : PatchProcessor<I>::ChooserContents
{
public:
    using Base = PatchProcessor<I>::ChooserContents;

    ChooserContents_RandomPoisson() :
        Base() {}

    //Parameters

    //

    //Generator

    void generate() override;
};

template<typename I>
void PatchProcessor<I>::ChooserContents_RandomPoisson::generate()
{
    assert(this &&
           "Base::generate: a patchProcessor must have been set (it can only be set by the PatchProcessor class)");
    assert(this->m_patches.size()>0 &&
           "PatchProcessor::contents_initRandom: a patch initializer must be called before being able to choose contents");
    srand(this->m_seed);
    unsigned i, j;
    int randomShiftX, randomShiftY;
    I shiftedTexture;
    shiftedTexture.initItk(this->m_texture.width(), this->m_texture.height());
    for(j=0; j<this->m_patches.size(); ++j)
    {
        Patch<I> &patch=this->m_patches[j];
        for(i=0; i<Base::m_nbContents-1; ++i)
        {
            randomShiftX = rand();
            randomShiftY = rand();
            shiftedTexture.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
            {
                pix = this->m_texture.pixelAbsolute((x + randomShiftX)%shiftedTexture.width(), (y + randomShiftY)%shiftedTexture.height());
            });
            Content<I> c(shiftedTexture, patch);
            patch.addContent(c);
        }
    }
}

//ChooserContents_PCTSEnhanced

template<typename I>
template<typename CC_Base>
class PatchProcessor<I>::ChooserContents_PCTSEnhanced : PatchProcessor<I>::ChooserContents
{
public:
    using Base = PatchProcessor<I>::ChooserContents;

    ChooserContents_PCTSEnhanced() :
        Base(), m_pctsArgFile() {}

    //Parameters

    void setPCTSArgFile(const std::string &pctsArgFile);
    void setChooserContentsBase(CC_Base *chooserContentsBase){m_ccb = chooserContentsBase;}

    //Generator

    void generate() override;

private:
    std::string m_pctsArgFile;
    CC_Base *m_ccb;
};

template<typename I>
template<typename CC_Base>
void PatchProcessor<I>::ChooserContents_PCTSEnhanced<CC_Base>::generate()
{
    assert(this &&
           "Base::generate: a patchProcessor must have been set (it can only be set by the PatchProcessor class)");
    assert(m_ccb && "ChooserContents_PCTSEnhanced: A base ChooserContents object must be set");
        m_ccb->generate();
    assert(this->m_patches.size()>0 && this->m_patches[0].nbContents()>1);
    size_t arraySize = sizeof(typename I::PixelType) / sizeof(typename I::DataType);
    typename I::DataType *pix_one = new typename I::DataType[arraySize];
    static typename I::PixelType pix_zero;
    ASTex::Pcts<I> pcts(this->m_texture);
    pcts.setWidth(this->m_texture.width());
    pcts.setHeight(this->m_texture.height());
    pcts.setPeriodicity(true);
    if(!m_pctsArgFile.empty())
        pcts.loadParametersFromFile(m_pctsArgFile);
    for(int i=0; i<this->m_patches.size(); ++i)
    {
        Patch<I> &p = this->m_patches[i];
        MipmapCEPatch::PixelPos origin = p.originAt(0, 0);
        I mask;
        mask.initItk(this->m_texture.width(), this->m_texture.height());
        mask.for_all_pixels([] (typename I::PixelType &pix)
        {
            pix = pix_zero;
        });
        for(int x=0; x<p.alphaMap().width(); ++x)
        {
            for(int y=0; y<p.alphaMap().height(); ++y)
            {
                if(p.alphaMap().pixelAbsolute(x, y)>0)
                {
                    mask.pixelAbsolute( (x+origin[0])%mask.width(),
                                        (y+origin[1])%mask.height())
                            = pix_one;
                }
            }
        }
        mask.save(std::string("/home/nlutz/ieee2019/_debug/mask_p") + std::to_string(i) + ".png");
        for(int j=1; j<p.nbContents(); ++j)
        {
            Content<I> &c = p.contentAt(j);
            I contentTexture;
            contentTexture.initItk(this->m_texture.width(), this->m_texture.height());
            contentTexture.copy_pixels(this->m_texture);

            for(int x=0; x<c.texture().width(); ++x)
            {
                for(int y=0; y<c.texture().height(); ++y)
                {
                    if(p.alphaMap().pixelAbsolute(x, y)>0)
                    {
                        contentTexture.pixelAbsolute(   (x+origin[0])%contentTexture.width(),
                                                        (y+origin[1])%contentTexture.height())
                                = c.texture().pixelAbsolute(x, y);
                    }
                }
            }
            contentTexture.save(std::string("/home/nlutz/ieee2019/_debug/inputWithContent_p") + std::to_string(i) + "_c" + std::to_string(j) + ".png");
            pcts.setMask(contentTexture, mask);
            contentTexture = pcts.generate();
            contentTexture.save(std::string("/home/nlutz/ieee2019/_debug/pctsd_p") + std::to_string(i) + "_c" + std::to_string(j) + ".png");
            Content<I> pctsdContent(contentTexture, p);
            c = pctsdContent;
        }
    }
    delete[](pix_one);
}

//Overrides of ChooserPatchesAndContents

//ChooserPatchesAndContents_OldMethod

//ChooserPatchesAndContents

//template<typename I>
//class PatchProcessor<I>::ChooserPatchesAndContents_OldMethod<typename std::enable_if<std::is_same<I, ASTex::ImageRGBu8>::value>::type*> :
//        public PatchProcessor<I>::ChooserPatchesAndContents
//{
//public:
//    using Base = PatchProcessor<I>::ChooserPatchesAndContents;

//    ChooserPatchesAndContents_OldMethod() :
//        ChooserPatchesAndContents<I>() {}

//    //Parameters

//    void setFragmentMinSize(unsigned fragmentMinSize) {m_fragmentMinSize = fragmentMinSize;}
//    void setFragmentMaxSize(unsigned fragmentMaxSize) {m_fragmentMaxSize = fragmentMaxSize;}
//    void setFragmentColorThreshold(unsigned fragmentColorThreshold) {m_fragmentColorThreshold = fragmentColorThreshold;}
//    void setDownsamplingMinSize(unsigned downsamplingMinSize) {m_downsamplingMinSize = downsamplingMinSize;}

//    //Generator

//    void generate() override;

//private:

//    unsigned m_fragmentMinSize;
//    unsigned m_fragmentMaxSize;
//    unsigned m_fragmentColorThreshold;
//    unsigned m_downsamplingMinSize;
//};

template<typename I>
void ChooserPatchesAndContents_OldMethod<typename std::enable_if<std::is_same<I, ASTex::ImageRGBu8>::value>::type*>::generate()
{
    srand(this->m_seed);

    ContentExchg::FragmentProcessor fProc( this->m_texture );
    fProc.createFragments( m_fragmentMaxSize, m_fragmentColorThreshold );
    fProc.cleanupSmallestFragments( m_fragmentMinSize );
    fProc.updateFragmentsAttributes();

    // Patches generation from old content exchange

    ContentExchg::PatchProcessor pProc( fProc );
    pProc.createPatches( int(Base::m_requestedNbPatches) );
    pProc.computePatchBoundaries();

    ASTex::ImageRGBu8::PixelType *idColors = new ASTex::ImageRGBu8::PixelType [ fProc.fragmentCount() ];
    for( int i=0; i<pProc.patchCount(); ++i )
    {
        idColors[i].SetRed  ( uint8_t((i * 255)/pProc.patchCount()-1) );
        idColors[i].SetGreen( uint8_t((i * 255)/pProc.patchCount()-1) );
        idColors[i].SetBlue ( uint8_t((i * 255)/pProc.patchCount()-1) );
    }

    ASTex::ImageRGBu8 *imgPatch = new ASTex::ImageRGBu8( this->m_texture.width(), this->m_texture.height() );
    for( auto &p : pProc.patches() )
        for( auto &frag : p.fragments )
            for( auto &pixel : fProc.fragmentById(int(frag)).pixels )
                imgPatch->pixelAbsolute(pixel) = idColors[p.id];

    //content generation from old content exchange

    std::vector<double> rotations;
    rotations.push_back( 0.0 );

    std::vector<double> scales;
    scales.push_back( 1.0000 );

    pProc.findAlternativeContents( rotations, scales, Base::m_nbContentsPerPatch-1, Base::m_downsamplingMinSize );

    //translation to new content exchange

    ChooserPatches_FromImage<I, ASTex::ImageRGBu8> cpfi(this);
    cpfi.setPatchImage(*imgPatch);
    cpfi.generate();

    this->contents_initDefault();
    int pId = 0;
    for(auto &p : pProc.patches())
    {
        for(auto &c : p.contents)
        {
            I shiftedTexture;
            shiftedTexture.initItk(this->m_texture.width(), this->m_texture.height());

            Patch<I> &patch=this->m_patches[pId];
            shiftedTexture.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
            {
                pix = this->m_texture.pixelAbsolute(  (x+c.offset[0])%shiftedTexture.width(),
                        (y+c.offset[1])%shiftedTexture.height());
            });
            Content<I> content(shiftedTexture, patch);
            patch.addContent(content);
            shiftedTexture.initItk(this->m_texture.width(), this->m_texture.height());
        }
        ++pId;
    }

    //Debug only
    for( int i=0; i<pProc.patchCount(); ++i )
    {
        idColors[i].SetRed  ( uint8_t((i * 255)/pProc.patchCount()-1) );
        idColors[i].SetGreen( uint8_t((i * 255)/pProc.patchCount()-1) );
        idColors[i].SetBlue ( uint8_t((i * 255)/pProc.patchCount()-1) );
    }
    for( auto &p : pProc.patches() )
        for( auto &frag : p.fragments )
            for( auto &pixel : fProc.fragmentById(int(frag)).pixels )
                imgPatch->pixelAbsolute(pixel) = idColors[p.id];
}

}

}

#endif
