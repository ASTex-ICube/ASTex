#ifndef __PCTS_H__
#define __PCTS_H__

#include "ASTex/mipmap.h"
#include "ASTex/easy_io.h"
#include "ASTex/Stamping/sampler.h"
#include "ASTex/utils.h"

#define PCTS_DEBUG_DIRECTORY "/home/nlutz/img/PCTS_debug/"

namespace ASTex
{

using ImageIndex2 = ASTex::ImageRGB32;

template<typename I>
class Pcts
{
public:

    typedef struct
    {
        unsigned    nbPasses = 2,
                    minimumSizeLog = 5,
                    lvl0BlocksSize = 4,
                    correctionNeighborhood = 2,
                    nbSamplesNNM = 60,
                    nbRefinementsNNM = 1,
                    radiusScaleNNM = 30;

        bool        periodicity = false,
                    useSynthesisForInit = false;
    } Parameters;

    Pcts();
    Pcts(const I& texture);

    void setTexture(const I& texture) {m_texture = texture; m_textureSet=true;}


	void setLabel(const I& label, double weight) { 
		assert(m_textureSet &&
			"Pcts::setLabel: a texture must be set (try using Pcts::setTexture()).");
        assert(m_texture.width()==label.width() && m_texture.height()==label.height() &&
			"Pcts::setLabel: the label map must match the texture.");
		m_label = label; m_labelSet = true; m_labelWeight = weight;
	}
	void setGuidance(const I& guid, const I& segmented, double weight, double strength) {
		assert(m_textureSet &&
			"Pcts::setGuidance: a texture must be set (try using Pcts::setTexture()).");
        assert(m_texture.width() == guid.width() && m_texture.height() == guid.height() &&
            "Pcts::setGuidance: the guidance map must match the texture's size.");
        assert(m_texture.width() == segmented.width() && m_texture.height() == segmented.height() &&
			"Pcts::setGuidance: the segmented map must match the texture.");
		m_guidance = guid; m_segmented = segmented;  m_guidanceSet = true;
		m_guidanceWeight = weight; m_strength = strength;
	}
	void setMask(const I& synthesis, const I& mask) {
        assert(m_textureSet &&
            "Pcts::setMask: a texture must be set (try using Pcts::setTexture()).");
        assert(m_texture.width() == synthesis.width() && m_texture.height() == synthesis.height() &&
            "Pcts::setMask: the synthesis must match the texture's size.");
        assert(m_texture.width() == mask.width() && m_texture.height() == mask.height() &&
            "Pcts::setMask: the mask must match the texture's size.");
		m_synthesis = synthesis; m_mask = mask;  m_maskSet = true;
	}
	void setStencil(const I& stencil, double w) {
		assert(m_textureSet &&
			"Pcts::setStencil: a texture must be set (try using Pcts::setTexture()).");
		m_stencil = stencil; m_stencilSet = true; m_stencilWeight = w;
	}
	void setWidth(int width) {m_width=width;}
    void setHeight(int height) {m_height=height;}
    void setSize(int width, int height) {m_width=width; m_height=height;}

    void setNbPasses(unsigned nbPasses) {m_nbPasses=nbPasses;}

    void setMinimumSizeLog(unsigned logValue) {m_minimumSizeLog=logValue;}

    void setLvl0BlockSize(unsigned blockSize) {m_lvl0BlockSize=blockSize;}

    void setCorrectionNeighborhood(unsigned neighborHood) {m_neighborhood=neighborHood;}

    void setNbSamplesNNM(unsigned nbSamples) {m_nbSamplesNNM=nbSamples;}
    void setRadiusScaleNNM(unsigned radiusScale) {m_radiusScaleNNM=radiusScale;}
    void setNbRefinementsNNM(unsigned nbRefinements) {m_nbRefinementsNNM=nbRefinements;}

    void setPeriodicity(bool periodicity) {m_periodicity = periodicity;}
    void setUseSynthesisForInit(bool b) {m_useSynthesisForInit = b;}

    void setParameters(const Parameters& p);
    bool loadParametersFromFile(const std::string &file);
    bool saveParametersIntoFile(const std::string &file);

    I generate();

private:

    I m_texture, m_label, m_guidance, m_segmented, m_synthesis, m_mask, m_stencil;
	ImageRGBf m_error;
    int m_width;
    int m_height;
	double m_labelWeight, m_stencilWeight;
	double m_guidanceWeight, m_strength;

    unsigned m_nbPasses;
    unsigned m_minimumSizeLog;
    unsigned m_lvl0BlockSize;

    unsigned m_neighborhood;

    unsigned m_nbSamplesNNM;
    unsigned m_radiusScaleNNM;
    unsigned m_nbRefinementsNNM;

    bool m_textureSet;
	bool m_labelSet;
	bool m_guidanceSet;
	bool m_maskSet;
	bool m_stencilSet;

    bool m_periodicity;
    bool m_useSynthesisForInit;
};

template<typename I>
Pcts<I>::Pcts() :
	m_texture(),
	m_width(800),
	m_height(800),
	m_nbPasses(2),
	m_minimumSizeLog(5),
	m_lvl0BlockSize(8),
	m_neighborhood(2),
	m_nbSamplesNNM(10),
	m_radiusScaleNNM(5),
	m_nbRefinementsNNM(2),
	m_textureSet(false),
	m_labelSet(false),
	m_labelWeight(0.0),
	m_guidanceSet(false),
	m_guidanceWeight(0.0),
	m_strength(1.0),
	m_maskSet(false),
    m_stencilSet(false),
    m_stencilWeight(0.0),
    m_periodicity(false),
    m_useSynthesisForInit(false)
{}

template<typename I>
Pcts<I>::Pcts(const I& texture) :
    m_texture(texture),
    m_width(800),
    m_height(800),
    m_nbPasses(2),
    m_minimumSizeLog(5),
    m_lvl0BlockSize(8),
    m_neighborhood(2),
    m_nbSamplesNNM(5),
    m_radiusScaleNNM(5),
    m_nbRefinementsNNM(2),
    m_textureSet(true),
	m_labelSet(false),
	m_labelWeight(0.0),
	m_guidanceSet(false),
	m_guidanceWeight(0.0),
	m_strength(1.0),
	m_maskSet(false),
    m_stencilSet(false),
    m_stencilWeight(0.0),
    m_periodicity(false),
    m_useSynthesisForInit(false)
{}

template<typename I>
I Pcts<I>::generate()
{
    assert(m_textureSet &&
           "Pcts::generate: a texture must be set (try using Pcts::setTexture()).");
    assert(!m_useSynthesisForInit || m_width == m_texture.width() || m_height == m_texture.height() ||
           m_width == m_synthesis.width() || m_height == m_synthesis.height());

    Mipmap<I> pyramidInput, pyramidLabel, pyramidGuidance, pyramidSegmented, pyramidSynthesis, pyramidMask, pyramidStencil;
    ImageIndex2 indexImageLevel0;

    //pyramid building

    unsigned maxReductionLevel = std::max(0.0, std::log2(std::min(m_texture.width(), m_texture.height())) - m_minimumSizeLog);
	std::cout << "levels:" << maxReductionLevel;

    pyramidInput.setTexture(m_texture);
    pyramidInput.setMode(ISOTROPIC);
    pyramidInput.setMaxPowReductionLevel(maxReductionLevel);
    pyramidInput.generate();
	if (m_labelSet)
	{
		pyramidLabel.setTexture(m_label);
		pyramidLabel.setMode(ISOTROPIC);
		pyramidLabel.setMaxPowReductionLevel(maxReductionLevel);
		pyramidLabel.generate();
	}
	if (m_guidanceSet)
	{
		pyramidGuidance.setTexture(m_guidance);
		pyramidGuidance.setMode(ISOTROPIC);
		pyramidGuidance.setMaxPowReductionLevel(maxReductionLevel);
		pyramidGuidance.generate();
		pyramidSegmented.setTexture(m_segmented);
		pyramidSegmented.setMode(ISOTROPIC);
		pyramidSegmented.setMaxPowReductionLevel(maxReductionLevel);
		pyramidSegmented.generate();
	}
	if (m_maskSet)
	{
		pyramidMask.setTexture(m_mask);
		pyramidMask.setMode(ISOTROPIC);
		pyramidMask.setMaxPowReductionLevel(maxReductionLevel);
		pyramidMask.generate();
		pyramidSynthesis.setTexture(m_synthesis);
		pyramidSynthesis.setMode(ISOTROPIC);
		pyramidSynthesis.setMaxPowReductionLevel(maxReductionLevel);
		pyramidSynthesis.generate();
	}
	if (m_stencilSet)
	{
		pyramidStencil.setTexture(m_stencil);
		pyramidStencil.setMode(ISOTROPIC);
		pyramidStencil.setMaxPowReductionLevel(maxReductionLevel);
		pyramidStencil.generate();
	}

    const I& lvl0mipmap = pyramidInput.mipmap(maxReductionLevel, maxReductionLevel); //mipmap of the input
    if (m_guidanceSet)
    {
        m_width = m_guidance.width();
        m_height = m_guidance.height();
    }

    int     widthBlockyImage = m_width/std::pow(2.0, maxReductionLevel), //width of the output initialization
            heightBlockyImage= m_height/std::pow(2.0, maxReductionLevel);
	int		ninit = 1;

    if (m_guidanceSet)
    {
		m_error.initItk(widthBlockyImage / m_lvl0BlockSize+1, heightBlockyImage / m_lvl0BlockSize+1);
		ninit = 100;
	}

    //output initialization
    indexImageLevel0.initItk(widthBlockyImage, heightBlockyImage);
    indexImageLevel0.parallel_for_all_pixels([&] (ImageIndex2::PixelType &pix, int x, int y)
    {
        pix[0]=x;
        pix[1]=y;
    });
    if(!m_useSynthesisForInit)
    {
        for (int kinit=0; kinit<ninit; kinit++) //1 without guidance
            for(int x=0; x<indexImageLevel0.width(); x+=m_lvl0BlockSize)
                for(int y=0; y<indexImageLevel0.height(); y+=m_lvl0BlockSize)
                {
                    int xT, yT;
                    //choice of random translation
                    if(m_periodicity)
                    {
                        xT=rand()%(std::max(1, int(lvl0mipmap.width())));
                        yT=rand()%(std::max(1, int(lvl0mipmap.height())));
                    }
                    else
                    {
                        xT=rand()%(std::max(1, int(lvl0mipmap.width()-m_lvl0BlockSize)));
                        yT=rand()%(std::max(1, int(lvl0mipmap.height()-m_lvl0BlockSize)));
                    }
                    bool do_init = true;
    //				if (m_guidanceSet)
    //				{
    //                    if (kinit == 0)
    //                    {
    //                        xT = x;
    //                        yT = y;
    //                    }
    //					float err = 0.0f;
    //					int count = 0;
    //					for (unsigned x2 = x; x2 < m_lvl0BlockSize + x && x2 < (unsigned)indexImageLevel0.width(); ++x2)
    //						for (unsigned y2 = y; y2 < m_lvl0BlockSize + y && y2 < (unsigned)indexImageLevel0.height(); ++y2) //TODO: preconditions
    //						{
    //							count++;
    //							for (int ii = 0; ii < 3; ++ii)
    //							{
    //								int indx = xT + x2 - x, indy = yT + y2 - y;
    //								err += ((pyramidSegmented.mipmap(maxReductionLevel, maxReductionLevel).pixelAbsolute(indx, indy)[ii]) - pyramidGuidance.mipmap(maxReductionLevel, maxReductionLevel).pixelAbsolute(x2, y2)[ii]) *
    //									((pyramidSegmented.mipmap(maxReductionLevel, maxReductionLevel).pixelAbsolute(indx, indy)[ii]) - pyramidGuidance.mipmap(maxReductionLevel, maxReductionLevel).pixelAbsolute(x2, y2)[ii]);
    //							}
    //						}
    //					err /= (float)count;
    //					std::cout << "Pass" << kinit << " x=" << x << " y=" << y << "error=" << err <<"\n";
    //					if (kinit == 0)
    //					{
    //                        m_error.pixelAbsolute(x / m_lvl0BlockSize, y / m_lvl0BlockSize)[0] = err;
    //					}
    //					else
    //					{
    //						float olderr = m_error.pixelAbsolute(x / m_lvl0BlockSize, y / m_lvl0BlockSize)[0];
    //                        if (err > olderr)
    //                            do_init = false;
    //						else m_error.pixelAbsolute(x / m_lvl0BlockSize, y / m_lvl0BlockSize)[0] = err;
    //					}
    //				}
                    if (do_init)
                        for(unsigned x2=x; x2<m_lvl0BlockSize+x && x2<(unsigned)indexImageLevel0.width(); ++x2)
                            for(unsigned y2=y; y2<m_lvl0BlockSize+y && y2<(unsigned)indexImageLevel0.height(); ++y2) //TODO: preconditions
                            {
                                if(!m_maskSet || !I::is_zero(pyramidMask.mipmap(maxReductionLevel,maxReductionLevel).pixelAbsolute(x2, y2)))
                                {
                                    //this formula takes periodicity enabled/disabled in account.
                                    indexImageLevel0.pixelAbsolute(x2, y2)[0] = (xT+x2-x + indexImageLevel0.width()) % indexImageLevel0.width();
                                    indexImageLevel0.pixelAbsolute(x2, y2)[1] = (yT+y2-y + indexImageLevel0.height()) % indexImageLevel0.height();
                                }
                                else
                                {
                                    //if mask is set, do not scramble the texture at x,y.
                                    indexImageLevel0.pixelAbsolute(x2, y2)[0] = x2;
                                    indexImageLevel0.pixelAbsolute(x2, y2)[1] = y2;
                                }
                            }
                }
    }

    I imageLevel0;
    I labelLevel0;
    I guidanceLevel0;
    I segmentedLevel0;

    /**
      * Lambda lmbd_lookupIndexIntoImage performs the lookup pyramid.mipmap(s) o indexImage and sets the result to image.
    **/
    auto lmbd_lookupIndexIntoImage = [&] (const ImageIndex2& indexImage, I& image, const Mipmap<I>& pyramid, int level)
    {
        image.initItk(indexImage.width(), indexImage.height());
        image.parallel_for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
        {
            const I& mipmap = pyramid.mipmap(level, level);
            const ImageIndex2::PixelType &pixIndirect = indexImage.pixelAbsolute(x, y);
            pix = mipmap.pixelAbsolute(pixIndirect[0], pixIndirect[1]);
        });
        return;
    };

    /**
      * Lambda lmbd_fakeColorIndex turns image into a vizualizable index map based on indexImage,
      * of colors (0,0)=red, (1,0)=black, (0,1)=yellow, (1,1)=green.
    **/
    auto lmbd_fakeColorIndex = [&] (const ImageIndex2& indexImage, ImageRGBd& image)
    {
        image.initItk(indexImage.width(), indexImage.height());
        ImageRGBd::PixelType color0, color1, color2;
        color0[0]=1; color0[1]=0; color0[2]=0;
        color1[0]=0; color1[1]=1; color1[2]=0;
        color2[0]=1; color2[1]=1; color2[2]=0;
        double s, t;
        image.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
        {
            s=indexImage.pixelAbsolute(x, y)[0]/double(image.width());
            t=indexImage.pixelAbsolute(x, y)[1]/double(image.height());
            pix = color0*(1.0-s)*(1.0-t) + color1*s*t + color2*(1-s)*t;
        });
        return;
    };

    if(m_useSynthesisForInit)
        imageLevel0 = pyramidSynthesis.mipmap(maxReductionLevel, maxReductionLevel);
    else
        lmbd_lookupIndexIntoImage(indexImageLevel0, imageLevel0, pyramidInput, maxReductionLevel);
    if (m_maskSet)
        imageLevel0.parallel_for_all_pixels([&](typename I::PixelType &pix, int x, int y)
        {
            if (I::is_zero(pyramidMask.mipmap(maxReductionLevel, maxReductionLevel).pixelAbsolute(x, y)))
                pix = pyramidSynthesis.mipmap(maxReductionLevel, maxReductionLevel).pixelAbsolute(x, y);
        });
	if (m_labelSet) lmbd_lookupIndexIntoImage(indexImageLevel0, labelLevel0, pyramidLabel, maxReductionLevel);
	if (m_guidanceSet) {
		lmbd_lookupIndexIntoImage(indexImageLevel0, guidanceLevel0, pyramidGuidance, maxReductionLevel);
		lmbd_lookupIndexIntoImage(indexImageLevel0, segmentedLevel0, pyramidSegmented, maxReductionLevel);
	}
    imageLevel0.save(std::string(PCTS_DEBUG_DIRECTORY) + "testBlock" + std::to_string(maxReductionLevel) + ".png");

    //Synthesis
    for(int s=maxReductionLevel; s>=0; --s)
    {
        float coeff = pow((float)s / (float)maxReductionLevel, m_strength);
        float gweight = m_guidanceWeight*coeff;
        //Correction
        for (int npass = 0; npass < m_nbPasses; npass++)
        {
			std::cout << "level:" << s << " pass:" << npass << "\n";
            if (m_guidanceSet)
                std::cout << "guidance weight:" << gweight << "\n";
            if(s!=maxReductionLevel) //(already done, and can break the useSynthesisInit option)
                lmbd_lookupIndexIntoImage(indexImageLevel0, imageLevel0, pyramidInput, s);
            if (m_maskSet) imageLevel0.parallel_for_all_pixels([&](typename I::PixelType &pix, int x, int y)
			{
                if (I::is_zero(pyramidMask.mipmap(s, s).pixelAbsolute(x, y)))
                    pix = pyramidSynthesis.mipmap(s, s).pixelAbsolute(x, y);
			});
            if (m_labelSet)
                lmbd_lookupIndexIntoImage(indexImageLevel0, labelLevel0, pyramidLabel, s);
            if (m_guidanceSet)
            {
				lmbd_lookupIndexIntoImage(indexImageLevel0, guidanceLevel0, pyramidGuidance, s);
				lmbd_lookupIndexIntoImage(indexImageLevel0, segmentedLevel0, pyramidSegmented, s);
			}
			for (unsigned nx = 0; nx < 2 * m_neighborhood; ++nx)
				for (unsigned ny = 0; ny < 2 * m_neighborhood; ++ny)
					for (unsigned x = nx; x<unsigned(indexImageLevel0.width()); x += 2 * m_neighborhood)
						for (unsigned y = ny; y<unsigned(indexImageLevel0.height()); y += 2 * m_neighborhood)
						{
                            if(m_maskSet && I::is_zero(pyramidMask.mipmap(s, s).pixelAbsolute(x, y)))
                            {
                                imageLevel0.pixelAbsolute(x, y) = pyramidSynthesis.mipmap(s, s).pixelAbsolute(x, y);
                                indexImageLevel0.pixelAbsolute(x, y)[0] = x;
                                indexImageLevel0.pixelAbsolute(x, y)[1] = y;
                            }
                            else
                            {
                                itk::Index<2> bestIdErrMin;
                                double bestErrMin;

                                bestIdErrMin[0] = indexImageLevel0.pixelAbsolute(x, y)[0];
                                bestIdErrMin[1] = indexImageLevel0.pixelAbsolute(x, y)[1];

                                itk::Index<2> idErrMin;
                                double errMin;

                                auto lmbd_updateError = [&] (itk::Index<2> checkedPos, double &checkedError)
                                {
                                    checkedError = mse(pyramidInput.mipmap(s, s), imageLevel0, checkedPos[0], checkedPos[1], x, y, m_neighborhood, m_periodicity);
                                    if (m_labelSet)
                                        checkedError = (1.0- m_labelWeight)*checkedError
                                                + m_labelWeight * mse(pyramidLabel.mipmap(s, s), labelLevel0, checkedPos[0], checkedPos[1], x, y, m_neighborhood, m_periodicity);
                                    if (m_guidanceSet)
                                        checkedError = (1.0 - gweight)*checkedError
                                                + gweight * mse(pyramidGuidance.mipmap(s, s), guidanceLevel0, checkedPos[0], checkedPos[1], x, y, m_neighborhood, m_periodicity);
                                    if (m_stencilSet)
                                        checkedError = checkedError*(1.0 + m_stencilWeight*100.0* (1.0 - pyramidStencil.mipmap(s, s).pixelAbsolute(checkedPos[0], checkedPos[1])[0]));
                                };

                                lmbd_updateError(bestIdErrMin, bestErrMin);

                                /**
                                  * Lambda nearestNeighborMatch (NMM) searches for the nnm in x and y.
                                  * dx and dy are used to easily expand the borders of each block.
                                **/
                                auto nearestNeighborMatch = [&](int dx, int dy)
                                {
                                    int xdx = m_periodicity ? (x+dx)%indexImageLevel0.width() : x+dx;
                                    int ydy = m_periodicity ? (y+dy)%indexImageLevel0.height() : y+dy;

                                    int px = indexImageLevel0.pixelAbsolute(xdx, ydy)[0];
                                    int py = indexImageLevel0.pixelAbsolute(xdx, ydy)[1];
                                    double err2 = mse(pyramidInput.mipmap(s, s), imageLevel0, px, py, x, y, m_neighborhood, m_periodicity);
                                    int radius;
                                    //std::cout << "init:=" << minerr << ",ind=" << minpx << "," << minpy << "\n";
                                    itk::Index<2> pxpy;
                                    pxpy[0] = px;
                                    pxpy[1] = py;
                                    errMin = err2;
                                    idErrMin = pxpy;

                                    for (unsigned n = 0; n < m_nbRefinementsNNM &&
                                        (radius = m_radiusScaleNNM*m_neighborhood / (int)pow(2.0, (double)n)) >= 1;
                                        ++n)
                                    {
                                        Stamping::SamplerPoisson sp;
                                        sp.setNbPoints(m_nbSamplesNNM);
                                        std::vector<Eigen::Vector2f> spResult = sp.generate();

                                        itk::Index<2> idPoisson;
                                        for (unsigned i = 0; i < m_nbSamplesNNM; ++i)
                                        {
                                            //std::cout << "selected:" << spResult[i][0] << "," << spResult[i][1] << "\n";
                                            idPoisson[0] = (int)((2.0*spResult[i][0] - 1.0)*(float)radius) + px;
                                            idPoisson[1] = (int)((2.0*spResult[i][1] - 1.0)*(float)radius) + py;
                                            int pyramidInputWidth =  pyramidInput.mipmap(s, s).width();
                                            int pyramidInputHeight = pyramidInput.mipmap(s, s).height();
                                            if(m_periodicity)
                                            {
                                                idPoisson[0] = (idPoisson[0] + pyramidInputWidth)%pyramidInputWidth;
                                                idPoisson[1] = (idPoisson[1] + pyramidInputHeight)%pyramidInputHeight;
                                            }
                                            if (idPoisson[0] >= 0 && idPoisson[0] < pyramidInputWidth &&
                                                idPoisson[1] >= 0 && idPoisson[1] < pyramidInputHeight)
                                            {
                                                lmbd_updateError(idPoisson, err2);
                                                if (err2 < errMin)
                                                {
                                                    errMin = err2;
                                                    idErrMin = idPoisson;
                                                }
                                                //else std::cout << "check(" << i << "," << n << "):= " << err2 << ",ind=" << idPoisson[0] << "," << idPoisson[1] << "not better\n";
                                            }
                                        }
                                        px = idErrMin[0];
                                        py = idErrMin[1];
                                    }
                                };

                                nearestNeighborMatch(0, 0);
                                if (errMin < bestErrMin)
                                {
                                    bestErrMin = errMin;
                                    bestIdErrMin = idErrMin;
                                }
                                if (m_periodicity || y + 1  < (unsigned)indexImageLevel0.height())
                                {
                                    nearestNeighborMatch(0, 1);
                                    if (errMin < bestErrMin)
                                    {
                                        bestErrMin = errMin;
                                        bestIdErrMin = idErrMin;
                                    }
                                }
                                if (m_periodicity || x + 1 < (unsigned)indexImageLevel0.width())
                                {
                                    nearestNeighborMatch(1, 0);
                                    if (errMin < bestErrMin)
                                    {
                                        bestErrMin = errMin;
                                        bestIdErrMin = idErrMin;
                                    }
                                }
                                if (m_periodicity || (x + 1 < (unsigned)indexImageLevel0.width() &&
                                                      y + 1 < (unsigned)indexImageLevel0.height()) )
                                {
                                    nearestNeighborMatch(1, 1);
                                    if (errMin < bestErrMin)
                                    {
                                        bestErrMin = errMin;
                                        bestIdErrMin = idErrMin;
                                    }
                                }

                                indexImageLevel0.pixelAbsolute(x, y)[0] = bestIdErrMin[0];
                                indexImageLevel0.pixelAbsolute(x, y)[1] = bestIdErrMin[1];
                                imageLevel0.pixelAbsolute(x, y) = pyramidInput.mipmap(s, s).pixelAbsolute(bestIdErrMin);
                                if(m_stencilSet && pyramidStencil.mipmap(s, s).pixelAbsolute(bestIdErrMin)[0]==0)
                                    imageLevel0.pixelAbsolute(x, y) = pyramidSynthesis.mipmap(s, s).pixelAbsolute(x, y);

                                //std::cout << "FINAL:=" << besterrMin << "ind=" << bestidErrMin[0] << "," << bestidErrMin[1] << "\n";
                            }
						}
		}
        ImageIndex2 indexImageLevel1;
        imageLevel0.save(std::string(PCTS_DEBUG_DIRECTORY) + "imageLevel0_" + std::to_string(s) + ".png");

        ImageRGBd fakeColorIndex0;
        lmbd_fakeColorIndex(indexImageLevel0, fakeColorIndex0);
        IO::save01_in_u8(fakeColorIndex0, std::string(PCTS_DEBUG_DIRECTORY) + "fakeColorIndex0_" + std::to_string(s) + ".png");

        if(s>0)
        {
            indexImageLevel1.initItk(indexImageLevel0.width()*2, indexImageLevel0.height()*2);
            indexImageLevel0.parallel_for_all_pixels([&] (const typename ImageIndex2::PixelType &pix, int x, int y)
            {
                indexImageLevel1.pixelAbsolute(2*x, 2*y) = pix*2;

                indexImageLevel1.pixelAbsolute(2*x+1, 2*y) = pix*2;
                indexImageLevel1.pixelAbsolute(2*x+1, 2*y)[0]++;

                indexImageLevel1.pixelAbsolute(2*x, 2*y+1) = pix*2;
                indexImageLevel1.pixelAbsolute(2*x, 2*y+1)[1]++;

                indexImageLevel1.pixelAbsolute(2*x+1, 2*y+1) = pix*2;
                indexImageLevel1.pixelAbsolute(2*x+1, 2*y+1)[0]++;
                indexImageLevel1.pixelAbsolute(2*x+1, 2*y+1)[1]++;
            });
        }
        else
            indexImageLevel1 = indexImageLevel0;
        I imageLevel1;
        imageLevel1.initItk(indexImageLevel1.width(), indexImageLevel1.height());
        lmbd_lookupIndexIntoImage(indexImageLevel1, imageLevel1, pyramidInput, std::max(0, s-1));
        imageLevel1.save(std::string(PCTS_DEBUG_DIRECTORY) + "imageLevel1_" + std::to_string(s) + ".png");

        ImageRGBd fakeColorIndexImageLevel0;
        lmbd_fakeColorIndex(indexImageLevel0, fakeColorIndexImageLevel0);

        indexImageLevel0 = indexImageLevel1;

        IO::save01_in_u8(fakeColorIndexImageLevel0, std::string(PCTS_DEBUG_DIRECTORY) + "fakeColorIndex1_" + std::to_string(s) + ".png");
    }

    return imageLevel0;
}

template<typename I>
void Pcts<I>::setParameters(const Parameters& p)
{
    setNbPasses(p.nbPasses);
    setMinimumSizeLog(p.minimumSizeLog);
    setLvl0BlockSize(p.lvl0BlocksSize);
    setCorrectionNeighborhood(p.correctionNeighborhood);
    setNbSamplesNNM(p.nbSamplesNNM);
    setNbRefinementsNNM(p.nbRefinementsNNM);
    setRadiusScaleNNM(p.radiusScaleNNM);
    setPeriodicity(p.periodicity);
    setUseSynthesisForInit(p.useSynthesisForInit);
}

template<typename I>
bool Pcts<I>::loadParametersFromFile(const std::string &file)
{
    Parameters p;
    bool o;

    std::ifstream ifs(file);
    if((o=ifs.is_open()))
    {
        ifs >> p.nbPasses;
        ifs >> p.minimumSizeLog;
        ifs >> p.lvl0BlocksSize;
        ifs >> p.correctionNeighborhood;
        ifs >> p.nbSamplesNNM;
        ifs >> p.nbRefinementsNNM;
        ifs >> p.radiusScaleNNM;
        ifs >> p.periodicity;
        ifs >> p.useSynthesisForInit;
        ifs.close();
        setParameters(p);
    }
    return o;
}

template<typename I>
bool Pcts<I>::saveParametersIntoFile(const std::string &file)
{
    std::ofstream ofs(file);
    bool o;

    if((o=ofs.is_open()))
    {
        Parameters p;
        p.nbPasses = m_nbPasses;
        p.minimumSizeLog = m_minimumSizeLog;
        p.lvl0BlocksSize = m_lvl0BlockSize;
        p.correctionNeighborhood = m_neighborhood;
        p.nbSamplesNNM = m_nbSamplesNNM;
        p.nbRefinementsNNM = m_nbRefinementsNNM;
        p.radiusScaleNNM = m_radiusScaleNNM;
        p.periodicity = m_periodicity;

        ofs << p.nbPasses << std::endl;
        ofs << p.minimumSizeLog << std::endl;
        ofs << p.lvl0BlocksSize << std::endl;
        ofs << p.correctionNeighborhood << std::endl;
        ofs << p.nbSamplesNNM << std::endl;
        ofs << p.nbRefinementsNNM << std::endl;
        ofs << p.radiusScaleNNM << std::endl;
        ofs << p.periodicity << std::endl;
        ofs << p.useSynthesisForInit << std::endl;
        ofs.close();
    }
    return o;
}

}

#endif
