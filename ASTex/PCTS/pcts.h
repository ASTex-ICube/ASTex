#ifndef __PCTS_H__
#define __PCTS_H__

#include "ASTex/mipmap.h"
#include "ASTex/easy_io.h"
#include "ASTex/Stamping/sampler.h"
#include "ASTex/utils.h"

#define PCTS_DEBUG_DIRECTORY "E:/developpement/AsTex/AsTex/Data/Debug/"

namespace ASTex
{

using ImageIndex2 = ASTex::ImageRGB32;

template<typename I>
class Pcts
{
public:

    Pcts();
    Pcts(const I& texture);

    void setTexture(const I& texture) {m_texture = texture; m_textureSet=true;}
	void setLabel(const I& label, double weight) { 
		assert(m_textureSet &&
			"Pcts::setLabel: a texture must be set (try using Pcts::setTexture()).");
		assert(m_texture.width()==label.width && m_texture.height()==label.height() &&
			"Pcts::setLabel: the label map must match the texture.");
		m_label = label; m_labelSet = true; m_labelWeight = weight;
	}
	void setGuidance(const I& guid, const I& segmented, double weight, double strength) {
		assert(m_textureSet &&
			"Pcts::setLabel: a texture must be set (try using Pcts::setTexture()).");
		assert(m_texture.width() == guid.width && m_texture.height() == guid.height() &&
			"Pcts::setGuidance: the guidance map must match the texture.");
		assert(m_texture.width() == segmented.width && m_texture.height() == segmented.height() &&
			"Pcts::setGuidance: the segmented map must match the texture.");
		m_guidance = guid; m_segmented = segmented;  m_guidanceSet = true;
		m_guidanceWeight = weight; m_strength = strength;
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

    I generate();

private:

    I m_texture, m_label, m_guidance, m_segmented;
	ImageRGBf m_error;
    int m_width;
    int m_height;
	double m_labelWeight;
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
	m_strength(1.0)
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
	m_strength(1.0)
{}

template<typename I>
I Pcts<I>::generate()
{
    assert(m_textureSet &&
           "Pcts::generate: a texture must be set (try using Pcts::setTexture()).");

    Mipmap<I> pyramidInput, pyramidLabel, pyramidGuidance, pyramidSegmented;
    ImageIndex2 indexImageLevel0;

    //pyramid building

    unsigned maxReductionLevel = std::log2(std::min(m_texture.width(), m_texture.height())) - m_minimumSizeLog;
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

    const I& lvl0mipmap = pyramidInput.mipmap(maxReductionLevel, maxReductionLevel);
	for(unsigned i=0; i<=maxReductionLevel; ++i)
        IO::save01_in_u8(pyramidGuidance.mipmap(i, i), std::string(PCTS_DEBUG_DIRECTORY) + "guid_pyramid" + std::to_string(i) + ".png");
	for (unsigned i = 0; i <= maxReductionLevel; ++i)
		IO::save01_in_u8(pyramidSegmented.mipmap(i, i), std::string(PCTS_DEBUG_DIRECTORY) + "seg_pyramid" + std::to_string(i) + ".png");

	if (m_guidanceSet) { m_width = m_guidance.width(); m_height = m_guidance.height();  }

    int     widthBlockyImage = m_width/std::pow(2.0, maxReductionLevel),
            heightBlockyImage= m_height/std::pow(2.0, maxReductionLevel);
	int		ninit = 1;

	if (m_guidanceSet) {
		m_error.initItk(widthBlockyImage / m_lvl0BlockSize+1, heightBlockyImage / m_lvl0BlockSize+1);
		ninit = 100;
	}

	indexImageLevel0.initItk(widthBlockyImage, heightBlockyImage);
	for (int kinit=0; kinit<ninit; kinit++)
		for(int x=0; x<indexImageLevel0.width(); x+=m_lvl0BlockSize)
			for(int y=0; y<indexImageLevel0.height(); y+=m_lvl0BlockSize)
			{
				int xT=(rand())%(std::max(1, int(lvl0mipmap.width()-m_lvl0BlockSize)));
				int yT=(rand())%(std::max(1, int(lvl0mipmap.height()-m_lvl0BlockSize)));
				bool do_init = true;
				if (m_guidanceSet)
				{
					if (kinit == 0) { xT = x; yT = y; }
					float err = 0.0f;
					int count = 0;
					for (unsigned x2 = x; x2 < m_lvl0BlockSize + x && x2 < (unsigned)indexImageLevel0.width(); ++x2)
						for (unsigned y2 = y; y2 < m_lvl0BlockSize + y && y2 < (unsigned)indexImageLevel0.height(); ++y2) //TODO: preconditions
						{
							count++;
							for (int ii = 0; ii < 3; ++ii)
							{
								int indx = xT + x2 - x, indy = yT + y2 - y;
								err += ((pyramidSegmented.mipmap(maxReductionLevel, maxReductionLevel).pixelAbsolute(indx, indy)[ii]) - pyramidGuidance.mipmap(maxReductionLevel, maxReductionLevel).pixelAbsolute(x2, y2)[ii]) *
									((pyramidSegmented.mipmap(maxReductionLevel, maxReductionLevel).pixelAbsolute(indx, indy)[ii]) - pyramidGuidance.mipmap(maxReductionLevel, maxReductionLevel).pixelAbsolute(x2, y2)[ii]);
							}
						}
					err /= (float)count;
					std::cout << "Pass" << kinit << " x=" << x << " y=" << y << "error=" << err <<"\n";
					if (kinit == 0)
					{
						m_error.pixelAbsolute(x / m_lvl0BlockSize, y / m_lvl0BlockSize)[0] = err;
					}
					else
					{
						float olderr = m_error.pixelAbsolute(x / m_lvl0BlockSize, y / m_lvl0BlockSize)[0];
						if (err > olderr) do_init = false;
						else m_error.pixelAbsolute(x / m_lvl0BlockSize, y / m_lvl0BlockSize)[0] = err;
					}
				}
				if (do_init) std::cout << "update value\n";
				if (do_init) for(unsigned x2=x; x2<m_lvl0BlockSize+x && x2<(unsigned)indexImageLevel0.width(); ++x2)
					for(unsigned y2=y; y2<m_lvl0BlockSize+y && y2<(unsigned)indexImageLevel0.height(); ++y2) //TODO: preconditions
					{
						indexImageLevel0.pixelAbsolute(x2, y2)[0] = xT+x2-x;
						indexImageLevel0.pixelAbsolute(x2, y2)[1] = yT+y2-y;
					}
			}
    I imageLevel0, labelLevel0, guidanceLevel0, segmentedLevel0;

    /**
      * Lambda lmbd_lookupIndexIntoImage performs the lookup pyramid.mipmap(s) o indexImage and copies the result to image.
    **/
    auto lmbd_lookupIndexIntoImage = [&] (const ImageIndex2& indexImage, I& image, const Mipmap<I>& pyramid, int s)
    {
        image.initItk(indexImage.width(), indexImage.height());
        image.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
        {
            const I& mipmap = pyramid.mipmap(s, s);
            const ImageIndex2::PixelType &pixIndirect = indexImage.pixelAbsolute(x, y);
            pix = mipmap.pixelAbsolute(pixIndirect[0], pixIndirect[1]);
        });
        return;
    };

    /**
      * Lambda lmbd_fakeColorIndex turns image into a vizualizable index map,
      * of colors (0,0)=red, (1,0)=black, (0,1)=yellow, (1,1)=green.
    **/
    auto lmbd_fakeColorIndex = [&] (const ImageIndex2& indexImage, ImageRGBd& image)
    {
        image.initItk(indexImage.width(), indexImage.height());
        ImageRGBd::PixelType color0, color1, color2;
        color0[0]=1; color0[1]=0; color0[2]=0;
        color1[0]=0; color1[1]=1; color1[2]=0;
        color2[0]=1; color2[1]=1; color2[2]=0;
        float s, t;
        image.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
        {
            s=indexImage.pixelAbsolute(x, y)[0]/float(image.width());
            t=indexImage.pixelAbsolute(x, y)[1]/float(image.height());
            pix = color0*(1-s)*(1-t) + color1*s*t + color2*(1-s)*t;
        });
        return;
    };

    lmbd_lookupIndexIntoImage(indexImageLevel0, imageLevel0, pyramidInput, maxReductionLevel);
	if (m_labelSet) lmbd_lookupIndexIntoImage(indexImageLevel0, labelLevel0, pyramidLabel, maxReductionLevel);
	if (m_guidanceSet) {
		lmbd_lookupIndexIntoImage(indexImageLevel0, guidanceLevel0, pyramidGuidance, maxReductionLevel);
		lmbd_lookupIndexIntoImage(indexImageLevel0, segmentedLevel0, pyramidSegmented, maxReductionLevel);
	}
	IO::save01_in_u8(imageLevel0, std::string(PCTS_DEBUG_DIRECTORY) + "testBlock" + std::to_string(maxReductionLevel) + ".png");

    //Synthesis
    for(int s=maxReductionLevel; s>=0; --s)
    {
		float coeff = pow((float)s / (float)maxReductionLevel, m_strength);
        float gweight = m_guidanceWeight*coeff;
        //Correction
		for (int npass = 0; npass < m_nbPasses; npass++)
		{
			std::cout << "level:" << s << " pass:" << npass << "\n";
			if (m_guidanceSet) std::cout << "guidance weight:" << gweight << "\n";
			lmbd_lookupIndexIntoImage(indexImageLevel0, imageLevel0, pyramidInput, s);
			if (m_labelSet) lmbd_lookupIndexIntoImage(indexImageLevel0, labelLevel0, pyramidLabel, s);
			if (m_guidanceSet) {
				lmbd_lookupIndexIntoImage(indexImageLevel0, guidanceLevel0, pyramidGuidance, s);
				lmbd_lookupIndexIntoImage(indexImageLevel0, segmentedLevel0, pyramidSegmented, s);
			}
			for (unsigned nx = 0; nx < 2 * m_neighborhood; ++nx)
				for (unsigned ny = 0; ny < 2 * m_neighborhood; ++ny)
					for (unsigned x = nx; x<unsigned(indexImageLevel0.width()); x += 2 * m_neighborhood)
						for (unsigned y = ny; y<unsigned(indexImageLevel0.height()); y += 2 * m_neighborhood)
						{

							itk::Index<2> bestidErrMin;
							double besterrMin;
							bestidErrMin[0] = indexImageLevel0.pixelAbsolute(x, y)[0];
							bestidErrMin[1] = indexImageLevel0.pixelAbsolute(x, y)[1];

							besterrMin = mse(pyramidInput.mipmap(s, s), imageLevel0, bestidErrMin[0], bestidErrMin[1], x, y, m_neighborhood);
							if (m_labelSet) besterrMin = (1.0- m_labelWeight)*besterrMin
								+ m_labelWeight * mse(pyramidLabel.mipmap(s, s), labelLevel0, bestidErrMin[0], bestidErrMin[1], x, y, m_neighborhood);
							if (m_guidanceSet) besterrMin = (1.0 - gweight)*besterrMin
								+ gweight * mse(pyramidGuidance.mipmap(s, s), guidanceLevel0, bestidErrMin[0], bestidErrMin[1], x, y, m_neighborhood);

							//std::cout << "\nAT pixel:" << x << "," << y << "=" << besterrMin << ",ind="<< bestidErrMin[0]<<","<< bestidErrMin[1]<< "\n";

							itk::Index<2> idErrMin;
							double errMin;

                            /**
							  * Lambda nearestNeighborMatch (NMM) searches for the nnm in x and y.
							  * dx and dy are used to easily expand the borders of each block.
							**/
							auto nearestNeighborMatch = [&](int dx, int dy)
							{
								int px = indexImageLevel0.pixelAbsolute(x + dx, y + dy)[0];
								int py = indexImageLevel0.pixelAbsolute(x + dx, y + dy)[1];
								double err2 = mse(pyramidInput.mipmap(s, s), imageLevel0, px, py, x, y, m_neighborhood);
								int radius;
								int minpx = px, minpy = py;
								double minerr = err2;
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
									Stamping::SamplerUniform sp;
									sp.setNbPoints(m_nbSamplesNNM);
									std::vector<Eigen::Vector2f> spResult = sp.generate();

									itk::Index<2> idPoisson;
									for (unsigned i = 0; i < m_nbSamplesNNM; ++i)
									{
										//std::cout << "selected:" << spResult[i][0] << "," << spResult[i][1] << "\n";
										idPoisson[0] = (int)((2.0*spResult[i][0] - 1.0)*(float)radius) + px;
										idPoisson[1] = (int)((2.0*spResult[i][1] - 1.0)*(float)radius) + py;
										if (idPoisson[0] >= 0 && idPoisson[0] < pyramidInput.mipmap(s, s).width() &&
											idPoisson[1] >= 0 && idPoisson[1] < pyramidInput.mipmap(s, s).height())
										{
											double err2 = mse(pyramidInput.mipmap(s, s), imageLevel0, idPoisson[0], idPoisson[1], x, y, m_neighborhood);
											if (m_labelSet) err2 = (1.0 - m_labelWeight)*err2
												+ m_labelWeight * mse(pyramidLabel.mipmap(s, s), labelLevel0, idPoisson[0], idPoisson[1], x, y, m_neighborhood);
											if (m_guidanceSet) err2 = (1.0 - gweight)*err2
												+ gweight * mse(pyramidGuidance.mipmap(s, s), guidanceLevel0, idPoisson[0], idPoisson[1], x, y, m_neighborhood);
											if (err2 < errMin)
											{
												errMin = err2; idErrMin = idPoisson;
											}
											if (err2 < minerr) {
												//std::cout << "check(" << i << "," << n << "):= " << err2 << ",rad="<< radius <<",ind=" << idPoisson[0] << "," << idPoisson[1] << "is better\n";
												minerr = err2; minpx = idPoisson[0]; minpy = idPoisson[1];
											}
											//else std::cout << "check(" << i << "," << n << "):= " << err2 << ",ind=" << idPoisson[0] << "," << idPoisson[1] << "not better\n";
										}
									}
									//std::cout << "best is (" << n << "):= " << minerr << " ind=" << minpx << "," << minpy << "\n";
									px = idErrMin[0]; py = idErrMin[1];
								}
								//std::cout << "last is := " << errMin << " ind=" << idErrMin[0] << "," << idErrMin[1] << "\n";
							};

							nearestNeighborMatch(0, 0);
							if (errMin < besterrMin) {
								besterrMin = errMin; bestidErrMin = idErrMin;
							}
							if (y + 1 < (unsigned)indexImageLevel0.height())
							{
								nearestNeighborMatch(0, 1);
								if (errMin < besterrMin) {
									besterrMin = errMin; bestidErrMin = idErrMin;
								}
							}
							if (x + 1 < (unsigned)indexImageLevel0.width())
							{
								nearestNeighborMatch(1, 0);
								if (errMin < besterrMin) {
									besterrMin = errMin; bestidErrMin = idErrMin;
								}
							}
							if (x + 1 < (unsigned)indexImageLevel0.width() && y + 1 < (unsigned)indexImageLevel0.height())
							{
								nearestNeighborMatch(1, 1);
								if (errMin < besterrMin) {
									besterrMin = errMin; bestidErrMin = idErrMin;
								}
							}

							indexImageLevel0.pixelAbsolute(x, y)[0] = bestidErrMin[0];
							indexImageLevel0.pixelAbsolute(x, y)[1] = bestidErrMin[1];
							imageLevel0.pixelAbsolute(x, y) = pyramidInput.mipmap(s, s).pixelAbsolute(bestidErrMin);
							//std::cout << "FINAL:=" << besterrMin << "ind=" << bestidErrMin[0] << "," << bestidErrMin[1] << "\n";
						}
		}
        if(s>0)
        {
            ImageIndex2 indexImageLevel1;
            indexImageLevel1.initItk(indexImageLevel0.width()*2, indexImageLevel0.height()*2);
            indexImageLevel0.for_all_pixels([&] (const typename I::PixelType &pix, int x, int y)
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
            I imageLevel1;
            imageLevel1.initItk(indexImageLevel1.width(), indexImageLevel1.height());
            lmbd_lookupIndexIntoImage(indexImageLevel1, imageLevel1, pyramidInput, s-1);
            IO::save01_in_u8(imageLevel1, std::string(PCTS_DEBUG_DIRECTORY) + "testBlock" + std::to_string(s) + ".png");

            ImageRGBd fakeColorIndexImageLevel0;
            lmbd_fakeColorIndex(indexImageLevel0, fakeColorIndexImageLevel0);

            indexImageLevel0 = indexImageLevel1;

            IO::save01_in_u8(fakeColorIndexImageLevel0, std::string(PCTS_DEBUG_DIRECTORY) + "fakeColorIndex" + std::to_string(s) + ".png");
        }

        ImageRGBd fakeColorIndexImageLevel0;
        lmbd_fakeColorIndex(indexImageLevel0, fakeColorIndexImageLevel0);

        IO::save01_in_u8(fakeColorIndexImageLevel0, std::string(PCTS_DEBUG_DIRECTORY) + "fakeColorIndex" + std::to_string(s) + ".png");
    }

    return imageLevel0;
}

}

#endif
