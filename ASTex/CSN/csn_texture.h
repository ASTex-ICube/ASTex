#ifndef CSN_TEXTURE_H_
#define CSN_TEXTURE_H_

#include <Eigen/Eigen>
#include "Algo/ProceduralNoiseFiltering/gaussian_transfer.h"
#include "ASTex/utils.h"
#include "ASTex/pca.h"
#include <cmath>
#include "ASTex/colorspace_filters.h"
#include "ASTex/easy_io.h"
#include "ASTex/histogram.h"
#include "ASTex/rpn_utils.h" //delete this idiot

//#define CSN_ENABLE_DEBUG
#ifdef CSN_ENABLE_DEBUG
#define print_debug(x) std::cout << x << std::endl
#define save_debug(x, y) IO::save01_in_u8(x, y)
#else
#define print_debug(x)
#define save_debug(x, y)
#endif

namespace ASTex
{

namespace CSN
{

template<typename I>
class CSN_Texture
{
public:

	using ImageType					= I;
	using PixelType					= typename ImageType::PixelType;
	using DataType					= typename ImageType::DataType;
	using PixelPosType				= itk::Index<2>;
	using PcaImageType				= ImageRGB<DataType>;
	using PcaPixelType				= typename PcaImageType::PixelType;
	using PcaType					= PCA<DataType>;
	using GaussianTransferType		= Gaussian_transfer<PcaImageType>;
	using LutType					= PcaImageType;
	using ProceduralBlendingType	= std::function<PcaImageType(const PcaImageType &)>;
	using PtrImageType				= ImageGrayu64;

	CSN_Texture();
	~CSN_Texture();
	void setTexture(const ImageType &texture);
	void setCycles(const Eigen::Vector2d &firstCycle, const Eigen::Vector2d &secondCycle);
	void setUseCycles(bool b);
    void setGamma(double gamma); //<represents the exponant in the blending weights
	void setProceduralBlendingSubstitute(ProceduralBlendingType substitute);
	void setUsePca(bool b);
	void setUseGaussianTransfer(bool b);
	void setUseYCbCr(bool b);
	void setUseCyclicTransfer(bool b);
	void setUVScale(double uvScale);
	void setCyclicTransferPolicy(unsigned radius=0, unsigned samples=1);

	ImageRGBd debug_cycleEvaluationMap(int width, int height, Eigen::Vector2d center, double radius) const;
	void estimateCycles(const Eigen::Vector2d &guidX, const Eigen::Vector2d &guidY, double searchRadius, bool periodicity, unsigned stochasticGradientSearchEnergy=16);

	ImageType synthesize(unsigned width, unsigned height);
	PixelType testCycles(const I& texture, const Eigen::Vector2d &tx, const Eigen::Vector2d &ty, bool useCyclicAverage=false) const;
	PixelType testCycles_v2(const I& texture, const Eigen::Vector2d &tx, const Eigen::Vector2d &ty, unsigned nbSamplesX, unsigned nbSamplesY) const;

	const Eigen::Vector2d &cycleX() const;
	const Eigen::Vector2d &cycleY() const;

private:

	PcaPixelType proceduralTilingAndBlending(const PcaImageType &image, Eigen::Vector2d uv) const;
	void TriangleGrid(Eigen::Vector2d uv, float &w1, float &w2, float &w3,
						Eigen::Vector2i &vertex1, Eigen::Vector2i &vertex2, Eigen::Vector2i &vertex3) const;
	Eigen::Vector2d hash(const Eigen::Vector2d &p) const;
	Eigen::Vector2d cyclicHash(const Eigen::Vector2d &p) const;
	Eigen::Vector2d floor(const Eigen::Vector2d &v) const;
	Eigen::Vector2d fract(const Eigen::Vector2d &v) const;

	void estimateCycles_periodic(const Eigen::Vector2d &guidX, const Eigen::Vector2d &guidY, int errorWindow);

	ImageType				m_exemplar;
	Eigen::Vector2d			m_cycles[2];
	bool					m_useCycles;
	double					m_gamma;
	unsigned				m_largestCycleProduct;
	ProceduralBlendingType	m_func_proceduralBlendingSubstitute;
	bool					m_usePca;
	bool					m_useGaussianTransfer;
	bool					m_useYCbCr;
	bool					m_useCyclicTransfer;
	bool					m_useCyclicPCA;
	PtrImageType			m_transferPtrImage;
	PtrImageType			m_pcaPtrImage;
	double					m_uvScale;
	double					m_cyclicTransferRadius;
	unsigned				m_cyclicTransferSamples;

};

template<typename I>
CSN_Texture<I>::CSN_Texture() :
	m_exemplar(),
	m_cycles(),
	m_useCycles(false),
	m_gamma(2),
	m_largestCycleProduct(512),
	m_func_proceduralBlendingSubstitute(nullptr),
	m_usePca(true),
	m_useGaussianTransfer(true),
	m_useYCbCr(false),
	m_useCyclicTransfer(false),
	m_uvScale(1.0),
	m_cyclicTransferRadius(0),
	m_cyclicTransferSamples(1)
{}

template<typename I>
CSN_Texture<I>::~CSN_Texture()
{
}

template<typename I>
void CSN_Texture<I>::setTexture(const ImageType &texture)
{
	m_exemplar = texture;
}

template<typename I>
void CSN_Texture<I>::setCycles(const Eigen::Vector2d &firstCycle, const Eigen::Vector2d &secondCycle)
{
	m_cycles[0] = firstCycle;
	m_cycles[1] = secondCycle;
}

template<typename I>
void CSN_Texture<I>::setUseCycles(bool b)
{
	m_useCycles = b;
}

template<typename I>
void CSN_Texture<I>::setGamma(double gamma)
{
	m_gamma = gamma;
}

template<typename I>
void CSN_Texture<I>::setProceduralBlendingSubstitute(ProceduralBlendingType substitute)
{
	m_func_proceduralBlendingSubstitute = substitute;
}

template<typename I>
void CSN_Texture<I>::setUsePca(bool b)
{
	m_usePca = b;
}

template<typename I>
void CSN_Texture<I>::setUseGaussianTransfer(bool b)
{
	m_useGaussianTransfer = b;
}

template<typename I>
void CSN_Texture<I>::setUseYCbCr(bool b)
{
	m_useYCbCr = b;
}

template<typename I>
void CSN_Texture<I>::setUseCyclicTransfer(bool b)
{
	m_useCyclicTransfer = b;
}

template<typename I>
void CSN_Texture<I>::setUVScale(double uvScale)
{
	m_uvScale = uvScale;
}

template<typename I>
void CSN_Texture<I>::setCyclicTransferPolicy(unsigned radius, unsigned samples)
{
	assert(samples%2!=0 && "There needs to be an odd number of samples in order to allow the middle pixel to be read back after transfer");
	m_cyclicTransferRadius = radius;
	m_cyclicTransferSamples = samples;
}






template<typename I>
typename CSN_Texture<I>::ImageType CSN_Texture<I>::synthesize(unsigned width, unsigned height)
{
	assert(m_exemplar.is_initialized());
	if(width==0)
		width=m_exemplar.width()*1;
	if(height==0)
		height=m_exemplar.height()*1;
	unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	assert(pixelSize <= 3 && "CSN_Texture::synthesize: Cannot use PCA with images of dimensions higher than 3!");

	//find the largest product of the cycles. This is to ensure the patterns are evenly distributed throughout the domain.
	//I'm not sure this works as intended yet. It's either that or each cycle gets its own divisor. Or something else.
	if(m_useCycles)
	{
		double cycleProduct = std::ceil(1.0/std::max(m_cycles[0][0], m_cycles[1][0])) * std::ceil(1.0/std::max(m_cycles[0][1], m_cycles[1][1]));
		if(std::isnan(cycleProduct) || std::isinf(cycleProduct))
		{
			cycleProduct = std::max(std::max(m_cycles[0][0], m_cycles[1][0]), std::max(m_cycles[0][1], m_cycles[1][1]));
			assert(cycleProduct != 0 && "Cycle product is always 0! (Cycles not set?)");
			cycleProduct = 1.0/cycleProduct;
		}
		m_largestCycleProduct = unsigned(cycleProduct);
		std::cout << "largest cycle product: " << m_largestCycleProduct << std::endl;
	}

	ImageType output;
	output.initItk(width, height, true);
	LutType lut;
	lut.initItk(128, 1);
	DataType *dataPix;
	const DataType *dataPixConst;

	auto toPcaImageType = [&] (const ImageType &texture) -> PcaImageType
	{
		PcaImageType pcaTexture;
		pcaTexture.initItk(texture.width(), texture.height());
		pcaTexture.for_all_pixels([&] (PcaPixelType &pcaPix, int x, int y)
		{
			const PixelType &pix = texture.pixelAbsolute(x, y);
			dataPixConst = reinterpret_cast<const DataType *>(&pix);
			if(pixelSize == 3)
			{ //image is rgb
				for(unsigned i=0; i<3; ++i)
				{
					pcaPix[i] = dataPixConst[i];
				}
			}
			else
			{ //image is most likely gray
				for(unsigned i=0; i<3; ++i)
				{
					unsigned j=i%pixelSize;
					pcaPix[i] = dataPixConst[j];
				}
			}
		});
		return pcaTexture;
	};

	auto fromPcaImageType = [&] (const PcaImageType &pcaTexture) -> ImageType
	{
		ImageType texture;
		texture.initItk(pcaTexture.width(), pcaTexture.height());
		pcaTexture.for_all_pixels([&] (const PcaPixelType &pcaPix, int x, int y)
		{
			PixelType &pix = texture.pixelAbsolute(x, y);
			dataPix = reinterpret_cast<DataType *>(&pix);
			for(unsigned i=0; i<pixelSize; ++i)
				dataPix[i] = pcaPix[i];
		});
		return texture;
	};

	PcaImageType pcaTexture = toPcaImageType(m_exemplar);
	if(m_useYCbCr)
	{
		pcaTexture.for_all_pixels([] (typename PcaImageType::PixelType &pix)
		{
			ColorSpace::fonctorRGBtoYCbCr<	typename PcaImageType::PixelType,
											typename PcaImageType::PixelType> functor;
			pix=functor(pix);
		});
	}
	Eigen::Vector3f colorSpaceVector1, colorSpaceVector2, colorSpaceVector3, colorSpaceOrigin;
	PcaImageType pcaGaussianTexture;
	GaussianTransferType gtt;

	//////////////////////////////INPUT////////////////////////////////////

	PcaType pca(pcaTexture);
	if(m_usePca)
	{
		if(m_useCycles && m_useCyclicPCA)
		{
			int pcaWidth, pcaHeight, GWidth, GHeight;
			pcaWidth = int(std::floor(1.0/m_cycles[0][0]));
			pcaHeight = int(std::floor(1.0/m_cycles[1][1]));
			GWidth = pcaWidth*m_cyclicTransferSamples;
			GHeight = pcaHeight*m_cyclicTransferSamples;
			PcaImageType pcaSubTexture, pcaGaussianSubTexture, pcaInputCopy;
			pcaInputCopy.initItk(pcaTexture.width(), pcaTexture.height());
			pcaInputCopy.copy_pixels(pcaTexture);
			pcaSubTexture.initItk(GWidth, GHeight);
			pcaGaussianSubTexture.initItk(GWidth, GHeight);
			m_pcaPtrImage.initItk(	std::round((m_cycles[0][0])*(m_exemplar.width()) +0.5),
										std::round((m_cycles[1][1])*(m_exemplar.height()) +0.5) );
			m_pcaPtrImage.for_all_pixels([&] (PtrImageType::PixelType &ptr, int inX, int inY)
			{
				auto getColor = [&] (int subX, int subY, int wx, int wy) -> PcaPixelType
				{
					double x, y;
					x = std::max(std::round(inX + (subX*m_cycles[0][0] + subY*m_cycles[1][0])*pcaInputCopy.width()) + (wx - m_cyclicTransferSamples/2.0 + 0.5)*m_cyclicTransferRadius, 0.0);
					y = std::max(std::round(inY + (subX*m_cycles[0][1] + subY*m_cycles[1][1])*pcaInputCopy.height()) + (wy - m_cyclicTransferSamples/2.0 + 0.5)*m_cyclicTransferRadius, 0.0);
					return bilinear_interpolation(pcaInputCopy, x, y, true);
				};
				for(int subX=0; subX<pcaWidth; ++subX)
					for(int subY=0; subY<pcaHeight; ++subY)
						for(unsigned wx=0; wx<m_cyclicTransferSamples; ++wx)
							for(unsigned wy=0; wy<m_cyclicTransferSamples; ++wy)
							{
								pcaSubTexture.pixelAbsolute(subX*m_cyclicTransferSamples+wx, subY*m_cyclicTransferSamples+wy) = getColor(subX, subY, int(wx), int(wy));
							}
				PcaType *subPca = new PcaType(pcaSubTexture);
				ptr = reinterpret_cast<PtrImageType::PixelType>(subPca);
				MaskBool mb_alwaysTrue(pcaSubTexture.width(), pcaSubTexture.height());
				mb_alwaysTrue |= [] (int, int) {return true;};
				subPca->computePCA(mb_alwaysTrue);
				subPca->project(pcaSubTexture);

				for(int subX=0; subX<pcaWidth; ++subX)
					for(int subY=0; subY<pcaHeight; ++subY)
					{
						PixelPosType inCoordinates;
						inCoordinates[0] = int(std::round(inX + (subX*m_cycles[0][0] + subY*m_cycles[1][0])*pcaInputCopy.width()))%pcaInputCopy.width();
						inCoordinates[1] = int(std::round(inY + (subX*m_cycles[0][1] + subY*m_cycles[1][1])*pcaInputCopy.height()))%pcaInputCopy.height();
						pcaTexture.pixelAbsolute(inCoordinates) = pcaSubTexture.pixelAbsolute(	subX*m_cyclicTransferSamples+int(m_cyclicTransferSamples/2.0-0.5),
																								subY*m_cyclicTransferSamples+int(m_cyclicTransferSamples/2.0-0.5));
					}
			});
		}
		else
		{
			MaskBool mb_alwaysTrue(pcaTexture.width(), pcaTexture.height());
			mb_alwaysTrue |= [] (int, int) {return true;};
			pca.computePCA(mb_alwaysTrue);
			pca.project(pcaTexture);
		}
	}
	pcaGaussianTexture.initItk(pcaTexture.width(), pcaTexture.height());
	if(m_useGaussianTransfer)
	{
		if(m_useCycles && m_useCyclicTransfer)
		{
			int transferWidth, transferHeight, GWidth, GHeight;
			transferWidth = int(std::floor(1.0/m_cycles[0][0]));
			transferHeight = int(std::floor(1.0/m_cycles[1][1]));
			GWidth = transferWidth*m_cyclicTransferSamples;
			GHeight = transferHeight*m_cyclicTransferSamples;
			PcaImageType pcaSubTexture, pcaGaussianSubTexture, pcaInputCopy;
			pcaInputCopy.initItk(pcaTexture.width(), pcaTexture.height());
			pcaInputCopy.copy_pixels(pcaTexture);
			pcaSubTexture.initItk(GWidth, GHeight);
			pcaGaussianSubTexture.initItk(GWidth, GHeight);
			m_transferPtrImage.initItk(	std::round((m_cycles[0][0])*(m_exemplar.width()) +0.5),
										std::round((m_cycles[1][1])*(m_exemplar.height()) +0.5) );
			m_transferPtrImage.for_all_pixels([&] (PtrImageType::PixelType &ptr, int inX, int inY)
			{
				LutType *subLut = new LutType;
				subLut->initItk(256, 1);
				ptr = reinterpret_cast<PtrImageType::PixelType>(subLut);
				auto getColor = [&] (int subX, int subY, int wx, int wy) -> PcaPixelType
				{
					double x, y;
					x = std::max(std::round(inX + (subX*m_cycles[0][0] + subY*m_cycles[1][0])*pcaInputCopy.width()) + (wx - m_cyclicTransferSamples/2.0 + 0.5)*m_cyclicTransferRadius, 0.0);
					y = std::max(std::round(inY + (subX*m_cycles[0][1] + subY*m_cycles[1][1])*pcaInputCopy.height()) + (wy - m_cyclicTransferSamples/2.0 + 0.5)*m_cyclicTransferRadius, 0.0);
					return bilinear_interpolation(pcaInputCopy, x, y, true);
				};
				for(int subX=0; subX<transferWidth; ++subX)
					for(int subY=0; subY<transferHeight; ++subY)
						for(unsigned wx=0; wx<m_cyclicTransferSamples; ++wx)
							for(unsigned wy=0; wy<m_cyclicTransferSamples; ++wy)
							{
								pcaSubTexture.pixelAbsolute(subX*m_cyclicTransferSamples+wx, subY*m_cyclicTransferSamples+wy) = getColor(subX, subY, wx, wy);
							}
				gtt.ComputeTinput(pcaSubTexture, pcaGaussianSubTexture);
				gtt.ComputeinvT(pcaSubTexture, *subLut);
				for(int subX=0; subX<transferWidth; ++subX)
					for(int subY=0; subY<transferHeight; ++subY)
					{
						PixelPosType inCoordinates;
						inCoordinates[0] = int(std::round(inX + (subX*m_cycles[0][0] + subY*m_cycles[1][0])*pcaInputCopy.width()))%pcaInputCopy.width();
						inCoordinates[1] = int(std::round(inY + (subX*m_cycles[0][1] + subY*m_cycles[1][1])*pcaInputCopy.height()))%pcaInputCopy.height();
						pcaGaussianTexture.pixelAbsolute(inCoordinates) = pcaGaussianSubTexture.pixelAbsolute(subX*m_cyclicTransferSamples+int(m_cyclicTransferSamples/2.0-0.5),
																											  subY*m_cyclicTransferSamples+int(m_cyclicTransferSamples/2.0-0.5));
					}
			});
		}
		else
		{
			gtt.ComputeTinput(pcaTexture, pcaGaussianTexture);
			gtt.ComputeinvT(pcaTexture, lut);
		}
	}
	else
	{
		pcaGaussianTexture.copy_pixels(pcaTexture);
	}

	//////////////////////////////OUTPUT///////////////////////////////////

	PcaImageType pcaOutput;
	pcaOutput.initItk(output.width(), output.height());
	if(m_func_proceduralBlendingSubstitute == nullptr)
	{
		pcaOutput.for_all_pixels([&] (PcaPixelType &pix, int x, int y)
		{
			Eigen::Vector2d uv;
			uv[0] = double(x)/m_exemplar.width();
			uv[1] = double(y)/m_exemplar.height();
			pix = proceduralTilingAndBlending(pcaGaussianTexture, uv);
		});
	}
	else
	{
		pcaOutput = m_func_proceduralBlendingSubstitute(pcaGaussianTexture);
	}
	if(m_useGaussianTransfer)
	{
		if(m_useCycles && m_useCyclicTransfer)
		{
			unsigned transferWidth, transferHeight, GWidth, GHeight;
			transferWidth = unsigned(std::ceil(double(width)/m_exemplar.width()) * std::floor(1.0/m_cycles[0][0]));
			transferHeight = unsigned(std::ceil(double(height)/m_exemplar.height()) * std::floor(1.0/m_cycles[1][1]));
			GWidth = transferWidth*m_cyclicTransferSamples;
			GHeight = transferHeight*m_cyclicTransferSamples;
			print_debug("GWidth: " << GWidth);
			print_debug("GHeight: " << GHeight);
			PcaImageType pcaSubTexture, pcaGaussianSubTexture, pcaOutputCopy;
			pcaOutputCopy.initItk(pcaOutput.width(), pcaOutput.height());
			pcaOutputCopy.copy_pixels(pcaOutput);
			pcaSubTexture.initItk(GWidth, GHeight);
			m_transferPtrImage.for_all_pixels([&] (PtrImageType::PixelType &ptr, int inX, int inY)
			{
				LutType *subLut = reinterpret_cast<LutType *>(ptr);
				pcaSubTexture.for_all_pixels([&] (typename PcaImageType::PixelType &pix, int subX, int subY)
				{
					PixelPosType inCoordinates;
					if(unsigned(std::ceil(double(width)/m_exemplar.width()) * std::floor(1.0/m_cycles[0][0])))
					inCoordinates[0] = int(std::round(inX + (subX*m_cycles[0][0] + subY*m_cycles[1][0])*m_exemplar.width()))%pcaOutputCopy.width();
					inCoordinates[1] = int(std::round(inY + (subX*m_cycles[0][1] + subY*m_cycles[1][1])*m_exemplar.height()))%pcaOutputCopy.height();
					pix = pcaOutputCopy.pixelAbsolute(inCoordinates);
				});
				pcaGaussianSubTexture = GaussianTransferType::invT(pcaSubTexture, *subLut);
				pcaGaussianSubTexture.for_all_pixels([&] (const typename PcaImageType::PixelType &pix, int subX, int subY)
				{
					PixelPosType inCoordinates;
					inCoordinates[0] = int(std::round(inX + (subX*m_cycles[0][0] + subY*m_cycles[1][0])*m_exemplar.width()))%pcaOutputCopy.width();
					inCoordinates[1] = int(std::round(inY + (subX*m_cycles[0][1] + subY*m_cycles[1][1])*m_exemplar.height()))%pcaOutputCopy.height();
					pcaOutput.pixelAbsolute(inCoordinates) = pix;
				});
				save_debug(pcaGaussianSubTexture, "/home/nlutz/pcaInput.png");
			});
		}
		else
		{
			pcaOutput = GaussianTransferType::invT(pcaOutput, lut);
		}
	}
	if(m_usePca)
	{
		if(m_useCycles && m_useCyclicPCA)
		{
			pcaTexture.initItk(pcaOutput.width(), pcaOutput.height());
			unsigned transferWidth, transferHeight, GWidth, GHeight;
			transferWidth = unsigned(std::ceil(double(width)/m_exemplar.width() * std::floor(1.0/m_cycles[0][0])));
			transferHeight = unsigned(std::ceil(double(height)/m_exemplar.height() * std::floor(1.0/m_cycles[1][1])));
			GWidth = transferWidth*m_cyclicTransferSamples;
			GHeight = transferHeight*m_cyclicTransferSamples;
			PcaImageType pcaSubTexture, pcaUnfoldedSubTexture, pcaOutputCopy;
			pcaOutputCopy.initItk(pcaOutput.width(), pcaOutput.height());
			pcaOutputCopy.copy_pixels(pcaOutput);
			pcaSubTexture.initItk(GWidth, GHeight);
			pcaUnfoldedSubTexture.initItk(pcaOutput.width(), pcaOutput.height());
			m_pcaPtrImage.for_all_pixels([&] (PtrImageType::PixelType &ptr, int inX, int inY)
			{
				PcaType *subPca = reinterpret_cast<PcaType *>(ptr);
				pcaSubTexture.for_all_pixels([&] (typename PcaImageType::PixelType &pix, int subX, int subY)
				{
					PixelPosType inCoordinates;
					//if(unsigned(std::ceil(double(width)/m_exemplar.width()) * std::floor(1.0/m_cycles[0][0])))
					inCoordinates[0] = int(std::round(inX + (subX*m_cycles[0][0] + subY*m_cycles[1][0])*m_exemplar.width()))%pcaOutputCopy.width();
					inCoordinates[1] = int(std::round(inY + (subX*m_cycles[0][1] + subY*m_cycles[1][1])*m_exemplar.height()))%pcaOutputCopy.height();
					pix = pcaOutputCopy.pixelAbsolute(inCoordinates);
				});

				subPca->back_project(pcaSubTexture, pcaUnfoldedSubTexture);

				pcaUnfoldedSubTexture.for_all_pixels([&] (const typename PcaImageType::PixelType &pix, int subX, int subY)
				{
					PixelPosType inCoordinates;
					inCoordinates[0] = int(std::round(inX + (subX*m_cycles[0][0] + subY*m_cycles[1][0])*m_exemplar.width()))%pcaOutputCopy.width();
					inCoordinates[1] = int(std::round(inY + (subX*m_cycles[0][1] + subY*m_cycles[1][1])*m_exemplar.height()))%pcaOutputCopy.height();
					pcaTexture.pixelAbsolute(inCoordinates) = pix;
				});
			});
		}
		else
		{
			pca.back_project(pcaOutput, pcaTexture);
		}
	}
	else
		pcaTexture = pcaOutput;
	if(m_useYCbCr)
	{
		pcaTexture.for_all_pixels([] (typename PcaImageType::PixelType &pix)
		{
			ColorSpace::fonctorYCbCrtoRGB<	typename PcaImageType::PixelType,
											typename PcaImageType::PixelType> functor;
			pix=functor(pix);
		});
	}
	if(m_useGaussianTransfer && m_useCycles && m_useCyclicTransfer)
	{
		m_transferPtrImage.for_all_pixels([&] (PtrImageType::PixelType &pix)
		{
			GaussianTransferType *gtt = reinterpret_cast<GaussianTransferType *>(pix);
			delete gtt;
		});
	}
	if(m_usePca && m_useCycles && m_useCyclicPCA)
	{
		m_pcaPtrImage.for_all_pixels([&] (PtrImageType::PixelType &pix)
		{
			PcaType *pca = reinterpret_cast<PcaType *>(pix);
			delete pca;
		});
	}
	output = fromPcaImageType(pcaTexture);
	return output;
}



///////////////////////////////////////////////////////////////////////////////


template<typename I>
void CSN_Texture<I>::estimateCycles(const Eigen::Vector2d &guidX, const Eigen::Vector2d &guidY, double searchRadius, bool periodicity, unsigned stochasticGradientSearchEnergy)
{
	assert(m_exemplar.is_initialized());
	assert(searchRadius < 0.25);
	unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	if(periodicity)
	{
		int errorWindow = std::max(	int(std::round(1.0/guidX[0])) - int(std::round(1.0/(guidX[0] - searchRadius))),
									int(std::round(1.0/guidY[1])) - int(std::round(1.0/(guidY[1] - searchRadius))));
		estimateCycles_periodic(guidX, guidY, errorWindow);
	}
	else
	{
		ImageType centeredTexture;
		centeredTexture.initItk(m_exemplar.width(), m_exemplar.height());
		PixelType mean;
		Histogram<ImageType> histogram(m_exemplar);
		mean = histogram.meanPixelType();
		centeredTexture.for_all_pixels([&] (PixelType &pix, int x, int y)
		{
			pix = m_exemplar.pixelAbsolute(x, y) - mean;
		});
		auto adjustSizeForEvaluation = [&] (const ImageType &texture, Eigen::Vector2d &cycle) -> ImageType //scale to match the cycle, and cut incomplete parallelograms.
		{
			print_debug("cycle: " << cycle);
			double xScale, yScale;
			double epsilonX = 1.0/texture.width();
			double epsilonY = 1.0/texture.height();
			Eigen::Vector2d adjustedCycle; //prevents the texture from having more than 4x the number of texels for a cycle too close to 0 anyway, and division by 0.
			adjustedCycle[0] = cycle[0] >= epsilonX ? cycle[0] : 1.0;
			adjustedCycle[1] = cycle[1] >= epsilonY ? cycle[1] : 1.0;
			xScale = ceil(texture.width()*adjustedCycle[0]) / (texture.width()*adjustedCycle[0]);
			yScale = ceil(texture.height()*adjustedCycle[1]) / (texture.height()*adjustedCycle[1]);
			print_debug("xScale: " << xScale);
			print_debug("yScale: " << yScale);
			ImageType scaledTexture;
			scaledTexture.initItk(int(std::floor(texture.width()*xScale)),
								  int(std::floor(texture.height()*yScale)), true);
			scaledTexture.for_all_pixels([&] (PixelType &pix, int x, int y)
			{
				double xd, yd;
				xd = double(x)/(scaledTexture.width()-1) * (texture.width()-1);
				yd = double(y)/(scaledTexture.height()-1) * (texture.height()-1);
				pix = bilinear_interpolation(texture, xd, yd, false);
			});
			ImageType cutTexture;
			itk::Size<2> rectangleSize;
			rectangleSize[0] = int(round(texture.width()*xScale*adjustedCycle[0]));
			rectangleSize[1] = int(round(texture.height()*yScale*adjustedCycle[1]));
			int nbRectanglesWidth =		int(std::floor(texture.width()*xScale / rectangleSize[0]));
			int nbRectanglesHeight =	int(std::floor(texture.height()*yScale / rectangleSize[1]));
			cutTexture.initItk(nbRectanglesWidth*rectangleSize[0], nbRectanglesHeight*rectangleSize[1]);
			cutTexture.for_all_pixels([&] (PixelType &pix, int x, int y)
			{
				pix = scaledTexture.pixelAbsolute(x, y);
			});
			cycle[0] = rectangleSize[0]/double(cutTexture.width());
			cycle[1] = rectangleSize[1]/double(cutTexture.height());
			return cutTexture;
		};

		double epsilonX = 1.0/m_exemplar.width();
		double epsilonY = 1.0/m_exemplar.height();

		auto projectionPixelTypeToDouble = [pixelSize] (PixelType projection) -> double
		{
			double projectionDouble = 0;
			DataType *dataProjection = reinterpret_cast<DataType *>(&projection);
			for(unsigned i=0; i<pixelSize; ++i)
				projectionDouble += dataProjection[i];
			print_debug("projectionDouble: " << projectionDouble);
			return projectionDouble*projectionDouble;
		};

		Eigen::Vector2d cycleX = guidX, cycleY = guidY;
		double epsilon = searchRadius/256.0;
		//first step: find out a good cycleX[0]
		Eigen::Vector2d alteredGuide = guidX;
		ImageType texture = adjustSizeForEvaluation(centeredTexture, alteredGuide);
		PixelType projection = testCycles(texture, alteredGuide, Eigen::Vector2d(0, 1));
		double projectionValue = projectionPixelTypeToDouble(projection);
		double maxProjectionValue = projectionValue;
		Eigen::Vector2d maxProjectionCycle = cycleX;
		double gradDistance = searchRadius/2.0;

		auto evaluateAndUpdateCycles = [&] (const Eigen::Vector2d &cycleX, const Eigen::Vector2d &cycleY, bool xDominates)
		{
			if(xDominates)
			{
				Eigen::Vector2d alteredCycle = cycleX;
				texture = adjustSizeForEvaluation(centeredTexture, alteredCycle);
				projection = testCycles(texture, alteredCycle, Eigen::Vector2d(0, 1));
			}
			else
			{
				Eigen::Vector2d alteredCycle = cycleY;
				texture = adjustSizeForEvaluation(centeredTexture, alteredCycle);
				projection = testCycles(texture, Eigen::Vector2d(1, 0), alteredCycle);
			}
			save_debug(texture, "/home/nlutz/interpolatedTexture.png");
			projectionValue = projectionPixelTypeToDouble(projection);
			if(projectionValue < maxProjectionValue)
			{
				maxProjectionValue = projectionValue;
				maxProjectionCycle = xDominates ? cycleX : cycleY;
			}
		};

		Eigen::Vector2d cycle;
		//stochastic gradient search
		for(unsigned i=0; i<stochasticGradientSearchEnergy; ++i)
		{
			Eigen::Vector2d point;
			point[0] = std::cos(std::rand()/double(RAND_MAX) * 2*M_PI);
			point[1] = std::sin(std::rand()/double(RAND_MAX) * 2*M_PI);
			cycle = cycleX + point*searchRadius;
			cycle[1] = cycle[1] > epsilonY ? cycle[1] : 0;
			if(cycle[0]>=0 && cycle[1]>=0 && cycle[0]<=0.5 && cycle[1]<=0.5)
			{
				evaluateAndUpdateCycles(cycle, Eigen::Vector2d(0, 1), true);
				cycleX = maxProjectionCycle;
			}
		}
		while(gradDistance > epsilon)
		{
			//+x
			cycle = cycleX;
			cycle[0] += gradDistance;
			evaluateAndUpdateCycles(cycle, Eigen::Vector2d(0, 1), true);

			//-x
			cycle = cycleX;
			cycle[0] -= gradDistance;

			if(cycle[0]>=0)
				evaluateAndUpdateCycles(cycle, Eigen::Vector2d(0, 1), true);

			//+y
			cycle = cycleX;
			if(cycle[1] > epsilonY)
			{
				cycle[1] += gradDistance;
				evaluateAndUpdateCycles(cycle, Eigen::Vector2d(0, 1), true);
			}

			//-y
			cycle = cycleX;
			cycle[1] -= gradDistance;
			if(cycle[1]>=0)
				evaluateAndUpdateCycles(cycle, Eigen::Vector2d(0, 1), true);

			cycleX = maxProjectionCycle;
			gradDistance /= 2.0;
		}
		m_cycles[0] = cycleX;

		alteredGuide = guidY;
		texture = adjustSizeForEvaluation(centeredTexture, alteredGuide);
		projection = testCycles(texture, Eigen::Vector2d(1, 0), alteredGuide);
		projectionValue = projectionPixelTypeToDouble(projection);
		maxProjectionValue = projectionValue;
		maxProjectionCycle = cycleY;
		gradDistance = searchRadius/2.0;

		for(unsigned i=0; i<stochasticGradientSearchEnergy; ++i)
		{
			Eigen::Vector2d point;
			double radius = std::sqrt(std::rand()/double(RAND_MAX));
			point[0] = radius * std::cos(std::rand()/double(RAND_MAX) * 2*M_PI);
			point[1] = radius * std::sin(std::rand()/double(RAND_MAX) * 2*M_PI);
			cycle = cycleY + point*searchRadius;
			cycle[0] = cycle[0] > epsilonX ? cycle[0] : 0;
			print_debug("cycle: " << cycle);
			if(cycle[0]>=0 && cycle[1]>=0 && cycle[0]<=0.5 && cycle[1]<=0.5)
			{
				evaluateAndUpdateCycles(Eigen::Vector2d(1, 0), cycle, false);
				cycleY = maxProjectionCycle;
			}
		}
		while(gradDistance > epsilon)
		{
			Eigen::Vector2d cycle;
			//+x
			cycle = cycleY;
			if(cycle[0] > epsilonX)
			{
				cycle[0] += gradDistance;
				evaluateAndUpdateCycles(Eigen::Vector2d(1, 0), cycle, false);
			}

			//-x
			cycle = cycleY;
			cycle[0] -= gradDistance;
			if(cycle[0]>=0)
				evaluateAndUpdateCycles(Eigen::Vector2d(1, 0), cycle, false);

			//+y
			cycle = cycleY;
			cycle[1] += gradDistance;
			evaluateAndUpdateCycles(Eigen::Vector2d(1, 0), cycle, false);

			//-y
			cycle = cycleY;
			cycle[1] -= gradDistance;
			if(cycle[1]>=0)
				evaluateAndUpdateCycles(Eigen::Vector2d(1, 0), cycle, false);

			cycleY = maxProjectionCycle;
			gradDistance /= 2.0;
		}
		m_cycles[1] = cycleY;
	}
}

template<typename I>
ImageRGBd CSN_Texture<I>::debug_cycleEvaluationMap(int width, int height, Eigen::Vector2d center, double radius) const
{
	assert(m_exemplar.is_initialized());
	unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	ImageType centeredTexture;
	centeredTexture.initItk(m_exemplar.width(), m_exemplar.height());
	centeredTexture.for_all_pixels([&] (PixelType &pix, int x, int y)
	{
		pix = m_exemplar.pixelAbsolute(x, y);
	});

	auto adjustSizeForEvaluation = [&] (const ImageType &texture, Eigen::Vector2d &cycle) -> ImageType //scale to match the cycle, and cut incomplete parallelograms.
	{
		double xScale, yScale;
		double epsilonX = 1.0/texture.width();
		double epsilonY = 1.0/texture.height();
		Eigen::Vector2d adjustedCycle; //prevents the texture from having more than 4x the number of texels for a cycle too close to 0 anyway, and division by 0.
		adjustedCycle[0] = cycle[0] > epsilonX ? cycle[0] : 1.0;
		adjustedCycle[1] = cycle[1] > epsilonY ? cycle[1] : 1.0;
		xScale = ceil(texture.width()*adjustedCycle[0]) / (texture.width()*adjustedCycle[0]);
		yScale = ceil(texture.height()*adjustedCycle[1]) / (texture.height()*adjustedCycle[1]);
		print_debug("xScale: " << xScale);
		print_debug("yScale: " << yScale);
		ImageType scaledTexture;
		scaledTexture.initItk(int(std::floor(texture.width()*xScale)),
							  int(std::floor(texture.height()*yScale)), true);
		print_debug("scaledTexture theorical size: " << texture.width()*xScale << ", " << texture.height()*yScale);
		print_debug("scaledTexture size: " << scaledTexture.width() << ", " << scaledTexture.height());
		print_debug("yScale: " << yScale);
		scaledTexture.for_all_pixels([&] (PixelType &pix, int x, int y)
		{
			double xd, yd;
			xd = double(x)/(scaledTexture.width()-1) * (texture.width()-1);
			yd = double(y)/(scaledTexture.height()-1) * (texture.height()-1);
			pix = bilinear_interpolation(texture, xd, yd, false);
		});
		ImageType cutTexture;
		itk::Size<2> rectangleSize;
		rectangleSize[0] = int(round(texture.width()*xScale*adjustedCycle[0]));
		rectangleSize[1] = int(round(texture.height()*yScale*adjustedCycle[1]));
		print_debug("rectangleSizeX (theorical): " << texture.width()*xScale*adjustedCycle[0]);
		print_debug("rectangleSize[0]: " << rectangleSize[0]);
		print_debug("rectangleSize[1]: " << rectangleSize[1]);
		int nbRectanglesWidth =		int(std::floor(texture.width()*xScale / rectangleSize[0]));
		int nbRectanglesHeight =	int(std::floor(texture.height()*yScale / rectangleSize[1]));
		print_debug("nbRectanglesWidth: " << nbRectanglesWidth);
		print_debug("nbRectanglesHeight: " << nbRectanglesHeight);
		cutTexture.initItk(nbRectanglesWidth*rectangleSize[0], nbRectanglesHeight*rectangleSize[1]);
		cutTexture.for_all_pixels([&] (PixelType &pix, int x, int y)
		{
			pix = scaledTexture.pixelAbsolute(x, y);
		});
		print_debug("cutTexture size: " << cutTexture.width() << ", " << cutTexture.height());
		cycle[0] = rectangleSize[0]/double(cutTexture.width());
		cycle[1] = rectangleSize[1]/double(cutTexture.height());
		return cutTexture;
	};

	double epsilonX = 1.0/m_exemplar.width();
	double epsilonY = 1.0/m_exemplar.height();

	Eigen::Vector2d cycleX, cycleY;
	cycleX[0]=center[0]-radius;
	cycleX[1] = 0;
	cycleY[1]=center[1]-radius;
	cycleY[0] = 0;
	double minimumProjection = std::numeric_limits<double>::infinity(), maximumProjection = -std::numeric_limits<double>::infinity();
	ImageGrayd d;
	d.initItk(width, height, true);
	int printableCount = 0;
	int percentageSave = -1;
	d.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
	{
		++printableCount;
		int percentage = (printableCount*100+1)/(d.width()*d.height());
		if(percentage != percentageSave)
		{
			if(percentageSave != -1)
			{
				if(percentageSave <10)
					std::cout << "\b\b";
				else
					std::cout << "\b\b\b";
			}
			percentageSave = percentage;
			std::cout << percentage << "%";
			std::flush(std::cout);
		}
		cycleX[0] = (center[0]-radius) + double(x)/(d.width()-1) * radius * 2;
		cycleY[1] = (center[1]-radius) + double(y)/(d.height()-1) * radius * 2;
		print_debug("cycleX[0]: " << cycleX[0]);
		print_debug("cycleY[1]: " << cycleY[1]);
		if(cycleX[0]>epsilonX && cycleY[1]>epsilonY && cycleX[0]<=0.5 && cycleY[1]<=0.5)
		{
			Eigen::Vector2d alteredCycle = Eigen::Vector2d(cycleX[0], cycleY[1]);
			ImageType texture = centeredTexture; //adjustSizeForEvaluation(centeredTexture, alteredCycle);
			texture.for_all_pixels([&] (PixelType &pix)
			{
				//I don't know why, I don't want to know why, I shouldn't have to wonder why, but without this line,
				//the red channel of "texture" fucks up the testCycles function. Something is fucked up with the core.
			});
			print_debug("alteredCycle: " << alteredCycle);
			print_debug("alteredCycle[0] * alteredWidth: " << alteredCycle[0]*texture.width());
			PixelType projection = testCycles_v2(texture, Eigen::Vector2d(alteredCycle[0], 0), Eigen::Vector2d(0, alteredCycle[1]), 50, 50);
			DataType *dataProjection = reinterpret_cast<DataType *>(&projection);
			for(unsigned i=0; i<pixelSize; ++i)
				pix += dataProjection[i];
			if(!std::isnan(pix))
			{
				minimumProjection = std::min(pix, minimumProjection);
				maximumProjection = std::max(pix, maximumProjection);
			}
			print_debug("d(cycle): " << pix);
		}
	});
	d.for_all_pixels([&] (ImageGrayd::PixelType &pix)
	{
		pix -= minimumProjection;
		pix /= (maximumProjection-minimumProjection);
	});

	ImageRGBd::PixelType node0, node1, node2;
	node0[0]=0.2;
	node0[1]=0.1;
	node0[2]=1.0;

	node1[0]=0.3;
	node1[1]=1.0;
	node1[2]=0.1;

	node2[0]=1.0;
	node2[1]=0.5;
	node2[2]=0.3;

	ImageRGBd map;
	map.initItk(d.width(), d.height());
	map.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
	{
		ImageGrayd::PixelType projection = d.pixelAbsolute(x, y);
		pix =	projection < 0.5 ?	node0*(1.0-projection*2) + node1*projection*2 :
									node1*(1.0-(projection-0.5)*2) + node2*(projection-0.5)*2;
	});
	return map;
}


template<typename I>
void CSN_Texture<I>::estimateCycles_periodic(const Eigen::Vector2d &guidX, const Eigen::Vector2d &guidY, int errorWindow)
{
	//Texel values must be floating and between 0 and 1 for this function to work.
	//If not, use a separate CSN_Texture class with a normalized input to estimate the cycles.
	assert(m_exemplar.is_initialized());
	unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	ImageType centeredTexture;
	centeredTexture.initItk(m_exemplar.width(), m_exemplar.height());
	PixelType mean;
	Histogram<ImageType> histogram(m_exemplar);
	mean = histogram.meanPixelType();
	centeredTexture.for_all_pixels([&] (PixelType &pix, int x, int y)
	{
		pix = m_exemplar.pixelAbsolute(x, y) - mean;
	});

	ImageType projectionImage, averageProjectionImage;
	//TODO: turn the two following steps in a single lambda.
	//first step: finding Cx
	PixelType projection;
	double projectionValue = 0;
	double maxProjectionValue = 0;
	int xProjection = 0;
	int yProjection = 0;
	int maxXProjection = 0;
	int maxYProjection = 0;
	xProjection = std::round(1.0/guidX[0]);
	if(guidX[1] == 0)
		yProjection = 0;
	else
		yProjection = std::round(1.0/guidX[1]);
	for(int cxInv=std::max(1, xProjection-errorWindow); cxInv<=xProjection+errorWindow; ++cxInv)
	{
		for(int cyInv=std::max(0, yProjection-errorWindow); cyInv<=yProjection+errorWindow; ++cyInv)
		{
			double cx = 1.0/cxInv, cy = cyInv>1 ? 1.0/cyInv : 0;
			projection = testCycles(centeredTexture, Eigen::Vector2d(cx, cy), Eigen::Vector2d(0, 1.0));
			DataType *dataProjection = reinterpret_cast<DataType *>(&projection);
			for(unsigned i=0; i<pixelSize; ++i)
			{
				projectionValue += dataProjection[i]*dataProjection[i];
			}
			if(maxProjectionValue < projectionValue)
			{
				maxProjectionValue = projectionValue;
				maxXProjection = cxInv;
				maxYProjection = cyInv > 1 ? cyInv : 0;
			}
		}
	}
	m_cycles[0][0] = 1.0/maxXProjection;
	m_cycles[0][1] = maxYProjection > 0 ? (maxYProjection == 1 ? 0 : 1.0/maxYProjection) : 0; //1 is th

	//second step: finding Cy
	projectionValue = 0;
	maxProjectionValue = 0;
	xProjection = 0;
	yProjection = 0;
	maxXProjection = 0;
	maxYProjection = 0;
	yProjection = std::round(1.0/guidY[1]);
	if(guidY[0] == 0)
		xProjection = 0;
	else
		xProjection = std::round(1.0/guidY[0]);
	for(int cyInv=std::max(1, yProjection-errorWindow); cyInv<=yProjection+errorWindow; ++cyInv)
	{
		for(int cxInv=std::max(0, xProjection-errorWindow); cxInv<=xProjection+errorWindow; ++cxInv)
		{
			double cx = cxInv>1 ? 1.0/cxInv : 0, cy = 1.0/cyInv;
			projection = testCycles(centeredTexture, Eigen::Vector2d(1.0, 0), Eigen::Vector2d(cx, cy));
			DataType *dataProjection = reinterpret_cast<DataType *>(&projection);
			for(unsigned i=0; i<pixelSize; ++i)
			{
				projectionValue += dataProjection[i]*dataProjection[i];
			}
			if(maxProjectionValue < projectionValue)
			{
				maxProjectionValue = projectionValue;
				maxXProjection = cxInv > 1 ? cxInv : 0;
				maxYProjection = cyInv;
			}
		}
	}
	m_cycles[1][0] = maxXProjection > 0 ? (maxXProjection == 1 ? 0 : 1.0/maxXProjection) : 0; //1 is th
	m_cycles[1][1] = 1.0/maxYProjection;
}

template<typename I>
typename CSN_Texture<I>::PixelType CSN_Texture<I>::testCycles(const I& texture, const Eigen::Vector2d &tx, const Eigen::Vector2d &ty, bool useCyclicAverage) const
{
	assert(tx[0] != 0 && ty[1] != 0);
	assert(texture.is_initialized());
	unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);

	auto pixelTypeProduct = [pixelSize] (PixelType p1, PixelType p2) -> PixelType //please add a PixelType*PixelType
	{
		PixelType result;
		DataType *dataResult = reinterpret_cast<DataType *>(&result);

		DataType *dataP1 = reinterpret_cast<DataType *>(&p1);
		DataType *dataP2 = reinterpret_cast<DataType *>(&p2);
		for(unsigned i=0; i<pixelSize; ++i)
		{
			dataResult[i] = dataP1[i]*dataP2[i];
		}
		return result;
	};

	auto pixelTypeCompare = [pixelSize] (PixelType p1, PixelType p2) -> PixelType //please add a PixelType*PixelType
	{
		PixelType result;
		DataType *dataResult = reinterpret_cast<DataType *>(&result);

		DataType *dataP1 = reinterpret_cast<DataType *>(&p1);
		DataType *dataP2 = reinterpret_cast<DataType *>(&p2);
		for(unsigned i=0; i<pixelSize; ++i)
		{
			dataResult[i] = (dataP1[i]-dataP2[i]) * (dataP1[i]-dataP2[i]); //make sure there's a square of that somewhere
		}
		return result;
	};

	auto pixelTypeDivision = [pixelSize] (PixelType p1, PixelType p2) -> PixelType //please
	{
		PixelType result;
		DataType *dataResult = reinterpret_cast<DataType *>(&result);

		DataType *dataP1 = reinterpret_cast<DataType *>(&p1);
		DataType *dataP2 = reinterpret_cast<DataType *>(&p2);
		for(unsigned i=0; i<pixelSize; ++i)
		{
			dataResult[i] = dataP1[i]/dataP2[i];
		}
		return result;
	};

	auto sqrtPixelType = [pixelSize] (PixelType p1) -> PixelType //please add a sqrt(PixelType)
	{
		PixelType result;
		DataType *dataResult = reinterpret_cast<DataType *>(&result);

		DataType *dataP1 = reinterpret_cast<DataType *>(&p1);
		for(unsigned i=0; i<pixelSize; ++i)
		{
			dataResult[i] = sqrt(dataP1[i]);
		}
		return result;
	};

	ImageType projectionImage, averageProjectionImage;
	ImageType parallelogram1, parallelogram2;
	PixelType projection = ImageType::zero();
	itk::Size<2> projectionSize;
	projectionSize[0] = int(texture.width() * tx[0]);
	projectionSize[1] = int(texture.height() * ty[1]);
	//Rectangle size is (tx, ty)
	//Number of different rectangles is (txInv, tyInv)
	int txInv = int(std::round(1.0/tx[0]));
	int tyInv = int(std::round(1.0/ty[1]));
	unsigned nudge=0;
	print_debug("txInv (double): " << 1.0/tx[0]);
	print_debug("tyInv (double): " << 1.0/ty[1]);

	ImageType average;
	unsigned energy = 0;
	average.initItk(projectionSize[0], projectionSize[1], true);
	if(useCyclicAverage)
	{
		for(int i=0; i<txInv; ++i)
			for(int j=0; j<tyInv; ++j)
			{
				if(1.0/tx[0] != double(txInv) || 1.0/ty[0] != double(tyInv) || tx[1] != 0 || ty[0] != 0)
					average.for_all_pixels([&] (PixelType &pix, int x, int y)
					{
						PixelType ijPix = bilinear_interpolation(texture,	((nudge+i)*tx[0] + (nudge+j)*ty[0])*texture.width()+x,
																			((nudge+i)*tx[1] + (nudge+j)*ty[1])*texture.height()+y, true);
						pix += ijPix;
					});
				else //if it doesn't significantly decrease execution time, please remove this conditional expr. and keep only bilinear interpolation result.
					average.for_all_pixels([&] (PixelType &pix, int x, int y)
					{
						PixelType ijPix = texture.pixelAbsolute(	unsigned(((nudge+i)*tx[0] + (nudge+j)*ty[0])*texture.width()+x)%texture.width(),
																	unsigned(((nudge+i)*tx[1] + (nudge+j)*ty[1])*texture.height()+y)%texture.height());
						pix += ijPix;
					});
				++energy;
			}
		average.for_all_pixels([&] (PixelType &pix)
		{
			pix = pix * (1.0/energy);
		});
	}

	averageProjectionImage.initItk(projectionSize[0], projectionSize[1], true);
	projectionImage.initItk(projectionSize[0], projectionSize[1], true);
	parallelogram1.initItk(projectionSize[0], projectionSize[1], true);
	parallelogram2.initItk(projectionSize[0], projectionSize[1], true);
	energy = 0;

	for(int i=0; i<txInv; ++i)
	{
		for(int j=0; j<tyInv; ++j)
		{
			for(int k=0; k<txInv; ++k)
			{
				for(int l=0; l<tyInv; ++l)
				{
					if(i != k || j != l)
					{
						PixelType p1Norm = ImageType::zero();
						PixelType p2Norm = ImageType::zero();
						projectionImage.for_all_pixels([&] (PixelType &pix, int x, int y)
						{
							parallelogram1.pixelAbsolute(x, y)= bilinear_interpolation(texture, ((nudge+i)*tx[0] + (nudge+j)*ty[0])*texture.width()+x,
																								((nudge+i)*tx[1] + (nudge+j)*ty[1])*texture.height()+y, true);
							if(useCyclicAverage)
								parallelogram1.pixelAbsolute(x, y) -= average.pixelAbsolute(x, y);
							p1Norm += pixelTypeProduct(parallelogram1.pixelAbsolute(x, y), parallelogram1.pixelAbsolute(x, y));

							parallelogram2.pixelAbsolute(x, y) =bilinear_interpolation(texture, ((nudge+k)*tx[0] + (nudge+l)*ty[0])*texture.width()+x,
																								((nudge+k)*tx[1] + (nudge+l)*ty[1])*texture.height()+y, true);
							if(useCyclicAverage)
								parallelogram2.pixelAbsolute(x, y) -= average.pixelAbsolute(x, y);
							p2Norm += pixelTypeProduct(parallelogram2.pixelAbsolute(x, y), parallelogram2.pixelAbsolute(x, y));
						});

						//save_debug(parallelogram1, "/home/nlutz/parallelogram1.png");

						p1Norm = sqrtPixelType(p1Norm);
						p2Norm = sqrtPixelType(p2Norm);

						projectionImage.for_all_pixels([&] (PixelType &pix, int x, int y)
						{
//							parallelogram1.pixelAbsolute(x, y) = pixelTypeDivision(parallelogram1.pixelAbsolute(x, y), p1Norm);
//							parallelogram2.pixelAbsolute(x, y) = pixelTypeDivision(parallelogram2.pixelAbsolute(x, y), p2Norm);

							pix = pixelTypeCompare(parallelogram1.pixelAbsolute(x, y), parallelogram2.pixelAbsolute(x, y));

							averageProjectionImage.pixelAbsolute(x, y) += pix;
						});
						++energy;
					}
				}
			}
		}
	}
	averageProjectionImage.for_all_pixels([&] (PixelType &pix)
	{
		projection += pix * (1.0/energy);
	});
	ImageRGBd normalizedAverageProjectionImage;
	normalizedAverageProjectionImage.initItk(averageProjectionImage.width(), averageProjectionImage.height());
	normalizedAverageProjectionImage.copy_pixels(averageProjectionImage);
	ImageRGBd::PixelType minAverageProjection, maxAverageProjection;
	for(unsigned i=0; i<3; ++i)
	{
		minAverageProjection[i] = std::numeric_limits<double>::infinity();
		maxAverageProjection[i] = -std::numeric_limits<double>::infinity();
	}
//	normalizedAverageProjectionImage.for_all_pixels([&] (ImageRGBd::PixelType &pix) //option 1
//	{
//		for(unsigned i=0; i<3; ++i)
//		{
//			minAverageProjection[i] = std::min(minAverageProjection[i], pix[i]);
//			maxAverageProjection[i] = std::max(maxAverageProjection[i], pix[i]);
//		}
//	});
	for(unsigned i=0; i<3; ++i) //option 2
	{
		minAverageProjection[i] = -2.0/energy;
		maxAverageProjection[i] = 2.0/energy;
	}
	normalizedAverageProjectionImage.for_all_pixels([&] (ImageRGBd::PixelType &pix)
	{
		for(unsigned i=0; i<3; ++i)
		{
			pix[i] = (pix[i] - minAverageProjection[i]) / (maxAverageProjection[i] - minAverageProjection[i]);
		}
	});
	save_debug(averageProjectionImage, "/home/nlutz/averageProjectionImage.png");
	projection = projection * (1.0/(projectionSize[0] * projectionSize[1]));
	return projection;
}

template<typename I>
typename CSN_Texture<I>::PixelType CSN_Texture<I>::testCycles_v2(const I& texture, const Eigen::Vector2d &tx, const Eigen::Vector2d &ty, unsigned nbSamplesX, unsigned nbSamplesY) const
{
	assert(tx[0] >= 0.02 && ty[1] >= 0.02);
	assert(texture.is_initialized());
	assert(nbSamplesX>1 && nbSamplesY>1);
	unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	Eigen::Vector2d txScale;
	txScale[0] = tx[0]*(texture.width());
	txScale[1] = tx[1]*(texture.height());
	Eigen::Vector2d tyScale;
	tyScale[0] = ty[0]*(texture.width());
	tyScale[1] = ty[1]*(texture.height());

	auto pixelTypeProduct = [pixelSize] (PixelType p1, PixelType p2) -> PixelType //please add a PixelType*PixelType
	{
		PixelType result;
		DataType *dataResult = reinterpret_cast<DataType *>(&result);

		DataType *dataP1 = reinterpret_cast<DataType *>(&p1);
		DataType *dataP2 = reinterpret_cast<DataType *>(&p2);
		for(unsigned i=0; i<pixelSize; ++i)
		{
			dataResult[i] = dataP1[i]*dataP2[i];
		}
		return result;
	};

	auto pixelTypeCompare = [pixelSize] (PixelType p1, PixelType p2) -> PixelType //please add a PixelType*PixelType
	{
		PixelType result;
		DataType *dataResult = reinterpret_cast<DataType *>(&result);

		DataType *dataP1 = reinterpret_cast<DataType *>(&p1);
		DataType *dataP2 = reinterpret_cast<DataType *>(&p2);
		for(unsigned i=0; i<pixelSize; ++i)
		{
			dataResult[i] = (dataP1[i]-dataP2[i]) * (dataP1[i]-dataP2[i]); //make sure there's a square of that somewhere
		}
		return result;
	};

	auto pixelTypeSqrt = [pixelSize] (PixelType p1) -> PixelType //please add a sqrt(PixelType)
	{
		PixelType result;
		DataType *dataResult = reinterpret_cast<DataType *>(&result);

		DataType *dataP1 = reinterpret_cast<DataType *>(&p1);
		for(unsigned i=0; i<pixelSize; ++i)
		{
			dataResult[i] = sqrt(dataP1[i]);
		}
		return result;
	};
	//requires tx or ty to be aligned with the axis x (resp. y).
	//To make this function more robust, sample the tx, ty parallelogram and extract the intersection between the lattice described by (tx, ty) and the defined texture space.
	int txInv = int(std::round(1.0/tx[0]))-1;
	int tyInv = int(std::round(1.0/ty[1]))-1;
	PixelType finalValue = ImageType::zero();
	for(unsigned i=0; i<nbSamplesX; ++i)
	{
		double dx = txScale[0]*(double(i)/(nbSamplesX-1));
		for(unsigned j=0; j<nbSamplesY; ++j)
		{
			double dy = tyScale[1]*(double(j)/(nbSamplesY-1));
			unsigned hitCount = 0;
			PixelType hitValue = ImageType::zero();
			for(int k0=0; k0<txInv; ++k0)
			{
				for(int l0=0; l0<tyInv; ++l0)
				{
					for(int k1=0; k1<txInv; ++k1)
					{
						for(int l1=0; l1<tyInv; ++l1)
						{
							if(k0 != k1 || l0 != l1)
							{

								PixelType value = bilinear_interpolation(texture, dx+k0*txScale[0]+l0*tyScale[0], dy+k0*txScale[1]+l0*tyScale[1], true)
												- bilinear_interpolation(texture, dx+k1*txScale[0]+l1*tyScale[0], dy+k1*txScale[1]+l1*tyScale[1], true);
								PixelType debugValue = ImageType::zero();
								if(false)
								{
									print_debug("coordinates - x0: " << dx+k0*txScale[0]+l0*tyScale[0] << ", y0: " << dy+k0*txScale[1]+l0*tyScale[1]);
									print_debug("coordinates - x1: " << dx+k1*txScale[0]+l1*tyScale[0] << ", y1: " << dy+k1*txScale[1]+l1*tyScale[1]);
									print_debug("dx: " << dx << ", dy: " << dy);
									print_debug("k0: " << k0 << ", l0: " << l0);
									print_debug("k1: " << k1 << ", l1: " << l1);
									print_debug("value: " << value << std::endl);
								}
								hitValue += pixelTypeSqrt(pixelTypeProduct(value, value));
								++hitCount;
							}
						}
					}
				}
			}
			hitValue = hitValue * (1.0/hitCount);
			finalValue += hitValue;
		}
	}
	finalValue = finalValue * (1.0/(nbSamplesX*nbSamplesY));
	return finalValue;
}


template<typename I>
const Eigen::Vector2d &CSN_Texture<I>::cycleX() const
{
	return m_cycles[0];
}

template<typename I>
const Eigen::Vector2d &CSN_Texture<I>::cycleY() const
{
	return m_cycles[1];
}

template<typename I>
typename CSN_Texture<I>::PcaPixelType CSN_Texture<I>::proceduralTilingAndBlending (const PcaImageType &image,
																				   Eigen::Vector2d uv) const
{
	// Get triangle info
	float w1, w2, w3;
	Eigen::Vector2i vertex1, vertex2, vertex3;
	TriangleGrid(uv, w1, w2, w3, vertex1, vertex2, vertex3);
	// Assign random offset to each triangle vertex

	auto lmbd_hashFunction = [&](const Eigen::Vector2i &vec) -> Eigen::Vector2d
	{
		if(m_useCycles)
			return cyclicHash(vec.cast<double>());
		else
			return hash(vec.cast<double>());
	};

	Eigen::Vector2d uv1 = fract(uv + lmbd_hashFunction(vertex1));
	Eigen::Vector2d uv2 = fract(uv + lmbd_hashFunction(vertex2));
	Eigen::Vector2d uv3 = fract(uv + lmbd_hashFunction(vertex3));
	// Fetch input
	auto lmbd_Vector2PixelPos = [&] (Eigen::Vector2d v) -> PixelPosType
	{
		PixelPosType pos;
		pos[0] = v[0]*(image.width()); //or width-1? Just make sure v can't be 1
		pos[1] = v[1]*(image.height());
		return pos;
	};
	PcaPixelType I1 = image.pixelAbsolute(lmbd_Vector2PixelPos(uv1));
	PcaPixelType I2 = image.pixelAbsolute(lmbd_Vector2PixelPos(uv2));
	PcaPixelType I3 = image.pixelAbsolute(lmbd_Vector2PixelPos(uv3));
	// Linear blending
	PcaPixelType color = I1 * w1 + I2 * w2 + I3 * w3;
	for(unsigned i=0; i<3; ++i)
	{
		color[i] -= 0.5;
		color[i] /= std::sqrt(w1*w1 + w2*w2 + w3*w3);
		color[i] += 0.5;
	}
	return color;
}

template<typename I>
void CSN_Texture<I>::TriangleGrid (	Eigen::Vector2d uv, float &w1, float &w2, float &w3,
									Eigen::Vector2i &vertex1, Eigen::Vector2i &vertex2, Eigen::Vector2i &vertex3) const
{
	// Scaling of the input
	uv *= m_uvScale*3.464; // 2 * sqrt (3)
	//uv[1] *= 0.4*3.464;
	// Skew input space into simplex triangle grid
	Eigen::Matrix2d gridToSkewedGrid;
	gridToSkewedGrid << 1.0, 0.0, -0.57735027, 1.15470054;
	Eigen::Vector2d skewedCoord = gridToSkewedGrid * uv;
	Eigen::Vector2d temp2d = floor(skewedCoord);
	Eigen::Vector2i baseId = temp2d.cast<int>();
	temp2d = fract(skewedCoord);
	Eigen::Vector3d temp;
	temp[0] = temp2d[0];
	temp[1] = temp2d[1];
	temp[2] = 1.0 - temp[0] - temp[1] ;
	if (temp[2] > 0.0)
	{
        w1 = std::pow(float(temp[2]), m_gamma);
        w2 = std::pow(float(temp[1]), m_gamma);
        w3 = std::pow(float(temp[0]), m_gamma);
		vertex1 = baseId ;
		vertex2 = baseId + Eigen::Vector2i(0 , 1);
		vertex3 = baseId + Eigen::Vector2i(1 , 0);
	}
	else
	{
        w1 = std::pow(float(-temp[2]), m_gamma);
        w2 = std::pow(float(1.0 - temp[1]), m_gamma);
        w3 = std::pow(float(1.0 - temp[0]), m_gamma);
		vertex1 = baseId + Eigen::Vector2i(1, 1);
		vertex2 = baseId + Eigen::Vector2i(1, 0);
		vertex3 = baseId + Eigen::Vector2i(0, 1);
	}
	double sumW = w1 + w2 + w3;
	w1 = w1/sumW;
	w2 = w2/sumW;
	w3 = w3/sumW;
}

template<typename I>
Eigen::Vector2d CSN_Texture<I>::hash(const Eigen::Vector2d &p) const
{
	Eigen::Matrix2d hashMat;
	hashMat << 127.1, 269.5, 311.7, 183.3;
	Eigen::Vector2d q = hashMat * p;
	q[0] = sin(q[0]);
	q[1] = sin(q[1]);
	return fract ( q * 43758.5453 );
}

template<typename I>
Eigen::Vector2d CSN_Texture<I>::cyclicHash(const Eigen::Vector2d &p) const
{
	int randMax = m_largestCycleProduct;
	Eigen::Matrix2d hashMat;
	hashMat << 127.1, 269.5, 311.7, 183.3;
	Eigen::Vector2d q = hashMat * p;
	q[0] = sin(q[0]);
	q[1] = sin(q[1]);
	Eigen::Vector2d h = fract ( q * 43758.5453 );
	int cycle1 = int(h[0]*randMax);
	int cycle2 = int(h[1]*randMax);
	return double(cycle1)*m_cycles[0] + double(cycle2)*m_cycles[1];
}

template<typename I>
Eigen::Vector2d CSN_Texture<I>::floor(const Eigen::Vector2d &v) const
{
	Eigen::Vector2d w;
	w[0] = std::floor(v[0]);
	w[1] = std::floor(v[1]);
	return w;
}

template<typename I>
Eigen::Vector2d CSN_Texture<I>::fract(const Eigen::Vector2d &v) const
{
	Eigen::Vector2d w;
	w[0] = v[0]-std::floor(v[0]);
	w[1] = v[1]-std::floor(v[1]);
	return w;
}

}

}

#endif
