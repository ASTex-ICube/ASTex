#ifndef __DICTIONARY_PROCESSOR__H
#define __DICTIONARY_PROCESSOR__H

#include "dictionary.h"
#include "ASTex/fourier.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkScaleTransform.h"

static bool firstTimeIdiot=false;
#define CHECK_FOR_NAN(x) do{if(!firstTimeIdiot && std::isnan(x)) {std::cout << "File " << __FILE__ << ", line " << __LINE__ << ": found an unexpected nan float (you may want to restart the routine)" << std::endl; firstTimeIdiot=true;} } while(0)

namespace ASTex
{

template<typename I>
class DictionaryProcessor
{
public:

	std::string external_imageMatcher_path; //< the alien attribute

	using ImageType = I;
	using DictionaryType = Dictionary<I>;
	using AtomType = Atom<I>;
	using PixelType = typename ImageType::PixelType;
	using DataType = typename ImageType::DataType;
	using SizeType = itk::Size<2>;
	using IndexType = itk::Index<2>;

	DictionaryProcessor();
	DictionaryProcessor(const DictionaryType &dictionary);

	//step 1

	void setInput(const I &input);

	void setPatchSize(unsigned width, unsigned height);
	void setSparsity(unsigned sparsity);
	void setMaxAtomsNb(unsigned nbAtoms);
	void setPatchOffset(unsigned offsetX, unsigned offsetY);
	void setEnableHistogramMatching(bool b);
	void setEnableFourierMatching(bool b);


	unsigned nbAtoms() const;
	unsigned sparsity() const;
	SizeType patchSize() const;
	SizeType patchOffset() const;
	bool histogramMatchingEnabled() const;
	bool fourierMatchingEnabled() const;

	void setShowDebugMessages(bool b);
	void setIntermediateOutputPath(const std::string &path);

	//step 2

	template<typename Compare = std::less<PixelType>>
	void dictionaryLearning(unsigned nbIterations);

	//step 3

	const DictionaryType &dictionary() const					{return m_dictionary;}
	const std::vector<std::vector<PixelType>> &weights() const	{return m_weights;}

	ImageType reconstructInput() const;

	template<typename Compare = std::less<PixelType>>
	ImageType synthesize(unsigned width, unsigned height, unsigned nbIterations, const I* imageInitialize = nullptr);

	void save(const std::string &directory) const;
	bool load(const std::string &directory);

	//step 4

	void savePatches(const std::string &directory);
	void saveReconstructedPatches(const std::string &directory);

	bool saveVizualisableAtoms(const std::string &directory);

	template<typename Compare = std::less<PixelType>>
	I amplify(const I& lowResMap, unsigned scale);

private:

	template<typename Compare = std::less<PixelType>>
	void orthogonalMatchingPursuit(const DictionaryType &dictionary, size_t patchIndex,
								   const std::vector<std::pair<IndexType, I>>&patches,
								   std::vector<std::vector<PixelType>> &weights);
	void updateDictionaryApproximateKSVD(size_t atomId);

	void readInput();

	void readImage(const I &image, std::vector<std::pair<IndexType, I>> &patches) const;

	I reconstructImage(const DictionaryType &dictionary, const std::vector<std::vector<PixelType>> &weights, unsigned width, unsigned height) const;

	DictionaryType m_dictionary;
	std::vector<std::vector<PixelType>> m_weights; //lost on the dimensionality? m_weights[patch][atom]
	std::vector<std::pair<IndexType, I>> m_patches;
	I m_input;

	unsigned m_patchOffsetX;
	unsigned m_patchOffsetY;

	unsigned m_sparsity;
	unsigned m_nbAtoms;
	SizeType m_patchSize;

	bool m_showDebugMessages;
	std::string m_intermediateOutputPath;

	bool m_histogramMatchingEnabled;
	bool m_fourierMatchingEnabled;
};

template<typename I>
DictionaryProcessor<I>::DictionaryProcessor() :
	external_imageMatcher_path(),
	m_dictionary(),
	m_weights(),
	m_patches(),
	m_input(),
	m_patchOffsetX(1),
	m_patchOffsetY(1),
	m_sparsity(8),
	m_nbAtoms(32),
	m_patchSize({{5, 5}}),
	m_showDebugMessages(false),
	m_intermediateOutputPath(""),
	m_histogramMatchingEnabled(false),
	m_fourierMatchingEnabled(false)
{
}

template<typename I>
DictionaryProcessor<I>::DictionaryProcessor(const DictionaryType &dictionary) :
	external_imageMatcher_path(),
	m_dictionary(dictionary),
	m_weights(),
	m_patches(),
	m_input(),
	m_patchOffsetX(1),
	m_patchOffsetY(1),
	m_sparsity(8),
	m_nbAtoms(32),
	m_patchSize({{5, 5}}),
	m_showDebugMessages(false),
	m_intermediateOutputPath(""),
	m_histogramMatchingEnabled(false),
	m_fourierMatchingEnabled(false)
{}

template<typename I>
void DictionaryProcessor<I>::setInput(const I &input)
{
	m_input = input;
}

template<typename I>
void DictionaryProcessor<I>::setSparsity(unsigned sparsity)
{
	m_sparsity = sparsity;
}

template<typename I>
void DictionaryProcessor<I>::setPatchSize(unsigned width, unsigned height)
{
	m_dictionary.setAtomSize(width, height);
	m_patchSize[0] = width;
	m_patchSize[1] = height;
}

template<typename I>
void DictionaryProcessor<I>::setMaxAtomsNb(unsigned nbAtoms)
{
	m_nbAtoms = nbAtoms;
}

template<typename I>
void DictionaryProcessor<I>::setPatchOffset(unsigned offsetX, unsigned offsetY)
{
	m_patchOffsetX = offsetX;
	m_patchOffsetY = offsetY;
}

template<typename I>
void DictionaryProcessor<I>::setEnableHistogramMatching(bool b)
{
	m_histogramMatchingEnabled = b;
}

template<typename I>
void DictionaryProcessor<I>::setEnableFourierMatching(bool b)
{
	m_fourierMatchingEnabled = b;
}

template<typename I>
unsigned DictionaryProcessor<I>::nbAtoms() const
{
	return m_nbAtoms;
}

template<typename I>
unsigned DictionaryProcessor<I>::sparsity() const
{
	return m_sparsity;
}

template<typename I>
typename DictionaryProcessor<I>::SizeType DictionaryProcessor<I>::patchSize() const
{
	return m_patchSize;
}

template<typename I>
typename DictionaryProcessor<I>::SizeType DictionaryProcessor<I>::patchOffset() const
{
	itk::Size<2> offsets;
	offsets[0] = m_patchOffsetX;
	offsets[1] = m_patchOffsetY;
	return offsets;
}

template<typename I>
bool DictionaryProcessor<I>::histogramMatchingEnabled() const
{
	return m_histogramMatchingEnabled;
}

template<typename I>
bool DictionaryProcessor<I>::fourierMatchingEnabled() const
{
	return m_fourierMatchingEnabled;
}

template<typename I>
void DictionaryProcessor<I>::setShowDebugMessages(bool b)
{
	m_showDebugMessages = b;
}

template<typename I>
void DictionaryProcessor<I>::setIntermediateOutputPath(const std::string &path)
{
	m_intermediateOutputPath = path;
}

//step 2

//private but the two following functions go well together
template<typename I>
void DictionaryProcessor<I>::readImage(const I& image, std::vector<std::pair<IndexType, I> > &patches) const
{
	assert(image.is_initialized() && image.width()>=m_patchSize[0] && image.height()>=m_patchSize[1] && "Dictionnary::read: exemplar uninitialized or too small compared to atoms size");
	unsigned int previousPatchSize=patches.size();
	patches.clear();
	patches.reserve(previousPatchSize);
	//modify region dynamically. Store region data into atoms.
	Region reg;
	int x, y=0;
	auto makePatch = [&] (int x, int y)
	{
		reg = gen_region(x, y, m_patchSize[0], m_patchSize[1]);
		ImageType patch;
		patch.initItk(m_patchSize[0], m_patchSize[1]);
		image.for_region_pixels(reg, [&] (const typename I::PixelType& pix, int w, int h)
		{
			patch.pixelAbsolute(w-x, h-y)=pix;
		});
		IndexType origin = {{x, y}};
		patches.push_back(std::make_pair(origin, patch));
	};
	for(x=0; x+m_patchSize[0]<=unsigned(image.width()); x+=m_patchOffsetX)
		for(y=0; y+m_patchSize[1]<=unsigned(image.height()); y+=m_patchOffsetY)
			makePatch(x, y);
	unsigned int corner=0;
	if(x-m_patchOffsetX+m_patchSize[0]!=unsigned(image.width())) //check previous patch
	{
		++corner;
		x=image.width()-m_patchSize[0];
		for(y=0; y+m_patchSize[1]<=unsigned(image.height()); y+=m_patchOffsetY)
			makePatch(x, y);
	}
	if(y-m_patchOffsetY+m_patchSize[1]!=unsigned(image.height())) //check previous patch
	{
		++corner;
		y=image.height()-m_patchSize[1];
		for(x=0; x+m_patchSize[0]<=unsigned(image.width()); x+=m_patchOffsetX)
			makePatch(x, y);
	}
	if(corner==2)
	{
		x=image.width()-m_patchSize[0];
		y=image.height()-m_patchSize[1];
		makePatch(x, y);
	}
}

template <typename I>
void DictionaryProcessor<I>::readInput()
{

	if(m_showDebugMessages)
		std::cout << "Reading input... ";

	return readImage(m_input, m_patches);

	if(m_showDebugMessages)
		std::cout << "done!" << std::endl;
}

template <typename I>
I DictionaryProcessor<I>::reconstructImage(const DictionaryType &dictionary, const std::vector<std::vector<PixelType> > &weights, unsigned width, unsigned height) const
{
	static PixelType pix_zero;
	unsigned pixSize = sizeof(PixelType)/sizeof(DataType);
	I reconstructedImage;
	reconstructedImage.initItk(width, height, true);
	reconstructedImage.parallel_for_all_pixels([&] (PixelType &pix) {pix = pix_zero;}); //TODO: remove when initItk is fixed

	//first compute the hitmap which finds the number of time each texel will be written
	ASTex::ImageGrayu32 hitmap;
	hitmap.initItk(int(width), int(height), true);

	for(auto &patch : m_patches)
	{
		ASTex::Region reg = ASTex::gen_region(patch.first[0], patch.first[1], patch.second.width(), patch.second.height());
		hitmap.for_region_pixels(reg, [&] (ASTex::ImageGrayu32::PixelType &pix)
		{
			++pix;
		});
	}
	unsigned w=0;
	for(auto &patch : m_patches)
	{
		I reconstructedPatch = dictionary * weights[w++];
		reconstructedPatch.for_all_pixels([&] (const PixelType &pix, int x, int y)
		{
			reconstructedImage.pixelAbsolute(patch.first[0]+x, patch.first[1]+y) += pix;
		});
	}
	reconstructedImage.for_all_pixels([&] (PixelType &pix, int x, int y)
	{
		pix = pix*(1.0/hitmap.pixelAbsolute(x, y));
		for(unsigned k=0; k<pixSize; ++k)
		{
			CHECK_FOR_NAN(reinterpret_cast<DataType *>(&pix)[k]);
			reinterpret_cast<DataType *>(&pix)[k] = std::min(DataType(1.0), std::max(DataType(0.0), reinterpret_cast<DataType *>(&pix)[k]));
		}
	});
	return reconstructedImage;
}

template<typename I>
typename DictionaryProcessor<I>::ImageType DictionaryProcessor<I>::reconstructInput() const
{
	size_t pixelSize = sizeof(PixelType)/sizeof(DataType);
	if(m_showDebugMessages)
		std::cout << "Reconstructing the input from dictionary and found weights...";
	I reconstructedInput = reconstructImage(m_dictionary, m_weights, m_input.width(), m_input.height());
	if(m_showDebugMessages)
	{
		PixelType l1Dif{};
		I imDif = m_input - reconstructedInput;
		imDif.for_all_pixels([&] (const typename I::PixelType &pix)
		{
			for(unsigned p=0; p<pixelSize; ++p)
				reinterpret_cast<DataType*>(&l1Dif)[p] += reinterpret_cast<const DataType*>(&pix)[p] * reinterpret_cast<const DataType*>(&pix)[p];
		});
		for(unsigned p=0; p<pixelSize; ++p)
			reinterpret_cast<DataType*>(&l1Dif)[p] = std::sqrt(reinterpret_cast<DataType*>(&l1Dif)[p])/(imDif.width()*imDif.height());
		std::cout << "L2 difference between input and reconstructed input: " << l1Dif << std::endl;
		std::cout << "done!" << std::endl;
	}
	return reconstructedInput;
}

//step 3

template<typename I>
template<typename Compare>
void DictionaryProcessor<I>::dictionaryLearning(unsigned nbIterations)
{
	size_t pixelSize = sizeof(PixelType)/sizeof(DataType);
	readInput();
	assert(m_input.is_initialized()
		   && m_dictionary.atomWidth()>0 && m_dictionary.atomHeight()>0
		   && m_patches.size()>0);
	m_weights.clear();
	m_weights.resize(m_patches.size());
	m_dictionary.initRandom(nbAtoms()); //1. Initialization

	//Alternating OMP (to update weights) and KSVD (to update atoms and weights)
	unsigned currentProgress=0;
	if(m_showDebugMessages)
	{
		std::cout << "Starting dictionary learning";
		std::cout.flush();
	}
	for(unsigned i=0; i<nbIterations; ++i)
	{
		if(m_showDebugMessages && 100*i/nbIterations > currentProgress)
		{
			std::cout << '.';
			std::cout.flush();
			++currentProgress;
		}
		for(size_t w=0; w<m_patches.size(); ++w) //2. Sparse coding (update weights)
			orthogonalMatchingPursuit<Compare>(m_dictionary, w, m_patches, m_weights);

		if(m_intermediateOutputPath != "")
			IO::save01_in_u8(reconstructInput(), m_intermediateOutputPath + "/dp_learning_" + std::to_string(i) + ".png");

		for(size_t a=0; a<m_dictionary.nbAtoms(); ++a)
			updateDictionaryApproximateKSVD(a);
	}
	if(m_showDebugMessages)
	{
		std::cout << ". Cleaning useless atoms...";
		std::cout.flush();
	}

//	//Deleting useless atoms (otherwise they may interfere with the synthesis)
	for(size_t a=0; a<m_dictionary.nbAtoms(); ++a)
	{
		bool used=false;
		for(size_t w=0; w<m_weights.size() && !used; ++w)
		{
			PixelType weight = m_weights[w][a];
			for(unsigned p=0; p<pixelSize; ++p)
				if(reinterpret_cast<const DataType *>(&weight)[p] != 0)
					used=true;
		}
		if(!used)
		{
			m_dictionary.erase(a);
			for(size_t w=0; w<m_weights.size(); ++w)
			{
				m_weights[w].erase(m_weights[w].begin()+a);
			}
			--a;
		}
	}
	if(m_showDebugMessages)
		std::cout << " Done!" << std::endl;
}

template<typename I>
template<typename Compare>
typename DictionaryProcessor<I>::ImageType DictionaryProcessor<I>::amplify(const I& lowResMap, unsigned scale)
{
	float invScale = 1.0f/scale;
	if(m_patchSize[0]<=1 || m_patchSize[1]<=1)
	{
		std::cerr << "DictionaryProcessor::matchLowResMap: scale must be >0 & <1" << std::endl;
		exit(EXIT_FAILURE);
	}
	DictionaryType dictionaryLowRes;
	dictionaryLowRes.setAtomSize(	std::min(std::max(1, int(m_patchSize[0] * invScale)), int(m_patchSize[0]-1)),
									std::min(std::max(1, int(m_patchSize[1] * invScale)), int(m_patchSize[1]-1)));
	dictionaryLowRes.initEmpty(m_dictionary.nbAtoms());
	SizeType atomSize;
	atomSize[0]=dictionaryLowRes.atomWidth();
	atomSize[1]=dictionaryLowRes.atomHeight();

	for(unsigned i=0; i<dictionaryLowRes.nbAtoms(); ++i)
	{
		using ImageType = typename I::ItkImg;

		typename ImageType::RegionType region = m_dictionary.atom(i).content().itk()->GetLargestPossibleRegion();
		typename ImageType::SizeType size = region.GetSize();
		typename ImageType::SpacingType spacing = m_dictionary.atom(i).content().itk()->GetSpacing();
		spacing[0] = spacing[0] * scale;
		spacing[1] = spacing[1] * scale;

		IndexType centralPixel;
		centralPixel[0] = size[0]/2;
		centralPixel[1] = size[0]/2;
		itk::Point< DataType, 2 > centralPoint;
		centralPoint[0] = centralPixel[0];
		centralPoint[1] = centralPixel[1];

		using ScaleTransformType = itk::ScaleTransform< typename I::DataType, 2 >;
		typename ScaleTransformType::Pointer scaleTransform = ScaleTransformType::New();

		typename ScaleTransformType::ParametersType parameters = scaleTransform->GetParameters();
		parameters[0] = 1;
		parameters[1] = 1;

		scaleTransform->SetParameters( parameters );
		scaleTransform->SetCenter( centralPoint );

		using LinearInterpolatorType = itk::LinearInterpolateImageFunction< ImageType, typename I::DataType>;
		typename LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();

		using ResampleFilterType = itk::ResampleImageFilter< ImageType, ImageType >;
		typename ResampleFilterType::Pointer rif = ResampleFilterType::New();

		rif->SetInput(m_dictionary.atom(i).content().itk());
		rif->SetOutputParametersFromImage(m_dictionary.atom(i).content().itk());
		rif->SetSize(atomSize);
		rif->SetTransform(scaleTransform);
		rif->SetInterpolator(interpolator);
		rif->SetOutputSpacing(spacing);
		rif->Update();
		I filterOutput(rif->GetOutput());
		dictionaryLowRes.atom(i).content() = filterOutput;
	}
	DictionaryProcessor dpLowRes(dictionaryLowRes);
	std::vector<std::pair<IndexType, I>> patches;
	std::vector<std::vector<PixelType>> weights;
	SizeType savePatchSize = m_patchSize;
	m_patchSize[0] = m_patchSize[0] * invScale;
	m_patchSize[1] = m_patchSize[1] * invScale;
	readImage(lowResMap, patches);
	weights.resize(patches.size());
	m_patches.resize(patches.size());
	for(size_t w=0; w<patches.size(); ++w) //2. Sparse coding (update weights)
	{
		orthogonalMatchingPursuit<Compare>(dictionaryLowRes, w, patches, weights);
		m_patches[w].first[0] = patches[w].first[0];
		m_patches[w].first[1] = patches[w].first[1];
		m_patches[w].second = patches[w].second;
	}
	I reconstructedLowResMap = reconstructImage(dictionaryLowRes, weights, lowResMap.width(), lowResMap.height());
	for(size_t w=0; w<patches.size(); ++w)
	{
		m_patches[w].first[0] = patches[w].first[0] * scale;
		m_patches[w].first[1] = patches[w].first[1] * scale;
		m_patches[w].second.initItk(savePatchSize[0], savePatchSize[1]);
	}
	m_patchSize = savePatchSize;
	I output=reconstructImage(m_dictionary, weights, lowResMap.width()*scale, lowResMap.height()*scale);
	return output;
}

//Private

template<typename I>
template<typename Compare>
void DictionaryProcessor<I>::orthogonalMatchingPursuit(const DictionaryType &dictionary, size_t patchIndex,
													   const std::vector<std::pair<IndexType, I>>&patches,
													   std::vector<std::vector<PixelType>> &weights)
{
	static PixelType s_zero;
	unsigned s = sparsity();
	weights[patchIndex].clear();
	weights[patchIndex].resize(dictionary.nbAtoms());
	for(unsigned i=0; i<dictionary.nbAtoms(); ++i)
		weights[patchIndex][i] = s_zero;
	const I& patch = patches[patchIndex].second;
	I residual;
	residual.initItk(patch.width(), patch.height());
	residual.copy_pixels(patch);
	assert(unsigned(residual.width()) == dictionary.atomWidth() && unsigned(residual.height()) == dictionary.atomHeight());
	unsigned i = 0;
	Compare comparatorLeast;
	size_t pixelSize = sizeof(PixelType)/sizeof(DataType);
	typename I::DataType *a_pix = new typename I::DataType[pixelSize];

	double residualEvaluation = 0;
	residual.for_all_pixels([&] (PixelType &pix)
	{
		for(unsigned k=0; k<pixelSize; ++k)
		{
			residualEvaluation += reinterpret_cast<DataType *>(&pix)[k] * reinterpret_cast<DataType *>(&pix)[k];
			CHECK_FOR_NAN(reinterpret_cast<DataType *>(&pix)[k]);
		}
	});

	while(i++ < s && residualEvaluation > std::numeric_limits<double>::epsilon()) //customizable stop condition (sparsity constraint)
	{
		residualEvaluation = 0;
		PixelType bestCorrelation = s_zero;
		PixelType correlation = s_zero;
		PixelType bestDotProduct = s_zero;
		PixelType dotProduct = s_zero;

		const DataType *operand1, *operand2;

		unsigned bestAtom=0;
		for(unsigned a=0; a<dictionary.nbAtoms(); ++a)
		{
			dotProduct = s_zero;
			const AtomType &atom = dictionary.atom(a);
			residual.for_all_pixels([&] (const PixelType &pix, int x, int y)
			{
				operand1 = reinterpret_cast<const DataType *>(&pix);
				PixelType atomAtxy = atom.content().pixelAbsolute(x, y);
				operand2 = reinterpret_cast<const DataType *>(&atomAtxy);
				for(unsigned k=0; k<pixelSize; ++k)
				{
					reinterpret_cast<DataType *>(&dotProduct)[k] += operand1[k] * operand2[k];
					CHECK_FOR_NAN(reinterpret_cast<DataType *>(&dotProduct)[k]);
				}
			});

			std::memcpy(a_pix, &dotProduct, sizeof(PixelType)); //absolution of dotproduct
			for(unsigned k=0; k<pixelSize; ++k)
			{
				a_pix[k] = std::abs(a_pix[k]);
				CHECK_FOR_NAN(a_pix[k]);
			}
			std::memcpy(&correlation, a_pix, sizeof(PixelType));

			//correlation = dotProduct;
			if(comparatorLeast(bestCorrelation, correlation)) //meaning that bestCorrelation < correlation
			{
				bestCorrelation = correlation;
				bestAtom = a;
				bestDotProduct = dotProduct;
			}
		} //best k found: represents the index of the atom that has the biggest absolute dot product with the residual
		const AtomType &atom = dictionary.atom(bestAtom);
		weights[patchIndex][bestAtom] = weights[patchIndex][bestAtom] + bestDotProduct;
		residual.for_all_pixels([&] (PixelType &pix, int x, int y)
		{
			operand1 = reinterpret_cast<const DataType *>(&bestDotProduct);
			PixelType atomAtxy = atom.content().pixelAbsolute(x, y);
			operand2 = reinterpret_cast<const DataType *>(&atomAtxy);
			for(unsigned k=0; k<pixelSize; ++k)
			{
				reinterpret_cast<DataType *>(&pix)[k] -= operand1[k] * operand2[k];
				residualEvaluation += reinterpret_cast<DataType *>(&pix)[k] * reinterpret_cast<DataType *>(&pix)[k];
				CHECK_FOR_NAN(reinterpret_cast<DataType *>(&pix)[k]);
			}
		});
	}
	delete[] a_pix;
}

template<typename I>
void DictionaryProcessor<I>::updateDictionaryApproximateKSVD(size_t atomId)
{
	//based on https://www.caam.rice.edu/~optimization/L1/optseminar/K-SVD_talk_lijun.pdf
	static PixelType pix_zero; //TODO: remove when initItk is fixed
	size_t pixelSize = sizeof(PixelType)/sizeof(DataType);
	//check if this atom is actually used
	bool used=false;
	for(unsigned w=0; w<m_weights.size() && !used; ++w)
	{
		PixelType weight = m_weights[w][atomId];
		for(unsigned p=0; p<pixelSize; ++p)
			if(reinterpret_cast<const DataType *>(&weight)[p] != 0)
				used=true;
	}
	if(!used)
		return; //sorry, dear reader
	//update atom
	AtomType &atom = m_dictionary.atom(atomId);
	atom.reset(); //Step 7: d_j := 0
	DictionaryType resultDictionaryByXi;
	std::vector<PixelType> G;
	G.reserve(sparsity());
	I resultPatchesByG, resultDictionaryXiByG;
	resultPatchesByG.initItk(m_patchSize[0], m_patchSize[1]);
	resultDictionaryXiByG.initItk(m_patchSize[0], m_patchSize[1]);
	resultPatchesByG.for_all_pixels([&] (PixelType &pix) {pix=pix_zero;});

	for(unsigned w=0; w<m_patches.size(); ++w)
	{
		PixelType weight = m_weights[w][atomId];
		if(!I::is_zero(weight)) //Step 8: I := {indices of the signals in Y (patches) whose representations use d_j}
		{
			G.push_back(weight); //Step 9: g := X^T_{j,I}

			const I& patch = m_patches[w].second;
			resultPatchesByG += patch * weight;
		}
	}
	//IO::save01_in_u8(resultPatchesByG, TEMPO_PATH+"debugResultYg.png");

	resultDictionaryByXi.setAtomSize(m_patchSize[0], m_patchSize[1]);
	resultDictionaryByXi.initEmpty(G.size());
	unsigned a2=0;
	for(unsigned w=0; w<m_weights.size(); ++w)
	{
		PixelType weight = m_weights[w][atomId];
		if(!I::is_zero(weight))
		{
			for(unsigned a=0; a<m_dictionary.nbAtoms(); ++a)
			{
				resultDictionaryByXi.atom(a2).content() += m_dictionary.atom(a).content() * m_weights[w][a];
			}
			++a2;
		}
	}
	resultDictionaryXiByG = resultDictionaryByXi * G;

	AtomType d;
	d.content() = resultPatchesByG - resultDictionaryXiByG; //Step 10: d := Y_Ig - DX_Ig
	// Histogram<I>::saveImageToCsv(resultPatchesByG, TEMPO_PATH+"resultPatchedByG.csv");
	d.forceUpdateStatistics();
	PixelType norm = d.norm();

	d.content().for_all_pixels([&] (PixelType &pix)
	{
		for(unsigned k=0; k<pixelSize; ++k)
		{
			DataType norm_component = reinterpret_cast<const DataType *>(&norm)[k];
			if(norm_component<std::numeric_limits<float>::epsilon())
			{
				std::cerr << "Dictionary::updateDictionaryApproximateKSVD: the norm of the atom d is 0" << std::endl;
				exit(EXIT_FAILURE);
			}
			if(norm_component > std::numeric_limits<float>::epsilon())
				reinterpret_cast<DataType *>(&pix)[k] /= norm_component; //Step 11: d := d/||d||_2
			CHECK_FOR_NAN(reinterpret_cast<DataType *>(&pix)[k]);
		}
	});

	a2=0;
	for(unsigned w=0; w<m_patches.size(); ++w)
	{
		PixelType &weight = m_weights[w][atomId];
		if(!I::is_zero(weight))
		{
			weight = pix_zero;
			d.content().for_all_pixels([&] (const PixelType &pix, int x, int y)
			{
				for(unsigned k=0; k<pixelSize; ++k)
				{
					reinterpret_cast<DataType *>(&weight)[k] +=
							reinterpret_cast<const DataType *>(&m_patches[w].second.pixelAbsolute(x, y))[k]
						*	reinterpret_cast<const DataType *>(&pix)[k]
					-		reinterpret_cast<const DataType *>(&resultDictionaryByXi.atom(a2).content().pixelAbsolute(x, y))[k]
						*	reinterpret_cast<const DataType *>(&pix)[k];
					CHECK_FOR_NAN(reinterpret_cast<DataType *>(&weight)[k]);
				}
			});
			++a2;
			//Step 12/14 : X_{j, I} := (Y^T_I d - (DX_I)^T d)^T

		}
	}
	atom.content() = d.content(); //Step 13: d_j := d
}

template<typename I>
template<typename Compare>
typename DictionaryProcessor<I>::ImageType DictionaryProcessor<I>::synthesize(unsigned width, unsigned height, unsigned nbIterations, const ImageType *imageInitializer)
{
	I output;
	std::vector<std::vector<PixelType>> weights;
	std::vector<std::pair<itk::Index<2>, I>> patches;
	const unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	output.initItk(width, height);
	if(imageInitializer!=nullptr)
	{
		ImageType scaledInitializer;
		scaledInitializer.initItk(width, height, true);
		scaledInitializer.for_all_pixels([&] (PixelType &pix, int x, int y)
		{
			double xd, yd;
			xd = double(x)/(scaledInitializer.width()-1) * (scaledInitializer.width()-1);
			yd = double(y)/(scaledInitializer.height()-1) * (scaledInitializer.height()-1);
			pix = bilinear_interpolation(*imageInitializer, xd, yd, false);
		});
		output.copy_pixels(*imageInitializer);
	}
	else //white noise
	{
		output.for_all_pixels([&] (PixelType &pix)
		{
			for(unsigned k=0; k<pixelSize; ++k)
			{
				reinterpret_cast<DataType *>(&pix)[k] = DataType(std::rand())/RAND_MAX;
				CHECK_FOR_NAN(reinterpret_cast<DataType *>(&pix)[k]);
			}
		});
	}
	readImage(output, patches);
	weights.resize(patches.size());
	unsigned currentProgress=0;
	if(m_showDebugMessages)
	{
		std::cout << "Starting synthesis";
		std::cout.flush();
	}
	for(unsigned i=0; i<nbIterations; ++i)
	{
		if(m_showDebugMessages && currentProgress > 100*i/nbIterations)
		{
			++currentProgress;
			std::cout << '.';
			std::cout.flush();
		}
		readImage(output, patches);

		for(size_t w=0; w<m_patches.size(); ++w) //2. Sparse coding (update weights)
		{
			orthogonalMatchingPursuit<Compare>(m_dictionary, w, patches, weights);
		}
		output=reconstructImage(m_dictionary, weights, output.width(), output.height());
		if(i!=i) //TODO : change this idiot
			matchImage(output, m_input, external_imageMatcher_path);
		if(i!=i)
		{
			ImageGrayd r_source, g_source, b_source, r_target, g_target, b_target;
			extract3Channels(output, r_source, g_source, b_source);
			extract3Channels(m_input, r_target, g_target, b_target);

			ImageSpectrald modulus_source, modulus_target, phase_source;
			modulus_source.initItk(output.width(), output.height());
			modulus_target.initItk(m_input.width(), m_input.height());
			phase_source.initItk(output.width(), output.height());

			Fourier::fftForwardModulusAndPhase(r_target, modulus_target, phase_source);
			Fourier::fftForwardModulusAndPhase(r_source, modulus_source, phase_source);
			Fourier::fftInverseModulusAndPhase(modulus_target, phase_source, r_source);

			Fourier::fftForwardModulusAndPhase(g_target, modulus_target, phase_source);
			Fourier::fftForwardModulusAndPhase(g_source, modulus_source, phase_source);
			Fourier::fftInverseModulusAndPhase(modulus_target, phase_source, g_source);

			Fourier::fftForwardModulusAndPhase(b_target, modulus_target, phase_source);
			Fourier::fftForwardModulusAndPhase(b_source, modulus_source, phase_source);
			Fourier::fftInverseModulusAndPhase(modulus_target, phase_source, b_source);

			fold3Channels(output, r_source, g_source, b_source);
		}
		if(m_intermediateOutputPath != "")
			IO::save01_in_u8(output, m_intermediateOutputPath + "/dp_synthesis_" + std::to_string(i) + ".png");
	}
	if(m_showDebugMessages)
		std::cout << " Done!" << std::endl;
	return output;
}

template<typename I>
void DictionaryProcessor<I>::save(const std::string &directory) const
{
	assert(m_weights.size() > 0);
	create_directory(directory);
	Histogram<I>::saveImageToCsv(m_input, directory + "/input.csv");

	std::ofstream ofs_data_out(directory + "/data.csv");
	ofs_data_out << m_patchOffsetX << std::endl;
	ofs_data_out << m_patchOffsetY << std::endl;
	ofs_data_out << m_sparsity << std::endl;
	ofs_data_out << m_patchSize[0] << ' ' << m_patchSize[1] << std::endl;
	ofs_data_out.close();

	m_dictionary.save(directory);

	std::ofstream ofs_weights_out(directory + "/weights.csv");
	ofs_weights_out << m_weights.size() << ' ' << (*m_weights.begin()).size() << std::endl;
	for(typename std::vector<std::vector<PixelType>>::const_iterator it=m_weights.begin(); it!=m_weights.end(); ++it)
		for(typename std::vector<PixelType>::const_iterator it2=(*it).begin(); it2!=(*it).end(); ++it2)
			ofs_weights_out << (*it2) << std::endl;
}

template<typename I>
bool DictionaryProcessor<I>::load(const std::string &directory)
{
	const unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);


	if(!Histogram<I>::loadImageFromCsv(m_input, directory + "/input.csv"))
		return false;

	std::ifstream ifs_data_in(directory + "/data.csv");
	if(ifs_data_in.is_open())
	{
		ifs_data_in >> m_patchOffsetX;
		ifs_data_in >> m_patchOffsetY;
		ifs_data_in >> m_sparsity;
		ifs_data_in >> m_patchSize[0] >> m_patchSize[1];
		ifs_data_in.close();
	}
	else
		return false;
	readInput();

	m_dictionary.load(directory);

	std::ifstream ifs_weights_in(directory + "/weights.csv");
	unsigned size1, size2;
	ifs_weights_in >> size1;
	ifs_weights_in >> size2;

	m_weights.resize(size1);
	for(typename std::vector<std::vector<PixelType>>::iterator it=m_weights.begin(); it!=m_weights.end(); ++it)
	{
		(*it).resize(size2);
		for(typename std::vector<PixelType>::iterator it2=(*it).begin(); it2!=(*it).end(); ++it2)
		{
			DataType tmpPixel[pixelSize];
			for(unsigned k=0; k<pixelSize; ++k)
				ifs_weights_in >> tmpPixel[k];
			memcpy(&(*it2), &tmpPixel, sizeof(PixelType));
		}
	}

	m_nbAtoms = m_dictionary.nbAtoms();
	return true;
}

template<typename I>
void DictionaryProcessor<I>::savePatches(const std::string &directory)
{
	for(unsigned w=0; w<m_patches.size(); ++w)
	{
		IO::save01_in_u8(m_patches[w].second,
						 directory + "/patch_"	+ std::to_string(m_patches[w].first[0]) + "_"
												+ std::to_string(m_patches[w].first[1]) + ".png");
	}
}

template<typename I>
void DictionaryProcessor<I>::saveReconstructedPatches(const std::string &directory)
{
	for(unsigned w=0; w<m_patches.size(); ++w)
	{
		I reconstructedPatch;
		reconstructedPatch = m_dictionary * m_weights[w];
		IO::save01_in_u8(reconstructedPatch,
						 directory + "/patch_reconstructed_"	+ std::to_string(m_patches[w].first[0]) + "_"
																+ std::to_string(m_patches[w].first[1]) + ".png");
	}
}

template<typename I>
bool DictionaryProcessor<I>::saveVizualisableAtoms(const std::string &directory)
{
	const unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	std::ofstream vizualisableWeightListOut(directory + "/viualisableWeightList.csv", std::ios_base::trunc);
	if(!vizualisableWeightListOut)
		return false;
	for(unsigned a=0; a<m_dictionary.nbAtoms(); ++a)
	{
		const AtomType &atom = m_dictionary.atom(a);
		typename I::PixelType pixMax = atom.content().pixelAbsolute(0, 0);
		typename I::PixelType pixMin = atom.content().pixelAbsolute(0, 0);
		atom.content().for_all_pixels([&] (const typename I::PixelType &pix)
		{
			for(unsigned k=0; k<pixelSize; ++k)
			{
				reinterpret_cast<DataType*>(&(pixMax))[k] = std::max(reinterpret_cast<const DataType*>(&(pix))[k],
																	 reinterpret_cast<DataType*>(&(pixMax))[k]);
				reinterpret_cast<DataType*>(&(pixMin))[k] = std::min(reinterpret_cast<const DataType*>(&(pix))[k],
																	 reinterpret_cast<DataType*>(&(pixMin))[k]);
			}
		});
		vizualisableWeightListOut << a << std::endl;
		vizualisableWeightListOut << pixMin << std::endl;
		vizualisableWeightListOut << pixMax << std::endl << std::endl;
		I vizualisableAtom;
		vizualisableAtom.initItk(atom.content().width(), atom.content().height());
		atom.content().for_all_pixels([&] (const typename I::PixelType &pix, int x, int y)
		{
			typename I::PixelType &vpix = vizualisableAtom.pixelAbsolute(x, y);
			for(unsigned k=0; k<pixelSize; ++k)
			{
				reinterpret_cast<DataType*>(&(vpix))[k] =
						(reinterpret_cast<const DataType*>(&pix)[k] - reinterpret_cast<DataType*>(&(pixMin))[k])
						/ (reinterpret_cast<DataType*>(&(pixMax))[k] - reinterpret_cast<DataType*>(&(pixMin))[k]);
			}
		});
		IO::save01_in_u8(vizualisableAtom, directory + "/atom_vizualisable_" + std::to_string(a) + ".png");
	}
	return true;
}

}

#endif
