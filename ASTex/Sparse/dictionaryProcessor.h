#ifndef __DICTIONARY_PROCESSOR__H
#define __DICTIONARY_PROCESSOR__H

#include "dictionary.h"

namespace ASTex
{

template<typename I>
class DictionaryProcessor
{
public:

	using DictionaryType = Dictionary<I>;
	using AtomType = Atom<I>;
	using PixelType = typename I::PixelType;
	using DataType = typename I::DataType;

	DictionaryProcessor();

	//step 1

	void setInput(const I &input);

	void setPatchSize(unsigned width, unsigned height);
	void setSparsity(unsigned sparsity);
	void setNbAtoms(unsigned nbAtoms);
	void setPatchOffset(unsigned offsetX, unsigned offsetY);


	unsigned nbAtoms() const;
	unsigned sparsity() const;
	itk::Size<2> patchSize() const;
	itk::Size<2> patchOffset() const;

	//step 2

	template<typename Compare = std::less<PixelType>>
	void dictionaryLearning(unsigned nbIterations);

	//step 3

	const DictionaryType &dictionary() const					{return m_dictionary;}
	const std::vector<std::vector<PixelType>> &weights() const	{return m_weights;}

	I reconstructInput() const;

	template<typename Compare = std::less<PixelType>>
	I synthesize(unsigned width, unsigned height, unsigned nbIterations);

	void save(const std::string &directory) const;
	void load(const std::string &directory);

private:

	template<typename Compare = std::less<PixelType>>
	void orthogonalMatchingPursuit(size_t patchIndex, const std::vector<std::pair<itk::Index<2>, I>>&patches,
								   std::vector<std::vector<PixelType>> &weights);
	void updateDictionaryApproximateKSVD(size_t atomId);

	void readInput();

	void readImage(I& image, std::vector<std::pair<itk::Index<2>, I>> &patches);

	I reconstructImage(const std::vector<std::vector<PixelType>> &weights, unsigned width, unsigned height);

	DictionaryType m_dictionary;
	std::vector<std::vector<PixelType>> m_weights;
	std::vector<std::pair<itk::Index<2>, I>> m_patches;
	I m_input;

	unsigned m_patchOffsetX;
	unsigned m_patchOffsetY;

	unsigned m_sparsity;
	unsigned m_nbAtoms;
	itk::Size<2> m_patchSize;
};

template<typename I>
DictionaryProcessor<I>::DictionaryProcessor() :
	m_dictionary(), m_weights(), m_patches(), m_input(),
	m_patchOffsetX(1), m_patchOffsetY(1), m_sparsity(8), m_nbAtoms(32), m_patchSize({{5, 5}})
{
}

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
void DictionaryProcessor<I>::setNbAtoms(unsigned nbAtoms)
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
itk::Size<2> DictionaryProcessor<I>::patchSize() const
{
	return m_patchSize;
}

template<typename I>
itk::Size<2> DictionaryProcessor<I>::patchOffset() const
{
	itk::Size<2> offsets;
	offsets[0] = m_patchOffsetX;
	offsets[1] = m_patchOffsetY;
	return offsets;
}

//step 2

//private but the two following functions go well together
template<typename I>
void DictionaryProcessor<I>::readImage(I& image, std::vector<std::pair<itk::Index<2>, I>> &patches)
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
		I patch;
		patch.initItk(m_patchSize[0], m_patchSize[1]);
		image.for_region_pixels(reg, [&] (typename I::PixelType& p, int w, int h)
		{
			patch.pixelAbsolute(w-x, h-y)=p;
		});
		itk::Index<2> origin = {{x, y}};
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
	return readImage(m_input, m_patches);
}

template <typename I>
I DictionaryProcessor<I>::reconstructImage(const std::vector<std::vector<PixelType> > &weights, unsigned width, unsigned height)
{
	static PixelType pix_zero;
	unsigned pixSize = sizeof(PixelType)/sizeof(DataType);
	I reconstructedImage;
	reconstructedImage.initItk(width, height, true);
	reconstructedImage.parallel_for_all_pixels([&] (PixelType &pix) {pix = pix_zero;}); //TODO: remove when initItk is fixed

	//first compute the hitmap which finds the number of time each texel will be written
	ASTex::ImageGrayu32 hitmap;
	hitmap.initItk(width, height, true);

	for(auto &patch : m_patches)
	{
		ASTex::Region reg = ASTex::gen_region(patch.first[0], patch.first[1], patch.second.width(), patch.second.height());
		hitmap.for_region_pixels(reg, [&] (ASTex::ImageGrayu32::PixelType &pix)
		{
			++pix;
		});
	}
	unsigned i=0;
	for(auto &patch : m_patches)
	{
		I reconstructedPatch = m_dictionary * weights[i++];
		reconstructedPatch.for_all_pixels([&] (const typename I::PixelType &pix, int x, int y)
		{
			reconstructedImage.pixelAbsolute(patch.first[0]+x, patch.first[1]+y) += pix;
		});
	}
	reconstructedImage.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
	{
		pix = pix*(1.0/hitmap.pixelAbsolute(x, y));
		for(unsigned k=0; k<pixSize; ++k)
		{
			reinterpret_cast<DataType *>(&pix)[k] = std::min(DataType(1.0), std::max(DataType(0.0), reinterpret_cast<DataType *>(&pix)[k]));
		}
	});
	return reconstructedImage;
}

template<typename I>
I DictionaryProcessor<I>::reconstructInput() const
{
	I reconstructedInput = reconstructImage(m_weights, m_input.width(), m_input.height());
	PixelType l1Dif{};
	I imDif = m_input - reconstructedInput;
	imDif.for_all_pixels([&] (const typename I::PixelType &pix)
	{
		l1Dif += pix;
	});
	for(unsigned j=0; j<m_weights[0].size(); ++j)
	{
		std::cout << "weight: " << m_weights[0][j] << std::endl;
	}
	std::cout << "L1 difference between input and reconstructed input: " << l1Dif << std::endl;
	return reconstructedInput;
}

//step 3

template<typename I>
template<typename Compare>
void DictionaryProcessor<I>::dictionaryLearning(unsigned nbIterations)
{
	readInput();
	assert(m_input.is_initialized()
		   && m_dictionary.atomWidth()>0 && m_dictionary.atomHeight()>0
		   && m_patches.size()>0);
	m_weights.resize(m_patches.size());
	m_dictionary.initRandom(nbAtoms()); //1. Initialization
	for(unsigned i=0; i<nbIterations; ++i)
	{
		for(size_t p=0; p<m_patches.size(); ++p) //2. Sparse coding (update weights)
			orthogonalMatchingPursuit<Compare>(p, m_patches, m_weights);

		for(size_t a=0; a<m_dictionary.nbAtoms(); ++a)
			updateDictionaryApproximateKSVD(a); //3. Dictionary update
	}
}

template<typename I>
template<typename Compare>
void DictionaryProcessor<I>::orthogonalMatchingPursuit(size_t patchIndex,
													   const std::vector<std::pair<itk::Index<2>, I>>&patches,
													   std::vector<std::vector<PixelType>> &weights)
{
	static PixelType s_zero;
	unsigned s = sparsity();
	weights[patchIndex].clear();
	weights[patchIndex].resize(m_dictionary.nbAtoms());
	for(unsigned i=0; i<m_dictionary.nbAtoms(); ++i)
		weights[patchIndex][i] = s_zero;
	const I& patch = patches[patchIndex].second;
	I residual;
	residual.initItk(patch.width(), patch.height());
	residual.copy_pixels(patch);
	assert(unsigned(residual.width()) == m_dictionary.atomWidth() && unsigned(residual.height()) == m_dictionary.atomWidth());
	unsigned i = 0;
	Compare comparatorLeast;
	size_t pixelSize = sizeof(PixelType)/sizeof(DataType);
	typename I::DataType *a_pix = new typename I::DataType[pixelSize];

	while(i++ < s) //customizable stop condition (sparsity constraint)
	{
		PixelType bestCorrelation = s_zero;
		PixelType correlation = s_zero;
		PixelType bestDotProduct = s_zero;
		PixelType dotProduct = s_zero;

		const DataType *operand1, *operand2;

		unsigned bestK=0;
		for(unsigned k=0; k<m_dictionary.nbAtoms(); ++k)
		{
			dotProduct = s_zero;
			const AtomType &atom = m_dictionary.atom(k);
			residual.for_all_pixels([&] (const PixelType &pix, int x, int y)
			{
				operand1 = reinterpret_cast<const DataType *>(&pix);
				PixelType atomAtxy = atom.content().pixelAbsolute(x, y);
				operand2 = reinterpret_cast<const DataType *>(&atomAtxy);
				for(unsigned j=0; j<pixelSize; ++j)
				{
					reinterpret_cast<DataType *>(&dotProduct)[j] += operand1[j] * operand2[j];
				}
			});

			std::memcpy(a_pix, &dotProduct, sizeof(PixelType)); //absolution of dotproduct
			for(unsigned i=0; i<pixelSize; ++i)
				a_pix[i] = std::abs(a_pix[i]);
			std::memcpy(&correlation, a_pix, sizeof(PixelType));

			correlation = dotProduct;
			if(comparatorLeast(bestCorrelation, correlation)) //meaning that bestCorrelation < correlation
			{
				bestCorrelation = correlation;
				bestK = k;
				bestDotProduct = dotProduct;
			}
		} //best k found: represents the index of the atom that has the biggest absolute dot product with the residual
		const AtomType &atom = m_dictionary.atom(bestK);
		weights[patchIndex][bestK] = weights[patchIndex][bestK] + bestDotProduct;
		residual.for_all_pixels([&] (PixelType &pix, int x, int y)
		{
			operand1 = reinterpret_cast<const DataType *>(&bestDotProduct);
			PixelType atomAtxy = atom.content().pixelAbsolute(x, y);
			operand2 = reinterpret_cast<const DataType *>(&atomAtxy);
			for(unsigned j=0; j<pixelSize; ++j)
			{
				reinterpret_cast<DataType *>(&pix)[j] -= operand1[j] * operand2[j];
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
	AtomType &atom = m_dictionary.atom(atomId);
	atom.reset(); //Step 7: d_j := 0
	DictionaryType resultDictionaryByXi;
	std::vector<PixelType> G;
	G.reserve(sparsity());
	I resultPatchesByG, resultDictionaryXiByG;
	resultPatchesByG.initItk(m_patchSize[0], m_patchSize[1]);
	resultDictionaryXiByG.initItk(m_patchSize[0], m_patchSize[1]);
	resultPatchesByG.for_all_pixels([&] (PixelType &pix) {pix=pix_zero;});

	for(unsigned i=0; i<m_patches.size(); ++i)
	{
		PixelType weight = m_weights[i][atomId];
		if(!I::is_zero(weight)) //Step 8: I := {indices of the signals in Y (patches) whose representations use d_j}
		{
			G.push_back(weight); //Step 9: g := X^T_{j,I}

			const I& patch = m_patches[i].second;
			resultPatchesByG += patch * weight;
		}
	}
	//IO::save01_in_u8(resultPatchesByG, "/home/nlutz/debugResultYg.png");

	resultDictionaryByXi.setAtomSize(m_patchSize[0], m_patchSize[1]);
	resultDictionaryByXi.initEmpty(G.size());
	unsigned k=0;
	for(unsigned i=0; i<m_weights.size(); ++i)
	{
		PixelType weight = m_weights[i][atomId];
		if(!I::is_zero(weight))
		{
			for(unsigned j=0; j<m_dictionary.nbAtoms(); ++j)
			{
				resultDictionaryByXi.atom(k).content() += m_dictionary.atom(j).content() * m_weights[i][j];
			}
			++k;
		}
	}
	resultDictionaryXiByG = resultDictionaryByXi * G;

	AtomType d;
	d.content() = resultPatchesByG - resultDictionaryXiByG; //Step 10: d := Y_Ig - DX_Ig
	d.forceUpdateStatistics();
	PixelType norm = d.norm();

	size_t pixelSize = sizeof(PixelType)/sizeof(DataType); //nothing to see here, get moving
	d.content().for_all_pixels([&] (PixelType &pix)
	{
		for(unsigned i=0; i<pixelSize; ++i)
			reinterpret_cast<DataType *>(&pix)[i] /= reinterpret_cast<const DataType *>(&norm)[i]; //Step 11: d := d/||d||_2
	});
//	IO::save01_in_u8(debugAtom, "/home/nlutz/debugAtom.png");

	k=0;
	for(unsigned i=0; i<m_patches.size(); ++i)
	{
		PixelType &weight = m_weights[i][atomId];
		if(!I::is_zero(weight))
		{
			weight = pix_zero;
			d.content().for_all_pixels([&] (const PixelType &pix, int x, int y)
			{
				for(unsigned l=0; l<pixelSize; ++l)
				{
					reinterpret_cast<DataType *>(&weight)[l] +=
							reinterpret_cast<const DataType *>(&m_patches[i].second.pixelAbsolute(x, y))[l]
						*	reinterpret_cast<const DataType *>(&pix)[l]
					-		reinterpret_cast<const DataType *>(&resultDictionaryByXi.atom(k).content().pixelAbsolute(x, y))[l]
						*	reinterpret_cast<const DataType *>(&pix)[l];
				}
			});
			++k;
			//Step 12/14 : X_{j, I} := (Y^T_I d - (DX_I)^T d)^T

		}
	}
	atom.content() = d.content(); //Step 13: d_j := d
}

template<typename I>
template<typename Compare>
I DictionaryProcessor<I>::synthesize(unsigned width, unsigned height, unsigned nbIterations)
{
	I output;
	std::vector<std::vector<PixelType>> weights;
	std::vector<std::pair<itk::Index<2>, I>> patches;
	const unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	output.initItk(width, height);
	output.for_all_pixels([&] (PixelType &pix)
	{
		for(unsigned i=0; i<pixelSize; ++i)
		{
			reinterpret_cast<DataType *>(&pix)[i] = DataType(std::rand())/RAND_MAX;
		}
	});
	readImage(output, patches);
	weights.resize(patches.size());
	for(unsigned i=0; i<nbIterations; ++i)
	{
		matchImage(output, m_input);
		readImage(output, patches);

		for(size_t p=0; p<m_patches.size(); ++p) //2. Sparse coding (update weights)
		{
			orthogonalMatchingPursuit<Compare>(p, patches, weights);
		}
		output=reconstructImage(weights, output.width(), output.height());
	}
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
void DictionaryProcessor<I>::load(const std::string &directory)
{
	const unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);


	Histogram<I>::loadImageFromCsv(m_input, directory + "/input.csv");

	std::ifstream ifs_data_in(directory + "/data.csv");
	ifs_data_in >> m_patchOffsetX;
	ifs_data_in >> m_patchOffsetY;
	ifs_data_in >> m_sparsity;
	ifs_data_in >> m_patchSize[0] >> m_patchSize[1];
	ifs_data_in.close();

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
			for(unsigned i=0; i<pixelSize; ++i)
				ifs_weights_in >> tmpPixel[i];
			memcpy(&(*it2), &tmpPixel, sizeof(PixelType));
		}
	}

	m_nbAtoms = m_dictionary.nbAtoms();
}

}

#endif
