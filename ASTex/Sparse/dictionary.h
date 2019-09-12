#ifndef _DICTIONNARY_H_
#define _DICTIONNARY_H_

#include "ASTex/fourier.h"
#include "ASTex/image_common.h"
#include "ASTex/easy_io.h" //TODO: remove when fully debugged

#include "ASTex/histogram.h"
#include "ASTex/rpn_utils.h"

namespace ASTex
{

template <typename I>
class Atom
{
public:

	Atom(typename std::enable_if<std::is_floating_point<typename I::DataType>::value>::type* = 0);

	typename I::PixelType norm();

	const I &content() const;
	I &content();

	void forceUpdateStatistics() {m_statisticsNeedUpdate = true;}

	void reset();

	Atom &operator=(const Atom& other);

	template <typename FUNC>
	typename std::enable_if<function_traits<FUNC>::arity == 1, I>::type visualize(const FUNC &f) const;

private:

	typename I::PixelType			m_norm;
	bool							m_statisticsNeedUpdate;
	I								m_content;
};

template <typename I>
class Dictionary
{
public:

	Dictionary();
	Dictionary(const Dictionary &other);

	size_t atomWidth() const {return m_atomWidth;}
	size_t atomHeight() const {return m_atomHeight;}

	void setAtomSize(size_t width, size_t height);

	Atom<I> &atom(size_t i);
	const Atom<I> &atom(size_t i) const;

	void erase(size_t i);

	size_t nbAtoms() const {return m_atoms.size();}

	void initRandom(unsigned nbAtoms);
	void initEmpty(unsigned nbAtoms);

	I operator*(const std::vector<typename I::PixelType> &other) const;

	bool save(const std::string &directory) const;
	bool load(const std::string &directory);

private:

	std::vector<Atom<I>> m_atoms;

	size_t m_atomWidth;
	size_t m_atomHeight;
};

////////////////////ATOM/////////////////////

template <typename I>
Atom<I>::Atom(typename std::enable_if<std::is_floating_point<typename I::DataType>::value>::type*) :
	m_norm(), m_statisticsNeedUpdate(false)
{
}

template <typename I>
typename I::PixelType Atom<I>::norm()
{
	using PixelType = typename I::PixelType;
	using DataType = typename I::DataType;
	const unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	if(m_statisticsNeedUpdate)
	{
		//encapsulate in a function if statistics are not only comprised of the norm.
		m_statisticsNeedUpdate = false;
		m_norm = PixelType{};
		m_content.for_all_pixels([&] (const typename I::PixelType &pix)
		{
			for(unsigned i=0; i<pixelSize; ++i)
				reinterpret_cast<DataType *>(&m_norm)[i] +=
					reinterpret_cast<const DataType *>(&pix)[i]
					*reinterpret_cast<const DataType *>(&pix)[i];
		});
	}

	size_t arraySize = sizeof(typename I::PixelType) / sizeof(typename I::DataType);
	typename I::DataType *normSquared = reinterpret_cast<typename I::DataType *>(&m_norm);
	for(unsigned i=0; i<arraySize; ++i)
	{
		normSquared[i] = std::sqrt(normSquared[i]);
	}
	return m_norm;
}

template <typename I>
const I &Atom<I>::content() const
{
	return m_content;
}

template <typename I>
I &Atom<I>::content()
{
	m_statisticsNeedUpdate = true;
	return const_cast<I&>(static_cast<const Atom<I> *>(this)->content());
}

template <typename I>
void Atom<I>::reset()
{
	static typename I::PixelType pix_zero;
	m_content.parallel_for_all_pixels([&] (typename I::PixelType &pix)
	{
		pix = pix_zero;
	});
}

template <typename I>
Atom<I> &Atom<I>::operator=(const Atom<I>& other)
{
	m_content.initItk(other.content().width(), other.content().height());
	m_content.copy_pixels(other.content());
	return (*this);
}

template <typename I>
template <typename FUNC>
typename std::enable_if<function_traits<FUNC>::arity == 1, I>::type Atom<I>::visualize(const FUNC &f) const
{
	I output;
	output.initItk(m_content.width(), m_content.height());
	output.copy_pixels(m_content);
	output.for_all_pixels(f);
	return output;
}

//////////////////DICTIONNARY///////////////////

template <typename I>
Dictionary<I>::Dictionary() :
	m_atoms(), m_atomWidth(0), m_atomHeight(0)
{}

template <typename I>
Dictionary<I>::Dictionary(const Dictionary &other) :
	m_atoms(other.nbAtoms()), m_atomWidth(other.atomWidth()), m_atomHeight(other.atomHeight())
{
	for(unsigned i=0; i<other.nbAtoms(); ++i)
	{
		m_atoms[i] = other.atom(i);
	}
}

template <typename I>
void Dictionary<I>::setAtomSize(size_t width, size_t height)
{
	m_atomWidth = width;
	m_atomHeight = height;
}

template <typename I>
const Atom<I> &Dictionary<I>::atom(size_t i) const
{
	assert((unsigned)i<m_atoms.size() && "Dictionnary::atom: index out of range");
	return m_atoms[i];
}

template <typename I>
Atom<I> &Dictionary<I>::atom(size_t i)
{
	return const_cast<Atom<I>&>(static_cast<const Dictionary<I> *>(this)->atom(i));
}

template <typename I>
void Dictionary<I>::erase(size_t i)
{
	m_atoms.erase(m_atoms.begin()+i);
}

template <typename I>
void Dictionary<I>::initRandom(unsigned nbAtoms)
{
	m_atoms.resize(nbAtoms);
	size_t pixelSize = sizeof(typename I::PixelType) / sizeof(typename I::DataType);
	typename I::DataType *pixd = new typename I::DataType[pixelSize]();
	for(Atom<I> &atom : m_atoms)
	{
		atom.content().initItk(m_atomWidth, m_atomHeight);
		atom.content().for_all_pixels([&] (typename I::PixelType &pix) //first iteration, random vector
		{
			for(unsigned k=0; k<pixelSize; ++k)
			{
				pixd[k] = typename I::DataType(std::rand())/RAND_MAX;
			}
			std::memcpy(&pix, pixd, sizeof(typename I::PixelType));
		});
		atom.forceUpdateStatistics();
		typename I::PixelType norm = atom.norm();
		typename I::DataType *normD = reinterpret_cast<typename I::DataType *>(&norm);
		atom.content().for_all_pixels([&] (typename I::PixelType &pix) //second iteration, random unit vector
		{
			std::memcpy(pixd, &pix, sizeof(typename I::PixelType));
			for(unsigned k=0; k<pixelSize; ++k)
			{
				if(normD[k]<std::numeric_limits<float>::epsilon())
				{
					std::cerr << "Dictionary::initRandom: one of the norms of the initialized atoms is 0" << std::endl;
					exit(EXIT_FAILURE);
				}
				pixd[k] /= normD[k];
			}
			std::memcpy(&pix, pixd, sizeof(typename I::PixelType));
		});
		atom.forceUpdateStatistics();
	}
	delete[] pixd;
}

template <typename I>
void Dictionary<I>::initEmpty(unsigned nbAtoms)
{
	static typename I::PixelType pix_zero; //TODO: remove when initItk is fixed
	m_atoms.resize(nbAtoms);
	for(Atom<I> &atom : m_atoms)
	{
		atom.content().initItk(m_atomWidth, m_atomHeight, true);
		atom.content().for_all_pixels([&] (typename I::PixelType &pix) //first iteration, random vector
		{
			pix = pix_zero;
		});
	}
}

template <typename I>
I Dictionary<I>::operator*(const std::vector<typename I::PixelType> &other) const
{
	assert(other.size() == m_atoms.size());
	using DataType = typename I::DataType;
	const unsigned pixelSize = sizeof(typename I::PixelType)/sizeof(DataType);
	static typename I::PixelType pix_zero; //TODO: remove if initItk fixed
	I output;
	output.initItk(m_atomWidth, m_atomHeight , false);
	output.for_all_pixels([&] (typename I::PixelType &pix)
	{
		pix = pix_zero;
	}); //TODO: remove if initItk fixed
	for(unsigned i = 0; i<m_atoms.size(); ++i)
	{
		const typename I::PixelType &weight=other[i];
		if(!I::is_zero(weight))
		{
			output.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
			{
				for(unsigned l=0; l<pixelSize; ++l)

				reinterpret_cast<DataType *>(&pix)[l] +=
						reinterpret_cast<const DataType *>(&m_atoms[i].content().pixelAbsolute(x, y))[l]
					*	reinterpret_cast<const DataType *>(&weight)[l];
			});
		}
	}
	return output;
}

template <typename I>
bool Dictionary<I>::save(const std::string &directory) const
{
	if(!create_directory(directory))
		return false;
	unsigned i=0;
	for(typename std::vector<Atom<I>>::const_iterator it=m_atoms.begin(); it!=m_atoms.end(); ++it, ++i)
	{
		if(!Histogram<I>::saveImageToCsv((*it).content(), std::string(directory) + "/atom" + std::to_string(i) + ".png"))
			return false;
	}
	std::ofstream ofs_data_out(directory + "/dictionaryData.csv");
	if(!ofs_data_out)
		return false;
	ofs_data_out << m_atoms.size() << std::endl;
	ofs_data_out.close();
	return true;
}

template <typename I>
bool Dictionary<I>::load(const std::string &directory)
{
	unsigned i=0, sizeAtoms=0;

	std::ifstream ifs_data_in(directory + "/dictionaryData.csv");
	if(ifs_data_in.is_open())
	{
		ifs_data_in >> sizeAtoms;
		ifs_data_in.close();
	}
	else
		return false;

	m_atoms.resize(sizeAtoms);
	for(typename std::vector<Atom<I>>::iterator it=m_atoms.begin(); it!=m_atoms.end(); ++it, ++i)
	{
		if(!Histogram<I>::loadImageFromCsv((*it).content(), std::string(directory) + "/atom" + std::to_string(i) + ".png"))
			return false;
	}
	m_atomWidth =  (*m_atoms.begin()).content().width();
	m_atomHeight = (*m_atoms.begin()).content().height();

	return true;
}

}

#endif
