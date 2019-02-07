#ifndef _DICTIONNARY_H_
#define _DICTIONNARY_H_

#include "ASTex/fourier.h"
#include "ASTex/image_common.h"
#include "ASTex/easy_io.h" //TODO: remove when fully debugged

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

	size_t atomWidth() const {return m_atomWidth;}
	size_t atomHeight() const {return m_atomHeight;}

	void setAtomSize(size_t width, size_t height);

	Atom<I> &atom(size_t i);
	const Atom<I> &atom(size_t i) const;

	size_t nbAtoms() const {return m_atoms.size();}

	void initRandom(unsigned nbAtoms);
	void initEmpty(unsigned nbAtoms);

	I operator*(const std::vector<typename I::PixelType> &other) const;

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
void Dictionary<I>::initRandom(unsigned nbAtoms)
{
	m_atoms.resize(nbAtoms);
	size_t arraySize = sizeof(typename I::PixelType) / sizeof(typename I::DataType);
	typename I::DataType *pixd = new typename I::DataType[arraySize]();
	for(Atom<I> &atom : m_atoms)
	{
		atom.content().initItk(m_atomWidth, m_atomHeight);
		atom.content().for_all_pixels([&] (typename I::PixelType &pix) //first iteration, random vector
		{
			for(unsigned i=0; i<arraySize; ++i)
			{
				pixd[i] = typename I::DataType(std::rand())/RAND_MAX;
			}
			std::memcpy(&pix, pixd, sizeof(typename I::PixelType));
		});
		atom.forceUpdateStatistics();
		typename I::PixelType norm = atom.norm();
		typename I::DataType *normD = reinterpret_cast<typename I::DataType *>(&norm);
		atom.content().for_all_pixels([&] (typename I::PixelType &pix) //second iteration, random unit vector
		{
			std::memcpy(pixd, &pix, sizeof(typename I::PixelType));
			for(unsigned i=0; i<arraySize; ++i)
			{
				pixd[i] /= normD[i];
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
	}); 	//TODO: remove if initItk fixed
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

}

#endif
