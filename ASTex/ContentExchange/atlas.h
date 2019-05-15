#ifndef __CTEXCH_ATLAS_H__
#define __CTEXCH_ATLAS_H__

#include "patchProcessor.h"
#include <queue>

namespace ASTex
{

namespace ContentExchange
{

using PixelPos = itk::Index<2>;

//TODO : atlas generation mode (default: vertical)
template<typename I>
/**
 * @brief The Atlas class generates an atlas texture from an input patch processor.
 * Usage is the following: call Atlas(patchProcessor), then generate() with any content ID you want
 * (usually one call for each of them). Use patchPositionAt to get the position of a patch in the atlas.
 * Use generatedImage() to get the last generated image (you have to copy it).
 * You can also use saveOrigins and loadOrigins to preserve both the atlas (with classical image save)
 * and the computed origins, so you can precompute the atlas and store it in a disk.
 *
 * Format is the following: Each atlas generated from the same patchProcessor has the same size as well as
 * patch locations, no matter the contentId. The main level is put on top of the atlas, and bellow are
 * the mipmapped contents, separated from each others. The algorithm tries to put in the largest levels first.
 * Atlas accepts periodic contents, as well as non-periodic contents, but periodic contents WILL take less place in
 * the atlas.
 */
class Atlas
{
public:

	Atlas(const PatchProcessor<I> &patchProcessor);
	~Atlas();

	void generate(int contentId);

	PixelPos patchPositionAt(int patchId, size_t mw, size_t mh) const;
	const I &generatedImage() {return m_generativeAtlas;}

	void saveOrigins(const std::string &path);
	void loadOrigins(const std::string &path);

	void generateAndSaveAtlas(std::string directory) const;

	void setStoreHighestLevel(bool b) {m_storeHighestLevel = b;}

private:
	class ComparableScorePatches
	{
	public:
		ComparableScorePatches(const Patch<I>* patch, unsigned k, unsigned l, unsigned int id) :
			m_patch(patch), m_k(k), m_l(l), m_id(id) {}

		const Patch<I> *patch() const {return m_patch;}

		unsigned k() const {return m_k;}
		unsigned l() const {return m_l;}
		unsigned id() const{return m_id;}

		bool operator<(const ComparableScorePatches &other) const
		{
			return m_patch->contentAt(0).mipmap(m_k, m_l).height() < other.patch()->contentAt(0).mipmap(m_k, m_l).height();
		}

		bool operator>(const ComparableScorePatches &other) const
		{
			return m_patch->contentAt(0).mipmap(m_k, m_l).height() > other.patch()->contentAt(0).mipmap(m_k, m_l).height();
		}

	private:
		const Patch<I> *m_patch;
		unsigned m_k, m_l;
		unsigned int m_id;

	};

	const PatchProcessor<I> &m_patchProcessor;

	void _generateFirstTime(int contentId); //function used when origins were not generated yet (hard)
	void _regenerate(int contentId); //function used when origins were generated (easy)

	void init_emplaceAlgorithm(int width, int height);
	bool check_emplace(const I& image, const ImageAlphad &alpha, unsigned x, unsigned y); //check if I can be placed on x, y
	void find_emplace(const I& image, const ImageAlphad &alpha, unsigned &x, unsigned &y); //finds x and y such that I can be placed
	void emplaceLevel0(const I& image, const ImageAlphad &alpha, unsigned x, unsigned y, unsigned height); //places I on x, y with a special periodicity
	void emplace(const I& image, const ImageAlphad &alpha, unsigned x, unsigned y); //places I on x, y if alpha(x, y) isn't 0
	void release_emplaceAlgorithm();

	I m_generativeAtlas;
	ImageGrayb m_occupationMap;
	std::vector<std::vector<std::vector<PixelPos>>> m_origins;
	bool m_generatedAtLeastOnce;
	bool m_storeHighestLevel;
};

template<typename I>
Atlas<I>::Atlas(const PatchProcessor<I> &patchProcessor) :
	m_patchProcessor(patchProcessor),
	m_generatedAtLeastOnce(false),
	m_storeHighestLevel(true)
{
}

template<typename I>
Atlas<I>::~Atlas()
{}

template<typename I>
void Atlas<I>::generate(int contentId)
{
	if(!m_generatedAtLeastOnce)
	{
		_generateFirstTime(contentId);
		m_generatedAtLeastOnce = true;
	}
	else
		_regenerate(contentId);
}

template<typename I>
void Atlas<I>::_generateFirstTime(int contentId)
{
	//resize origins vector
	m_origins.resize(m_patchProcessor.nbPatches());
	for(std::vector<std::vector<std::vector<PixelPos>>>::iterator it=m_origins.begin(); it!=m_origins.end(); ++it)
	{
		(*it).resize(m_patchProcessor.patchMapMipmap().numberMipmapsWidth());
		for(std::vector<std::vector<PixelPos>>::iterator it2=(*it).begin(); it2!=(*it).end(); ++it2)
		{
			(*it2).clear();
			(*it2).resize(m_patchProcessor.patchMapMipmap().numberMipmapsHeight());
		}
	}
	int width = m_patchProcessor.patchMapMipmap().texture().width();
	int height= m_patchProcessor.patchMapMipmap().texture().height();
	//resize image
	for(unsigned n=0; n<m_patchProcessor.nbPatches(); ++n)
	{
		const Patch<I> &patch = m_patchProcessor.patchAt(n);
		for(unsigned k=0; k<patch.alphaMipmap().numberMipmapsWidth(); ++k)
		{
			if(patch.alphaMipmap().mode() == ANISOTROPIC)
			{
				for(unsigned l=0; l<patch.alphaMipmap().numberMipmapsHeight(); ++l)
				{
					if(k!=0 || l!=0 /*|| !m_patchesDoNotOverlap*/) //In case of emergency, break this comment boundary
					{
						height += patch.contentAt(contentId).mipmap(k, l).height();
					}
				}
			}
			else
			{
				if(k!=0)
					height += patch.contentAt(contentId).mipmap(k, k).height();
			}
		}
	}

	init_emplaceAlgorithm(width, height);
	unsigned x=0, y=0;
	unsigned yFindEmplace = 0;
	//first emplace the main image (if patches do not overlap)
	if(!m_patchProcessor.patchesOverlap() && m_storeHighestLevel)
	{
		for(unsigned n=0; n<m_patchProcessor.nbPatches(); ++n)
		{
			const Patch<I> &patch = m_patchProcessor.patchAt(n);
			const ImageAlphad &alpha = patch.mipmap(0, 0);
			emplaceLevel0(patch.contentAt(contentId).texture(), alpha, patch.originAt(0, 0)[0], patch.originAt(0, 0)[1],
			m_patchProcessor.texture().height());
		}
		y=m_patchProcessor.texture().height();
		yFindEmplace = y;
	}

	//then emplace the rest

	//v nested class used to compare contents' heights so as to sort them by highest height first
	for(unsigned k=0; k<m_patchProcessor.patchAt(0).alphaMipmap().numberMipmapsWidth(); ++k)
	{
		if(m_patchProcessor.patchAt(0).alphaMipmap().mode() == ANISOTROPIC)
		{
			for(unsigned l=0; l<m_patchProcessor.patchAt(0).alphaMipmap().numberMipmapsHeight(); ++l)
			{
				if(k!=0 || l!=0 || (m_patchProcessor.patchesOverlap() && m_storeHighestLevel))
				{
					std::priority_queue<ComparableScorePatches> pq;
					for(unsigned n=0; n<m_patchProcessor.nbPatches(); ++n)
					{
						ComparableScorePatches csp(&m_patchProcessor.patchAt(n), k, l, n);
						pq.push(csp);
					}

					for(unsigned n=0; n<m_patchProcessor.nbPatches(); ++n)
					{
						const ComparableScorePatches &csp = pq.top();
						const Patch<I> *patch = csp.patch();
						x=0;
						y=yFindEmplace;
						const ImageAlphad &alpha = patch->mipmap(k, l);
						find_emplace(patch->contentAt(contentId).mipmap(k, l), alpha, x, y);
						m_origins[csp.id()][k][l][0] = x;
						m_origins[csp.id()][k][l][1] = y;
						emplace(patch->contentAt(contentId).mipmap(k, l), alpha, x, y);
						pq.pop();
					}
				}
			}
		}
		else
		{
			if(k!=0 || (m_patchProcessor.patchesOverlap() && m_storeHighestLevel) )
			{
				std::priority_queue<ComparableScorePatches> pq;
				for(unsigned n=0; n<m_patchProcessor.nbPatches(); ++n)
					pq.push(ComparableScorePatches(&m_patchProcessor.patchAt(n), k, k, n));

				for(unsigned n=0; n<m_patchProcessor.nbPatches(); ++n)
				{
					const ComparableScorePatches &csp = pq.top();
					const Patch<I> *patch = csp.patch();
					x=0;
					y=yFindEmplace;
					const ImageAlphad &alpha = patch->mipmap(k, k);
					find_emplace(patch->contentAt(contentId).mipmap(k, k), alpha, x, y);
					m_origins[csp.id()][k][k][0] = x;
					m_origins[csp.id()][k][k][1] = y;
					emplace(patch->contentAt(contentId).mipmap(k, k), alpha, x, y);
					pq.pop();
				}
			}
		}
	}
	release_emplaceAlgorithm();
	m_generatedAtLeastOnce = true;
}

template<typename I>
void Atlas<I>::_regenerate(int contentId)
{
	int width = m_generativeAtlas.width();
	int height= m_generativeAtlas.height();
	//resize image

	init_emplaceAlgorithm(width, height);
	unsigned x=0, y=0;
	//first emplace the main image, unless patches overlap
	if(!m_patchProcessor.patchesOverlap() && m_storeHighestLevel)
		for(unsigned n=0; n<m_patchProcessor.nbPatches(); ++n)
		{
			const Patch<I> &patch = m_patchProcessor.patchAt(n);
			const ImageAlphad &alpha = patch.mipmap(0, 0);
			emplaceLevel0(patch.contentAt(contentId).texture(), alpha, patch.originAt(0, 0)[0], patch.originAt(0, 0)[1],
			m_patchProcessor.texture().height());
		}
		y+=m_patchProcessor.texture().height();
	//then emplace the rest

	//v nested class used to compare contents' heights so as to sort them by highest height first
	for(unsigned k=0; k<m_patchProcessor.patchAt(0).alphaMipmap().numberMipmapsWidth(); ++k)
	{
		if(m_patchProcessor.patchAt(0).alphaMipmap().mode() == ANISOTROPIC)
		{
			for(unsigned l=0; l<m_patchProcessor.patchAt(0).alphaMipmap().numberMipmapsHeight(); ++l)
			{
				if(k!=0 || l!=0 || (m_patchProcessor.patchesOverlap() && m_storeHighestLevel))
				{
					for(unsigned n=0; n<m_patchProcessor.nbPatches(); ++n)
					{
						const Patch<I> &patch = m_patchProcessor.patchAt(n);
						const ImageAlphad &alpha = patch.mipmap(k, l);
						int x = m_origins[n][k][l][0];
						int y = m_origins[n][k][l][1];
						emplace(patch.contentAt(contentId).mipmap(k, l), alpha, x, y);
					}
				}
			}
		}
		else
		{
			if(k!=0 || (m_patchProcessor.patchesOverlap() && m_storeHighestLevel))
			{
				for(unsigned n=0; n<m_patchProcessor.nbPatches(); ++n)
				{
					const Patch<I> &patch = m_patchProcessor.patchAt(n);
					const ImageAlphad &alpha = patch.mipmap(k, k);
					int x = m_origins[n][k][k][0];
					int y = m_origins[n][k][k][1];
					emplace(patch.contentAt(contentId).mipmap(k, k), alpha, x, y);
				}
			}
		}
	}
	release_emplaceAlgorithm();
}

template<typename I>
PixelPos Atlas<I>::patchPositionAt(int patchId, size_t mw, size_t mh) const
{
	assert(patchId < (int)m_patchProcessor.nbPatches() &&
		   "ContentExchange::Atlas::patchPositionAt: patch id is out of range (try bounding it by patchProcessor.nbPatches())");
	mw = std::min(mw, m_origins[patchId].size()-1);
	mh = std::min(mh, m_origins[patchId][mw].size()-1);
	return m_origins[patchId][mw][mh];
}

template<typename I>
void Atlas<I>::init_emplaceAlgorithm(int width, int height)
{
	m_occupationMap.initItk(width, height);
	m_generativeAtlas.initItk(width, height);
	m_occupationMap.for_all_pixels([] (bool &pix)
	{
		pix=false;
	});
}

template<typename I>
bool Atlas<I>::check_emplace(const I& image, const ImageAlphad& alpha, unsigned x, unsigned y)
{
	if(x + (unsigned)image.width()>(unsigned)m_occupationMap.width())
		return false;
	if(y + (unsigned)image.height()>(unsigned)m_occupationMap.height())
		return false;

//    //The following is optimal only if not using alpha
//    if(m_occupationMap.pixelAbsolute(x, y))
//        return false;
//    if(m_occupationMap.pixelAbsolute(x+image.width()-1, y))
//        return false;
//    if(m_occupationMap.pixelAbsolute(x+image.width()-1, y+image.height()-1))
//        return false;
//    if(m_occupationMap.pixelAbsolute(x, y+image.height()-1))
//        return false;

	unsigned i, j;
	bool occupied=false;
	for(i=0; i<(unsigned)image.width() && !occupied; ++i)
		for(j=0; j<(unsigned)image.height() && !occupied; ++j)
			occupied = m_occupationMap.pixelAbsolute(x+i, y+j) && alpha.pixelAbsolute(i, j)>0.0;
	return !occupied;
}

template<typename I>
void Atlas<I>::find_emplace(const I& image, const ImageAlphad& alpha, unsigned &x, unsigned &y)
{
	//precondition: width and height of images are large enough
	while(!check_emplace(image, alpha, x, y))
		if(++x >=(unsigned)m_occupationMap.width())
		{
			x=0;
			++y;
		}
}

template<typename I>
void Atlas<I>::emplaceLevel0(const I& image, const ImageAlphad& alpha, unsigned x, unsigned y, unsigned height)
{
	image.for_all_pixels([&] (const typename I::PixelType &pix, int i, int j)
	{
		if(alpha.pixelAbsolute(i, j)>0.0)
		{
			m_generativeAtlas.pixelAbsolute((x+i)%m_generativeAtlas.width(), (y+j)%height) = pix;
			m_occupationMap.pixelAbsolute((x+i)%m_occupationMap.width(), (y+j)%height) = true;
		}
	});
}


template<typename I>
void Atlas<I>::emplace(const I& image, const ImageAlphad& alpha, unsigned x, unsigned y)
{
	image.for_all_pixels([&] (const typename I::PixelType &pix, int i, int j)
	{
		if(alpha.pixelAbsolute(i, j)>0.0)
		{
			m_generativeAtlas.pixelAbsolute((x+i)%m_generativeAtlas.width(), (y+j)%m_generativeAtlas.height()) = pix;
			m_occupationMap.pixelAbsolute((x+i)%m_occupationMap.width(), (y+j)%m_occupationMap.height()) = true;
		}
	});
}

template<typename I>
void Atlas<I>::release_emplaceAlgorithm()
{
	int y=m_occupationMap.height()-1;
	bool fullLineOrLimit = false;
	int i;
	while(!fullLineOrLimit)
	{
		for(i=0; i<m_generativeAtlas.width() && !m_occupationMap.pixelAbsolute(i, y); ++i);
		fullLineOrLimit = i!=m_generativeAtlas.width() || --y==0 /*debug only TODO*//* || y==1023*/;
		//again, debug TODO
//        if(y<1023)
//        {
//            y=1023;
//            fullLineOrLimit = true;
//        }
	}
	I temporaryCrop;
	temporaryCrop.initItk(m_generativeAtlas.width(), y+1);
	temporaryCrop.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
	{
		pix = m_generativeAtlas.pixelAbsolute(x, y);
	});
	m_generativeAtlas.initItk(temporaryCrop.width(), temporaryCrop.height(), true);
	m_generativeAtlas.copy_pixels(temporaryCrop);
	m_occupationMap.initItk(1, 1);
}

template<typename I>
void Atlas<I>::saveOrigins(const std::string &path)
{
	assert(m_origins.size() > 0 && "ContentExchange::Atlas::saveOrigins: atlas not generated (try using generate())");
	std::ofstream ofs(path);
	ofs << m_origins.size() << std::endl;
	ofs << m_origins[0].size() << std::endl;
	ofs << m_origins[0][0].size() << std::endl;
	for(unsigned i=0; i<m_origins.size(); ++i)
		for(unsigned j=0; j<m_origins[i].size(); ++j)
			for(unsigned k=0; k<m_origins[i][j].size(); ++k)
				ofs << m_origins[i][j][k][0] << ' ' << m_origins[i][j][k][1] << std::endl;
	ofs.close();
}

template<typename I>
void Atlas<I>::loadOrigins(const std::string &path)
{
	std::ifstream ifs(path);
	unsigned size1, size2, size3;
	ifs >> size1;
	ifs >> size2;
	ifs >> size3;
	m_origins.resize(size1);
	for(std::vector<std::vector<std::vector<PixelPos>>>::iterator it=m_origins.begin(); it!=m_origins.end(); ++it)
	{
		(*it).resize(size2);
		for(std::vector<std::vector<PixelPos>>::iterator it2=(*it).begin(); it2!=(*it).end(); ++it2)
		{
			(*it2).clear();
			(*it2).resize(size3);
		}
	}

	for(unsigned i=0; i<m_origins.size(); ++i)
		for(unsigned j=0; j<m_origins[i].size(); ++j)
			for(unsigned k=0; k<m_origins[i][j].size(); ++k)
			{
				int s, t;
				ifs >> s;
				ifs >> t;
				m_origins[i][j][k][0]=s;
				m_origins[i][j][k][1]=t;
			}

	ifs.close();
}

template<typename I>
void Atlas<I>::generateAndSaveAtlas(std::string directory) const
{
    for(size_t i=0; i<m_patchProcessor.nbContents(); ++i)
    {
        generate(i);
        m_generativeAtlas.save(directory + "/contentAtlas_" + std::to_string(i) + ".png");
		Histogram<I>::saveImageToCsv(m_generativeAtlas, directory + "/contentAtlas_" + std::to_string(i) + ".csv");
    }
    saveOrigins(directory + "/atlasOrigins.csv");
}

}

}

#endif
