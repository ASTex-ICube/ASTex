#ifndef __CTEXCH_AlternativeAtlas_H__
#define __CTEXCH_AlternativeAtlas_H__

#include "patchProcessor.h"
#include <queue>

namespace ASTex
{

namespace ContentExchange
{

using PixelPos = itk::Index<2>;

//TODO : AlternativeAtlas generation mode (default: vertical)
template<typename I>
/**
 * @brief The AlternativeAtlas class generates an AlternativeAtlas texture from an input patch processor.
 * Usage is the following: call AlternativeAtlas(patchProcessor), then generate() with any content ID you want
 * (usually one call for each of them). Use patchPositionAt to get the position of a patch in the AlternativeAtlas.
 * Use generatedImage() to get the last generated image (you have to copy it).
 * You can also use saveOrigins and loadOrigins to preserve both the AlternativeAtlas (with classical image save)
 * and the computed origins, so you can precompute the AlternativeAtlas and store it in a disk.
 *
 * Format is the following: Each AlternativeAtlas generated from the same patchProcessor has the same size as well as
 * patch locations, no matter the contentId. The main level is put on top of the AlternativeAtlas, and bellow are
 * the mipmapped contents, separated from each others. The algorithm tries to put in the largest levels first.
 * AlternativeAtlas accepts periodic contents, as well as non-periodic contents, but periodic contents WILL take less place in
 * the AlternativeAlternativeAtlas.
 */
class AlternativeAtlas
{
public:

    AlternativeAtlas(const PatchProcessor<I> &patchProcessor);
    ~AlternativeAtlas();

    void generate(int contentId);

    PixelPos patchPositionAt(int patchId, size_t mw, size_t mh) const;
    const I &generatedImage() {return m_generativeAlternativeAtlas;}

    void saveOrigins(const std::string &path);
    void loadOrigins(const std::string &path);

    void generateAndSaveAlternativeAtlas(std::string directory) const;

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

    I m_generativeAlternativeAtlas;
    ImageGrayb m_occupationMap;
    std::vector<std::vector<std::vector<PixelPos>>> m_origins;
    bool m_generatedAtLeastOnce;
    bool m_storeHighestLevel;
};

template<typename I>
AlternativeAtlas<I>::AlternativeAtlas(const PatchProcessor<I> &patchProcessor) :
    m_patchProcessor(patchProcessor),
    m_generatedAtLeastOnce(false),
    m_storeHighestLevel(true)
{
}

template<typename I>
AlternativeAtlas<I>::~AlternativeAtlas()
{}

template<typename I>
void AlternativeAtlas<I>::generate(int contentId)
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
void AlternativeAtlas<I>::_generateFirstTime(int contentId)
{
    static typename I::PixelType zero;
    std::vector<I> atlasMipmaps;
    assert(m_patchProcessor.numberMipmapsWidth()>0);
    atlasMipmaps.resize(std::max(m_patchProcessor.numberMipmapsWidth(), m_patchProcessor.numberMipmapsHeight()) -1);
    for(unsigned k=0; k<atlasMipmaps.size(); ++k)
    {
        I& atlasMipmap = atlasMipmaps[k];
        atlasMipmap.initItk(m_patchProcessor.nbPatches()*(std::pow(2.0, k+1)), std::pow(2.0, k+1), true);
        atlasMipmap.for_all_pixels([&] (typename I::PixelType &pix)
        {
            pix = zero;
        });
        for(unsigned i=0; i<m_patchProcessor.nbPatches(); ++i)
        {
            itk::Index<2> origin;
            origin[0] = (std::pow(2.0, k+1))*i;
            origin[1] = 0;
            const I& contentMipmap = m_patchProcessor.patchAt(i).contentAt(contentId).mipmap(atlasMipmaps.size()-k,
                                                                                             atlasMipmaps.size()-k);
            contentMipmap.for_all_pixels([&] (const typename I::PixelType &pix, int x, int y)
            {
                atlasMipmap.pixelAbsolute(origin[0]+x, origin[1]+y) = pix;
            });
        }
        atlasMipmap.save(std::string("/home/nlutz/help") + std::to_string(k) + ".png");
    }
    return;

}

template<typename I>
void AlternativeAtlas<I>::_regenerate(int contentId)
{
}

template<typename I>
PixelPos AlternativeAtlas<I>::patchPositionAt(int patchId, size_t mw, size_t mh) const
{
    assert(patchId < (int)m_patchProcessor.nbPatches() &&
           "ContentExchange::AlternativeAtlas::patchPositionAt: patch id is out of range (try bounding it by patchProcessor.nbPatches())");
    mw = std::min(mw, m_origins[patchId].size()-1);
    mh = std::min(mh, m_origins[patchId][mw].size()-1);
    return m_origins[patchId][mw][mh];
}

template<typename I>
void AlternativeAtlas<I>::saveOrigins(const std::string &path)
{
    assert(m_origins.size() > 0 && "ContentExchange::AlternativeAtlas::saveOrigins: AlternativeAtlas not generated (try using generate())");
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
void AlternativeAtlas<I>::loadOrigins(const std::string &path)
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
void AlternativeAtlas<I>::generateAndSaveAlternativeAtlas(std::string directory) const
{
    for(size_t i=0; i<m_patchProcessor.nbContents(); ++i)
    {
        generate(i);
        m_generativeAlternativeAtlas.save(directory + "/contentAlternativeAtlas_" + std::to_string(i) + ".png");
    }
    saveOrigins(directory + "/AlternativeAtlasOrigins.csv");
}

}

}

#endif
