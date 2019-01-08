#ifndef __DICTIONARY_H__
#define __DICTIONARY_H__

#include "residual.h"

namespace ASTex
{

namespace Sparse
{

class DictionaryProcessor
{
public:

    DictionaryProcessor();

    //step 0

    void setInput(const std::vector<ImageGrayd> &input);

    //step 1

    const std::vector<ImageGrayd> &input() const;
    void generateKSVD();

    //step 2

    const Residual &residual(unsigned index) const;

private:

    std::string *m_step_details = { "Define an input atom (a texture) by using setInputAtom().",
                                    "Use generateXXX() to generate a dictionary."};

    void _step_set(unsigned step);
    bool _step_assert(unsigned expectedStep);

    unsigned m_step;
    const std::vector<ImageGrayd> &m_input;
    std::vector<Residual> m_residuals;
    std::vector<Atom> m_dictionary;
};

class Atom
{
public:

    Atom(size_t size);

    //types

    using RealType = double;

    //operators

    RealType &operator[](unsigned index);
    const RealType &operator[](unsigned index) const;

private:
    itk::Array<RealType> m_array;
};

}

}

#endif
