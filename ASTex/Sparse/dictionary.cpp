#include "dictionary.h"

namespace ASTex
{

namespace Sparse
{

DictionaryProcessor::DictionaryProcessor() :
    m_step(0),
    m_residuals()
{}

//step 0

void DictionaryProcessor::setInput(const std::vector<ImageGrayd> &input)
{
    _step_set(1);
    m_input = input;
}

const std::vector<ImageGrayd> &DictionaryProcessor::input() const
{
    _step_assert(1);
    return m_input;
}

void DictionaryProcessor::generateKSVD()
{
    _step_assert(1);
    _step_set(2);
}

//step 1

const Residual & DictionaryProcessor::residual(unsigned index) const
{
    assert(index < m_residuals.size());
    _step_assert(1);
    return m_residuals[index];
}

//steo handling

void DictionaryProcessor::_step_set(unsigned step)
{
    m_step = step;
}

bool DictionaryProcessor::_step_assert(unsigned expectedStep)
{
    bool isCallTooEarly;
    if(isCallTooEarly=(m_step < expectedStep))
    {
        std::cerr << "Unexpected call detected. Step " << m_step << ": " << m_step_details[m_step];
    }
    assert(!isCallTooEarly);
    return !isCallTooEarly;
}

//ATOM

Atom::Atom(size_t size) :
    m_array(size)
{}

//operators

Atom::RealType &Atom::operator[](unsigned index)
{
    return const_cast<Atom::RealType&>(static_cast<const Atom*>(this)->m_array[index]);
}

const Atom::RealType &Atom::operator[](unsigned index) const
{
    return m_array[index];
}

}

}
