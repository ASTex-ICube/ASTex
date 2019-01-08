#include "residual.h"

namespace ASTex
{

namespace Sparse
{

//constructor

Residual::Residual() :
    m_layers(),
    m_size()
{}

//public

void Residual::setSize(SizeType size)
{
    m_size = size;
    for(ImageType &layer : m_layers)
    {
        layer.initItk(m_size[0], m_size[1], true);
    }
}

void Residual::setNbLayers(unsigned nbLayers)
{
    m_layers.resize(nbLayers);
}

Residual::ImageType & Residual::layer(unsigned layerId)
{
    return m_layers[layerId];
}

//private

void Residual::assert_integrity()
{
    for(ImageType &layer : m_layers)
    {
        assert(layer.size() == m_size);
    }
}

}

}
