#include "sampler.h"

namespace ASTex
{

namespace Stamping
{

std::vector<Eigen::Vector2f> RegularSampler::generate()
{
    std::vector<Eigen::Vector2f> SamplePoints;
    SamplePoints.reserve(m_nbPoints);
    float step = float(1.0/float(m_nbPoints));
    for (unsigned i = 0; i < m_nbPoints; ++i)
    {
        for (unsigned j = 0; j < m_nbPoints; ++j)
        {
            Eigen::Vector2f tmp_point(float(i)*step,float(j)*step);
            SamplePoints.push_back(tmp_point);
        }
    }

    return SamplePoints;
}

std::vector<Eigen::Vector2f> UniformSampler::generate()
{
    std::vector<Eigen::Vector2f> SamplePoints;
    SamplePoints.reserve(m_nbPoints);
    srand(time(NULL));
    for (unsigned i = 0; i < m_nbPoints; ++i)
    {
        Eigen::Vector2f tmp_point(double(rand())/RAND_MAX, double(rand())/RAND_MAX);
        SamplePoints.push_back(tmp_point);
    }
    return SamplePoints;
}

} //namsepace Stamping

} //namespace ASTex
