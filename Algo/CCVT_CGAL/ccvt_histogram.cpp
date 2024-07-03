//
// Created by grenier on 12/10/23.
//

#include "ccvt.h"


void CCVT::compute_capacity_histogram(std::vector<unsigned>& histogram) const
{
    unsigned nsites = count_visible_sites();
    if (nsites == 0) return;

    unsigned nbins = histogram.size();
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;

        FT Ci = m_capacities[i];
        FT Vi = vi->compute_area();
        FT value = 0.5 * (Vi / Ci);
        if (value > 1.0) value = 1.0;
        if (value < 0.0) value = 0.0;

        unsigned bin = unsigned(ceil(value*(nbins-1)));
        histogram[bin]++;
    }
}


void CCVT::compute_weight_histogram(const double,
                                     std::vector<unsigned>& histogram) const
{
    unsigned nsites = count_visible_sites();
    if (nsites == 0) return;

    std::vector<FT> weights;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        FT value = vi->get_weight();
        weights.push_back(value);
    }

    FT sum = 0.0;
    FT min_w = 0.0;
    FT max_w = 0.0;
    for (unsigned i = 0; i < weights.size(); ++i)
    {
        min_w = std::min(min_w, weights[i]);
        max_w = std::max(max_w, weights[i]);
        sum += weights[i];
    }

    FT range_w = max_w - min_w;
    if (range_w == 0.0) range_w = 1.0;
    FT avg_len2 = std::pow(m_rt.get_average_length(), 2);

    /*
    std::cout << "AvgLen2: " << avg_len2 << std::endl;
    std::cout << "WeightsRange: " << min_w << " ; " << max_w << std::endl;
    std::cout << "Weights\%: " << min_w/avg_len2 << " ; " << max_w/avg_len2 << std::endl;
    */

    unsigned nbins = histogram.size();
    for (unsigned i = 0; i < weights.size(); ++i)
    {
        FT value = (weights[i] - min_w) / range_w;
        if (value < 0.0) value = 0.0;
        value = std::pow(value, 0.3);
        unsigned bin = unsigned(ceil(value*(nbins-1)));
        histogram[bin]++;
    }
}
