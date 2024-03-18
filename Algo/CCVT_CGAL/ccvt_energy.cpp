//
// Created by grenier on 12/10/23.
//

#include "ccvt.h"
#include "timer.h"

FT CCVT::compute_weight_threshold(FT epsilon) const
{
    // reference: 1e-4 for 1000 sites
    FT A = compute_value_integral();
    unsigned N = count_visible_sites();
    return  (0.1*epsilon) * (A) / FT(N);
}

FT CCVT::compute_position_threshold(FT epsilon) const
{
    // reference: 1e-4 for 1000 sites
    FT A = compute_value_integral();
    unsigned N = count_visible_sites();
    return (0.1*epsilon) * (std::sqrt(A*A*A)) / FT(N);
}




FT CCVT::compute_wcvt_energy()
{
    if (m_timer_on) Timer::start_timer(m_timer, COLOR_BLUE, "Energy");

    FT cvt = 0.0;
    FT w_dot_V = 0.0;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;

        FT Vi = vi->compute_area();
        FT wi = vi->get_weight();
        FT Ci = m_capacities[i];
        w_dot_V += wi * (Vi - Ci);
        cvt += vi->compute_variance();
    }

    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);

    return (w_dot_V - cvt); // équation 2 ?
}




void CCVT::compute_weight_gradient(std::vector<FT>& gradient, FT coef)
{
    if (m_timer_on) Timer::start_timer(m_timer, COLOR_BLUE, "WGrad");

    gradient.clear();
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;

        FT Ci = m_capacities[i]; // m (objectif)
        FT Vi = vi->compute_area(); // mi (volume actuel)
        FT Gi = Vi - Ci; // m-mi is coef -1; mi-m si coef 1
        gradient.push_back(coef * Gi);
    }

    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);
}





void CCVT::compute_position_gradient(std::vector<Vector>& gradient, FT coef)
{
    if (m_timer_on) Timer::start_timer(m_timer, COLOR_BLUE, "XGrad");

    gradient.clear();
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;

        FT Vi = vi->compute_area(); // mi
        Point xi = vi->get_position(); // xi
        Point ci = vi->compute_centroid(); // bi
        Vector gi = -2.0*Vi*(xi - ci); // 2mi(xi-bi)
        gradient.push_back(coef * gi);
    }

    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);
}






void CCVT::compute_neightbour_gradient(std::vector<Vector>& gradient, FT coef)
{
    if (m_timer_on) Timer::start_timer(m_timer, COLOR_BLUE, "NGrad");

    gradient.clear();
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;

        FT Vi = vi->compute_area(); // mi
        Point xi = vi->get_position(); // xi
        Point ci = vi->compute_centroid(); // bi
        Vector gi = -2.0*Vi*(xi - ci); // 2mi(xi-bi) // TODO
        gradient.push_back(coef * gi);
    }

    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);
}

