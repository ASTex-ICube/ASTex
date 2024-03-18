//
// Created by grenier on 12/10/23.
//

#include "ccvt.h"
#include "timer.h"

bool CCVT::is_valid() const
{
    if (!m_domain.is_valid()) return false;
    if (m_vertices.empty()) return false;
    return true;
}

unsigned CCVT::count_visible_sites() const
{
    unsigned nb = 0;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        nb++;
    }
    return nb;
}

void CCVT::collect_visible_points(std::vector<Point>& points) const
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        points.push_back(vi->get_position());
    }
}

void CCVT::collect_visible_weights(std::vector<FT>& weights) const
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        weights.push_back(vi->get_weight());
    }
}

void CCVT::collect_sites(std::vector<Point>& points,
                          std::vector<FT>& weights) const
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        Point pi = vi->get_position();
        points.push_back(pi);

        FT wi = 0.0;
        wi = vi->get_weight();
        weights.push_back(wi);
    }
}

void CCVT::clear_triangulation()
{
    m_ratio.clear();
    m_vertices.clear();
    m_capacities.clear();
    m_rt.clear();
}

bool CCVT::update_triangulation(bool skip)
{
    std::vector<FT> weights;
    std::vector<Point> points;
    collect_sites(points, weights);
    return construct_triangulation(points, weights, skip);
}

bool CCVT::construct_triangulation(const std::vector<Point>& points,
                                    const std::vector<FT>& weights,
                                    bool skip)
{
    if (m_timer_on)
    {
        Timer::start_timer(m_timer, COLOR_BLUE, "Triangulation");
        std::cout << std::endl;
    }

    clear_triangulation();
    bool ok = populate_vertices(points, weights);
    if (ok || !skip)
    {
        pre_build_dual_cells();
        assign_pixels();
        pre_compute_area();
        compute_capacities(m_capacities);
    }

    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);
    return (ok || !skip);
}



bool CCVT::populate_vertices(const std::vector<Point>& points,
                              const std::vector<FT>& weights)
{
    if (m_timer_on) Timer::start_timer(m_timer, COLOR_YELLOW, "Populate");

    unsigned nb = 0;
    unsigned nsites = points.size();
    for (unsigned i = 0; i < nsites; ++i)
    {
        Vertex_handle vertex = insert_vertex(points[i], weights[i], nb);
        if (vertex == Vertex_handle()) continue;
        m_vertices.push_back(vertex);
        nb++;
    }

    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_YELLOW);

    bool none_hidden = true;
    if (count_visible_sites() != m_vertices.size())
        none_hidden = false;

    return none_hidden;
}

Vertex_handle CCVT::insert_vertex(const Point& point,
                                   const FT weight,
                                   const unsigned index)
{
    Weighted_point wp(point, weight);
    Vertex_handle vertex = m_rt.insert(wp);

    if (vertex->get_index() != -1)
        return Vertex_handle();

    vertex->set_index(index);
    return vertex;
}

FT CCVT::compute_mean_capacity() const
{
    FT domain_area = compute_value_integral();
    unsigned nb = count_visible_sites();
    return (domain_area / nb);
}

void CCVT::compute_capacities(std::vector<FT>& capacities) const
{
    if(!m_custom_proportions){
        FT C = compute_mean_capacity();
        for (unsigned i = 0; i < m_vertices.size(); ++i)
        {
            FT Ci = 0.0;
            Vertex_handle vi = m_vertices[i];
            if (!vi->is_hidden()) Ci = C;
            capacities.push_back(Ci);
        }
    }
    else{
        FT domain_area = compute_value_integral();
        for (unsigned i = 0; i < m_vertices.size(); ++i)
        {
            FT Ci = 0.0;
            Vertex_handle vi = m_vertices[i];
            if (!vi->is_hidden()) Ci = domain_area*m_proportions.at(i);
            capacities.push_back(Ci);
        }
    }
}


void CCVT::update_positions(const std::vector<Point>& points, bool clamp, bool hidden)
{
    unsigned j = 0;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (hidden && vi->is_hidden()) continue;

        Point pi = points[j++];
        if (clamp) pi = m_domain.clamp(pi);
        vi->set_position(pi);
    }
}

void CCVT::update_weights(const std::vector<FT>& weights, bool hidden)
{
    unsigned j = 0;
    FT mean = compute_mean(weights);
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (hidden && vi->is_hidden()) continue;
        vi->set_weight(weights[j++] - mean);
    }
}

void CCVT::reset_weights()
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        vertex->set_weight(0.0);
    }
    update_triangulation();
}

FT CCVT::compute_value_integral() const
{
    return m_domain.integrate_intensity();
}

void CCVT::pre_build_dual_cells()
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        if (vertex->is_hidden()) continue;

        bool ok = m_rt.pre_build_polygon(vertex, vertex->dual().points());
        /*
        if (!ok)
            std::cout << "Vertex " << vertex->get_index()
            << ": pre_build_dual_cell failed" << std::endl;
        */
    }
}

void CCVT::pre_compute_area()
{
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vertex = m_vertices[i];
        vertex->pre_compute_area();
    }
}
