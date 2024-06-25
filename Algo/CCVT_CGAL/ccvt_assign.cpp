//
// Created by grenier on 12/10/23.
//

#include "ccvt.h"
#include "timer.h"


Pixel CCVT::build_pixel(unsigned i, unsigned j) const
{
    FT hx = 0.5*m_domain.get_px();
    FT hy = 0.5*m_domain.get_py();
    FT x  = m_domain.compute_x(i);
    FT y  = m_domain.compute_y(j);

    ConvexPolygon shape;
    FT value = m_domain.get_value(i, j);
    shape.init_rectangle(hx, hy, Point(x+hx,y-hy));
    return Pixel(shape, value);
}

void CCVT::set_ratio(Edge edge, FT value)
{
    Vertex_handle vs = m_rt.get_source(edge);
    Vertex_handle vt = m_rt.get_target(edge);
    if (vs->get_index() > vt->get_index()) edge = m_rt.get_twin(edge);
    m_ratio[edge] = value;
}

FT CCVT::get_ratio(Edge edge) const
{
//    Vertex_handle vs = m_rt.get_source(edge);
//    Vertex_handle vt = m_rt.get_target(edge);
//    if (vs->get_index() > vt->get_index()) edge = m_rt.get_twin(edge);
//    std::map<Edge, FT>::const_iterator it = m_ratio.find(edge);
//    if (it == m_ratio.end()) return 0.0;
//    return it->second;

    // edge  =  e_ij
    Segment edge_star = m_rt.build_bounded_dual_edge(edge); // e*_ij
    double n_eij = m_rt.get_length(edge);
    double n_eij_star = m_rt.get_length(edge_star);
    if(n_eij <= 0 or n_eij_star <= 0){ return 0.; }
    double ratio = n_eij_star / n_eij; // |e*ij|/|eij|

    // intégral de rho sur un segment de 0 à 1 couvrant e*ij (produit des gaussiennes marginales)
    Point c_ijl = edge_star.target();
    Point c_ijk = edge_star.source();

    double a = c_ijl.x()-c_ijk.x();
    double b = c_ijl.y()-c_ijk.y();
    double mu_1 = m_domain.get_mu_x()-m_domain.get_dx() - c_ijk.x();
    double mu_2 = m_domain.get_mu_y()-m_domain.get_dy() - c_ijk.y();

    // constantes du produit des gaussiennes
    double A = product_gaussian_amplitude(a, b, mu_1, mu_2, m_domain.get_sig_x(), m_domain.get_sig_y());
    double mu = product_gaussian_mean(a, b, mu_1, mu_2, m_domain.get_sig_x(), m_domain.get_sig_y());
    double var = product_gaussian_variance(a, b, mu_1, mu_2, m_domain.get_sig_x(), m_domain.get_sig_y());

    double int_rho = A * compute_int01_gauss_t(mu, sqrt(var));
//    std::cout<<ratio*int_rho<<" "<<it->second<<std::endl;
    return ratio*int_rho*m_domain.get_max_value();
}

void CCVT::clean_pixels()
{
    for (Finite_vertices_iterator
                 vit = m_rt.finite_vertices_begin();
         vit != m_rt.finite_vertices_end();
         ++vit)
    {
        Vertex_handle vertex = vit;
        vertex->clear_pixels();
    }
}

void CCVT::assign_pixels()
{
    unsigned width  = m_domain.get_width();
    unsigned height = m_domain.get_height();

    if (m_timer_on) Timer::start_timer(m_timer, COLOR_GREEN, "Rasterize");

    m_ratio.clear();
    Grid grid(width, height);
    for (Finite_edges_iterator
                 eit = m_rt.finite_edges_begin();
         eit != m_rt.finite_edges_end();
         ++eit)
    {
        FT value = 0.0;
        Edge edge = *eit;

        Segment segment = m_rt.build_bounded_dual_edge(edge);
        if (std::sqrt(segment.squared_length()) > EPS)
        {
            // see convention in Rt
            EnrichedSegment enriched_segment(segment,
                                             m_rt.get_source(edge),
                                             m_rt.get_target(edge));
            value = rasterize(enriched_segment, grid);
        }

        FT ratio = value;
        FT len = m_rt.get_length(edge);
        if (len == 0.0) ratio = 0.0;
        else            ratio /= len; // len = |e_ij|
        set_ratio(edge, ratio);
    }

    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_GREEN);
    if (m_timer_on)
    {
        Timer::start_timer(m_timer, COLOR_GREEN, "Assign");
        std::cout << std::endl;
    }


    Vertex_handle vertex = Vertex_handle();
    if (m_timer_on) Timer::start_timer(m_timer, COLOR_YELLOW, "boucle");
    for (unsigned i = 0; i < width; ++i)
    {
        for (unsigned j = 0; j < height; ++j)
        {
            Pixel pixel = build_pixel(i, j);
            if (grid.is_empty(i, j)) // with no cuts
            {
                Point q = pixel.compute_centroid();
                vertex = m_rt.find_nearest_vertex(q, vertex);
                vertex->append_pixel(pixel);
            }
            else // with cuts
            {
                // tag corners
                std::vector<Vertex_handle> corner(pixel.nb_corners());
                for (unsigned k = 0; k < corner.size(); ++k)
                {
                    Point q = pixel.get_corner(k);
                    vertex = m_rt.find_nearest_vertex(q, vertex);
                    corner[k] = vertex;
                }
                split_pixel(pixel, corner, grid.get_enriched_segments(i, j));
            }
        }
    }
    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_YELLOW);

    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_GREEN);
}

FT CCVT::rasterize(const EnrichedSegment& enriched_segment, Grid& grid)
{
    const Segment& segment = enriched_segment.segment();
    Vertex_handle left  = enriched_segment.left();
    Vertex_handle right = enriched_segment.right();

    const Point& p0 = segment.source();
    const Point& p1 = segment.target();
    Vector velocity = p1 - p0;

    unsigned i0, j0, i1, j1;
    m_domain.locate(p0, i0, j0);
    m_domain.locate(p1, i1, j1);

    FT sum = 0.0;
    Point p = p0;
    unsigned i = i0;
    unsigned j = j0;
    while (i != i1 || j != j1)
    {
        // run one step of Bresenham
        Point q;
        unsigned u;
        unsigned v;
        bool ok = move(i, j, p, velocity, u, v, q);
        if (!ok)
        {
            // never come here
            std::cout << "move_pixel failed" << std::endl;
            return 0.0;
        }

        // process pixel
        EnrichedSegment es(Segment(p, q), left, right);
        grid.append(i, j, es);
        FT dist = std::sqrt(CGAL::squared_distance(p, q));
        sum += m_domain.get_value(i, j)*dist;

        // update for next iter
        p = q;
        i = u;
        j = v;
    }
    if (p == p1) return sum;

    // process last pixel
    EnrichedSegment es(Segment(p, p1), left, right);
    grid.append(i, j, es);
    FT dist = std::sqrt(CGAL::squared_distance(p, p1));
    sum += m_domain.get_value(i, j)*dist;
    return sum;
}

bool CCVT::move(const unsigned i, const unsigned j, const Point& source,
                const Vector& velocity,
                unsigned& u, unsigned& v, Point& target)
{
    FT xdelta = velocity.x();
    FT ydelta = velocity.y();
    if (std::abs(xdelta) < EPS && std::abs(ydelta) < EPS) return false;

    if (std::abs(ydelta) < EPS) return move_horizontal(i, j, source, velocity, u, v, target);
    if (std::abs(xdelta) < EPS) return move_vertical  (i, j, source, velocity, u, v, target);

    unsigned iprev = (i == 0)? 0 : i - 1;
    unsigned jprev = (j == 0)? 0 : j - 1;
    unsigned inext = std::min(i + 1, m_domain.get_width());
    unsigned jnext = std::min(j + 1, m_domain.get_height());

    FT xwall = m_domain.compute_x(i);
    if (xdelta > 0.0) xwall = m_domain.compute_x(inext);

    FT ywall = m_domain.compute_y(j);
    if (ydelta < 0.0) ywall = m_domain.compute_y(jnext);

    FT xtime = (xwall - source.x()) / xdelta;
    FT ytime = (ywall - source.y()) / ydelta;
    FT time = std::min(xtime, ytime);
    target = source + time*velocity;

    u = i;
    if (std::abs(xtime - time) < EPS)
    {
        if (xdelta > 0.0) u = inext;
        else              u = iprev;
    }
    v = j;
    if (std::abs(ytime - time) < EPS)
    {
        if (ydelta > 0.0) v = jprev;
        else              v = jnext;
    }

    if (i == u && j == v) return false;
    return true;
}

bool CCVT::move_horizontal(const unsigned i, const unsigned j,
                            const Point& source, const Vector& velocity,
                            unsigned& u, unsigned& v, Point& target)
{
    unsigned iprev = (i == 0)? 0 : i - 1;
    unsigned inext = std::min(i + 1, m_domain.get_width ());

    FT wall = m_domain.compute_x(i);
    if (velocity.x() > 0.0) wall = m_domain.compute_x(inext);

    FT time = (wall - source.x()) / velocity.x();
    target = source + time*velocity;

    v = j;
    if (velocity.x() > 0.0) u = inext;
    else                    u = iprev;

    if (i == u && j == v) return false;
    return true;
}

bool CCVT::move_vertical(const unsigned i, const unsigned j,
                          const Point& source, const Vector& velocity,
                          unsigned& u, unsigned& v, Point& target)
{
    unsigned jprev = (j == 0)? 0 : j - 1;
    unsigned jnext = std::min(j + 1, m_domain.get_height());

    FT wall = m_domain.compute_y(j);
    if (velocity.y() < 0.0) wall = m_domain.compute_y(jnext);

    FT time = (wall - source.y()) / velocity.y();
    target = source + time*velocity;

    u = i;
    if (velocity.y() > 0.0) v = jprev;
    else                    v = jnext;

    if (i == u && j == v) return false;
    return true;
}

void CCVT::split_pixel(const Pixel& original_pixel,
                        const std::vector<Vertex_handle>& corner_tags,
                        const std::vector<EnrichedSegment>& enriched_segments)
{
    std::map< Vertex_handle, std::vector<Point> > table;
    for (unsigned i = 0; i < original_pixel.nb_corners(); ++i)
    {
        const Point&  pi = original_pixel.get_corner(i);
        Vertex_handle vi = corner_tags[i];
        append_point_to_vertex(table, pi, vi);
    }

    for (unsigned i = 0; i < enriched_segments.size(); ++i)
    {
        const EnrichedSegment& es = enriched_segments[i];
        const Segment& segment = es.segment();
        const Point& ps = segment.source();
        const Point& pt = segment.target();
        Vertex_handle left  = es.left ();
        Vertex_handle right = es.right();

        append_point_to_vertex(table, ps, left );
        append_point_to_vertex(table, ps, right);

        append_point_to_vertex(table, pt, left );
        append_point_to_vertex(table, pt, right);
    }

    FT intensity = original_pixel.get_value();
    for (std::map< Vertex_handle, std::vector<Point> >::const_iterator
                 it = table.begin(); it != table.end(); ++it)
    {
        FT intensity = original_pixel.get_value();
        Vertex_handle vertex = it->first;
        const std::vector<Point>& points = it->second;

        Pixel pixel(intensity);
        pixel.compute_shape(points.begin(), points.end());
        vertex->append_pixel(pixel);
    }

    /*
    std::set<Vertex_handle> vertices;
    vertices.insert(corner_tags.begin(), corner_tags.end());
    for (unsigned i = 0; i < enriched_segments.size(); ++i)
    {
        const EnrichedSegment& es = enriched_segments[i];
        Vertex_handle left  = es.left ();
        Vertex_handle right = es.right();
        vertices.insert(left);
        vertices.insert(right);
    }

    FT value = original_pixel.get_value();
    for (std::set<Vertex_handle>::const_iterator
         vi = vertices.begin(); vi != vertices.end(); ++vi)
    {
        Vertex_handle vertex = *vi;

        ConvexPolygon polygon;
        m_rt.build_polygon(vertex, polygon.points());
        polygon.convex_clip(original_pixel.get_shape());

        Pixel pixel(polygon, value);
        vertex->append_pixel(pixel);
    }
    */
}

void CCVT::append_point_to_vertex(std::map< Vertex_handle, std::vector<Point> >& table,
                                   const Point& point, Vertex_handle vertex) const
{
    std::map< Vertex_handle, std::vector<Point> >::iterator it = table.find(vertex);
    if (it == table.end())
    {
        std::vector<Point> points;
        points.push_back(point);
        table[vertex] = points;
        return;
    }
    it->second.push_back(point);
}
