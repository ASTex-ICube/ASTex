//
// Created by grenier on 09/10/23.
//

#ifndef CCVT_TEST_CHARLINE_CCVT_H
#define CCVT_TEST_CHARLINE_CCVT_H


// STL
#include <map>
#include <vector>

// local
#include "matrix/sparse_matrix.h"
#include "types.h"


class CCVT
{
public:

private:
    RT m_rt;
    Domain m_domain;
    std::vector<FT> m_capacities; // vecteur des capcité des cellules (objectif de répartition de la densité totale du domaine)
    std::vector<FT> m_proportions; // proportions d'utilisation de la densité pour chaque cellule
    bool m_custom_proportions;

    double m_tau;
    std::map<Edge, FT> m_ratio; // les ratio des arrête du dual ?
    std::vector<double> m_r, m_g, m_b;
    std::vector<Vertex_handle> m_vertices; // vecteur des position des graines et les poids associés
    // pour chaque Vertex_handle on a
    // m_ptr : informations sur la graine
    //     _f : face ?
    //          V : vertex
    //          N : neighbors
    //     _p : position ?
    //     m_dual : les extrémité des arrêtes de la cellule ?
    //     m_pixels : les valeurs des pixel de la densité couvert par la cellule ?
    //     m_area : le volume de la densité présente dans la cellule ?

    bool m_timer_on;
    std::vector<double> m_timer;
    bool m_fixed_connectivity;



public:
    CCVT()
    {
//        srand(24);
        m_tau = 1.0;
        m_timer_on = false;
        m_fixed_connectivity = false;
        m_custom_proportions = false;
    }

    CCVT(unsigned seed)
    {
        srand(seed);
        m_tau = 1.0;
        m_timer_on = false;
        m_fixed_connectivity = false;
        m_custom_proportions = false;
    }

    ~CCVT()
    {
        clear();
    }

    double get_tau() const { return m_tau; }
    void set_tau(double tau) { m_tau = tau; }

    void set_custom_proportions(std::vector<FT>& proportions){
        m_proportions = proportions;
        m_custom_proportions = true;
    }

    std::vector<FT> get_capacities(){ return m_capacities; } // TODO : à retirer
    std::vector<FT> get_area(){
        std::vector<FT> areas;
        for(auto vi=m_vertices.begin(); vi<m_vertices.end(); vi++){
            areas.push_back((*vi)->compute_area());
        }
        return areas;
    }

    void toggle_invert() { m_domain.toggle_invert(); }

    void toggle_timer() { m_timer_on = !m_timer_on; }

    void toggle_connectivity() { m_fixed_connectivity = !m_fixed_connectivity; }

    bool connectivity_fixed() const { return m_fixed_connectivity; }

    void clear()
    {
        clear_triangulation();
    }


    void set_domain(double moy_x, double moy_y, double variance_x, double variance_y){
        bool ok = m_domain.set(moy_x, moy_y, variance_x, variance_y);
        if (!ok) return;

        m_rt.set_boundary(m_domain.get_dx(),
                          m_domain.get_dy());
//    std::cout << "Dx vs Dy: " << m_domain.get_dx() << " ; " << m_domain.get_dy() << std::endl;
    }




    // IO //

    void load_image(const std::string& filename);
//
//    void load_points(const QString& filename);
//
//    void load_dat(const QString& filename, std::vector<Point>& points) const;
//
//    void save_points(const QString& filename) const;
//
//    void save_dat(const QString& filename, const std::vector<Point>& points) const;
//
//    void save_txt(const QString& filename, const std::vector<Point>& points) const;
//
    void save_eps(const std::string& filename) const;




    // SITES //

    void generate_random_sites(const unsigned nb);

    void set_initial_sites(std::vector<Point> sites);

    void generate_random_sites_based_on_image(const unsigned nb);

    void generate_regular_grid(const unsigned nx, const unsigned ny);

    void init_colors(const unsigned nb);




    // RENDER //

//    void draw_point(const Point& a) const;
//
//    void draw_segment(const Point& a, const Point& b) const;
//
//    void draw_triangle(const Point& a, const Point& b, const Point& c) const;
//
//    void draw_polygon(const std::vector<Point>& polygon) const;
//
//    void draw_circle(const Point& center,
//                     const FT scale,
//                     const std::vector<Point>& pts) const;
//
//    void draw_image() const;
//
//    void draw_image_grid() const;
//
//    void draw_domain(const float line_width,
//                     const float red,
//                     const float green,
//                     const float blue) const;
//
//    void draw_sites(const float point_size,
//                    const float red,
//                    const float green,
//                    const float blue) const;
//
//    void draw_vertices(const float point_size) const;
//
//    void draw_faces(const float red,
//                    const float green,
//                    const float blue) const;
//
//    void draw_primal(const float line_width,
//                     const float red,
//                     const float green,
//                     const float blue) const;
//
//    void draw_dual(const float line_width,
//                   const float red,
//                   const float green,
//                   const float blue) const;
//
//    void draw_bounded_dual(const float line_width,
//                           const float red,
//                           const float green,
//                           const float blue) const;
//
//    void draw_weights() const;
//
//    void draw_pixels() const;
//
//    void draw_capacity() const;
//
//    void draw_variance() const;
//
//    void draw_regularity() const;
//
//    void draw_regular_sites() const;
//
//    void draw_cell(Vertex_handle vertex, bool filled) const;
//
//    void draw_barycenter(const float point_size,
//                         const float red,
//                         const float green,
//                         const float blue) const;
//
//    void draw_capacity_histogram(const unsigned nbins,
//                                 const double xmin,
//                                 const double xmax,
//                                 const double ymin,
//                                 const double ymax) const;
//
//    void draw_weight_histogram(const double range,
//                               const unsigned nbins,
//                               const double xmin,
//                               const double xmax,
//                               const double ymin,
//                               const double ymax) const;
//
//    void draw_histogram(const std::vector<unsigned>& histogram,
//                        const double xmin,
//                        const double xmax,
//                        const double ymin,
//                        const double ymax) const;
//
//    void setViewer(GlViewer* vi) { viewer = vi; }




    // HISTOGRAM //

    void compute_capacity_histogram(std::vector<unsigned>& histogram) const;

    void compute_weight_histogram(const double range, std::vector<unsigned>& histogram) const;





    // INIT //

    bool is_valid() const;

    FT compute_mean_capacity() const;

    unsigned count_visible_sites() const;

    void collect_visible_points(std::vector<Point>& points) const;

    void collect_visible_weights(std::vector<FT>& weights) const;

    void collect_sites(std::vector<Point>& points,
                       std::vector<FT>& weights) const;

    void clear_triangulation();

    bool update_triangulation(bool skip = false);

    bool construct_triangulation(const std::vector<Point>& points,
                                 const std::vector<FT>& weights,
                                 bool skip = false);

    bool populate_vertices(const std::vector<Point>& points,
                           const std::vector<FT>& weights);

    Vertex_handle insert_vertex(const Point& point,
                                const FT weight,
                                const unsigned index);

    void compute_capacities(std::vector<FT>& capacities) const;

    void update_positions(const std::vector<Point>& points,
                          bool clamp = true,
                          bool hidden = true);

    void update_weights(const std::vector<FT>& weights,
                        bool hidden = true);

    void reset_weights();

    FT compute_value_integral() const;

    void pre_build_dual_cells();

    void pre_compute_area();





    // ENERGY //

    FT compute_wcvt_energy();

    void compute_weight_gradient(std::vector<FT>& gradient, FT coef = 1.0);

    void compute_position_gradient(std::vector<Vector>& gradient, FT coef = 1.0);

    void compute_neightbour_gradient(std::vector<Vector>& gradient, FT coef = 1.0);

    FT compute_weight_threshold(FT epsilon) const;

    FT compute_position_threshold(FT epsilon) const;





    // OPTIMIZER //

    FT optimize_positions_via_lloyd(bool update);

    FT optimize_positions_via_gradient_ascent(FT& timestep, bool update);

    FT optimize_neightbour_via_gradient_descent(FT& timestep, bool update);

    FT optimize_weights_via_newton(FT& timestep, bool update);

    unsigned optimize_neightbour_via_gradient_descent_until_converge(FT& timestep,
                                                                  FT threshold,
                                                                  unsigned update,
                                                                  unsigned max_iters);

    unsigned optimize_weights_via_newton_until_converge(FT& timestep,
                                                        FT epsilon,
                                                        unsigned update,
                                                        unsigned max_iters);


    unsigned optimize_all(FT& wstep, FT& xstep, unsigned max_newton_iters,
                          FT epsilon, unsigned max_iters,
                          std::ostream& out);

    bool solve_newton_step(const std::vector<FT>& b,
                           std::vector<FT>& x);

    void build_laplacian(const FT scale,
                         const std::map<unsigned, unsigned>& indices,
                         SparseMatrix& A) const;

    bool solve_linear_system(const SparseMatrix& A,
                             std::vector<double>& x,
                             const std::vector<double>& b) const;

    unsigned optimize_H(FT& wstep, FT& xstep, unsigned max_newton_iters,
                        FT epsilon, unsigned max_iters);




    // ASSIGN //

    Pixel build_pixel(unsigned i, unsigned j) const;

    void set_ratio(Edge edge, FT value);

    FT get_ratio(Edge edge) const;

    void clean_pixels();

    void assign_pixels();

    FT rasterize(const EnrichedSegment& enriched_segment, Grid& grid);

    bool move(const unsigned i, const unsigned j,
              const Point& source, const Vector& velocity,
              unsigned& u, unsigned& v, Point& target);

    bool move_horizontal(const unsigned i, const unsigned j,
                         const Point& source, const Vector& velocity,
                         unsigned& u, unsigned& v, Point& target);

    bool move_vertical(const unsigned i, const unsigned j,
                       const Point& source, const Vector& velocity,
                       unsigned& u, unsigned& v, Point& target);

    void split_pixel(const Pixel& original_pixel,
                     const std::vector<Vertex_handle>& corner_tags,
                     const std::vector<EnrichedSegment>& enriched_segments);

    void append_point_to_vertex(std::map< Vertex_handle, std::vector<Point> >& table,
                                const Point& point, Vertex_handle vertex) const;





    // REGULARITY //

//    void detect_and_break_regularity(FT max_radius, unsigned max_teleport);
//
//    FT compute_regularity_threshold() const;
//
//    FT compute_max_regularity() const;
//
//    void compute_variance_vector(std::vector<FT>& variance) const;
//
//    void compute_regularity_vector(const std::vector<FT>& variance,
//                                   std::vector<FT>& regularity) const;
//
//    FT compute_regularity(Vertex_handle vi, const std::vector<FT>& variance) const;
//
//    void jitter_vertices(const std::set<Vertex_handle>& vertices, const FT max_radius);
//
//    Point jitter_point(const Point& p, const FT max_radius) const;
//
//    void count_sites_per_bin(unsigned N) const;


};




#endif //CCVT_TEST_CHARLINE_CCVT_H
