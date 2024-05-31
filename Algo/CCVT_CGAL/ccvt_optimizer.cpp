//
// Created by grenier on 12/10/23.
//

#include "ccvt.h"
#include "timer.h"

#include "pw_line_search.h"

#include "matrix/sparse_array.h"
#include "matrix/suite_sparse_qr.h"

typedef CLSWeights<CCVT, FT> LSWeights;
typedef CLSPositions<CCVT, Point, Vector> LSPositions;


FT CCVT::optimize_positions_via_lloyd(bool update)
{
    if(m_verbose){std::cout<<"optimizing positions..."<<std::endl;}
    if (m_timer_on) Timer::start_timer(m_timer, COLOR_BLUE, "Centroid");
    std::vector<Point> points;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        Point ci = vi->compute_centroid(); // barycentre des cellules
        points.push_back(ci);
    }
    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);

    update_positions(points);
    if (update) update_triangulation(); // on place les graines au centre des cellules

    std::vector<Vector> gradient;
    compute_position_gradient(gradient);
    return compute_norm(gradient);
}

FT CCVT::optimize_positions_via_gradient_ascent(FT& timestep, bool update)
{
    if(m_verbose){std::cout<<"optimizing positions..."<<std::endl;}
    std::vector<Point> points;
    collect_visible_points(points);

    std::vector<Vector> gradient;
    compute_position_gradient(gradient);

    if (timestep <= 0.0) // TODO : ça marche ça ?
    {
        double mean_capacity = compute_mean(m_capacities);
        double max_alpha = 1.0 / mean_capacity;
        LSPositions line_search(this, 10, max_alpha);
        timestep = line_search.run_bt(points, gradient);
    }
    else {
        for (unsigned i = 0; i < points.size(); ++i)
        {
            Point  pi = points[i];
            Vector gi = gradient[i];
            points[i] = pi + timestep*gi; // + ou - ??
        }
        update_positions(points);
        if (update) update_triangulation();
    }

    compute_position_gradient(gradient);
    return compute_norm(gradient);
}




// pas utilisé
FT CCVT::optimize_neightbour_via_gradient_descent(FT& timestep, bool update) // TODO
{
    if(m_verbose){std::cout<<"optimizing positions..."<<std::endl;}
    std::vector<Point> points;
    collect_visible_points(points);

    std::vector<Vector> gradient;
    compute_neightbour_gradient(gradient);

    if (timestep <= 0.0)  // TODO n'a pas l'aire de marcher ...
    {
        double mean_capacity = compute_mean(m_capacities);
        double max_alpha = 1.0 / mean_capacity;
        LSPositions line_search(this, 10, max_alpha);
        timestep = line_search.run_bt(points, gradient);
    }
    else {
        for (unsigned i = 0; i < points.size(); ++i)
        {
            Point  pi = points[i];
            Vector gi = gradient[i];
            points[i] = pi + timestep*gi; // + ou - ??
        }
        update_positions(points);
        if (update) update_triangulation();
    }


    compute_neightbour_gradient(gradient);
    return compute_norm(gradient);
}





FT CCVT::optimize_neightbour(FT& timestep, bool update) // TODO
{
    if(m_verbose){std::cout<<"optimizing positions..."<<std::endl;}
    std::vector<Point> points;

    std::vector<std::vector<FT>> current_neightbour = get_neightbour_val();

    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        double norm  =  0.;
        double facteur_grad = 0.;
        Vector grad{0., 0.};

        Vertex_handle vi = m_vertices[i]; // x_i
        if (vi->is_hidden()) continue;

        Edge_circulator ecirc = m_rt.incident_edges(vi); // liste des eij
        Edge_circulator eend  = ecirc;

        CGAL_For_all(ecirc, eend)   // for j in Omega_i
        {
            Edge edge = *ecirc; // e_ij
            if (!m_rt.is_inside(edge)) continue;

            // position graine x_j
            Vertex_handle vj = m_rt.get_source(edge); // x_j
            if (vj == vi) vj = m_rt.get_target(edge);
            unsigned j = vj->get_index();

            Segment dual = m_rt.build_bounded_dual_edge(edge); // extrémités de e*ij

            double prop_m_ij_obj = m_neightbour_proportions.at(i).at(j);//*compute_value_integral();
            double sum_m_ij_cur = std::accumulate(current_neightbour.at(i).begin(), current_neightbour.at(i).end(), 0.);// current_neightbour.at(i).at(j);//*m_domain.get_max_value();

            double m_ij_obj = current_neightbour.at(i).at(j);// prop_m_ij_obj*sum_m_ij_cur; // TODO utiliser les m_ij obj
            double n_eij = m_rt.get_length(edge); // |e_ij|
            double n_eij_star = m_rt.get_length(dual); // |e*_ij|
            double wi = vi->get_weight();
            double wj = vj->get_weight();
            double rho_xi = m_domain.get_value(vi->get_position(),true);
            double d_ij = (n_eij*n_eij + wi - wj)/(2.*n_eij);

//            double grad_dij = (wi-wj-n_eij*n_eij)/(n_eij*n_eij*n_eij);
//            double facteur = (2.*m_ij_cur*(n_eij*n_eij+wi-wj) + m_ij_obj*(n_eij*n_eij+wj-wi))/(4.*n_eij*n_eij*n_eij);
//            double facteur =  m_ij_obj/n_eij;
            double facteur = (m_ij_obj + rho_xi*n_eij_star)*(wi - wj - n_eij*n_eij)/(2.*n_eij*n_eij*n_eij);

            facteur_grad += n_eij_star*d_ij;
            norm += facteur;
            grad += facteur*Vector{Point{0,0}, vj->get_position()};
        }
        double vol = vi->compute_area()/m_domain.integrate_intensity();
        double rho_xi = m_domain.get_value(vi->get_position(),true);
        Vector grad_rho = rho_xi*Vector{-(vi->get_position().x()-m_domain.get_mu_x())/(m_domain.get_var_x()),
                                        -(vi->get_position().y()-m_domain.get_mu_y())/(m_domain.get_var_y())};

        Point ci = Point{0,0} + (1./norm)*(grad + facteur_grad*grad_rho);
//        Point ci = Point{0,0} + (1./norm)*(grad - vol*grad_rho);
        if(i==5){
            points.push_back(ci);
        } else{
            points.push_back(vi->get_position());
        }
//        points.push_back(ci);



    }
    update_positions(points);
    if (update) update_triangulation();

    return 1.;
}




// contenue de la subroutine Enforce-Capacity-Constraints
FT CCVT::optimize_weights_via_newton(FT& timestep, bool update)
{
    std::vector<FT> gradient;
    compute_weight_gradient(gradient, -1.0); // vecteur des m - mi (capacité obj - aire actuelle)

    std::vector<FT> direction;
    bool ok = solve_newton_step(gradient, direction); // solve for delta in eq 4
    if (!ok) return 0.0;

    std::vector<FT> weights;
    collect_visible_weights(weights);

    if (timestep <= 0.0) // TODO n'a pas l'aire de marcher ...
    {
        // ça marcherait avec un timestep = 1/||direction||^2 ?
        LSWeights line_search(this, 20, 2.0); // (ccvt, max_iter, max_alpha)
        timestep = line_search.run_bt(weights, direction); // finding alpha (timestep) satisfying Armijo condition
//        std::cout<<timestep<<std::endl;
    }
    else {
        for (unsigned i = 0; i < weights.size(); ++i)
        {
            FT wi = weights[i];
            FT gi = direction[i]; // delta
            weights[i] = wi + timestep*gi; // W <- W + alpha * delta
        }
        update_weights(weights);
        if (update) update_triangulation();
    }

    compute_weight_gradient(gradient); // nabla_w F
    return compute_norm(gradient); // ||nabla_w F||
}




bool CCVT::solve_newton_step(const std::vector<FT>& b, std::vector<FT>& x) // solve for delta in eq 4
{
    if (m_timer_on) Timer::start_timer(m_timer, COLOR_BLUE, "LinearSolver");

    unsigned nb = 0;
    std::map<unsigned, unsigned> indices;
    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;
        indices[vi->get_index()] = nb++;
    }

    SparseMatrix L(nb, nb);
    build_laplacian(0.5, indices, L);

    bool ok = solve_linear_system(L, x, b); // x : vecteur delta
    if (!ok)
    {
        std::cout << red << "linear solver failed" << white << std::endl;
        return false;
    }

    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);
    return true;
}

void CCVT::build_laplacian(const FT scale, // build matrix Delta^{w,rho}
                            const std::map<unsigned, unsigned>& indices,
                            SparseMatrix& A) const
{
    unsigned nb = A.numRows();
    for (unsigned k = 0; k < m_vertices.size(); ++k)
    {
        Vertex_handle vi = m_vertices[k];
        if (vi->is_hidden()) continue;
        unsigned i = indices.find(vi->get_index())->second;

        double diagi = 0.0;
        SparseArray rowi(nb);
        Edge_circulator ecirc = m_rt.incident_edges(vi); // liste des eij ?
        Edge_circulator eend  = ecirc;
        CGAL_For_all(ecirc, eend)
        {
            Edge edge = *ecirc; // e_ij
            if (!m_rt.is_inside(edge)) continue;

            // position graine vj
            Vertex_handle vj = m_rt.get_source(edge);
            if (vj == vi) vj = m_rt.get_target(edge);

            unsigned j = vj->get_index();
            j = indices.find(j)->second;

            double coef = scale * get_ratio(edge); // get_ratio = rho*|e*ij|/|eij|, scale = 0.5 // TODO : calcul de int de rho
            if (std::abs(coef) < EPS) continue;

            rowi.setValue(j, -coef);
            diagi += coef;
        }

        rowi.setValue(i, diagi);
        A.setRow(i, rowi);
    }
}

bool CCVT::solve_linear_system(const SparseMatrix& A,
                                std::vector<double>& x,
                                const std::vector<double>& b) const
{
    SuiteSparseQRFactorizer solver;
    bool ok = solver.factorize(A);
    if (!ok) return false;

    ok = solver.solve(b, x);
    return ok;
}

unsigned CCVT::optimize_neightbour_via_gradient_descent_until_converge(FT& timestep, // TODO ?
                                                                     FT threshold,
                                                                     unsigned update,
                                                                     unsigned max_iters)
{
    for (unsigned i = 0; i < max_iters; ++i)
    {
        bool flag = (update == 0 || (i+1) % update == 0);
        FT norm = optimize_neightbour_via_gradient_descent(timestep, flag);
        if (norm < threshold) return i;
    }
    return max_iters;
}

// subroutine Enforce-Capacity-Constraints
unsigned CCVT::optimize_weights_via_newton_until_converge(FT& timestep,
                                                           FT threshold,
                                                           unsigned update,
                                                           unsigned max_iters)
{
    if(m_verbose){std::cout<<"optimizing weights..."<<std::endl;}
    for (unsigned i = 0; i < max_iters; ++i) // boucle 28 à 33
    {
        bool flag = (update == 0 || (i+1) % update == 0);
        FT norm = optimize_weights_via_newton(timestep, flag); // ||nabla_w F||
        if (norm < threshold) return i;
    }
//    std::cout<<"newton_max_iter reached"<<std::endl;
    return max_iters;
}









unsigned CCVT::optimize_all(FT& wstep, FT& xstep, unsigned max_newton_iters,
                             FT epsilon, unsigned max_iters,
                             std::ostream& out)
{
    bool global_connectivity = m_fixed_connectivity;
    unsigned nb0 = count_visible_sites();

    FT xthreshold = compute_position_threshold(epsilon);
    FT wthreshold = compute_weight_threshold(epsilon);

    out << "NbSites = " << nb0 << std::endl;
    out << "Threshold: " << xthreshold << " ; " << wthreshold << std::endl;

    m_fixed_connectivity = false;
    FT coarse_xthreshold = 2.0*xthreshold;
    FT coarse_wthreshold = 2.0*wthreshold;

    unsigned iters = 0;
    unsigned nb_assign = 0;




    // optimisation "grossière"
    while (iters < max_iters)
    {
        iters++;
        reset_weights();

        // on ajuste les poids pour coller au capacités (Newton method for W)
        nb_assign += optimize_weights_via_newton_until_converge(wstep, coarse_wthreshold, 0, max_newton_iters);

        // on replace les graine au milieu des cellules (Lloyd step for X)
        FT norm = optimize_positions_via_lloyd(true);

        nb_assign++;
        out << "(Coarse) Norm: " << norm << std::endl;
        if (norm <= coarse_xthreshold) break;
    }

    out << "Partial: " << iters << " iters" << std::endl;
    m_fixed_connectivity = global_connectivity;
    if (iters == max_iters) return iters;

    m_fixed_connectivity = false;
    FT fine_xthreshold = xthreshold;
    FT fine_wthreshold = wthreshold;





    // optimisation "fine"
    while (iters < max_iters)
    {
        iters++;
        unsigned nb1 = count_visible_sites();
        if (nb1 != nb0) reset_weights();

        // on ajuste les poids pour coller au capacités (Newton method for W)
        nb_assign += optimize_weights_via_newton_until_converge(wstep, fine_wthreshold, 0, max_newton_iters);

        // on replace les graine au milieu des cellules (gradient descent for X)
        FT norm = optimize_positions_via_gradient_ascent(xstep, true);

        nb_assign++;
        out << "(Fine) Norm: " << norm << std::endl;
        if (norm <= fine_xthreshold) break;
    }



    // dernière optimisation des volumes
    optimize_weights_via_newton_until_converge(wstep, 0.1*fine_wthreshold, 0, max_newton_iters);

    std::cout << "NbAssign: " << nb_assign << std::endl;

    m_fixed_connectivity = global_connectivity;
    return nb_assign;//iters;
}









unsigned CCVT::optimize_H(FT& wstep, FT& xstep, unsigned max_newton_iters, FT epsilon, unsigned max_iters){
    bool global_connectivity = m_fixed_connectivity;
    unsigned nb0 = count_visible_sites();

    FT xthreshold = compute_position_threshold(epsilon);
    FT wthreshold = compute_weight_threshold(epsilon);

    m_fixed_connectivity = false;
    unsigned iters = 0;
    unsigned nb_assign = 0;


    while (iters < max_iters)
    {
        FT norm = 0.;
        iters++;
        unsigned nb1 = count_visible_sites();
        if (nb1 != nb0) reset_weights();

        // on ajuste les poids pour coller au capacités (Newton method for W)
        nb_assign += optimize_weights_via_newton_until_converge(wstep, wthreshold, 0, max_newton_iters);
        verbose();

        // on replace les graine au milieu des cellules (Lloyd step for X)
        norm = optimize_positions_via_lloyd(true);
        verbose();


        nb_assign++;
        if (norm <= xthreshold) break;
    }

    // dernière optimisation des volumes
    optimize_weights_via_newton_until_converge(wstep, wthreshold, 0, max_newton_iters);
    verbose();


    // on ajuste les positions pour coller au voisinages
//    optimize_neightbour(xstep, true); // TODO
//    verbose();


    m_fixed_connectivity = global_connectivity;
    return nb_assign*iters;

}