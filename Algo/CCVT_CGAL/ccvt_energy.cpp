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
        Vector grad_mij{0., 0.};

        Vertex_handle vi = m_vertices[i]; // x_i
        if (vi->is_hidden()) continue;

        Edge_circulator ecirc = m_rt.incident_edges(vi); // liste des eij
        Edge_circulator eend  = ecirc;

        CGAL_For_all(ecirc, eend)
        {
            Vector int_grad_rho{0.,0.};
            Vector rho_grad_cijl{0.,0.};
            Vector rho_grad_cijk{0.,0.};

            Edge edge = *ecirc; // e_ij
            if (!m_rt.is_inside(edge)) continue;

            // position graine x_j
            Vertex_handle vj = m_rt.get_source(edge); // x_j
            if (vj == vi) vj = m_rt.get_target(edge);
            unsigned j = vj->get_index();

            // segment dual
            Segment dual = m_rt.build_bounded_dual_edge(edge); // extrémités de e*ij
            Point left_cw = dual.target(); // c_ijl
            Point right_cw = dual.source(); // c_ijk

            // triangle jik et jil
            Edge twin = m_rt.get_twin(edge); // e_ji
            Face_handle left_face  = edge.first; // T_ijl
            Face_handle right_face = twin.first; // T_ijk
            bool left_inside  = m_rt.is_inside( left_face);
            bool right_inside = m_rt.is_inside(right_face);

            if (left_inside) // si x_l existe
            {
                Vertex_handle vl = left_face->vertex(edge.second); // x_l
                double rho_cijl = m_domain.get_value(left_cw); // rho(c_ijl)

                double T = m_rt.get_area(left_face); // aire du triangle ijl
                Edge eil = m_rt.get_next(edge); // e_il
                Edge ejl = m_rt.get_twin(m_rt.get_next(eil)); // e_jl
                Vector p_ejl = m_rt.get_orthogonal_vector(ejl); // e_jl^perp

                double n_eil = m_rt.get_length(eil);  // |e_il|
                double n_eij = m_rt.get_length(edge); // |e_ij|
                double n_ejl = m_rt.get_length(ejl);  // |e_jl|

                double theta_j = std::acos((n_eij*n_eij + n_ejl*n_ejl - n_eil*n_eil)/(2.*n_ejl*n_eij)); // par Al-Kashi
                double theta_l = std::acos((n_ejl*n_ejl + n_eil*n_eil - n_eij*n_eij)/(2.*n_ejl*n_eil));

                Vector grad_hi = (- n_eil*n_eil*(cos(theta_j)/ sin(theta_j))*Vector{vi->get_position(),vj->get_position()}
                                  - n_eij*n_eij*(cos(theta_l)/ sin(theta_l))*Vector{vi->get_position(),vl->get_position()}
                                  + vj->get_weight()*Vector{vl->get_position(),vi->get_position()}.perpendicular(CGAL::CLOCKWISE)
                                  + vl->get_weight()*Vector{vi->get_position(),vj->get_position()}.perpendicular(CGAL::CLOCKWISE)
                                  + vi->get_weight()*Vector{vj->get_position(),vl->get_position()}.perpendicular(CGAL::CLOCKWISE))
                                        * n_ejl/(8.*T*T);


                // produit tensoriel et évaluation en xi
                Vector grad_cijl =  (1./n_ejl)* Vector{grad_hi.x()*p_ejl.x()*vi->get_position().x() + grad_hi.x()*p_ejl.y()*vi->get_position().y(),
                                                       grad_hi.y()*p_ejl.x()*vi->get_position().x() + grad_hi.y()*p_ejl.y()*vi->get_position().y()};

                rho_grad_cijl = rho_cijl*grad_cijl;
            }

            if (right_inside) // si x_k existe
            {
                Vertex_handle vk = right_face->vertex(twin.second); // x_k
                double rho_cijk = m_domain.get_value(right_cw); // rho(c_ijk)

                double T = m_rt.get_area(right_face); // aire du triangle ijk
                Edge ejk = m_rt.get_next(twin); // e_jk
                Edge eik = m_rt.get_twin(m_rt.get_next(ejk)); // e_ik
                Vector p_ejk = m_rt.get_orthogonal_vector(ejk); // e_jk^perp

                double n_eik = m_rt.get_length(eik);  // |e_ik|
                double n_eij = m_rt.get_length(edge); // |e_ij|
                double n_ejk = m_rt.get_length(ejk);  // |e_jk|

                double theta_j = std::acos((n_eij*n_eij + n_ejk*n_ejk - n_eik*n_eik)/(2.*n_ejk*n_eij)); // par Al-Kashi
                double theta_k = std::acos((n_ejk*n_ejk + n_eik*n_eik - n_eij*n_eij)/(2.*n_ejk*n_eik));

                Vector grad_hi = (- n_eik*n_eik*(cos(theta_j)/ sin(theta_j))*Vector{vi->get_position(),vj->get_position()}
                                  - n_eij*n_eij*(cos(theta_k)/ sin(theta_k))*Vector{vi->get_position(),vk->get_position()}
                                  + vj->get_weight()*Vector{vk->get_position(),vi->get_position()}.perpendicular(CGAL::CLOCKWISE)
                                  + vk->get_weight()*Vector{vi->get_position(),vj->get_position()}.perpendicular(CGAL::CLOCKWISE)
                                  + vi->get_weight()*Vector{vj->get_position(),vk->get_position()}.perpendicular(CGAL::CLOCKWISE))
                                 * n_ejk/(8.*T*T);

                // produit tensoriel et évaluation en xi
                Vector grad_cijk =  (1./n_ejk)* Vector{grad_hi.x()*p_ejk.x()*vi->get_position().x() + grad_hi.x()*p_ejk.y()*vi->get_position().y(),
                                                       grad_hi.y()*p_ejk.x()*vi->get_position().x() + grad_hi.y()*p_ejk.y()*vi->get_position().y()};

                rho_grad_cijk = rho_cijk*grad_cijk;
            }

            // calculer int_grad_rho
            double n_eij_star = m_rt.get_length(dual); // |e*_ij|

            // constantes du produit des gaussiennes
            double a = left_cw.x()-right_cw.x();
            double b = left_cw.y()-right_cw.y();
            double mu_1 = m_domain.get_mu_x()-m_domain.get_dx() - right_cw.x();
            double mu_2 = m_domain.get_mu_y()-m_domain.get_dy() - right_cw.y();

            double A = product_gaussian_amplitude(a, b, mu_1, mu_2, m_domain.get_sig_x(), m_domain.get_sig_y());
            double mu = product_gaussian_mean(a, b, mu_1, mu_2, m_domain.get_sig_x(), m_domain.get_sig_y());
            double var = product_gaussian_variance(a, b, mu_1, mu_2, m_domain.get_sig_x(), m_domain.get_sig_y());

            // intégralles produit
            double int01_rho = compute_int01_gauss_t(mu, sqrt(var));
            double int01_t_rho = compute_int01_t_gauss_t(mu, sqrt(var));

            //  constante des dérivations partielles
            double a_x = m_domain.get_mu_x() - right_cw.x();
            double b_x = left_cw.x()-right_cw.x();
            double a_y = m_domain.get_mu_y() - right_cw.y();
            double b_y = left_cw.y()-right_cw.y();

            // intégralles dérivées
            double int_grad_rho_dx = (A*n_eij_star/m_domain.get_var_x())*(a_x*int01_rho - b_x*int01_t_rho);
            double int_grad_rho_dy = (A*n_eij_star/m_domain.get_var_y())*(a_y*int01_rho - b_y*int01_t_rho);

            int_grad_rho = Vector{int_grad_rho_dx, int_grad_rho_dy};


            // somme des gradient des m_ij // TODO : pb  au calcul des m_ij et obj
            double m_ij = n_eij_star*A*int01_rho; // m_ij actuel
            double m_ij_obj = compute_value_integral()*m_neightbour_proportions.at(i).at(j);// vi->compute_area() * m_neightbour_proportions.at(i).at(j)/m_proportions.at(i);//m_i * p_ij/p_i

            grad_mij += (m_ij - m_ij_obj) * (int_grad_rho + rho_grad_cijl - rho_grad_cijk);

        }

        Vector gi = -grad_mij; // gi = -(nabla mij)
        gradient.push_back(coef * gi);
    }

    if (m_timer_on) Timer::stop_timer(m_timer, COLOR_BLUE);
}

