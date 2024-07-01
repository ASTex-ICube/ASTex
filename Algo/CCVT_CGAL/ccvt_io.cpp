//
// Created by grenier on 12/10/23.
//

#include "ccvt.h"
#include <iostream>


void CCVT::load_image(const std::string& filename)
{
    bool ok = m_domain.load(filename);
    if (!ok) return;

    m_rt.set_boundary(m_domain.get_dx(),
                      m_domain.get_dy());
//    std::cout << "Dx vs Dy: " << m_domain.get_dx() << " ; " << m_domain.get_dy() << std::endl;
}


void CCVT::save_point_eps(const std::string& filename) const
{
    double dx = m_domain.get_dx();
    double dy = m_domain.get_dy();

    double scale = 512.0;
    double radius = 0.004;

    double wx = dx * scale;
    double wy = dy * scale;

    double min_x = 0;
    double max_x = 2.0 * wx;
    double min_y = 0;
    double max_y = 2.0 * wy;

    std::ofstream ofs(filename);
    ofs.precision(20);

    ofs << "%!PS-Adobe-3.1 EPSF-3.0\n";
    ofs << "%%HiResBoundingBox: "
        << min_x << " " << min_y << " " << max_x << " " << max_y << std::endl;
    ofs << "%%BoundingBox: "
        << min_x << " " << min_y << " " << max_x << " " << max_y << std::endl;
    ofs << "%%CropBox: "
        << min_x << " " << min_y << " " << max_x << " " << max_y << "\n";

    ofs << "/radius { " << radius << " } def\n";
    ofs << "/p { radius 0 360 arc closepath fill stroke } def\n";
    ofs << "gsave " << scale << " " << scale << " scale\n";
    ofs << "0 0 0 setrgbcolor" << std::endl;

    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;

        const Point& pi = vi->get_position();
        ofs << pi.x() + dx << " " << pi.y() + dy << " p" << std::endl;
    }
    ofs << "grestore" << std::endl;
    ofs.close();
}




void CCVT::save_cell_eps(const std::string& filename) const
{
    double dx = m_domain.get_dx();
    double dy = m_domain.get_dy();

    double scale = 512.0;
    double radius = 0.004;
    double linewidth = 0.002;

    double wx = dx * scale;
    double wy = dy * scale;

    double min_x = 0;
    double max_x = 2.0 * wx;
    double min_y = 0;
    double max_y = 2.0 * wy;

    std::ofstream ofs(filename);
    ofs.precision(20);

    ofs << "%!PS-Adobe-3.1 EPSF-3.0\n";
    ofs << "%%HiResBoundingBox: "
        << min_x << " " << min_y << " " << max_x << " " << max_y << std::endl;
    ofs << "%%BoundingBox: "
        << min_x << " " << min_y << " " << max_x << " " << max_y << std::endl;
    ofs << "%%CropBox: "
        << min_x << " " << min_y << " " << max_x << " " << max_y << "\n";

    ofs << "/radius { " << radius << " } def\n";
    ofs << "/p { radius 0 360 arc closepath setrgbcolor fill stroke } def\n";
    ofs << "/l { newpath moveto lineto " << linewidth << " setlinewidth stroke } def\n";
    ofs << "gsave " << scale << " " << scale << " scale\n";

    for (unsigned i = 0; i < m_vertices.size(); ++i)
    {
        Vertex_handle vi = m_vertices[i];
        if (vi->is_hidden()) continue;

        const Point& pi = vi->get_position();
        ofs << (i)%2 << " " << 1.-double(i)/m_vertices.size() << " " << double(i)/m_vertices.size() << " " << pi.x() + dx << " " << pi.y() + dy << " p" << std::endl;
    }

    ofs << "0 0 0 setrgbcolor" << std::endl;
    for (Finite_edges_iterator
                 eit = m_rt.finite_edges_begin();
         eit != m_rt.finite_edges_end();
         ++eit)
    {
        Edge edge = *eit;
        Segment segment = m_rt.build_bounded_dual_edge(edge);
        ofs << segment.source().x() + dx << " "<< segment.source().y() + dy << " " << segment.target().x() + dx << " "<< segment.target().y() + dy << " l" << std::endl;
    }
    ofs << "grestore" << std::endl;
    ofs.close();

}

void CCVT::verbose() const
{
//   save_cell_eps(TEMPO_PATH+"results/CCVT_CGAL/cells.eps");

    save_cell_eps("cells.eps");

    if(m_verbose){
        std::vector<FT> weights;
        std::vector<Point> points;
        std::vector<FT> capacities;
        std::vector<FT> areas;
        std::vector<std::vector<FT>> neightbour_proportion;
        std::vector<std::vector<FT>> neightbour_val;

        collect_sites(points, weights);
        areas = get_area();
        capacities = get_capacities();
        neightbour_proportion = get_neightbour_proportion();
        neightbour_val = get_neightbour_val();

        std::cout<<"position (x,y), poids, aires (volumes objectifs)"<<std::endl;
        for(int i=0; i<points.size(); i++){
            std::cout<<"("<<points.at(i).x()<<", "<<points.at(i).y()<<"), "<<weights.at(i)<<", "<<areas.at(i)<<" ("<<capacities.at(i)<<")"<<std::endl;
        }
        std::cout<<std::endl;


        std::cout<<"voisinages proportion"<<std::endl;
        for(auto it=neightbour_proportion.begin(); it<neightbour_proportion.end(); it++){
            for(auto it2=(*it).begin(); it2<(*it).end(); it2++){
                std::cout<<(*it2)<<"; ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

        std::cout<<"voisinages valeurs"<<std::endl;
        for(auto it=neightbour_val.begin(); it<neightbour_val.end(); it++){
            for(auto it2=(*it).begin(); it2<(*it).end(); it2++){
                std::cout<<(*it2)<<"; ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;

//        std::cout<<"c_ij"<<std::endl;
//        for(auto vi=m_vertices.begin(); vi<m_vertices.end(); vi++){
//            Edge_circulator ecirc = m_rt.incident_edges(*vi); // liste des eij
//            Edge_circulator eend  = ecirc;
//
//            CGAL_For_all(ecirc, eend)   // for j in Omega_i
//            {
//                Edge edge = *ecirc; // e_ij
//                if (!m_rt.is_inside(edge)) continue;
//
//                Point c_ij = m_rt.get_edge_cw(edge);
//                std::cout<<c_ij.x()+0.5<<", "<<c_ij.y()+0.5<<std::endl;
//            }
//        }
        std::cout<<" ------------------------";
    }


    if(m_step_by_step){
        std::cout << " ENTER key to continue";
        if (std::cin.get() == '\n') {
            std::cout<<std::endl;
        }
    }
}