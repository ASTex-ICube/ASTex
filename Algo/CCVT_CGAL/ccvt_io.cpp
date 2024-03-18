//
// Created by grenier on 12/10/23.
//

#include "ccvt.h"

void CCVT::load_image(const std::string& filename)
{
    bool ok = m_domain.load(filename);
    if (!ok) return;

    m_rt.set_boundary(m_domain.get_dx(),
                      m_domain.get_dy());
//    std::cout << "Dx vs Dy: " << m_domain.get_dx() << " ; " << m_domain.get_dy() << std::endl;
}


void CCVT::save_eps(const std::string& filename) const
{
    double dx = m_domain.get_dx();
    double dy = m_domain.get_dy();

    double scale = 512.0;
    double radius = 0.002;

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