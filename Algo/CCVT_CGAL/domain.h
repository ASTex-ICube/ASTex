//
// Created by grenier on 12/10/23.
//

#ifndef CCVT_TEST_CHARLINE_DOMAIN_H
#define CCVT_TEST_CHARLINE_DOMAIN_H


// local
#include "util.h"
#include "pgm.h"

#define PIXEL_EPS 0.0 //1.0e-6

//struct DomainImage{
//    double mu_x;
//    double mu_y;
//    double var_x;
//    double var_y;
//
//    double m_max = 1.;
//    unsigned m_w = 256;
//    unsigned m_h = 256;
//
//    DomainImage(){
//        var_x = -0.1;
//        var_y = -0.1;
//    }
//
//    void set_data(double moy_x, double moy_y, double variance_x, double variance_y){
//        mu_x = moy_x;
//        mu_y = moy_y;
//        var_x = variance_x;
//        var_y = variance_y;
//    }
//    bool isNull() const {
//        return (var_x <= 0. or var_y <= 0.);
//    }
//    unsigned width() const {
//        return m_w;
//    }
//    unsigned height() const {
//        return m_h;
//    }
//    double max_value() const {
//        return m_max;
//    }
//    double pixel(unsigned i, unsigned j) const {
//        double Gx = std::exp(-(i-mu_x)*(i-mu_x)/(2.*var_x));
//        double Gy = std::exp(-(j-mu_y)*(j-mu_y)/(2.*var_y));
//        return Gx*Gy;
//    }
//
//    bool load(const std::string& filename){}
//};




template <class Kernel>
class CDomain
{
public:
    typedef typename Kernel::FT        FT;
    typedef typename Kernel::Point_2   Point;
    typedef typename Kernel::Segment_2 Segment;

private:
    double m_dx;
    double m_dy;
    double m_mu_x;
    double m_mu_y;
    double m_var_x;
    double m_var_y;
    PGMImage m_image;
//    DomainImage m_image;
    bool m_invert;

public:
    CDomain()
    {
        m_dx = 0.0;
        m_dy = 0.0;
        m_invert = false;
    }


    void toggle_invert() { m_invert = !m_invert; }

    double compute_area() const { return 4.0*std::abs(m_dx)*std::abs(m_dy); }

    double get_dx() const { return m_dx; }

    double get_dy() const { return m_dy; }

    double get_px() const { return 2.0*get_dx()/get_width(); } // pixel width

    double get_py() const { return 2.0*get_dy()/get_height(); } // pixel height

    double get_mu_x() const { return m_mu_x; }

    double get_mu_y() const { return m_mu_y; }

    double get_var_x() const { return m_var_x; } // variance
    double get_sig_x() const { return sqrt(m_var_x); }  // écart type

    double get_var_y() const { return m_var_y; }
    double get_sig_y() const { return sqrt(m_var_y); }

    void tonemap(double /*key*/) { /*m_image.tonemap(key);*/ }

    double compute_x(unsigned i) const
    {
        unsigned w = get_width();
        return m_dx*(2.0*double(i)/double(w) - 1.0);
    }

    double compute_y(unsigned j) const
    {
        unsigned h = get_height();
        return m_dy*(1.0 - 2.0*double(j)/double(h));
    }

    void locate(const Point& p,
                unsigned& i,
                unsigned& j,
                bool pixelid = true) const
    {
        double x = 0.5 * (p.x() + m_dx) / m_dx;
        double y = 0.5 * (m_dy - p.y()) / m_dy;
        x = std::max(0.0, x);
        x = std::min(1.0, x);
        y = std::max(0.0, y);
        y = std::min(1.0, y);
        unsigned w = get_width() ;
        unsigned h = get_height();
        i = (unsigned) std::floor( x * double(w) + EPS );
        j = (unsigned) std::floor( y * double(h) + EPS );
        if (pixelid && i == w) i--;
        if (pixelid && j == h) j--;
    }

    bool is_grid_line(const Segment& segment) const
    {
        const Point& p0 = segment.source();
        const Point& p1 = segment.target();
        if (p0 == p1) return true;

        unsigned i0, j0;
        locate(p0, i0, j0, false);
        double x0 = compute_x(i0);
        double y0 = compute_y(j0);
        if (std::abs(x0 - p0.x()) > EPS ||
            std::abs(y0 - p0.y()) > EPS) return false;

        unsigned i1, j1;
        locate(p1, i1, j1, false);
        double x1 = compute_x(i1);
        double y1 = compute_y(j1);
        if (std::abs(x1 - p1.x()) > EPS ||
            std::abs(y1 - p1.y()) > EPS) return false;

        if (i0 == i1) return true;
        if (j0 == j1) return true;
        return false;
    }

    bool is_valid() const
    {
        return !m_image.isNull();
    }

    unsigned get_width() const
    {
        return m_image.width();
    }

    unsigned get_height() const
    {
        return m_image.height();
    }

    double get_max_value() const
    {
        return m_image.max_value();
//        return 1.0; // TODO : pu*** les gars sérieux ??!
    }

    bool load(const std::string& filename)
    {
        bool ok = m_image.load(filename);
        //bool ok = m_image.load(filename.toStdString());
        if (!ok) return false;

        unsigned w = get_width();
        unsigned h = get_height();
        m_dx = 0.5;
        m_dy = m_dx * double(h) / double(w);
        return true;
    }

    bool set(double moy_x, double moy_y, double variance_x, double variance_y, unsigned size_x, unsigned size_y, double max_val){
        m_image.set_data(moy_x, moy_y, variance_x, variance_y, size_x, size_y, max_val);
        m_mu_x = moy_x;
        m_mu_y = moy_y;
        m_var_x = variance_x;
        m_var_y = variance_y;

        unsigned w = get_width();
        unsigned h = get_height();
        m_dx = 0.5;
        m_dy = m_dx * double(h) / double(w);
//        m_image.save("/home/grenier/Documents/ASTex_fork/results/CCVT_CGAL/density.pgm");
        return true;
    }

    double get_value(unsigned i, unsigned j, const bool normalized = false) const
    {
        double value = m_image.pixel(i, j);
//        double value = double(qGray(m_image.pixel(i, j))) / 255.0;
        if (m_invert) value = get_max_value() - value;
        if (normalized) value /= get_max_value();
        return value + PIXEL_EPS;
    }

    double get_value(const Point& p, const bool normalized = false) const
    {
        unsigned i, j;
        locate(p, i, j);
        return get_value(i, j, normalized);
    }

    double integrate_intensity() const
    {
        double sum = 0.0;
        double pixel_area = get_px()*get_py();
        for (unsigned i = 0; i < get_width(); ++i)
        {
            for (unsigned j = 0; j < get_height(); ++j)
            {
                double value = get_value(i, j);
                sum += value*pixel_area;
            }
        }
        return sum;
    }

    bool is_outside(const Point& p) const
    {
        if (std::abs(p.x()) > m_dx) return true;
        if (std::abs(p.y()) > m_dy) return true;
        return false;
    }

    Point clamp(const Point& p) const
    {
        double x = p.x();
        if (x < -m_dx) x = -m_dx;
        if (x >  m_dx) x =  m_dx;

        double y = p.y();
        if (y < -m_dy) y = -m_dy;
        if (y >  m_dy) y =  m_dy;

        return Point(x, y);
    }

    // Draw

//    void draw_boundary() const
//    {
//        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//        glBegin(GL_QUADS);
//        glVertex2d(-m_dx, -m_dy);
//        glVertex2d( m_dx, -m_dy);
//        glVertex2d( m_dx,  m_dy);
//        glVertex2d(-m_dx,  m_dy);
//        glEnd();
//    }
//
//    void draw_image() const
//    {
//        double dx = 2.0 * m_dx / get_width();
//        double dy = 2.0 * m_dy / get_height();
//        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//        for (unsigned i = 0; i < get_width(); ++i)
//        {
//            double x = - m_dx + i*dx;
//            for (unsigned j = 0; j < get_height(); ++j)
//            {
//                double y = m_dy - j*dy;
//                double value = 0.9*get_value(i, j, true);
//                glColor3d(value, value, value);
//                draw_quad(x, y, dx, dy);
//            }
//        }
//    }
//
//    void draw_grid() const
//    {
//        double dx = 2.0 * m_dx / get_width();
//        double dy = 2.0 * m_dy / get_height();
//        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//        for (unsigned i = 0; i < get_width(); ++i)
//        {
//            double x = - m_dx + i*dx;
//            for (unsigned j = 0; j < get_height(); ++j)
//            {
//                double y = m_dy - j*dy;
//                glColor3d(0.0, 0.0, 0.0);
//                draw_quad(x, y, dx, dy);
//            }
//        }
//    }
//
//    void draw_quad(double x, double y, double dx, double dy) const
//    {
//        glBegin(GL_QUADS);
//        glVertex2d(x     , y - dy);
//        glVertex2d(x + dx, y - dy);
//        glVertex2d(x + dx, y     );
//        glVertex2d(x     , y     );
//        glEnd();
//    }

};

#endif //CCVT_TEST_CHARLINE_DOMAIN_H
