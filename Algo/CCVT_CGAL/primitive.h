//
// Created by grenier on 12/10/23.
//

#ifndef CCVT_TEST_CHARLINE_PRIMITIVE_H
#define CCVT_TEST_CHARLINE_PRIMITIVE_H


// STL
#include <vector>

// local
//#include "convex_polygon.h"
#include "pixel.h"

template <class Kernel, class Vbb>
class My_vertex_base : public Vbb
{
public:
    typedef typename Kernel::FT FT;
//    typedef typename Kernel::Weight   Weight;
    typedef typename Kernel::Point_2  Point2;
    typedef typename Kernel::Weighted_point_2  Point;
    typedef typename Kernel::Vector_2 Vector;
    typedef typename Kernel::Weighted_point_2 Weighted_point;

    typedef typename Vbb::Triangulation_data_structure TDS;
    typedef typename TDS::Face_handle   Face_handle;
    typedef typename TDS::Vertex_handle Vertex_handle;

    template < typename TDS2 >
    struct Rebind_TDS {
        typedef typename Vbb::template Rebind_TDS<TDS2>::Other Vb2;
        typedef My_vertex_base<Kernel, Vb2> Other;
    };

    typedef CConvexPolygon<Kernel> ConvexPolygon;
    typedef CPixel<Kernel> Pixel;

private:
    int m_index;
    ConvexPolygon m_dual;
    std::vector<Pixel> m_pixels;
    FT m_area;

public:
    My_vertex_base() : Vbb()
    {
        m_index = -1;
        m_area = 0.0;
    }

    My_vertex_base(const Weighted_point& p) : Vbb(p)
    {
        m_index = -1;
        m_area = 0.0;
    }

    My_vertex_base(const Weighted_point& p, Face_handle f) : Vbb(p, f)
    {
        m_index = -1;
        m_area = 0.0;
    }

    ~My_vertex_base()
    {
        m_pixels.clear();
    }

    void set_index(const int x) { m_index = x; }
    const int get_index() const { return m_index; }

    // POSITION / WEIGHT //

    const Point2& get_position() const
    {
        return this->point().point();
    }

    void set_position(const Point2& p)
    {
        FT w = get_weight();
        Weighted_point wp(p, w);
        this->set_point(wp);
    }

    const FT get_weight() const
    {
        return this->point().weight();
    }

    void set_weight(const FT w)
    {
        const Point2& p = get_position();
        Weighted_point wp(p, w);
        this->set_point(wp);
    }

    // PIXELS //

    void clear_pixels()
    {
        m_pixels.clear();
    }

    unsigned nb_pixels() const
    {
        return m_pixels.size();
    }

    const Pixel& get_pixel(unsigned i) const
    {
        return m_pixels[i];
    }

    void append_pixel(const Pixel& pixel)
    {
        m_pixels.push_back(pixel);
    }

    // ATTRIBUTES //

    FT compute_area() const
    {
        return m_area;
    }

    void pre_compute_area()
    {
        m_area = 0.0;
        for (unsigned i = 0; i < nb_pixels(); ++i)
        {
            const Pixel& pixel = get_pixel(i);
            m_area += pixel.compute_area();
        }
    }

    Point2 compute_centroid() const
    {
        FT sum_area = 0.0;
        Vector sum_vector = CGAL::NULL_VECTOR;
        for (unsigned i = 0; i < nb_pixels(); ++i)
        {
            const Pixel& pixel = get_pixel(i);
            FT area = pixel.compute_area();
            Point2 centroid = pixel.compute_centroid();

            sum_area += area;
            sum_vector = sum_vector + area*(centroid - CGAL::ORIGIN);
        }
        if (sum_area == 0.0) return get_position();
        return CGAL::ORIGIN + (sum_vector / sum_area);
    }

    FT compute_variance() const
    {
        FT variance = 0.0;
        const Point2& q = get_position();
        for (unsigned i = 0; i < nb_pixels(); ++i)
        {
            const Pixel& pixel = get_pixel(i);
            variance += pixel.compute_variance(q);
        }
        return variance;
    }

    // DUAL CELL //

    ConvexPolygon& dual()
    { return m_dual; }

    const ConvexPolygon& get_dual() const
    { return m_dual; }

    void set_dual(const ConvexPolygon& dual)
    { m_dual = dual; }
};

template <class Kernel, class Fbb>
class My_face_base : public Fbb
{
public:
    typedef typename Fbb::Triangulation_data_structure  TDS;
    typedef typename TDS::Vertex_handle Vertex_handle;
    typedef typename TDS::Face_handle   Face_handle;

    template < typename TDS2 >
    struct Rebind_TDS {
        typedef typename Fbb::template Rebind_TDS<TDS2>::Other Fb2;
        typedef My_face_base<Kernel, Fb2> Other;
    };

public:
    My_face_base()
            : Fbb()
    {
    }

    My_face_base(Vertex_handle v1,
                 Vertex_handle v2,
                 Vertex_handle v3)
            : Fbb(v1, v2, v3)
    {
    }

    My_face_base(Vertex_handle v1,
                 Vertex_handle v2,
                 Vertex_handle v3,
                 Face_handle f1,
                 Face_handle f2,
                 Face_handle f3)
            : Fbb(v1,v2,v3,f1,f2,f3)
    {
    }
};


#endif //CCVT_TEST_CHARLINE_PRIMITIVE_H
