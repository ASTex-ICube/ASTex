//
// Created by grenier on 26/06/23.
//

#ifndef ASTEX_POINT_H
#define ASTEX_POINT_H

#include <vector>
#include <list>
#include <assert.h>


struct Point2 {

    typedef std::list<Point2>     List;
    typedef std::vector<Point2>   Vector;

    Point2()
            : x(0), y(0) {
    }

    Point2(const double x, const double y)
            :x(x), y(y) {
    }

    Point2(const Point2& p)
            : x(p.x), y(p.y) {
    }

    inline const double operator[](const int i) const {
        assert(i == 0 || i == 1);
        if (i == 0) {
            return x;
        }
        return y;
    }

    inline double operator[](const int i) {
        assert(i == 0 || i == 1);
        if (i == 0) {
            return x;
        }
        return y;
    }

    double  x;
    double  y;

    static const int D = 2;

};

struct Point3 {

    typedef std::list<Point3>   List;
    typedef std::vector<Point3> Vector;

    Point3()
            : x(0), y(0), z(0) {
    }

    Point3(const double x, const double y, const double z)
            : x(x), y(y), z(z) {
    }

    Point3(const Point3& p)
            : x(p.x), y(p.y), z(p.z) {
    }

    inline const double operator[](const int i) const {
        assert(i == 0 || i == 1 || i == 2);
        if (i == 0) {
            return x;
        }
        if (i == 1) {
            return y;
        }
        return z;
    }

    inline double operator[](const int i) {
        assert(i == 0 || i == 1 || i == 2);
        if (i == 0) {
            return x;
        }
        if (i == 1) {
            return y;
        }
        return z;
    }

    double  x;
    double  y;
    double  z;

    static const int D = 3;

};

template<int d>
struct Point {

    typedef std::list<Point>   List;
    typedef std::vector<Point> Vector;

    Point() {
        for (int i = 0; i < d; ++i) {
            coordinates[i] = 0.0;
        }
    }

    Point(const double coordinates[d]) {
        for (int i = 0; i < d; ++i) {
            this->coordinates[i] = coordinates[i];
        }
    }

    Point(const Point& p) {
        for (int i = 0; i < d; ++i) {
            coordinates[i] = p.coordinates[i];
        }
    }

    inline const double operator[](const int i) const {
        assert(i >= 0 && i < d);
        return coordinates[i];
    }

    inline double operator[](const int i) {
        assert(i >= 0 && i < d);
        return coordinates[i];
    }

    static const int D = d;

private:

    double coordinates[d];

};


#endif //ASTEX_POINT_H
