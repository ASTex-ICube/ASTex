//
// Created by grenier on 26/06/23.
//

#ifndef ASTEX_SITES_H
#define ASTEX_SITES_H

#include <vector>
#include <list>

template<class Point>
class Site{
public:
    typedef std::list<Site>     List;
    typedef std::vector<Site>   Vector;
    typedef std::vector<Site*>	VectorPtr;

    Site()
            : id(-1),
              capacity(0) {
    }

    Site(const int id, const int capacity, const Point& location)
            : id(id),
              capacity(capacity),
              location(location) {
    }

    int	  id;
    int   capacity;
    Point location;
};

#endif //ASTEX_SITES_H
