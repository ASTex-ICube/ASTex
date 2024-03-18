//
// Created by grenier on 12/10/23.
//

#ifndef CCVT_TEST_CHARLINE_RANDOM_H
#define CCVT_TEST_CHARLINE_RANDOM_H

#include <cstdlib>

inline
double random_double(const double min, const double max)
{
    double range = max - min;
    return min + (double(std::rand()) / double(RAND_MAX)) * range;
}

inline
int random_int(const int min, const int max)
{
    int range = max - min;
    return min + int((double(std::rand())/double(RAND_MAX)) * range);
}

template <class Vector>
Vector random_vec(const double scale)
{
    double dx = random_double(-scale, scale);
    double dy = random_double(-scale, scale);
    return Vector(dx, dy);
}

#endif //CCVT_TEST_CHARLINE_RANDOM_H
