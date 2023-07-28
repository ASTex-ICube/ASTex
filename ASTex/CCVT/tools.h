//
// Created by grenier on 26/06/23.
//

#ifndef ASTEX_TOOLS_H
#define ASTEX_TOOLS_H

#include <ASTex/easy_io.h>
#include <ASTex/image_rgb.h>

#include "ASTex/CCVT/point.h"
#include "ASTex/CCVT/sites.h"
#include "ASTex/CCVT/metric.h"


using namespace ASTex;

template<class Point>
void nonconstant_density(std::list<Point>& points, const int numberOfPoints, const double torusSize, float var1, float var2) {
    const double E = 2.718281828459;
    const double PI = 3.141592653590;
    while (points.size() < static_cast<unsigned int>(numberOfPoints)) {
        double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * 2 - 1;
        double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * 2 - 1;

//    double p = pow(E, -20.0 * x * x - 20.0 * y * y) + 0.2 * sin(PI * x) * sin(PI * x) * sin(PI * y) * sin(PI * y);
//        double p = pow(E, -4.0 * x * x - 4.0 * y * y);
        double p = std::exp(-0.5 * (x*x)/var1) * std::exp(-0.5 * (y*y)/var2);

        double r = static_cast<double>(rand() % RAND_MAX) / RAND_MAX;
        if (p >= r) {
            points.push_back(Point((x + 1) / 2 * torusSize, (y + 1) / 2 * torusSize));
        }
    }
}


template<class Point>
void constant_regular_density(std::list<Point>& points, const int numberOfPoints, const double torusSize) {
    double n = sqrt(static_cast<double>(numberOfPoints));
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < n; ++y) {
            double dx = x / n * torusSize;
            double dy = y / n * torusSize;
            points.push_back(Point(dx, dy));
        }
    }
}

template<class Site, class Metric>
bool save_res_cell(const std::vector<Site>& result, Metric metric, int img_size, std::string img_name){
    ImageGrayu8 image_(img_size, img_size);

    image_.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y) // cm
                                   {
                                       double dist = 12400.;
                                       for (unsigned int i = 0; i < result.size(); ++i)
                                       {
                                           double new_dist = metric.distance(Point2(x,y), result[i].location);
                                           dist = std::min(dist, new_dist);
                                       }

                                       dist = std::clamp(255.*(dist/500.), 0., 255.);
                                       P = ImageGrayu8::PixelType(dist);
                                   });

    if(image_.save(img_name)){
        return true;
    }
    return false;
}

template<class Site, class Metric>
bool save_res_point(const std::vector<Site>& result, Metric metric, float radius, int img_size, std::string img_name){
    ImageGrayu8 image_(img_size, img_size);

    image_.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y) // cm
                                   {
                                       double dist = 12400.;
                                       for (unsigned int i = 0; i < result.size(); ++i)
                                       {
                                           double new_dist = metric.distance(Point2(x,y), result[i].location);
                                           dist = std::min(dist, new_dist);
                                       }

                                       dist = std::clamp(255.*(dist/500.), 0., 255.) - radius;
                                       double point = dist < 0. ? 0. : 255.;

                                       P = ImageGrayu8::PixelType(point);
                                   });

    if(image_.save(img_name)){
        return true;
    }
    return false;
}

template<class Site, class Metric>
bool save_res_zone(const std::vector<Site>& result, Metric metric, std::vector<ImageRGB8::PixelType> colors, int img_size, std::string img_name){
    ImageGrayu8 image_(img_size, img_size);

    image_.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
                                   {
                                       double dist = 12400.;
                                       int color_id = -1;
                                       for (unsigned int i = 0; i < result.size(); ++i)
                                       {
                                           double new_dist = metric.distance(Point2(x,y), result[i].location);
                                           if(new_dist < dist){
                                               dist = new_dist;
                                               color_id = i;
                                           }
                                       }

//                                       dist = std::clamp(255.*(dist/500.), 0., 255.) - radius;
                                       double color = (color_id+1)*(255./4.);//colors.at(color_id);

                                       P = ImageGrayu8::PixelType(color);
                                   });

    if(image_.save(img_name)){
        return true;
    }
    return false;
}

#endif //ASTEX_TOOLS_H
