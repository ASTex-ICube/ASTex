//
// Created by grenier on 26/06/23.
//

#include <iostream>
#include <numeric>

#include <ASTex/easy_io.h>
#include <ASTex/image_rgb.h>

#include "ASTex/CCVT/point.h"
#include "ASTex/CCVT/sites.h"
#include "ASTex/CCVT/metric.h"
#include "ASTex/CCVT/optimizer.h"
#include "ASTex/CCVT/tools.h"

using namespace ASTex;



int main(){
    const int     NUMBER_SITES      = 4; // nombre de graine de cellule
    const int     NUMBER_POINTS     = 4096 * NUMBER_SITES; // nombre de points pour la contrainte de densité
    const double  TORUS_SIZE        = 256; // taille de l'image
    std::string directory = "/home/grenier/Documents/ASTex_fork/results/ccvt/";

    typedef MetricEuclidean2 Metric;
    typedef Point2 Point;

    // contrainte de densité
    Point::List points;
//    constant_regular_density(points, NUMBER_POINTS, TORUS_SIZE); // densité constante
//    nonconstant_density(points, NUMBER_POINTS, TORUS_SIZE, .16, .16); // densité non constante in
    nonconstant_density(points, NUMBER_POINTS, TORUS_SIZE, 0.081327, 0.0816067); // densité non constante


    // graine voronoi
    unsigned int overallCapacity = static_cast<int>(points.size());
    Site<Point>::List sites;

    // donnée initiales

//    // initialisation position et capacité controllées
//    Point::Vector init_loc{Point(600., 500.), Point(530., 595.), Point(419., 558.), Point(419., 441.), Point(530., 404.), Point(500., 500.)}; // 1 point au centre et 5 points en cercle (pour une image 1000x1000
//    std::vector<float> init_vol{1., 1., 3., 1., 1., 1.};
//
//    for (int i = 0; i < NUMBER_SITES; ++i) {
//        Point position = init_loc.at(i);
//        int capacity = std::round(init_vol.at(i) * overallCapacity / float(NUMBER_SITES - i));
//        overallCapacity -= capacity;
//        sites.push_back(Site<Point2>(i, capacity, position));
//    }


    // initialisation partiellement controllée : proportion de volume (pourcentage), position pour une distance en centre donnée (angle aléatoire)
//    std::vector<float> init_rad{100., 150., 200., 250., 300., 350.}; // entre 0 et taille image/2
//    std::vector<float> init_cap{0.131405, 0.356746, 0.109477, 0.116753, 0.157846, 0.127772}; // proportion de chaque couleur (somme à 1)
//    unsigned int globalCapacity = overallCapacity;
//
//    for (int i=0; i<NUMBER_SITES -1; i++){
//        double init_angle = (static_cast<double>(rand() % RAND_MAX) / RAND_MAX) * 2. * M_PI; // angle aléatoire sur 0, 2pi
//        double x = init_rad.at(i)* cos(init_angle) + TORUS_SIZE/2.;
//        double y = init_rad.at(i)* sin(init_angle) + TORUS_SIZE/2.;
//        Point position = Point{x, y};
//
//        int capacity = init_cap.at(i) * overallCapacity;
//        globalCapacity -= capacity;
//        sites.push_back(Site<Point2>(i, capacity, position));
//    }
//    // le dernier point récupère les miette de capacité des autres (due aux imprecisions numériques)
//    double init_angle = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * 2. * M_PI; // angle aléatoire sur 0, 2pi
//    double x = init_rad.at(NUMBER_SITES-1)* cos(init_angle) + TORUS_SIZE/2.;
//    double y = init_rad.at(NUMBER_SITES-1)* sin(init_angle) + TORUS_SIZE/2;
//    Point position = Point{x, y};
//
//    int capacity = globalCapacity;
//    sites.push_back(Site<Point2>(NUMBER_SITES-1, capacity, position));



//    // initialisation position aléatoire, capacité toutes égales
//    for (int i = 0; i < NUMBER_SITES; ++i) {
//        double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//        double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//        int capacity = overallCapacity / (NUMBER_SITES - i);
//        overallCapacity -= capacity;
//        sites.push_back(Site<Point2>(i, capacity, Point2(x, y)));
//    }


    // initialisation position aléatoire et capacité controllées (ish)
    unsigned int globalCapacity = overallCapacity;
//    std::vector<float> init_vol{1./NUMBER_SITES, 1.5/NUMBER_SITES, .5/NUMBER_SITES, 1./NUMBER_SITES}; // input
    std::vector<float> init_vol{0.368759,
                                0.270779,
                                0.121647,
                                0.238815}; // donnée de la composition d'entrée

    for (int i = 0; i < NUMBER_SITES-1; ++i) {
        double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
        double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;

        int capacity = std::round(init_vol.at(i) * overallCapacity);
        globalCapacity -= capacity;

        sites.push_back(Site<Point2>(i, capacity, Point2(x, y)));
    }
    // le dernier point récupère les miette de capacité des autres (due aux imprecisions numériques)
    double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
    double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;

    int capacity = globalCapacity;
    sites.push_back(Site<Point2>(NUMBER_SITES-1, capacity, Point2(x, y)));




    // initialisation du ccvt
    Metric metric(Point(TORUS_SIZE, TORUS_SIZE));
    Optimizer<Site<Point>, Point, Metric> optimizer;
    optimizer.initialize(sites, points, metric);

    // état initial
    save_res_point(optimizer.sites(), metric, 1, TORUS_SIZE, directory+"ccvt_point_initialisation.png");




    // optimization
    int iteration = 0;
    bool stable;
    do {
        printf("iteration %d...", ++iteration);
        stable = optimizer.optimize(true);
        printf("done\n");
    } while (!stable);






    // résultat
    const Site<Point>::Vector& result = optimizer.sites();
    printf("\nfinal positions :\n");
    for (unsigned int i = 0; i < result.size(); ++i) {
        printf("site %d: %f, %f\n", result[i].id, result[i].location.x, result[i].location.y);
    }

    std::vector<ImageRGB8::PixelType> colors {ImageRGB8::itkPixel(128, 12, 24),
                                              ImageRGB8::itkPixel(24, 128, 12),
                                              ImageRGB8::itkPixel(12, 24, 128),
                                              ImageRGB8::itkPixel(24, 124, 128)};

    // écriture dans des images
    save_res_point(result, metric, 1, TORUS_SIZE, directory+"ccvt_point_result.png");
    save_res_cell(result, metric, TORUS_SIZE, directory+"ccvt_cell_result.png");
    save_res_zone(result, metric, colors, TORUS_SIZE, directory+"ccvt_zone_result.png");


    return 0;
}