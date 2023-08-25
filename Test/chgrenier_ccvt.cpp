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
    const int     NUMBER_SITES      = 6; // nombre de graine de cellule
    const int     NUMBER_POINTS     = 4096 * NUMBER_SITES; // nombre de points pour la contrainte de densité
    const double  TORUS_SIZE        = 256; // taille de l'image
    std::string directory = "/home/grenier/Documents/ASTex_fork/results/ccvt/";

    typedef MetricEuclidean2 Metric;
    typedef Point2 Point;

    // contrainte de densité
    Point::List points;
//    constant_regular_density(points, NUMBER_POINTS, TORUS_SIZE); // densité constante
    nonconstant_density(points, NUMBER_POINTS, TORUS_SIZE, .16, .16); // densité non constante


    // graine voronoi
    int overallCapacity = static_cast<int>(points.size());
    Site<Point>::List sites;

    // donnée initiales
    // initialisation position aléatoire, capacité toutes égales
    for (int i = 0; i < NUMBER_SITES; ++i) {
        double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
        double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
        int capacity = overallCapacity / (NUMBER_SITES - i);
        overallCapacity -= capacity;
        sites.push_back(Site<Point2>(i, capacity, Point2(x, y)));
    }


//    // initialisation position et capacité controllés
//    int globalCapacity = overallCapacity;
//    std::vector<double> init_vol{0.25/NUMBER_SITES, 1./NUMBER_SITES, 1./NUMBER_SITES, 1./NUMBER_SITES, 1.75/NUMBER_SITES, 1./NUMBER_SITES};
//    std::vector<std::array<double, 2>> init_pos{{0.8, 0.75}, {0.4, 0.4}, {0.8, 0.3}, {0.4, 0.2}, {0.1, 0.5}, {0.35, 0.8}};
//    std::vector<double> colors {24., 124., 255., 190., 24., 124.};
//
//    for (int i = 0; i < NUMBER_SITES-1; ++i) {
//        double x = init_pos.at(i).at(0) * TORUS_SIZE;  //static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//        double y = init_pos.at(i).at(1) * TORUS_SIZE;  //static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//
//        int capacity = std::round(init_vol.at(i) * overallCapacity);
//        globalCapacity -= capacity;
//
//        sites.push_back(Site<Point2>(i, capacity, Point2(x, y)));
//    }
//
//    // le dernier point récupère les miette de capacité des autres (due aux imprecisions numériques)
//    double x = init_pos.at(NUMBER_SITES-1).at(0) * TORUS_SIZE;  //static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//    double y = init_pos.at(NUMBER_SITES-1).at(1) * TORUS_SIZE;  //static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//
//    int capacity = globalCapacity;
//    sites.push_back(Site<Point2>(NUMBER_SITES-1, capacity, Point2(x, y)));




















    // initialisation du ccvt
    Metric metric(Point(TORUS_SIZE, TORUS_SIZE));
    Optimizer<Site<Point>, Point, Metric> optimizer;
    optimizer.initialize(sites, points, metric);

    // écriture état initial
    save_res_point(optimizer.sites(), metric, 1, TORUS_SIZE, directory+"ccvt_point_initialisation.png");
//    save_res_zone(optimizer.sites(), metric, colors, TORUS_SIZE, directory+"ccvt_zone_initalisation.png");




    // optimization
    int iteration = 0;
    bool stable;
    do {
        printf("iteration %d...", ++iteration);
        stable = optimizer.optimize(true);
        printf("done\n");
    } while (!stable);




    // écriture du résultat
    const Site<Point>::Vector& result = optimizer.sites();
    printf("\nfinal positions :\n");
    for (unsigned int i = 0; i < result.size(); ++i) {
        printf("site %d: %f, %f\n", result[i].id, result[i].location.x, result[i].location.y);
    }



    // écriture dans des images
    save_res_point(result, metric, 1, TORUS_SIZE, directory+"ccvt_point_result.png");
    save_res_cell(result, metric, TORUS_SIZE, directory+"ccvt_cell_result.png");
//    save_res_zone(result, metric, colors, TORUS_SIZE, directory+"ccvt_zone_result.png");


    return 0;
}