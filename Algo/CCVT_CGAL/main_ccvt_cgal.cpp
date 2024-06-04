//
// Created by grenier on 09/10/23.
//

#include <iostream>

#include "ccvt.h"
#include "ccvt_application.h"
//#include "types.h"

//#include <GLFW/glfw3.h>

struct comput_err{ // foncteur de calcul d'erreur
    double operator()(const double &ref, const double &est) const{
        return 100*std::abs(ref-est)/ref;
    }
};


int main() {
    // description de la densité
    double mu_x = 0.5;
    double mu_y = 0.5;
    double var_x = 0.024025;
    double var_y = 0.024025;

    // position initial des graines
    std::vector<Point> init_sites{Point(0.81, 0.65), Point(0.15, 0.25), Point(0.52, 0.91), Point(0.65, 0.28), Point(0.12, 0.72), Point(0.45, 0.48)};
//    std::vector<FT> custom_capacities{0.0947948, 0.162811, 0.0396879, 0.138455, 0.149857, 0.414395};
    std::vector<FT> custom_capacities{0.166666667, 0.166666667, 0.166666667, 0.166666667, 0.166666667, 0.166666667};

    CCVT main_ccvt;
    main_ccvt.set_domain(mu_x, mu_y, var_x, var_y);
    main_ccvt.set_custom_proportions(custom_capacities);
    main_ccvt.set_initial_sites(init_sites);

    ccvt_application app;
    app.attacheCCVT(main_ccvt);

    if (!app.onInit()) return 1;

    while (app.isRunning()) {
        app.onFrame();
    }

    app.onFinish();
    return 0;
}


//int main()
//{
//
//    std::string working_directory = "/home/grenier/Documents/ASTex_fork/results/CCVT_CGAL/";
//
//    FT stepX = 0.01; // pour les positions
//    FT stepW = 0.1; // pour les poids
//    FT epsilon = 1.;
//    unsigned max_newton_iters = 500;
//    unsigned max_iters = 500;
//    unsigned nb_site = 24;//6;
//    unsigned seed = 624;
//
//    // description de la densité
//    double mu_x = 0.500233;
//    double mu_y = 0.499979;
//    double var_x = 0.0151757;
//    double var_y = 0.0198063;
//
//    // proportion de présence des couleurs
////    std::vector<FT> custom_capacities{0.0947948, 0.162811, 0.0396879, 0.138455, 0.149857, 0.414395}; // proportion de capacité objectif
//    std::vector<FT> custom_capacities{0.095978584176086,
//                                        0.031243307555027,
//                                        0.035431290898275,
//                                        0.05437239738251,
//                                        0.032659131469364,
//                                        0.041719214753123,
//                                        0.033491969066032,
//                                        0.0338251041047,
//                                        0.038953004164188,
//                                        0.049535990481856,
//                                        0.028941106484236,
//                                        0.034610350981559,
//                                        0.025074360499703,
//                                        0.059922665080309,
//                                        0.050826888756693,
//                                        0.032343842950625,
//                                        0.052325996430696,
//                                        0.028102320047591,
//                                        0.042028554431886,
//                                        0.042784057108864,
//                                        0.043468173706127,
//                                        0.039750148720999,
//                                        0.035776323616895,
//                                        0.036835217132659};
//    assert(custom_capacities.size() == nb_site); // une proportion par graine
//    assert(std::abs(std::accumulate(custom_capacities.begin(), custom_capacities.end(), 0.0) - 1.) < 0.000001); // somme des proportion = 1
//
////    // proportion de voisinage entre les couleurs
////    std::vector<std::vector<FT>> custom_neightbour_capacities{{0.0, 0.0, 0.178821415058645, 0.185253499810821, 0.0, 0.635925085130533},
////                                                              {0.0, 0.0, 0.0, 0.172562466051059, 0.186345736013036, 0.641091797935904},
////                                                              {0.308099739243807, 0.0, 0.0, 0.0, 0.360250977835724, 0.331649282920469},
////                                                              {0.145117071724956, 0.188315056312982, 0.0, 0.0, 0.0, 0.666567871962063},
////                                                              {0.0, 0.191301014184644, 0.154079392186248, 0.0, 0.0, 0.654619593629108},
////                                                              {0.183751281175265, 0.258066279467031, 0.055620088828152, 0.245876323881107, 0.256686026648445, 0.0}};
////    assert(custom_neightbour_capacities.size() == nb_site);
////    for(int cap=0; cap<custom_neightbour_capacities.size(); cap++){
////        assert(custom_neightbour_capacities.at(cap).size() == nb_site);
////        assert(std::abs(std::accumulate(custom_neightbour_capacities.at(cap).begin(), custom_neightbour_capacities.at(cap).end(), 0.0) - 1.) < 0.01);
////    }
////
////    // position initial des graines
////    std::vector<Point> init_sites{Point(0.81, 0.65),
////                                  Point(0.15, 0.25),
////                                  Point(0.52, 0.91),
////                                  Point(0.65, 0.28),
////                                  Point(0.12, 0.72),
////                                  Point(0.45, 0.48)};
////
////    assert(init_sites.size() == nb_site); // une position initial par graine
//
//
//
//
//
//
//
//
////    // ------------------------------------------------------------------------------------------------
////    CCVT main_ccvt{seed};
////
////    main_ccvt.set_domain(mu_x, mu_y, var_x, var_y);
////
////    main_ccvt.set_custom_proportions(custom_capacities);
////    main_ccvt.set_neightbour_proportions(custom_neightbour_capacities);
////
////    main_ccvt.set_initial_sites(init_sites);
//////    main_ccvt.generate_random_sites(nb_site);
////
////    main_ccvt.toggle_verbose(); // affichage des positions, poids, volumes et voisinages à chaque étapes
////    main_ccvt.toggle_step_by_step();  // appuyer sur "entrée" pour passer à l'étape suivante
////
////
////    main_ccvt.save_cell_eps(working_directory + "initial_cells.eps");
////
////
////
////    main_ccvt.verbose();
////
////
////    unsigned iter = main_ccvt.optimize_H(stepW, stepX, max_newton_iters, epsilon, max_iters);
////    std::cout<< "Total: "<<iter<<" iters "<<std::endl;
////
////    main_ccvt.verbose();
//
//
//
//
//
//
//
//
//
//
//    // ------------------------------------------------------------------------------------------------
//    for(unsigned int s=0; s<1; s++){
//        std::cout<<"seed : "<<s+seed<<std::endl;
//
//        CCVT main_ccvt{s+seed};
//        main_ccvt.set_domain(mu_x, mu_y, var_x, var_y);
//
//        main_ccvt.set_custom_proportions(custom_capacities);
////        main_ccvt.set_neightbour_proportions(custom_neightbour_capacities);
//
////        main_ccvt.toggle_verbose();
//
//        main_ccvt.generate_random_sites(nb_site);
//
//
//        main_ccvt.verbose();
//
//        unsigned iter = main_ccvt.optimize_all(stepW, stepX, max_newton_iters, epsilon, max_iters, std::cout);
//        seed = iter;
//
//        main_ccvt.verbose();
//
//
//
//        std::vector<FT> weights;
//        std::vector<Point> points;
//        std::vector<FT> capacities;
//        std::vector<FT> areas;
//
//        main_ccvt.collect_sites(points, weights);
//        areas = main_ccvt.get_area();
//        capacities = main_ccvt.get_capacities();
//
//        std::cout<<"position (x,y), poids, aires (volumes objectifs)"<<std::endl;
//        for(int i=0; i<points.size(); i++){
//            std::cout<<"("<<points.at(i).x()<<", "<<points.at(i).y()<<"), "<<weights.at(i)<<", "<<areas.at(i)<<" ("<<capacities.at(i)<<")"<<std::endl;
//        }
//        std::cout<<std::endl;
//
//        std::cout<<"graines"<<std::endl;
//        for(int i=0; i<points.size(); i++){
//            std::cout<<"Graine("<<points.at(i).x()+0.5<<", "<<points.at(i).y()+0.5<<", "<<weights.at(i)<<"),"<<std::endl;
//        }
//        std::cout<<std::endl;
//        std::cout<<"================================================================================="<<std::endl;
//    }
//
//
//
//
//
//
//
//
//    return 0;
//}