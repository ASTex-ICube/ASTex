//
// Created by grenier on 09/10/23.
//

#include <iostream>

#include "ccvt.h"
#include "ccvt_application.h"


int main() {
    // description de la densité
    double mu_x = 0.5;
    double mu_y = 0.5;
    double var_x = 0.02;//4025; // 0.0151757;//
    double var_y = 0.02;//4025; // 0.0151757;//

    unsigned size_x = 32;// 256;
    unsigned size_y = 32;// 256;
    double max_val = 255.;

    // position initial des graines
    std::vector<Point> init_sites{Point(0.81, 0.65), Point(0.15, 0.25), Point(0.52, 0.91), Point(0.65, 0.28), Point(0.12, 0.72), Point(0.45, 0.48)};
    std::vector<FT> custom_capacities{0.1313, 0.0534, 0.0546, 0.1644, 0.0493, 0.5469};

    CCVT main_ccvt;
    main_ccvt.set_domain(mu_x, mu_y, var_x, var_y, size_x, size_y, max_val);
    main_ccvt.set_custom_proportions(custom_capacities);
    main_ccvt.set_initial_sites(init_sites);
    // main_ccvt.toggle_timer();

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
//    std::string working_directory = TEMPO_PATH+"results/CCVT_CGAL/";
//
//    FT stepX = 0.01; // pour les positions
//    FT stepW = 0.1; // pour les poids
//    FT epsilon = 1.;
//    unsigned max_newton_iters = 500;
//    unsigned max_iters = 500;
//    unsigned nb_site = 6;
//    unsigned seed = 624;
//
//    // description de la densité
//    double mu_x = 0.5;// 0.500233;
//    double mu_y = 0.5;//0.499979;
//    double var_x = 0.024025;// 0.0151757;
//    double var_y = 0.024025;// 0.0198063;
//
//    // proportion de présence des couleurs
//    std::vector<FT> custom_capacities{0.0947948, 0.162811, 0.0396879, 0.138455, 0.149857, 0.414395}; // proportion de capacité objectif
//    assert(custom_capacities.size() == nb_site); // une proportion par graine
//    assert(std::abs(std::accumulate(custom_capacities.begin(), custom_capacities.end(), 0.0) - 1.) < 0.000001); // somme des proportion = 1
//
//    // position initial des graines
//    std::vector<Point> init_sites{Point(0.81, 0.65),
//                                  Point(0.15, 0.25),
//                                  Point(0.52, 0.91),
//                                  Point(0.65, 0.28),
//                                  Point(0.12, 0.72),
//                                  Point(0.45, 0.48)};
//
//    assert(init_sites.size() == nb_site); // une position initial par graine
//
//
//
//
//
//
//
//
//    // ------------------------------------------------------------------------------------------------
//    CCVT main_ccvt{seed};
//
//    main_ccvt.set_domain(mu_x, mu_y, var_x, var_y);
//
//    main_ccvt.set_custom_proportions(custom_capacities);
//
//    main_ccvt.set_initial_sites(init_sites);
//
//    main_ccvt.toggle_verbose(); // affichage des positions, poids, volumes et voisinages à chaque étapes
//    main_ccvt.toggle_step_by_step();  // appuyer sur "entrée" pour passer à l'étape suivante
//
//
//    main_ccvt.save_cell_eps(working_directory + "initial_cells.eps");
//
//
//
//    main_ccvt.verbose();
//
//
//    unsigned iter = main_ccvt.optimize_H(stepW, stepX, max_newton_iters, epsilon, max_iters);
//    std::cout<< "Total: "<<iter<<" iters "<<std::endl;
//
//    main_ccvt.verbose();
//
//
//
//
//    return 0;
//}