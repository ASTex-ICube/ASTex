//
// Created by grenier on 09/10/23.
//

#include <iostream>

#include "ccvt.h"
//#include "types.h"

struct comput_err{ // foncteur de calcul d'erreur
    double operator()(const double &ref, const double &est) const{
        return 100*std::abs(ref-est)/ref;
    }
};



int main()
{

    std::string working_directory = "/home/grenier/Documents/ASTex_fork/results/CCVT_CGAL/";

    FT stepX = 0.01; // pour les positions
    FT stepW = 0.1; // pour les poids
    FT epsilon = 1.;
    unsigned max_newton_iters = 500;
    unsigned max_iters = 500;
    unsigned nb_site = 6;
    unsigned seed = 24;

    // description de la densité
    double mu_x = 0.500052;
    double mu_y = 0.500039;
    double var_x = 0.015029;
    double var_y = 0.0156462;

    // proportion de présence des couleurs
    std::vector<FT> custom_capacities{0.0947948, 0.162811, 0.0396879, 0.138455, 0.149857, 0.414395}; // proportion de capacité objectif
    assert(custom_capacities.size() == nb_site); // une proportion par graine
    assert(std::abs(std::accumulate(custom_capacities.begin(), custom_capacities.end(), 0.0) - 1.) < 0.000001); // somme des proportion = 1

    // proportion de voisinage entre les couleurs
    std::vector<std::vector<FT>> custom_neightbour_capacities{{0.0933725, 0.0, 0.00025523, 0.000252651, 0.0, 0.000913634},
                                                              {0.0, 0.160857, 0.0, 0.000311421, 0.000335418, 0.00130924},
                                                              {0.00025523, 0.0, 0.0388644, 0.0, 0.000313686, 0.00025593},
                                                              {0.000252651, 0.000311421, 0.0, 0.136678, 0.0, 0.00121814},
                                                              {0.0, 0.000335418, 0.000313686, 0.0, 0.147925, 0.00128302},
                                                              {0.000913634, 0.00130924, 0.00025593, 0.00121814, 0.00128302, 0.409406}};
    assert(custom_neightbour_capacities.size() == nb_site);
    for(int cap=0; cap<custom_neightbour_capacities.size(); cap++){
        assert(custom_neightbour_capacities.at(cap).size() == nb_site);
        assert(std::abs(std::accumulate(custom_neightbour_capacities.at(cap).begin(), custom_neightbour_capacities.at(cap).end(), 0.0) - custom_capacities.at(cap)) < 0.01);
    }

    // position initial des graines
    std::vector<Point> init_sites{Point(0.81, 0.65),
                                  Point(0.15, 0.25),
                                  Point(0.52, 0.91),
                                  Point(0.65, 0.28),
                                  Point(0.12, 0.72),
                                  Point(0.45, 0.48)};
    assert(init_sites.size() == nb_site); // une position initial par graine



//    std::cout<<"graine (position, poids) & erreur moyenne & erreur min & erreur max"<<std::endl;
//    for(unsigned loop=0; loop<100; loop++){
//        CCVT main_ccvt{};
//        main_ccvt.load_image(working_directory + "density_theo.pgm");
//        main_ccvt.set_custom_proportions(custom_capacities);
//
//        main_ccvt.generate_random_sites(nb_site);
//        std::vector<FT> capacities = main_ccvt.get_capacities(); // capacité objectif
//        main_ccvt.optimize_all(stepW, stepX, max_newton_iters, epsilon, max_iters, std::cout);
//        std::vector<FT> areas_f= main_ccvt.get_area(); // aires finales
//
//
//        std::vector<FT> weights_f;
//        std::vector<Point> points_f;
//        main_ccvt.collect_sites(points_f, weights_f);
//
//        std::cout<<"std::vector<Graine>{";
//        for(int i=0; i<points_f.size(); i++){
//            std::cout<<"Graine("<<points_f.at(i).x()+0.5<<", "<<points_f.at(i).y()+0.5<<", "<<weights_f.at(i)<<"), ";
//        }
//        std::cout<<"}, & ";
//
//        // erreur
//        std::transform (capacities.begin(), capacities.end(), areas_f.begin(), capacities.begin(), comput_err());
//        double moy = std::accumulate(capacities.begin(), capacities.end(), 0.0)/capacities.size();
//        auto [min, max] = std::minmax_element(capacities.begin(), capacities.end());
//        std::cout<<moy<<" & "<<*min<<" & "<<*max<<std::endl;
//
//    }





//    std::cout<<"graine (position, poids) & erreur moyenne & erreur min & erreur max"<<std::endl;
//    for(unsigned loop=0; loop<100; loop++){
//        CCVT main_ccvt{};
//        main_ccvt.load_image(working_directory + "density_theo.pgm");
//        main_ccvt.set_custom_proportions(custom_capacities);
//
//        main_ccvt.generate_random_sites(nb_site);
//
//        std::vector<FT> weights;
//        std::vector<Point> points;
//        main_ccvt.collect_sites(points, weights);
//
//        std::vector<FT> capacities = main_ccvt.get_capacities(); // capacité objectif
//        main_ccvt.optimize_all(stepW, stepX, max_newton_iters, epsilon, max_iters, std::cout);
//        std::vector<FT> areas_f= main_ccvt.get_area(); // aires finales
//
//
//        std::vector<FT> weights_f;
//        std::vector<Point> points_f;
//        main_ccvt.collect_sites(points_f, weights_f);
//
//        // erreur
//        std::transform (capacities.begin(), capacities.end(), areas_f.begin(), capacities.begin(), comput_err());
//        double moy = std::accumulate(capacities.begin(), capacities.end(), 0.0)/capacities.size();
//        auto [min, max] = std::minmax_element(capacities.begin(), capacities.end());
//
//        if(moy<2.){
//            std::cout<<"result : {";
//            for(int i=0; i<points_f.size(); i++){
//                std::cout<<"Graine("<<points_f.at(i).x()+0.5<<", "<<points_f.at(i).y()+0.5<<", "<<weights_f.at(i)<<"), ";
//            }
//            std::cout<<"}, & ";
//
//        }
//        else{
//            std::cout<<"fail : {";
//            for(int i=0; i<points.size(); i++){
//                std::cout<<"Graine("<<points.at(i).x()<<", "<<points.at(i).y()<<", "<<weights.at(i)<<"), ";
//            }
//            std::cout<<"}, & ";
//        }
//        std::cout<<moy<<" & "<<*min<<" & "<<*max<<std::endl;
//
//    }











    CCVT main_ccvt{seed};

//    main_ccvt.load_image(working_directory + "density_theo.pgm");
    main_ccvt.set_domain(mu_x, mu_y, var_x, var_y);

    main_ccvt.set_custom_proportions(custom_capacities);
    main_ccvt.set_neightbour_proportions(custom_neightbour_capacities);
    main_ccvt.set_initial_sites(init_sites);

    main_ccvt.toggle_verbose(); // affichage des positions, poids, volumes et voisinages à chaque étapes
    main_ccvt.toggle_step_by_step();  // appuyer sur "entrée" pour passer à l'étape suivante

//    main_ccvt.generate_random_sites(nb_site);
//    main_ccvt.save_point_eps(working_directory + "generated_random_sites.eps");
    main_ccvt.save_cell_eps(working_directory + "initial_cells.eps");


    // ------------------------------------------------------------------------------------------------
//    std::cout<<std::endl;
//
//    std::vector<FT> weights;
//    std::vector<Point> points;
//    std::vector<FT> capacities;
//    std::vector<FT> areas;
//    std::vector<std::vector<FT>> neightbour_proportion;
//
//    main_ccvt.collect_sites(points, weights);
//    capacities = main_ccvt.get_capacities();
//    areas = main_ccvt.get_area();
//    neightbour_proportion = main_ccvt.get_neightbour_proportion();
//
//    std::cout<<"entrée (position, poids)"<<std::endl;
//    for(int i=0; i<points.size(); i++){
//        std::cout<<"Graine("<<points.at(i).x()<<", "<<points.at(i).y()<<", "<<weights.at(i)<<"),"<<std::endl;
//    }
//    std::cout<<std::endl;
//
//    std::cout<<"entrée (capacities)"<<std::endl;
//    for(auto it=capacities.begin(); it<capacities.end(); it++){
//        std::cout<<(*it)<<"; ";
//    }
//    std::cout<<std::endl;
//
//    std::cout<<"entrée (area)"<<std::endl;
//    for(auto it=areas.begin(); it<areas.end(); it++){
//        std::cout<<(*it)<<"; ";
//    }
//    std::cout<<std::endl;
//    std::cout<<std::endl;
//
//    std::cout<<"entrée (obj voisinage)"<<std::endl;
//    double rho = main_ccvt.compute_value_integral();
//    for(auto it=custom_neightbour_capacities.begin(); it<custom_neightbour_capacities.end(); it++){
//        for(auto it2=(*it).begin(); it2<(*it).end(); it2++){
//            std::cout<<rho*(*it2)<<"; ";
//        }
//        std::cout<<std::endl;
//    }
//    std::cout<<std::endl;
//
//    std::cout<<"entrée (voisinage actuel)"<<std::endl;
//    for(auto it=neightbour_proportion.begin(); it<neightbour_proportion.end(); it++){
//        for(auto it2=(*it).begin(); it2<(*it).end(); it2++){
//            std::cout<<(*it2)<<"; ";
//        }
//        std::cout<<std::endl;
//    }
//
//
//    std::cout<<std::endl;
//    std::cout<<std::endl;

    main_ccvt.verbose();


    // ------------------------------------------------------------------------------------------------
//    unsigned iter = main_ccvt.optimize_all(stepW, stepX, max_newton_iters, epsilon, max_iters, std::cout);
    unsigned iter = main_ccvt.optimize_H(stepW, stepX, max_newton_iters, epsilon, max_iters);
    std::cout<< "Total: "<<iter<<" iters "<<std::endl;

//    main_ccvt.save_point_eps(working_directory + "optimized_sites.eps");
    // ------------------------------------------------------------------------------------------------

    main_ccvt.verbose();

//    std::cout<<std::endl;
//
//    std::vector<FT> weights_f;
//    std::vector<Point> points_f;
//    std::vector<FT> capacities_f;
//    std::vector<FT> areas_f;
//    std::vector<std::vector<FT>> neightbour_proportion_f;
//
//    main_ccvt.collect_sites(points_f, weights_f);
//    capacities_f = main_ccvt.get_capacities();
//    areas_f = main_ccvt.get_area();
//    neightbour_proportion_f = main_ccvt.get_neightbour_proportion();
//
//
//    std::cout<<"sortie (position, poids)"<<std::endl;
//    for(int i=0; i<points_f.size(); i++){
//        std::cout<<"Graine("<<points_f.at(i).x()+0.5<<", "<<points_f.at(i).y()+0.5<<", "<<weights_f.at(i)<<"),"<<std::endl;
//    }
//    std::cout<<std::endl;
//
//    std::cout<<"sortie (capacities)"<<std::endl;
//    for(auto it=capacities_f.begin(); it<capacities_f.end(); it++){
//        std::cout<<(*it)<<"; ";
//    }
//    std::cout<<std::endl;
//
//    std::cout<<"sortie (area)"<<std::endl;
//    for(auto it=areas_f.begin(); it<areas_f.end(); it++){
//        std::cout<<(*it)<<"; ";
//    }
//    std::cout<<std::endl;
//    std::cout<<std::endl;
//
//    std::cout<<"sortie (obj voisinage)"<<std::endl;
//    double rho_f = main_ccvt.compute_value_integral();
//    for(auto it=custom_neightbour_capacities.begin(); it<custom_neightbour_capacities.end(); it++){
//        for(auto it2=(*it).begin(); it2<(*it).end(); it2++){
//            std::cout<<rho_f*(*it2)<<"; ";
//        }
//        std::cout<<std::endl;
//    }
//    std::cout<<std::endl;
//
//    std::cout<<"sortie (voisinage actuel)"<<std::endl;
//    for(auto it=neightbour_proportion_f.begin(); it<neightbour_proportion_f.end(); it++){
//        for(auto it2=(*it).begin(); it2<(*it).end(); it2++){
//            std::cout<<(*it2)<<"; ";
//        }
//        std::cout<<std::endl;
//    }
//
//    std::cout<<std::endl;
//
//
//
//    std::transform (capacities.begin(), capacities.end(), areas_f.begin(), capacities.begin(), comput_err());
//    double moy = std::accumulate(capacities.begin(), capacities.end(), 0.0)/capacities.size();
//    auto [min, max] = std::minmax_element(capacities.begin(), capacities.end());
//    std::cout<<moy<<", "<<*min<<", "<<*max<<std::endl;





    return 0;
}