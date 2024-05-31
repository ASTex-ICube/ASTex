//
// Created by grenier on 22/08/23.
//
#include <iostream>

#include "ASTex/image_gray.h"
#include "ASTex/image_common.h"
#include "ASTex/Noises/Gabor.h"

#include "ASTex/CCVT_tests/mesure_statistiques.h"
#include "ASTex/CCVT_tests/numerical_Tcontent.h"
#include "ASTex/CCVT_tests/numerical_Tvois.h"
#include "ASTex/CCVT_tests/large_numerical_test.h"
#include "ASTex/CCVT_tests/tools.h"
//#include "CCVT_CGAL/ccvt.h"
//#include "CCVT_CGAL/types.h"


using namespace ASTex;








int main(){
    std::string working_directory = "/home/grenier/Documents/ASTex_fork/results/CCVT_CGAL/";

    // ---------------------------------------------------------------------------
    int resolution = 256;
    int img_size = 4096;//4096;//2048; // nombre de pixel dans l'image
    int cm_size = 256;//256;//512;

    float number_of_impulses_per_kernel = 64.;// 64.0;
    unsigned period = 128; // non utilisé

    float K_ = 1.0; //laisser à 1
    float a_ = 0.01; // taille des noyaux (a*a = 1/variance)0.02

    unsigned random_offset_ = 954248632;
    unsigned seed_ = 4;


    // ---------------------------------------------------------------------------
    noise noise_1(K_, // isotrope
                  a_,
                  0.06,//0.012, //F_0_,
                  0.785, //omega_0_,
                  number_of_impulses_per_kernel,
                  period,
                  random_offset_,
                  seed_);
    ImageGrayd image_1 = storing_noise_d(resolution, img_size, noise_1);
    save_pgm(working_directory + "gabor1.pgm", image_1);


    noise noise_2(K_, // anisotrope
                  a_,
                  0.06,//0.02, //F_0_,
                  3., //omega_0_,
                  number_of_impulses_per_kernel,
                  period,
                  random_offset_-4,
                  seed_+12);
    ImageGrayd image_2 = storing_noise_d(resolution, img_size, noise_2);
    save_pgm(working_directory + "gabor2.pgm", image_2);


    // ---------------------------------------------------------------------------
    // statistiques
    double mu_1 = moyenne(image_1); // mu
    double mu_2 = moyenne(image_2); // mu'

    double mu_1_carre = moyenne_carre(image_1);
    double mu_2_carre = moyenne_carre(image_2);

    double var_1 = mu_1_carre - mu_1*mu_1; // sigma
    double var_2 = mu_2_carre - mu_2*mu_2; // sigma'

    // histogrammes
    ImageGrayd histo_N1_N2 = histo_2D(image_1, image_2, cm_size); // réel
    save_pgm(working_directory + "density_real.pgm", histo_N1_N2);

    ImageGrayd histo_theo = histo_2D_theo(mu_1, mu_2, var_1, var_2, cm_size); // théorique
    save_pgm(working_directory + "density_theo.pgm", histo_theo);




    // ---------------------------------------------------------------------------
    std::vector<Graine> H_seeds{Graine(0.81, 0.65, 0.),
                                Graine(0.15, 0.25, 0.1),
                                Graine(0.52, 0.91 ,0.),
                                Graine(0.65, 0.28, 0.),
                                Graine(0.12, 0.72, 0.1),
                                Graine(0.45, 0.48, 0.)};
//    std::vector<Graine> H_seeds{Graine(0.81, 0.65, -0.0344047),
//                                Graine(0.15, 0.25, 0.0645013),
//                                Graine(0.52, 0.91 ,0.019120),
//                                Graine(0.65, 0.28, -0.036392),
//                                Graine(0.12, 0.72, 0.06444),
//                                Graine(0.45, 0.48, -0.0390251)};
    std::vector<double> H_color{0.0, 0.2, 0.4, 0.6, 0.8, 1.0};


    ImageGrayd cm_(cm_size, cm_size);
    cm_.parallel_for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                                {
                                    double dist = 100.;
                                    int id_seed = -1;

                                    for(int i=0; i<H_seeds.size(); i++){
                                        double X = double(x)/(cm_size-1);
                                        double Y = double(y)/(cm_size-1);

                                        double new_dist = (H_seeds.at(i).x_ - X)*(H_seeds.at(i).x_ - X) +
                                                          (H_seeds.at(i).y_ - Y)*(H_seeds.at(i).y_ - Y) -
                                                          H_seeds.at(i).weight_;
                                        if(new_dist<dist){
                                            dist = new_dist;
                                            id_seed = i;
                                        }
                                        P = H_color.at(id_seed);
                                    }
                                });
    save_pgm(working_directory + "color_map_ref.pgm", cm_);


    // ---------------------------------------------------------------------------
    ImageGrayd res_composition(img_size, img_size);
    res_composition.parallel_for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                                            {
                                                // valeur sur [0, 1]
                                                double n1 = image_1.pixelAbsolute(x,y);
                                                double n2 = image_2.pixelAbsolute(x,y);

                                                int id_n1 = int(std::round(n1*(cm_size-1)));
                                                int id_n2 = int(std::round(n2*(cm_size-1)));

                                                P = cm_.pixelAbsolute(id_n1, id_n2);
                                            });
    save_pgm(working_directory + "res_composition_ref.pgm", res_composition);





    // ---------------------------------------------------------------------------

//    // proportion présence dans E
//    std::cout<<"mesure présence dans E"<<std::endl;
//    std::vector<color_info> couleurs_E = Tcontent(res_composition);
//    std::sort(couleurs_E.begin(), couleurs_E.end());
//
//    double tot_pixel_E = 0.;
//    for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
//    {
//        tot_pixel_E += (*it).compteur_;
//    }
//
//
////
////    // capacité présence dans H
////    std::cout<<"estimation présence dans H"<<std::endl;
//////    std::vector<color_info> couleurs_H = cell_capacity_real_histo(cm_, histo_N1_N2);
////    std::vector<color_info> couleurs_H = cell_capacity(cm_, mu_1, mu_2, var_1, var_2);
////    std::sort(couleurs_H.begin(), couleurs_H.end());
////
////    double tot_pixel_H = 0.;
////    for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
////    {
////        tot_pixel_H += (*it).compteurD_;
////    }
//
//
//
//
//    // affichage
//    std::cout<<"erreurs présence référence (%)"<<std::endl;
////    std::cout<<"couleurs & mesure E & estimation H & erreurs \\\\"<<std::endl;
////    std::cout<<"\\hline"<<std::endl;
//    for(int i=0; i<couleurs_E.size(); i++){
////        double couleur = couleurs_E.at(i).couleur_;
//        double proportion_E = couleurs_E.at(i).compteur_/tot_pixel_E;
////        double proportion_H = couleurs_H.at(i).compteurD_/tot_pixel_H;
////        double err = std::abs(proportion_E - proportion_H)/proportion_E;
//
////        std::cout<<couleur<<" & "<< proportion_E<<" & "<< proportion_H<<" & "<< 100.*err<<" \\\\"<<std::endl;
//        std::cout<<proportion_E<<", ";
//    }
//    std::cout<<std::endl;







    std::vector<color_info> couleurs_E;
    std::vector<color_vois> couleur_vois_E;
    int nb_voisinage = 0;
    int nb_voisinage_diff = 0;

    mesure_E(res_composition, couleurs_E, couleur_vois_E, true, nb_voisinage, nb_voisinage_diff);

    std::sort(couleurs_E.begin(), couleurs_E.end());
    std::sort(couleur_vois_E.begin(), couleur_vois_E.end());

    double tot_pixel_E = 0.;
    for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
    {
        tot_pixel_E += (*it).compteur_;
    }


    // affichage
    std::cout<<"présence"<<std::endl;
    for(int i=0; i<couleurs_E.size(); i++){
        double couleur = couleurs_E.at(i).couleur_;
        double proportion_E = couleurs_E.at(i).compteur_/tot_pixel_E;

        std::cout<<couleur<<" & "<< proportion_E<<std::endl;
    }
    std::cout<<std::endl;


    std::cout<<"voisinages"<<std::endl;
    for(int i=0; i<couleur_vois_E.size(); i++){
        double couleur1 = couleur_vois_E.at(i).couleur1_;
        double couleur2 = couleur_vois_E.at(i).couleur2_;
        double proportion_vois_E;


        if(couleur1 == couleur2){
            proportion_vois_E = couleur_vois_E.at(i).compteur_/double(nb_voisinage);
        }
        else{
            proportion_vois_E = couleur_vois_E.at(i).compteur_/(2.*double(nb_voisinage));
        }

        std::cout<<couleur1<<", "<<couleur2<<" & "<< proportion_vois_E<<std::endl;
    }
    std::cout<<std::endl;






}