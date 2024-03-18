//
// Created by grenier on 22/08/23.
//
#include <iostream>

#include "ASTex/image_gray.h"
#include "ASTex/image_common.h"
#include "ASTex/Noises/Gabor.h"

#include "ASTex/CCVT_tests/numerical_Tstat.h"
#include "ASTex/CCVT_tests/numerical_Tcontent.h"
#include "ASTex/CCVT_tests/tools.h"
#include "ASTex/CCVT_tests/data.h"
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
    ImageGrayd cm_(cm_size, cm_size);
    ImageGrayd res_composition(img_size, img_size);

    // ---------------------------------------------------------------------------
    for(auto s=CCVT_seeds.begin(); s<CCVT_seeds.end(); s++){
        std::vector<Graine> H_seeds = *s;
        std::vector<double> H_color{0.0, 0.2, 0.4, 0.6, 0.8, 1.0};


        cm_.parallel_for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                                    {
                                        double dist = 100.;
                                        int id_seed = -1;

                                        for(int i=0; i<H_seeds.size(); i++){
                                            double X = double(x)/(cm_size-1);
                                            double Y = double(y)/(cm_size-1);

//                                            double new_dist = (H_seeds.at(i).x_ - X)*(H_seeds.at(i).x_ - X) +
//                                                              (H_seeds.at(i).y_ - Y)*(H_seeds.at(i).y_ - Y) -
//                                                              H_seeds.at(i).weight_;

                                            double new_dist = (X - H_seeds.at(i).x_)*(X - H_seeds.at(i).x_) +
                                                              (Y - H_seeds.at(i).y_)*(Y - H_seeds.at(i).y_) -
                                                              H_seeds.at(i).weight_;

                                            if(new_dist<dist){
                                                dist = new_dist;
                                                id_seed = i;
                                            }
                                            P = H_color.at(id_seed);
                                        }
                                    });
//        save_pgm(working_directory + "color_map_ccvt.pgm", cm_);


        // ---------------------------------------------------------------------------

        res_composition.parallel_for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                                                {
                                                    // valeur sur [0, 1]
                                                    double n1 = image_1.pixelAbsolute(x,y);
                                                    double n2 = image_2.pixelAbsolute(x,y);

                                                    int id_n1 = int(std::round(n1*(cm_size-1)));
                                                    int id_n2 = int(std::round(n2*(cm_size-1)));

                                                    P = cm_.pixelAbsolute(id_n1, id_n2);
                                                });
//        save_pgm(working_directory + "res_composition_ccvt.pgm", res_composition);





        // ---------------------------------------------------------------------------

        // proportion présence dans E
//        std::cout<<"mesure présence dans E"<<std::endl;
        std::vector<color_info> couleurs_E = Tcontent(res_composition);
        std::sort(couleurs_E.begin(), couleurs_E.end());

        double tot_pixel_E = 0.;
        for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
        {
            tot_pixel_E += (*it).compteur_;
        }



        // capacité présence dans H
//        std::cout<<"estimation présence dans H"<<std::endl;
//    std::vector<color_info> couleurs_H = cell_capacity_real_histo(cm_, histo_N1_N2);
        std::vector<color_info> couleurs_H = cell_capacity(cm_, mu_1, mu_2, var_1, var_2);
        std::sort(couleurs_H.begin(), couleurs_H.end());

        double tot_pixel_H = 0.;
        for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
        {
            tot_pixel_H += (*it).compteurD_;
        }




        // affichage
        if(couleurs_E.size() != H_color.size()){
            std::cout<<"fail : ";
            for(auto g=H_seeds.begin(); g<H_seeds.end(); g++){
                std::cout<<"("<<(*g).x_<<", "<<(*g).y_<<", "<<(*g).weight_<<")";
            }
            std::cout<<std::endl;
        }
        else{
//            std::cout<<"erreurs présence ccvt (%)"<<std::endl;
            std::cout<<"couleurs & mesure E & estimation H & erreurs \\\\"<<std::endl;
//            std::cout<<"\\hline"<<std::endl;
            for(int i=0; i<couleurs_E.size(); i++){
                double couleur = couleurs_E.at(i).couleur_;
                double proportion_E = couleurs_E.at(i).compteur_/tot_pixel_E;
                double proportion_H = couleurs_H.at(i).compteurD_/tot_pixel_H;
                double err = std::abs(proportion_E - proportion_H)/proportion_E;

                std::cout<<couleur<<" & "<< proportion_E<<" & "<< proportion_H<<" & "<< 100.*err<<" \\\\"<<std::endl;
            }
            std::cout<<std::endl;
        }



    }


}