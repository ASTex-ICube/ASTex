//
// Created by grenier on 17/11/23.
//
//
// Created by grenier on 18/09/23.
//

#include <chrono>
#include <iostream>
#include <string>
#include <fstream>

#include <ASTex/image_gray.h>
#include "ASTex/Noises/Gabor.h"
#include "ASTex/fourier.h"

#include "ASTex/CCVT_tests/numerical_Tcontent.h"
#include "ASTex/CCVT_tests/numerical_Tvois.h"
#include "ASTex/CCVT_tests/numerical_Tstat.h"
#include "ASTex/CCVT_tests/numerical_Treynold.h"
#include "ASTex/CCVT_tests/large_numerical_test.h"

using namespace ASTex;
typedef ImageGrayd ImgType;

// Gabor.h est le code exemple du papier d'origine dont il est un peut touffu (aussi à cause des structure accélératrices)
// les principales partie de l'algo sont :
// fonction gabor : définition d'un noyaux de gabor (gaussienne * cosinus)
// fonction cell : boucle pour placer et sommer les noyaux




//int main()
//{
//    auto start = std::chrono::high_resolution_clock::now();
//    large_test(1000);
//    auto stop = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
//    std::cout << duration.count()<<" secondes" << std::endl;
//    return EXIT_SUCCESS;
//}







int main()
{
// ---------------------------------------------------------------------------
    int resolution = 256;
    int img_size = 256;//4096;//2048; // nombre de pixel dans l'image
    int cm_size = 256;//256;//512;
    double scale = double(resolution)/double(img_size);

    float F_0_ = 0.06;//0.04; // fréquence
    float omega_0_ = 0.;//M_PI/4.; // orientation (seulement dans le cas anisotrope, cf code gabor ligne 160)

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
                  4.7, //omega_0_,
                  number_of_impulses_per_kernel,
                  period,
                  random_offset_,
                  seed_);
// la valeur du bruit noise_ en x,y peut être récupérée par noise_(x,y)
    ImgType image_1 = storing_noise_d(resolution, img_size, noise_1); // pour écrire le bruit dans une image
    IO::save(image_1, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_1.png");


    noise noise_2(K_, // anisotrope
                  a_,
                  0.06,//0.02, //F_0_,
                  3.8, //omega_0_,
                  number_of_impulses_per_kernel,
                  period,
                  random_offset_-4,
                  seed_+12);
// la valeur du bruit noise_ en x,y peut être récupérée par noise_(x,y)
    ImgType image_2 = storing_noise_d(resolution, img_size, noise_2); // pour écrire le bruit dans une image
    IO::save(image_2, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_2.png");


    // histogrammes
    ImgType histo_N1_N2 = histo_2D(image_1, image_2, cm_size); // réel








    // ---------------------------------------------------------------------------
    std::vector<Graine> H_seeds{Graine(0.15, 0.25, 0.1),
                                Graine(0.65, 0.28, 0.),
                                Graine(0.81, 0.65, 0.),
                                Graine(0.52, 0.91 ,0.),
                                Graine(0.12, 0.72, 0.1),
                                Graine(0.45, 0.48, 0.)};
    std::vector<double> H_color{0.2, 0.6, 0.0, 0.4, 0.8, 1.0};


//    std::vector<Graine> H_seeds{Graine(0.14, 0.69, 0.),
//                                Graine(0.43, 0.83, 0.),
//                                Graine(0.81, 0.72, 0.),
//                                Graine(0.12, 0.24 ,0.),
//                                Graine(0.46, 0.36, 0.),
//                                Graine(0.80, 0.20, 0.)};
//    std::vector<double> H_color{0.2, 0.6, 0.0, 0.4, 0.8, 1.0};


//    std::vector<Graine> H_seeds{Graine(0.15, 0.25, 0.1),
//                                Graine(0.81, 0.65, 0.)};
//    std::vector<double> H_color{0.2, 0.0};


//    std::vector<Graine> H_seeds{Graine(0.1, 0.1, 0.1),
//                                Graine(0.1, 0.9, 0.1),
//                                Graine(0.9, 0.5, 0.1)};
//    std::vector<double> H_color{0.2, 0.6, 0.0};

//    std::vector<Graine> H_seeds{Graine(0.9, 0.1, 0.1),
//                                Graine(0.9, 0.9, 0.1),
//                                Graine(0.1, 0.5, 0.1)};
//    std::vector<double> H_color{0.2, 0.6, 0.0};

//    std::vector<Graine> H_seeds{Graine(0.81, 0.65, 0.),
//                                Graine(0.15, 0.25, 0.1),
//                                Graine(0.52, 0.91 ,0.),
//                                Graine(0.65, 0.28, 0.),
//                                Graine(0.12, 0.72, 0.1),
//                                Graine(0.45, 0.48, 0.)};
//    std::vector<double> H_color{0.0, 0.2, 0.4, 0.6, 0.8, 1.0};


    ImgType cm_(cm_size, cm_size);

    cm_.parallel_for_all_pixels([&] (typename ImgType::PixelType& P, int x, int y)
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
    IO::save(cm_, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/cm.png");


    // ---------------------------------------------------------------------------
    ImgType res_composition(img_size, img_size);

//    res_composition.parallel_for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
//                                            {
//                                                double dist = 100.;
//                                                int id_seed = -1;
//
//                                                for(int i=0; i<H_seeds.size(); i++){
//                                                    double X = image_1.pixelAbsolute(x,y);
//                                                    double Y = image_2.pixelAbsolute(x,y);
//
////                                                    double X = std::clamp(image_1.pixelAbsolute(x,y), 0., 1.);
////                                                    double Y = std::clamp(image_2.pixelAbsolute(x,y), 0., 1.);
//
//                                                    double new_dist = (H_seeds.at(i).x_ - X)*(H_seeds.at(i).x_ - X) +
//                                                                      (H_seeds.at(i).y_ - Y)*(H_seeds.at(i).y_ - Y) -
//                                                                      H_seeds.at(i).weight_;
//                                                    if(new_dist<dist){
//                                                        dist = new_dist;
//                                                        id_seed = i;
//                                                    }
//                                                    P = H_color.at(id_seed);
//                                                }
//                                            });
//    IO::save(res_composition, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/res_composition.png");


    //  composition en nutilisant lle cmm enn mémoire (auucunne apprroximmatitonn lors duu calcul ded l'histo après)
    res_composition.parallel_for_all_pixels([&] (typename ImgType::PixelType& P, int x, int y)
                                            {
                                                // valeur sur [0, 1]
                                                double n1 = image_1.pixelAbsolute(x,y);
                                                double n2 = image_2.pixelAbsolute(x,y);

                                                int id_n1 = int(std::round(n1*(cm_size-1)));
                                                int id_n2 = int(std::round(n2*(cm_size-1)));

                                                P = cm_.pixelAbsolute(id_n1, id_n2);
                                            });
    IO::save(res_composition, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/res_composition.png");







    // ---------------------------------------------------------------------------

// statistiques
    double mu_1 = moyenne(image_1); // mu
    double mu_2 = moyenne(image_2); // mu'

//    double mu_1_carre = moyenne_carre(image_1);
//    double mu_2_carre = moyenne_carre(image_2);
//
//    double var_1 = mu_1_carre - mu_1*mu_1; // sigma
//    double var_2 = mu_2_carre - mu_2*mu_2; // sigma'
//
//    double AC1 = covariance(image_1, 1, 0); // sigma_uv
//    double AC2 = covariance(image_2, 1, 0); // sigma_uv'

    ImgType Fourier_1 = covariance_Fourier(image_1);
    ImgType Fourier_2 = covariance_Fourier(image_2);

    double var_1 = Fourier_1.pixelAbsolute(0, 0); // sigma
    double var_2 = Fourier_2.pixelAbsolute(0, 0); // sigma'

    double AC1 = Fourier_1.pixelAbsolute(0, 1); // sigma_uv
    double AC2 = Fourier_2.pixelAbsolute(0, 1); // sigma_uv'


    // affichage
    std::cout<<"statistiques bruits"<<std::endl;

    std::cout<<"bruit N  : "<<std::endl;
    std::cout<<"moyenne = "<<mu_1<<std::endl;
    std::cout<<"variance = "<<var_1<<std::endl;
//    std::cout<<"variance sans biais = "<<var_1*(double(img_size)*double(img_size))/(double(img_size)*double(img_size)-1.)<<std::endl;
    std::cout<<"autocovariance = "<<AC1<<std::endl;
//    std::cout<<"variance Fourier = "<<Fourier_1.pixelAbsolute(0, 0)<<std::endl;
//    std::cout<<"autocovariance Fourier = "<<Fourier_1.pixelAbsolute(1, 0)<<std::endl;

    std::cout<<std::endl;

    std::cout<<"bruit N'  : "<<std::endl;
    std::cout<<"moyenne = "<<mu_2<<std::endl;
    std::cout<<"variance = "<<var_2<<std::endl;
//    std::cout<<"variance sans biais = "<<var_2*(double(img_size)*double(img_size))/(double(img_size)*double(img_size)-1.)<<std::endl;
    std::cout<<"autocovariance = "<<AC2<<std::endl;
//    std::cout<<"variance Fourier = "<<Fourier_2.pixelAbsolute(0, 0)<<std::endl;
//    std::cout<<"autocovariance Fourier = "<<Fourier_2.pixelAbsolute(1, 0)<<std::endl;

    std::cout<<std::endl;


// histogrammes
    // présences
    ImgType histo_theo = histo_2D_theo(mu_1, mu_2, var_1, var_2, cm_size); // théorique

        // voisinages
//    std::cout<<"histogramme voisinage"<<std::endl;
//    ImageGrayd histo_AC_N1 = histo_2D_AC(image_1, 1, 0, 1); // réel
//    ImageGrayd histo_AC_N2 = histo_2D_AC(image_2, 1, 0, 2);



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
////    for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
////    {
////        std::cout<<"("<<(*it).couleur_<<") : "<<(*it).compteur_/tot_pixel_E<<std::endl;
//////        std::cout<<"("<<(*it).couleur_<<") : "<<(*it).compteur_<<std::endl;
////    }
////    std::cout<<std::endl;
//
//
//
//
//    // capacité présence dans H
//    std::cout<<"estimation présence dans H"<<std::endl;
//    std::vector<color_info> couleurs_H = cell_capacity_real_histo(cm_, histo_N1_N2);
////    std::vector<color_info> couleurs_H = cell_capacity(cm_, mu_1, mu_2, var_1, var_2);
//    std::sort(couleurs_H.begin(), couleurs_H.end());
//
//    double tot_pixel_H = 0.;
//    for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
//    {
//        tot_pixel_H += (*it).compteurD_;
//    }
//
//
////    for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
////    {
////        std::cout<<"("<<(*it).couleur_<<") : "<<(*it).compteurD_/tot_pixel_H<<std::endl;
//////        std::cout<<"("<<(*it).couleur_<<") : "<<(*it).compteurD_<<std::endl;
////    }
////    std::cout<<std::endl;
//
//
//
//
//
//    // affichage
//    std::cout<<"erreurs présence (%)"<<std::endl;
//    std::cout<<"couleurs & mesure E & estimation H & erreurs \\\\"<<std::endl;
//    std::cout<<"\\hline"<<std::endl;
//    for(int i=0; i<couleurs_E.size(); i++){
//        double couleur = couleurs_E.at(i).couleur_;
//        double proportion_E = couleurs_E.at(i).compteur_/tot_pixel_E;
//        double proportion_H = couleurs_H.at(i).compteurD_/tot_pixel_H;
//        double err = std::abs(proportion_E - proportion_H)/proportion_E;
//
//        std::cout<<couleur<<" & "<< proportion_E<<" & "<< proportion_H<<" & "<< 100.*err<<" \\\\"<<std::endl;
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



    for(int i=0; i<couleur_vois_E.size(); i++){
        double couleur1 = couleur_vois_E.at(i).couleur1_;
        double couleur2 = couleur_vois_E.at(i).couleur2_;
        double proportion_E;


        if(couleur1 == couleur2){
            proportion_E = couleur_vois_E.at(i).compteur_/double(nb_voisinage);
        }
        else{
            proportion_E = couleur_vois_E.at(i).compteur_/(2.*double(nb_voisinage));
        }

        std::cout<<couleur1<<", "<<couleur2<<" & "<< proportion_E<<std::endl;
    }
    std::cout<<std::endl;





//    // ---------------------------------------------------------------------------
//    bool keep_ieqj = false;
//
//    // proportion voisinage dans E
//    std::cout<<"mesure voisinage dans E"<<std::endl;
//    std::vector<color_vois> couleur_vois_E = Tcontent_vois_test(res_composition, false);
//    std::sort(couleur_vois_E.begin(), couleur_vois_E.end());
//
//
//    // total tous voisinage
//    double tot_E = 0.;
//    for(auto it = couleur_vois_E.begin(); it != couleur_vois_E.end(); it++)
//    {
//        tot_E += (*it).compteur_;
//    }
////    std::cout<<"tot E : "<<tot_E<<std::endl;
////
//    for(auto it = couleur_vois_E.begin(); it != couleur_vois_E.end(); it++)
//    {
////        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") & "<<(*it).compteur_/tot_E<<std::endl;
//        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") & "<<(*it).compteur_<<std::endl;
//
//    }
//    std::cout<<std::endl;
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
//
//    // capacité voisinage dans H
//    std::cout<<"estimation voisinage dans H"<<std::endl;
//    std::vector<color_vois> couleur_vois_H = Hcontent_vois_test(cm_, mu_1, mu_2, var_1, var_2, histo_N1_N2, false);
//    std::sort(couleur_vois_H.begin(), couleur_vois_H.end());
//
//
//    // total tous voisinage
//    double tot_H = 0.;
//    for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
//    {
//        tot_H += (*it).compteurD_;
//    }
////    std::cout<<"tot H : "<<tot_H<<std::endl;
//
//
//    for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
//    {
////        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") & "<<(*it).compteurD_/tot_H<<std::endl;
//        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") & "<<(*it).compteurD_<<std::endl;
//    }
//    std::cout<<std::endl;





//    // affichage
//    std::cout<<"erreurs voisinages (%)"<<std::endl;
//    std::cout<<"couleurs & mesure E & estimation H & erreurs \\\\"<<std::endl;
//    std::cout<<"\\hline"<<std::endl;
//    for(int i=0; i<couleur_vois_E.size(); i++){
//        double couleur_1 = couleur_vois_E.at(i).couleur1_;
//        double couleur_2 = couleur_vois_E.at(i).couleur2_;
//
//        double proportion_E = couleur_vois_E.at(i).compteur_/tot_E;
//        double proportion_H = couleur_vois_H.at(i).compteurD_/tot_H;
//
//        double err = std::abs(proportion_E - proportion_H)/proportion_E;
//
//        std::cout<<"("<<couleur_1<<", "<<couleur_2<<") & "<<proportion_E<<" & "<<proportion_H<<" & "<<100.*err<<" \\\\"<<std::endl;
//
//    }
//    std::cout<<std::endl;










//    // capacité conjointe dans H
//    std::cout<<"estimation voisinage dans H intégrale 4D"<<std::endl;
////    std::vector<color_vois> couleur_vois_H_int = capacityH_vois_histo_reel(cm_, histo_AC_N1, histo_AC_N2);
//    std::vector<color_vois> couleur_vois_H_int = capacityH_vois(cm_, mu_1, mu_2, var_1, var_2, AC1, AC2);
//    std::sort(couleur_vois_H_int.begin(), couleur_vois_H_int.end());
//
//    double tot_H_int = 0.;
//    for(auto it = couleur_vois_H_int.begin(); it != couleur_vois_H_int.end(); it++)
//    {
////        if((*it).couleur1_ != (*it).couleur2_){
//            tot_H_int += (*it).compteurD_;
////        }
//    }
//
//    for(auto it = couleur_vois_H_int.begin(); it != couleur_vois_H_int.end(); it++)
//    {
//        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") & "<<(*it).compteurD_/tot_H_int<<std::endl;
////        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteurD_<<std::endl;
//    }


    return EXIT_SUCCESS;
}

