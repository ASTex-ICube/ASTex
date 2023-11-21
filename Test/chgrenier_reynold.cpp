//
// Created by grenier on 17/11/23.
//
//
// Created by grenier on 18/09/23.
//

#include <chrono>
#include <ASTex/image_gray.h>
#include "ASTex/Noises/Gabor.h"
#include "ASTex/fourier.h"

#include "ASTex/CCVT/numerical_Tcontent.h"
#include "ASTex/CCVT/numerical_Tvois.h"
#include "ASTex/CCVT/numerical_Tstat.h"
#include "ASTex/CCVT/reynold.h"

using namespace ASTex;

// Gabor.h est le code exemple du papier d'origine dont il est un peut touffu (aussi à cause des structure accélératrices)
// les principales partie de l'algo sont :
// fonction gabor : définition d'un noyaux de gabor (gaussienne * cosinus)
// fonction cell : boucle pour placer et sommer les noyaux












int main()
{
// ---------------------------------------------------------------------------
    int resolution = 512;
    int img_size = 2048; // nombre de pixel dans l'image
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
                  0.785, //omega_0_,
                  number_of_impulses_per_kernel,
                  period,
                  random_offset_,
                  seed_);
// la valeur du bruit noise_ en x,y peut être récupérée par noise_(x,y)
    ImageGrayd image_1 = storing_noise_d(resolution, img_size, noise_1); // pour écrire le bruit dans une image
    IO::save(image_1, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_1.png");


    noise noise_2(K_, // anisotrope
                  a_,
                  0.06,//0.02, //F_0_,
                  3., //omega_0_,
                  number_of_impulses_per_kernel,
                  period,
                  random_offset_-4,
                  seed_+12);
// la valeur du bruit noise_ en x,y peut être récupérée par noise_(x,y)
    ImageGrayd image_2 = storing_noise_d(resolution, img_size, noise_2); // pour écrire le bruit dans une image
    IO::save(image_2, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_2.png");











    // ---------------------------------------------------------------------------
    std::vector<Graine> H_seeds{Graine(0.15, 0.25, 0.1),
                                Graine(0.65, 0.28, 0.),
                                Graine(0.81, 0.65, 0.),
                                Graine(0.52, 0.91 ,0.),
                                Graine(0.12, 0.72, 0.1),
                                Graine(0.45, 0.48, 0.)};
    std::vector<double> H_color{0.2, 0.6, 0.0, 0.4, 0.8, 1.0};

    int cm_size = 512;
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
    IO::save(cm_, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/cm.png");


    // ---------------------------------------------------------------------------
    ImageGrayd res_composition(img_size, img_size);

    res_composition.parallel_for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                                            {
                                                double dist = 100.;
                                                int id_seed = -1;

                                                for(int i=0; i<H_seeds.size(); i++){
                                                    double X = image_1.pixelAbsolute(x,y);
                                                    double Y = image_2.pixelAbsolute(x,y);

//                                                    double X = std::clamp(image_1.pixelAbsolute(x,y), 0., 1.);
//                                                    double Y = std::clamp(image_2.pixelAbsolute(x,y), 0., 1.);

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
    IO::save(res_composition, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/res_composition.png");



    // ---------------------------------------------------------------------------

// statistiques
    double mu_1 = moyenne(image_1); // mu
    double mu_2 = moyenne(image_2); // mu'

    double mu_1_carre = moyenne_carre(image_1);
    double mu_2_carre = moyenne_carre(image_2);

    double var_1 = mu_1_carre - mu_1*mu_1; // sigma
    double var_2 = mu_2_carre - mu_2*mu_2; // sigma'

//    double AC1 = covariance(image_1, 1, 0); // sigma_uv
//    double AC2 = covariance(image_2, 1, 0); // sigma_uv'

    std::cout<<"statistiques bruits"<<std::endl;
    std::cout<<"bruit N  : moyenne = "<<mu_1<<", variance = "<<var_1<<std::endl;
    std::cout<<"bruit N' : moyenne = "<<mu_2<<", variance = "<<var_2<<std::endl;
    std::cout<<std::endl;


// histogrammes
    // présences
    ImageGrayd histo_N1_N2 = histo_2D(image_1, image_2, cm_size); // réel
    ImageGrayd histo_theo = histo_2D_theo(mu_1, mu_2, var_1, var_2, cm_size); // théorique



    // ---------------------------------------------------------------------------

    // proportion présence dans E
    std::cout<<"mesure présence dans E"<<std::endl;
    std::vector<color_info> couleurs_E = Tcontent(res_composition);
    std::sort(couleurs_E.begin(), couleurs_E.end());

    double tot_pixel_E = 0.;
    for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
    {
        tot_pixel_E += (*it).compteur_;
    }

    for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
    {
        std::cout<<"("<<(*it).couleur_<<") : "<<(*it).compteur_/tot_pixel_E<<std::endl;
//        std::cout<<"("<<(*it).couleur_<<") : "<<(*it).compteur_<<std::endl;
    }
    std::cout<<std::endl;



    // capacité présence dans H
    std::cout<<"estimation présence dans H"<<std::endl;
//    std::vector<color_info> couleurs_H = cell_capacity_real_histo(cm_, histo_N1_N2);
    std::vector<color_info> couleurs_H = cell_capacity(cm_, mu_1, mu_2, var_1, var_2);
    std::sort(couleurs_H.begin(), couleurs_H.end());

    double tot_pixel_H = 0.;
    for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
    {
        tot_pixel_H += (*it).compteurD_;
    }

    for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
    {
        std::cout<<"("<<(*it).couleur_<<") : "<<(*it).compteurD_/tot_pixel_H<<std::endl;
    }
    std::cout<<std::endl;





    // ---------------------------------------------------------------------------


    // proportion voisinage dans E
    std::cout<<"mesure voisinage dans E"<<std::endl;
    std::vector<color_vois> couleur_vois_E = Tcontent_vois_test(res_composition);
    std::sort(couleur_vois_E.begin(), couleur_vois_E.end());


    // total par couleur1_
    double last_col = couleur_vois_E.at(0).couleur1_;
    std::vector<std::array<double, 2>> tot_vois_E{std::array{last_col, 0.}}; // couleur, total

    for(auto it = couleur_vois_E.begin(); it != couleur_vois_E.end(); it++)
    {
        if((*it).couleur1_ == last_col){
            tot_vois_E.at(tot_vois_E.size()-1)[1] += (*it).compteur_;
        }
        else{
            tot_vois_E.push_back(std::array{(*it).couleur1_, double((*it).compteur_)});
            last_col = (*it).couleur1_;
        }
    }

//    for(auto it = tot_vois_E.begin(); it != tot_vois_E.end(); it++)
//    {
//        std::cout<<(*it)[0]<<" "<<(*it)[1]<<std::endl;
//    }

    int last_col_id = 0;
    for(auto it = couleur_vois_E.begin(); it != couleur_vois_E.end(); it++)
    {
        if((*it).couleur1_ == tot_vois_E.at(last_col_id)[0]){
            std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteur_/tot_vois_E.at(last_col_id)[1]<<std::endl;
        }
        else{
            last_col_id += 1;
            std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteur_/tot_vois_E.at(last_col_id)[1]<<std::endl;
        }
    }


//    // total tous voisinage
//    double tot = 0.;
//    for(auto it = couleur_vois_E.begin(); it != couleur_vois_E.end(); it++)
//    {
//        tot += (*it).compteur_;
//    }
//
//    for(auto it = couleur_vois_E.begin(); it != couleur_vois_E.end(); it++)
//    {
////        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteur_/tot<<std::endl;
//        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteur_<<std::endl;
//
//
//    }
    std::cout<<std::endl;











    // capacité voisinage dans H
    std::cout<<"estimation voisinage dans H"<<std::endl;
    std::vector<color_vois> couleur_vois_H = Hcontent_vois_test(cm_, mu_1, mu_2, var_1, var_2, histo_N1_N2);
    std::sort(couleur_vois_H.begin(), couleur_vois_H.end());


    // total par couleur1_
    last_col = couleur_vois_H.at(0).couleur1_;
    std::vector<std::array<double, 2>> tot_vois_H{std::array{last_col, 0.}}; // couleur, total

    for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
    {
        if((*it).couleur1_ == last_col){
            tot_vois_H.at(tot_vois_H.size()-1)[1] += (*it).compteurD_;
        }
        else{
            tot_vois_H.push_back(std::array{(*it).couleur1_, double((*it).compteurD_)});
            last_col = (*it).couleur1_;
        }
    }

//    for(auto it = tot_vois_E.begin(); it != tot_vois_E.end(); it++)
//    {
//        std::cout<<(*it)[0]<<" "<<(*it)[1]<<std::endl;
//    }

    last_col_id = 0;
    for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
    {
        if((*it).couleur1_ == tot_vois_H.at(last_col_id)[0]){
            std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteurD_/tot_vois_H.at(last_col_id)[1]<<std::endl;
        }
        else{
            last_col_id += 1;
            std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteurD_/tot_vois_H.at(last_col_id)[1]<<std::endl;
        }
    }



//    // total tous voisinage
//    double tot_H = 0.;
//    for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
//    {
//        tot_H += (*it).compteurD_;
//    }
//
//    for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
//    {
//        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteurD_/tot_H<<std::endl;
////        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteurD_<<std::endl;
//
////        double dist = grain_dist((*it).couleur1_, (*it).couleur2_, H_color, H_seeds);
////        std::cout<<2.*dist*cm_size<<std::endl;
////        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteurD_/(2.*dist*cm_size)<<std::endl;
//    }






    return EXIT_SUCCESS;
}
