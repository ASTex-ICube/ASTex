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

using namespace ASTex;

// Gabor.h est le code exemple du papier d'origine dont il est un peut touffu (aussi à cause des structure accélératrices)
// les principales partie de l'algo sont :
// fonction gabor : définition d'un noyaux de gabor (gaussienne * cosinus)
// fonction cell : boucle pour placer et sommer les noyaux







// ---------------------------------------------------------------------------

void capacity_large_test(int loop)
{
    int resolution = 128;
    int img_size = 512;//8192; // nombre de pixel dans l'image

    float F_0_ = 0.2;//0.04; // fréquence
    float omega_0_ = 0.;//M_PI/4.; // orientation (seulement dans le cas anisotrope, cf code gabor ligne 160)

    float number_of_impulses_per_kernel = 64.0;
    unsigned period = 128; // non utilisé

    float K_ = 1.0; //laisser à 1
    float a_ = 0.02; // taille des noyaux (a*a = 1/variance)


    // ---------------------------------------------------------------------------
    int cm_size = 256;

    ImageGrayd cm_(cm_size, cm_size);
    IO::load(cm_, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/cm.png");

    std::vector<color_info> couleur_cm = Tcontent(cm_);

    std::vector<comput_mean> presence_E;
    std::vector<comput_mean> capacity_H;
    std::vector<comput_mean> erreur_EH;

    for(auto it = couleur_cm.begin(); it != couleur_cm.end(); it++)
    {
        presence_E.push_back(comput_mean{(*it).couleur_});
        capacity_H.push_back(comput_mean{(*it).couleur_});
        erreur_EH.push_back(comput_mean{(*it).couleur_});
    }



    // ---------------------------------------------------------------------------
    // faire une boucle à partir d'ici, avec un tireage de seed à chaque itération
    for(int l=0; l<loop; l++)
    {
//        if(l%10==0){std::cout<<l<<std::endl;}
        unsigned random_offset_ = std::rand();//954248632;
        unsigned seed_ = 1 + std::rand() / ((RAND_MAX + 1u) / (l+1));//4;


        noise noise_1(K_, a_, F_0_, omega_0_, number_of_impulses_per_kernel, period, random_offset_, seed_);
        noise noise_2(K_, a_, F_0_, omega_0_, number_of_impulses_per_kernel, period, random_offset_, seed_+24);

        ImageGrayd image_1 = storing_noise_d(resolution, img_size, noise_1);
        ImageGrayd image_2 = storing_noise_d(resolution, img_size, noise_2);


        // ---------------------------------------------------------------------------
        ImageGrayd res_composition(img_size, img_size);
        res_composition.parallel_for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                                                {
                                                    double n1 = image_1.pixelAbsolute(x,y);
                                                    double n2 = image_2.pixelAbsolute(x,y);

                                                    int id_n1 = std::clamp(int(n1*255.), 0, cm_size);
                                                    int id_n2 = std::clamp(int(n2*255.), 0, cm_size);

                                                    P = cm_.pixelAbsolute(id_n1, id_n2);
                                                });

        double mu_1 = moyenne(image_1);
        double mu_2 = moyenne(image_2);

        double mu_1_carre = moyenne_carre(image_1);
        double mu_2_carre = moyenne_carre(image_2);

        double var_1 = mu_1_carre - mu_1*mu_1;
        double var_2 = mu_2_carre - mu_2*mu_2;


        // ---------------------------------------------------------------------------
        std::vector<color_info> couleurs_E = Tcontent(res_composition);
        double tot_pixel_E = 0.;
        for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
        {
            tot_pixel_E += (*it).compteur_;
        }

        std::vector<color_info> couleurs_H = cell_capacity(cm_, mu_1, mu_2, var_1, var_2);
        double tot_pixel_H = 0.;
        for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
        {
            tot_pixel_H += (*it).compteurD_;
        }


        // ---------------------------------------------------------------------------
        for(int col_id=0; col_id<presence_E.size(); col_id++)
        {
            int id_E;
            int id_H;

            is_in(couleurs_E, presence_E.at(col_id).couleur_, id_E);
            is_in(couleurs_H, presence_E.at(col_id).couleur_, id_H);

            // ARP : les compteur sont des int, toute les proportion vont tomber à 0 -> faire une autre structure
            presence_E.at(col_id).incr(couleurs_E.at(id_E).compteur_/tot_pixel_E);
            capacity_H.at(col_id).incr(couleurs_H.at(id_H).compteurD_/tot_pixel_H);
            erreur_EH.at(col_id).incr(couleurs_E.at(id_E).compteur_/tot_pixel_E - couleurs_H.at(id_H).compteurD_/tot_pixel_H);
        }

    }


    // ---------------------------------------------------------------------------
    for(int col_id=0; col_id<presence_E.size(); col_id++)
    {
        float tot = presence_E.at(col_id).values_.size();
        std::cout<<presence_E.at(col_id).couleur_<<" : "<<presence_E.at(col_id).sum_/tot<<", "<<capacity_H.at(col_id).sum_/tot<<", "<<erreur_EH.at(col_id).sum_/tot<<std::endl;
    }

}









int main()
{
// ---------------------------------------------------------------------------
    int resolution = 512;
    int img_size = 512; // nombre de pixel dans l'image
    double scale = double(resolution)/double(img_size);

    float F_0_ = 0.06;//0.04; // fréquence
    float omega_0_ = 0.;//M_PI/4.; // orientation (seulement dans le cas anisotrope, cf code gabor ligne 160)

    float number_of_impulses_per_kernel = 64.;// 64.0;
    unsigned period = 128; // non utilisé

    float K_ = 1.0; //laisser à 1
    float a_ = 0.02; // taille des noyaux (a*a = 1/variance)0.02

    unsigned random_offset_ = 954248632;
    unsigned seed_ = 4;


    // ---------------------------------------------------------------------------
    noise noise_1(K_,
                 a_,
                 0.016, //F_0_,
                 0.57, //omega_0_,
                 number_of_impulses_per_kernel,
                 period,
                 random_offset_,
                 seed_);
// la valeur du bruit noise_ en x,y peut être récupérée par noise_(x,y)
    ImageGrayd image_1 = storing_noise_d(resolution, img_size, noise_1); // pour écrire le bruit dans une image
    IO::save(image_1, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_1.png");


    noise noise_2(K_,
                 a_,
                 0.016, //F_0_,
                 -0.57, //omega_0_,
                 number_of_impulses_per_kernel,
                 period,
                 random_offset_+4,
                 seed_+4);
// la valeur du bruit noise_ en x,y peut être récupérée par noise_(x,y)
    ImageGrayd image_2 = storing_noise_d(resolution, img_size, noise_2); // pour écrire le bruit dans une image
    IO::save(image_2, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_2.png");











    // ---------------------------------------------------------------------------
    int cm_size = 256;

    ImageGrayd cm_(cm_size, cm_size);
    IO::load(cm_, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/cm.png");


    // ---------------------------------------------------------------------------
    ImageGrayd res_composition(img_size, img_size);
    double scale_1 = 3.6 * std::sqrt(noise_1.variance());
    double scale_2 = 3.6 * std::sqrt(noise_2.variance());

    res_composition.parallel_for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                                   {
//                                        float X = x * scale;
//                                        float Y = y * scale;
//
//                                        // valeur sur [??, ??]
//                                       double n1 = 0.5 + (0.5 * (noise_1(X, Y) / scale_1)); // sur  0, 1
//                                       double n2 = 0.5 + (0.5 * (noise_2(X, Y) / scale_2));
//
//                                       int id_n1 = std::clamp(int(std::round(n1*(cm_size-1))), 0, (cm_size-1));
//                                       int id_n2 = std::clamp(int(std::round(n2*(cm_size-1))), 0, (cm_size-1));


                                       // valeur sur [0, 1]
                                       double n1 = image_1.pixelAbsolute(x,y);
                                       double n2 = image_2.pixelAbsolute(x,y);

                                       int id_n1 = int(std::round(n1*(cm_size-1)));
                                       int id_n2 = int(std::round(n2*(cm_size-1)));

                                       P = cm_.pixelAbsolute(id_n1, id_n2);
                                   });
    IO::save(res_composition, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/res_composition.png");










    // ---------------------------------------------------------------------------


// histogrammes réels

//    ImageGrayd histo_N1_N2 = histo_2D(image_1,image_2);
//    std::cout<<std::endl;






// proportion présence dans E
    std::cout<<"mesure présence dans E"<<std::endl;
    std::vector<color_info> couleurs_E = Tcontent(res_composition);

    double tot_E = 0.;
    for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
    {
        tot_E += (*it).compteur_;
    }

    for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
    {
        std::cout<<"("<<(*it).couleur_<<") : "<<(*it).compteur_/tot_E<<std::endl;
//        std::cout<<"("<<(*it).couleur_<<") : "<<(*it).compteur_<<std::endl;
    }
    std::cout<<std::endl;





// proportion voisinage dans E
    std::cout<<"mesure voisinage dans E 1 0"<<std::endl;
    std::vector<color_vois> couleur_vois = Tcontent_vois(res_composition, 1, 0);

    double tot = 0.;
    for(auto it = couleur_vois.begin(); it != couleur_vois.end(); it++)
    {
        tot += (*it).compteur_;
    }

    for(auto it = couleur_vois.begin(); it != couleur_vois.end(); it++)
    {
        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteur_/tot<<std::endl;
//        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteur_<<std::endl;
    }
    std::cout<<std::endl;





    std::cout<<"mesure voisinage dans E 0 1"<<std::endl;
    couleur_vois = Tcontent_vois(res_composition, 0, 1);

    tot = 0.;
    for(auto it = couleur_vois.begin(); it != couleur_vois.end(); it++)
    {
        tot += (*it).compteur_;
    }

    for(auto it = couleur_vois.begin(); it != couleur_vois.end(); it++)
    {
        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteur_/tot<<std::endl;
//        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteur_<<std::endl;
    }
    std::cout<<std::endl;







// statistiques
    double mu_1 = moyenne(image_1); // mu
    double mu_2 = moyenne(image_2); // mu'

//    double mu_1_p = moyenne(noise_1, img_size, resolution); // mu
//    double mu_2_p = moyenne(noise_2, img_size, resolution); // mu'

    double mu_1_carre = moyenne_carre(image_1);
    double mu_2_carre = moyenne_carre(image_2);

    double var_1 = mu_1_carre - mu_1*mu_1; // sigma
    double var_2 = mu_2_carre - mu_2*mu_2; // sigma'

    double AC1 = covariance(image_1, 1, 0); // sigma_uv
    double AC2 = covariance(image_2, 1, 0); // sigma_uv'

    std::cout<<"bruit N  : moyenne = "<<mu_1<<", variance = "<<var_1<<", autocovariance = "<<AC1<<std::endl;
    std::cout<<"bruit N' : moyenne = "<<mu_2<<", variance = "<<var_2<<", autocovariance = "<<AC2<<std::endl;
    std::cout<<std::endl;
//    std::cout<<"bruit N  : moyenne = "<<mu_1_p<<", variance = "<<noise_1.variance()<<", autocovariance = "<<AC1<<std::endl;
//    std::cout<<"bruit N' : moyenne = "<<mu_2_p<<", variance = "<<noise_2.variance()<<", autocovariance = "<<AC2<<std::endl;
//    std::cout<<std::endl;



// histogrammes
    // présences
    std::cout<<"histogramme présence"<<std::endl;
    ImageGrayd histo_N1_N2 = histo_2D(image_1,image_2); // réel
    ImageGrayd histo_theo = histo_2D_theo(mu_1, mu_2, var_1, var_2); // théorique
    ImageGrayd histo_dist = histo_2D_dist(histo_N1_N2, mu_1, mu_2, var_1, var_2); // disance
    std::cout<<std::endl;

//    // voisinages
//    std::cout<<"histogramme voisinage"<<std::endl;
//    ImageGrayd histo_AC_N1 = histo_2D_AC(image_1, 1, 0, 1); // réel
//    ImageGrayd histo_AC_N2 = histo_2D_AC(image_2, 1, 0, 2);
//
//    ImageGrayd histo_theo_AC_N1 = histo_2D_theo_vois(mu_1, var_1, AC1, 1); // théorique
//    ImageGrayd histo_theo_AC_N2 = histo_2D_theo_vois(mu_2, var_2, AC2, 2);
//
//    ImageGrayd histo_dist_AC_N1 = histo_2D_dist_vois(histo_AC_N1, mu_1, var_1, AC1, 1); // disance
//    ImageGrayd histo_dist_AC_N2 = histo_2D_dist_vois(histo_AC_N2, mu_2, var_2, AC2, 2);
//    std::cout<<std::endl;






//// capacité présence dans H
//    std::cout<<"estimation présence dans H"<<std::endl;
////    std::vector<color_info> couleurs_H = cell_capacity_real_histo(cm_, histo_N1_N2);
//    std::vector<color_info> couleurs_H = cell_capacity(cm_, mu_1, mu_2, var_1, var_2);
////    std::vector<color_info> couleurs_H = cell_dist_real_histo(cm_, histo_dist, mu_1, mu_2, var_1, var_2);
//
//    double tot_pixel_H = 0.;
//    for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
//    {
//        tot_pixel_H += (*it).compteurD_;
//    }
//
//    for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
//    {
//        std::cout<<"("<<(*it).couleur_<<") : "<<(*it).compteurD_/tot_pixel_H<<std::endl;
////        std::cout<<"("<<(*it).couleur_<<") : "<<(*it).compteurD_<<std::endl;
//    }
//    std::cout<<std::endl;
//
//
//
//// capacité conjointe dans H
//    std::cout<<"estimation voisinage dans H"<<std::endl;
//    std::vector<color_vois> couleur_vois_H = capacityH_vois_histo_reel(cm_, histo_AC_N1, histo_AC_N2);
////    std::vector<color_vois> couleur_vois_H = capacityH_vois(cm_, mu_1, mu_2, var_1, var_2, AC1, AC2);
////    std::vector<color_vois> couleur_vois_H = capacityH_vois_histo_dist(cm_, histo_dist_AC_N1, histo_dist_AC_N2, mu_1, mu_2, var_1, var_2, AC1, AC2);
//
//    double tot_H = 0.;
//    for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
//    {
//        tot_H += (*it).compteurD_;
//    }
//
//    for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
//    {
//        std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteurD_/tot_H<<std::endl;
//    }













//// test proportion présence sur plusieurs itération de bruit
//        auto start = std::chrono::high_resolution_clock::now();
//        capacity_large_test(10);
//        auto stop = std::chrono::high_resolution_clock::now();
//        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
//        std::cout << duration.count() << std::endl;





    // tracer pour le manuscrit
//    ImageGrayd image_3(512, 512);
//    ImageGrayd image_4(512, 512);
//    IO::load(image_3, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/manuscrit/cos_bilobe.png");
//    IO::load(image_4, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/manuscrit/sin_bilobe.png");
//    ImageGrayd image_3 = dF_dx(image_1);
//    ImageGrayd image_4 = dF_dy(image_1);
//    histo_2D(image_3, image_4);







    return EXIT_SUCCESS;
}
