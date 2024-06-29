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
#include <ASTex/image_rgb.h>
#include "ASTex/Noises/Gabor.h"
#include "ASTex/fourier.h"

//#include "ASTex/CCVT_tests/numerical_Tcontent.h"
//#include "ASTex/CCVT_tests/numerical_Tvois.h"
#include "ASTex/CCVT_tests/mesure_statistiques.h"
//#include "ASTex/CCVT_tests/numerical_Treynold.h"
//#include "ASTex/CCVT_tests/large_numerical_test.h"

#include "ASTex/CCVT_tests/mesure_full.h"



using namespace ASTex;
//typedef ImageGrayd ImgType; // ImageRGBu8

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
    int img_size = 1024;//4096;//2048; // nombre de pixel dans l'image
    int cm_size = 256;//256;//512;
    double scale = double(resolution)/double(img_size);


// ---------------------------------------------------------------------------
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
    ImageGrayd image_1 = storing_noise_d(resolution, img_size, noise_1);
    IO::save(image_1, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_1.png");


    noise noise_2(K_, // anisotrope
                  a_,
                  0.06,//0.02, //F_0_,
                  3., //omega_0_,
                  number_of_impulses_per_kernel,
                  period,
                  random_offset_-4,
                  seed_+12);
    ImageGrayd image_2 = storing_noise_d(resolution, img_size, noise_2);
    IO::save(image_2, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_2.png");






// ---------------------------------------------------------------------------
// statistiques
    double mu_1 = moyenne(image_1); // mu
    double mu_2 = moyenne(image_2); // mu'

    double mu_1_carre = moyenne_carre(image_1);
    double mu_2_carre = moyenne_carre(image_2);

    double var_1 = mu_1_carre - mu_1*mu_1; // sigma
    double var_2 = mu_2_carre - mu_2*mu_2; // sigma'
//
//    double AC1 = covariance(image_1, 1, 0); // sigma_uv
//    double AC2 = covariance(image_2, 1, 0); // sigma_uv'

//    ImageGrayd Fourier_1 = covariance_Fourier(image_1);
//    ImageGrayd Fourier_2 = covariance_Fourier(image_2);
//
//    double var_1 = Fourier_1.pixelAbsolute(0, 0); // sigma
//    double var_2 = Fourier_2.pixelAbsolute(0, 0); // sigma'
//
//    double AC1 = Fourier_1.pixelAbsolute(0, 1); // sigma_uv
//    double AC2 = Fourier_2.pixelAbsolute(0, 1); // sigma_uv'


    // affichage
    std::cout<<"statistiques bruits"<<std::endl;

    std::cout<<"bruit N  : "<<std::endl;
    std::cout<<"moyenne = "<<mu_1<<std::endl;
    std::cout<<"variance = "<<var_1<<std::endl;
//    std::cout<<"variance sans biais = "<<var_1*(double(img_size)*double(img_size))/(double(img_size)*double(img_size)-1.)<<std::endl;
//    std::cout<<"autocovariance = "<<AC1<<std::endl;
//    std::cout<<"variance Fourier = "<<Fourier_1.pixelAbsolute(0, 0)<<std::endl;
//    std::cout<<"autocovariance Fourier = "<<Fourier_1.pixelAbsolute(1, 0)<<std::endl;

    std::cout<<std::endl;

    std::cout<<"bruit N'  : "<<std::endl;
    std::cout<<"moyenne = "<<mu_2<<std::endl;
    std::cout<<"variance = "<<var_2<<std::endl;
//    std::cout<<"variance sans biais = "<<var_2*(double(img_size)*double(img_size))/(double(img_size)*double(img_size)-1.)<<std::endl;
//    std::cout<<"autocovariance = "<<AC2<<std::endl;
//    std::cout<<"variance Fourier = "<<Fourier_2.pixelAbsolute(0, 0)<<std::endl;
//    std::cout<<"autocovariance Fourier = "<<Fourier_2.pixelAbsolute(1, 0)<<std::endl;

    std::cout<<std::endl;



// histogrammes
    // présences
    ImageGrayd histo_theo = histo_2D_theo(mu_1, mu_2, var_1, var_2, cm_size); // théorique
    IO::save(histo_theo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_theo.png");

//    double tot_pixel_t = 0;
//    histo_theo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y){
//        tot_pixel_t += P;
//    });
//
//    std::cout<<"histo théo : "<<tot_pixel_t<<std::endl;




    ImageGrayd histo_N1_N2 = histo_2D(image_1, image_2, cm_size); // réel
    IO::save(histo_N1_N2, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_reel.png");

//    double tot_pixel_r = 0;
//    histo_N1_N2.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y){
//        tot_pixel_r += P;
//    });
//
//    std::cout<<"histo réel : "<<tot_pixel_r<<std::endl;







// ---------------------------------------------------------------------------
    std::vector<Graine> H_seeds{Graine(0.81, 0.65, 0.),
                                Graine(0.15, 0.25, 0.1),
                                Graine(0.52, 0.91 ,0.),
                                Graine(0.65, 0.28, 0.),
                                Graine(0.12, 0.72, 0.1),
                                Graine(0.45, 0.48, 0.)};
    std::vector<ImageRGBu8::PixelType> H_color{itkRGBPixel(240,  120,  120),
                                                itkRGBPixel(120,  240,  120),
                                                itkRGBPixel(120,  120,  240),
                                                itkRGBPixel(240,  240,  120),
                                                itkRGBPixel(240,  120,  240),
                                                itkRGBPixel(120,  240,  240)};

//    std::vector<Graine> H_seeds{Graine(0.530276, 0.599851, -0.00779434),
//                                Graine(0.736087, 0.662678, 0.00726154),
//                                Graine(0.75289, 0.339547, 0.0102616),
//                                Graine(0.607775, 0.213809, 0.00904789),
//                                Graine(0.219568, 0.342634, 0.00887349),
//                                Graine(0.430012, 0.838025, 0.0173274),
//                                Graine(0.408312, 0.309818, -0.00889121),
//                                Graine(0.670424, 0.362264, 0.000172283),
//                                Graine(0.355614, 0.718189, 0.000232448),
//                                Graine(0.548694, 0.463016, -0.00989114),
//                                Graine(0.215937, 0.417019, 0.00265394),
//                                Graine(0.509193, 0.539335, -0.0108552),
//                                Graine(0.341272, 0.149167, 0.013753),
//                                Graine(0.493025, 0.475341, -0.0103001),
//                                Graine(0.599388, 0.353583, -0.00344984),
//                                Graine(0.423117, 0.378514, -0.0127338),
//                                Graine(0.363846, 0.541002, -0.012969),
//                                Graine(0.177101, 0.72001, 0.0242443),
//                                Graine(0.515779, 0.71665, -0.00264227),
//                                Graine(0.766039, 0.561647, 0.00575405),
//                                Graine(0.446871, 0.496666, -0.0124909),
//                                Graine(0.613095, 0.512305, -0.010515),
//                                Graine(0.513674, 0.377436, -0.0104248),
//                                Graine(0.690721, 0.766485, 0.0133756)};
//    std::vector<ImageRGBu8::PixelType> H_color{itkRGBPixel(29,  14,  5),
//                                               itkRGBPixel(50,  13,  2),
//                                               itkRGBPixel(54,  46,  30),
//                                               itkRGBPixel(66,  39,  23),
//                                               itkRGBPixel(80,  32,  14),
//                                               itkRGBPixel(85,  73,  62),
//                                                itkRGBPixel(92,  71,  47),
//                                                itkRGBPixel(107,  61,  43),
//                                                itkRGBPixel(111,  105,  96),
//                                                itkRGBPixel(120,  102,  78),
//                                                itkRGBPixel(127,  89,  47),
//                                                itkRGBPixel(140,  75,  37),
//                                                itkRGBPixel(140,  91,  71),
//                                                itkRGBPixel(141,  136,  125),
//                                                itkRGBPixel(151,  133,  101),
//                                                itkRGBPixel(170,  129,  76),
//                                                itkRGBPixel(172,  168,  158),
//                                                itkRGBPixel(183,  135,  107),
//                                                itkRGBPixel(185,  166,  132),
//                                                itkRGBPixel(198,  125,  73),
//                                                itkRGBPixel(200,  161,  104),
//                                                itkRGBPixel(213,  208,  199),
//                                                itkRGBPixel(226,  203,  167),
//                                                itkRGBPixel(233,  196,  135)};




    ImageRGBu8 cm_(cm_size, cm_size);

    cm_.for_all_pixels([&] (typename ImageRGBu8::PixelType& P, int x, int y)
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
    cm_.save("/home/grenier/Documents/ASTex_fork/results/equ_CCVT/cm.png");
//    IO::save(cm_, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/cm.png");




// ---------------------------------------------------------------------------
    ImageRGBu8 res_composition(img_size, img_size);


    // composition des bruits par la carte de couleurs
    res_composition.for_all_pixels([&] (typename ImageRGBu8::PixelType& P, int x, int y)
                                            {
                                                // valeur sur [0, 1]
                                                double n1 = image_1.pixelAbsolute(x,y);
                                                double n2 = image_2.pixelAbsolute(x,y);

                                                int id_n1 = int(std::round(n1*(cm_size-1)));
                                                int id_n2 = int(std::round(n2*(cm_size-1)));

                                                P = cm_.pixelAbsolute(id_n1, id_n2);
                                            });
    res_composition.save("/home/grenier/Documents/ASTex_fork/results/equ_CCVT/res_composition.png");



//     chargement d'une texture à partir d'une image
//    res_composition.load("/home/grenier/Documents/ASTex_fork/results/equ_CCVT/S_test_24.png");







// ---------------------------------------------------------------------------
    std::vector<color_info> couleurs_E;
    std::vector<color_vois> couleur_vois_E;
    std::vector<color_triple> couleur_tri_E;
    int nb_voisinage = 0;
    int nb_voisinage_diff = 0;

    mesure_E(res_composition,
             couleurs_E, couleur_vois_E, couleur_tri_E,
             false, nb_voisinage, nb_voisinage_diff);

    // tri
    std::sort(couleurs_E.begin(), couleurs_E.end());
    std::sort(couleur_vois_E.begin(), couleur_vois_E.end());
    std::sort(couleur_tri_E.begin(), couleur_tri_E.end());

    // compte du total
    double tot_pixel_E = 0.;
    for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
    {
        tot_pixel_E += (*it).compteur_;
    }


    // affichage
    std::cout<<"présence texture"<<std::endl;
    for(auto col=couleurs_E.begin(); col<couleurs_E.end(); col++){
        ImageRGBu8::PixelType couleur = (*col).couleur_;
        double proportion_E = (*col).compteur_;//tot_pixel_E;

        std::cout<<couleur<<" & "<< proportion_E<<std::endl;
    }
    std::cout<<std::endl;


//    std::cout<<"voisinage texture"<<std::endl;
//    for(auto vois=couleur_vois_E.begin(); vois<couleur_vois_E.end(); vois++){
//        ImageRGBu8::PixelType couleur1 = (*vois).couleur1_;
//        ImageRGBu8::PixelType couleur2 = (*vois).couleur2_;
//        double proportion_E = (*vois).compteur_;
//
//        std::cout<<couleur1<<", "<<couleur2<<" & "<< proportion_E<<std::endl;
//    }
//    std::cout<<std::endl;
//
//
//
//    std::cout<<"triplet texture"<<std::endl;
//    for(auto tri=couleur_tri_E.begin(); tri<couleur_tri_E.end(); tri++){
//        ImageRGBu8::PixelType couleur1 = (*tri).couleur1_;
//        ImageRGBu8::PixelType couleur2 = (*tri).couleur2_;
//        ImageRGBu8::PixelType couleur3 = (*tri).couleur3_;
//        double proportion_E = (*tri).compteur_;
//
//        std::cout<<couleur1<<", "<<couleur2<<", "<<couleur3<<" & "<< proportion_E<<std::endl;
//    }
//    std::cout<<std::endl;








//// ---------------------------------------------------------------------------
//    std::vector<color_info> couleurs_H;
//    std::vector<color_info> couleurs_H_r;
//
//    std::vector<color_vois> couleur_vois_H;
//    std::vector<color_vois> couleur_vois_H_r;
//
//    std::vector<color_triple> couleur_tripl_H;
//    std::vector<color_triple> couleur_tripl_H_r;
//
//    mesure_H(cm_, couleurs_H, couleurs_H_r,
//             couleur_vois_H, couleur_vois_H_r,
//             couleur_tripl_H, couleur_tripl_H_r,
//             mu_1, mu_2, var_1, var_2, histo_N1_N2, true);
//
//    // tri
//    std::sort(couleurs_H.begin(), couleurs_H.end());
//    std::sort(couleurs_H_r.begin(), couleurs_H_r.end());
//
//    std::sort(couleur_vois_H.begin(), couleur_vois_H.end());
//    std::sort(couleur_vois_H_r.begin(), couleur_vois_H_r.end());
//
//    std::sort(couleur_tripl_H.begin(), couleur_tripl_H.end());
//    std::sort(couleur_tripl_H_r.begin(), couleur_tripl_H_r.end());
//
//
//    // compte total
//    double tot_pixel_H = 0.;
//    for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
//    {
//        tot_pixel_H += (*it).compteurD_;
//    }
//
//    double tot_pixel_H_r = 0.;
//    for(auto it = couleurs_H_r.begin(); it != couleurs_H_r.end(); it++)
//    {
//        tot_pixel_H_r += (*it).compteurD_;
//    }
//
//
//
//    // affichage
//    std::cout<<"présence carte H histo théo"<<std::endl;
//    for(auto col=couleurs_H.begin(); col<couleurs_H.end(); col++){
//        ImageRGBu8::PixelType couleur = (*col).couleur_;
//        double proportion_E = (*col).compteurD_;//tot_pixel_E;
//
//        std::cout<<couleur<<" & "<< proportion_E<<std::endl;
//    }
//    std::cout<<std::endl;
//
//
//
//    std::cout<<"présence carte H histo réel"<<std::endl;
//    for(auto col=couleurs_H_r.begin(); col<couleurs_H_r.end(); col++){
//        ImageRGBu8::PixelType couleur = (*col).couleur_;
//        double proportion_E = (*col).compteur_;//tot_pixel_E;
//
//        std::cout<<couleur<<" & "<< proportion_E<<std::endl;
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
////    std::cout<<"voisinage carte H histo théo"<<std::endl;
////    for(int i=0; i<couleur_vois_H.size(); i++){
////        ImageRGBu8::PixelType couleur1 = couleur_vois_H.at(i).couleur1_;
////        ImageRGBu8::PixelType couleur2 = couleur_vois_H.at(i).couleur2_;
////        double proportion_H;
////
////
////        if(couleur1 == couleur2){
////            proportion_H = couleur_vois_H.at(i).compteurD_/double(nb_voisinage);
////        }
////        else{
////            proportion_H = couleur_vois_H.at(i).compteurD_/(2.*double(nb_voisinage));
////        }
////
////        std::cout<<couleur1<<", "<<couleur2<<" & "<< proportion_H<<std::endl;
////    }
////    std::cout<<std::endl;
////
////
////
//    std::cout<<"voisinage carte H histo théo"<<std::endl;
//    for(auto vois=couleur_vois_H.begin(); vois<couleur_vois_H.end(); vois++){
//        ImageRGBu8::PixelType couleur1 = (*vois).couleur1_;
//        ImageRGBu8::PixelType couleur2 = (*vois).couleur2_;
//        double proportion_E = (*vois).compteurD_;
//
//        std::cout<<couleur1<<", "<<couleur2<<" & "<< proportion_E<<std::endl;
//    }
//    std::cout<<std::endl;
//
//
//    std::cout<<"voisinage carte H histo réel"<<std::endl;
//    for(auto vois=couleur_vois_H_r.begin(); vois<couleur_vois_H_r.end(); vois++){
//        ImageRGBu8::PixelType couleur1 = (*vois).couleur1_;
//        ImageRGBu8::PixelType couleur2 = (*vois).couleur2_;
//        double proportion_E = (*vois).compteur_;
//
//        std::cout<<couleur1<<", "<<couleur2<<" & "<< proportion_E<<std::endl;
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
//    std::cout<<"triplet carte H histo théo"<<std::endl;
//    for(auto tri=couleur_tripl_H.begin(); tri<couleur_tripl_H.end(); tri++){
//        ImageRGBu8::PixelType couleur1 = (*tri).couleur1_;
//        ImageRGBu8::PixelType couleur2 = (*tri).couleur2_;
//        ImageRGBu8::PixelType couleur3 = (*tri).couleur3_;
//        double proportion_E = (*tri).compteurD_;
//
//        std::cout<<couleur1<<", "<<couleur2<<", "<<couleur3<<" & "<< proportion_E<<std::endl;
//    }
//    std::cout<<std::endl;
//
//
//    std::cout<<"triplet carte H histo réel"<<std::endl;
//    for(auto tri=couleur_tripl_H_r.begin(); tri<couleur_tripl_H_r.end(); tri++){
//        ImageRGBu8::PixelType couleur1 = (*tri).couleur1_;
//        ImageRGBu8::PixelType couleur2 = (*tri).couleur2_;
//        ImageRGBu8::PixelType couleur3 = (*tri).couleur3_;
//        double proportion_E = (*tri).compteur_;
//
//        std::cout<<couleur1<<", "<<couleur2<<", "<<couleur3<<" & "<< proportion_E<<std::endl;
//    }
//    std::cout<<std::endl;






    return EXIT_SUCCESS;
}

