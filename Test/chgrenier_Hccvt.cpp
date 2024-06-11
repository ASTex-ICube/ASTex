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

#include "ASTex/CCVT_tests/mesure_statistiques.h"
#include "ASTex/CCVT_tests/mesure_full.h"



using namespace ASTex;



int main()
{
// ---------------------------------------------------------------------------
    int img_size = 800;//4096;//2048; // nombre de pixel dans l'image
    int cm_size = 256;//256;//512;


// ---------------------------------------------------------------------------
    ImageGrayd noise_1(img_size, img_size);
    noise_1.load("/home/grenier/Documents/ASTex_fork/results/CCVT_CGAL/application/noise_1.png");

    ImageGrayd noise_2(img_size, img_size);
    noise_2.load("/home/grenier/Documents/ASTex_fork/results/CCVT_CGAL/application/noise_2.png");



// ---------------------------------------------------------------------------
// statistiques
    double mu_1 = moyenne_255(noise_1); // mu
    double mu_2 = moyenne_255(noise_2); // mu'

    double mu_1_carre = moyenne_carre_255(noise_1);
    double mu_2_carre = moyenne_carre_255(noise_2);

    double var_1 = mu_1_carre - mu_1*mu_1; // sigma
    double var_2 = mu_2_carre - mu_2*mu_2; // sigma'



    // affichage
    std::cout<<"statistiques bruits"<<std::endl;

    std::cout<<"bruit N  : "<<std::endl;
    std::cout<<"moyenne = "<<mu_1<<std::endl;
    std::cout<<"variance = "<<var_1<<std::endl;
    std::cout<<std::endl;

    std::cout<<"bruit N'  : "<<std::endl;
    std::cout<<"moyenne = "<<mu_2<<std::endl;
    std::cout<<"variance = "<<var_2<<std::endl;
    std::cout<<std::endl;


    // histogrammes
    ImageGrayd histo_theo = histo_2D_theo(mu_1, mu_2, var_1, var_2, cm_size); // théorique
    IO::save(histo_theo, "/home/grenier/Documents/ASTex_fork/results/CCVT_CGAL/application/histo_theo.png");

    ImageGrayd histo_theo_app = histo_2D_theo(0.5, 0.5, 0.02, 0.02, cm_size); // théorique
    IO::save(histo_theo_app, "/home/grenier/Documents/ASTex_fork/results/CCVT_CGAL/application/histo_theo_app.png");



    ImageGrayd histo_N1_N2 = histo_2D(noise_1, noise_2, cm_size); // réel
    IO::save(histo_N1_N2, "/home/grenier/Documents/ASTex_fork/results/CCVT_CGAL/application/histo_reel.png");



//// ---------------------------------------------------------------------------
//    ImageRGBu8 cm_(cm_size, cm_size);
//    cm_.load("file");




// ---------------------------------------------------------------------------
    ImageRGBu8 res_composition(img_size, img_size);
    res_composition.load("/home/grenier/Documents/ASTex_fork/results/CCVT_CGAL/application/composition.png");







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
        double proportion_E = (*col).compteur_/tot_pixel_E;

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

