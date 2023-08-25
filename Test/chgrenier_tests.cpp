//
// Created by grenier on 28/03/23.
//
#include <ASTex/image_gray.h>
#include "ASTex/Noises/Gabor.h"

using namespace ASTex;

// Gabor.h est le code exemple du papier d'origine dont il est un peut touffu (aussi à cause des structure accélératrices)
// les principales partie de l'algo sont :
// fonction gabor : définition d'un noyaux de gabor (gaussienne * cosinus)
// fonction cell : boucle pour placer et sommer les noyaux

int main()
{
// ---------------------------------------------------------------------------
    int resolution = 128;
    int img_size = 512; // nombre de pixel dans l'image

    float F_0_ = 0.1;//0.04; // fréquence
    float omega_0_ = 0.; // orientation (seulement dans le cas anisotrope, cf code gabor ligne 160)

    float number_of_impulses_per_kernel = 64.0;
    unsigned period = 128; // non utilisé

    float K_ = 1.0; //laisser à 1
    float a_ = 0.02; // taille des noyaux (a*a = 1/variance)

    unsigned random_offset_ = 954248632;
    unsigned seed_ = 1;


    // ---------------------------------------------------------------------------
    noise noise_(K_,
                 a_,
                 F_0_,
                 omega_0_,
                 number_of_impulses_per_kernel,
                 period,
                 random_offset_,
                 seed_);
// la valeur du bruit noise_ en x,y peut être récupérée par noise_(x,y)
    ImageGrayu8 image_ = storing_noise(resolution, img_size, noise_); // pour écrire le bruit dans une image
    image_.save("/home/grenier/Documents/ASTex_fork/results/T_analysis/gabor_test.png");

    return EXIT_SUCCESS;
}


