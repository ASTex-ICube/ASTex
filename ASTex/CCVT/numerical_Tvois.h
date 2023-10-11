//
// Created by grenier on 11/10/23.
//

#ifndef ASTEX_NUMERICAL_TVOIS_H
#define ASTEX_NUMERICAL_TVOIS_H




// ---------------------------------------------------------------------------
struct color_vois{
    ImageGrayd::PixelType couleur1_;
    ImageGrayd::PixelType couleur2_;
    int compteur_;
    double compteurD_;

    color_vois(ImageGrayd::PixelType couleur1, ImageGrayd::PixelType couleur2){
        couleur1_ = couleur1;
        couleur2_ = couleur2;
        compteur_ = 1;
    }

    color_vois(ImageGrayd::PixelType couleur1, ImageGrayd::PixelType couleur2, int i){
        couleur1_ = couleur1;
        couleur2_ = couleur2;
        compteur_ = i;
    }

    color_vois(ImageGrayd::PixelType couleur1, ImageGrayd::PixelType couleur2, double i){
        couleur1_ = couleur1;
        couleur2_ = couleur2;
        compteurD_ = i;
    }

    void incr()
    {
        compteur_ += 1;
    }

    void incr(int i)
    {
        compteur_ += i;
    }

    void incr(double i)
    {
        compteurD_ += i;
    }

    bool operator==(const color_vois& couleur2)
    {
        return (couleur1_ == couleur2.couleur1_) and (couleur2_ == couleur2.couleur2_);
    }
};



bool is_in(std::vector<color_vois>& vecteur, ImageGrayd::PixelType& couleur1, ImageGrayd::PixelType& couleur2, int& id)
{
    id = -1;
//    for(ImageRGB8::PixelType col : vecteur)
    for(int i=0; i<vecteur.size(); i++)
    {
        if(vecteur.at(i) == color_vois{couleur1, couleur2})
        {
            id = i;
            return true;
        }
    }
    return false;
}



std::vector<color_vois> Tcontent_vois(ImageGrayd image_, int tauX, int tauY)
{
    std::vector<color_vois> couleurs = {};

    int id;
    bool presence;

    // voisinage gauche-droite
    for(int x=0; x<image_.width()-tauX; x++)
    {
        for(int y=0; y<image_.height()-tauY; y++)
        {
//            int id;
//            if(couleurs.empty())
//            {
//                couleurs.push_back(color_vois{image_.pixelAbsolute(x,y), image_.pixelAbsolute(x+1, y)});
//            }
//            else
//            {
//                bool presence = is_in(couleurs, image_.pixelAbsolute(x,y), image_.pixelAbsolute(x+1, y), id);
//                if(not presence)
//                {
//                    couleurs.push_back(color_vois{image_.pixelAbsolute(x,y), image_.pixelAbsolute(x+1, y)});
//                }
//                else
//                {
//                    couleurs.at(id).incr();
//                }
//            }
//            if(couleurs.empty())
//            {
//                couleurs.push_back(color_vois{image_.pixelAbsolute(x,y), image_.pixelAbsolute(x+1, y)});
//            }
//            else
//            {
            presence = is_in(couleurs, image_.pixelAbsolute(x,y), image_.pixelAbsolute(x+tauX, y+tauY), id);
            if(not presence)
            {
                couleurs.push_back(color_vois{image_.pixelAbsolute(x,y), image_.pixelAbsolute(x+tauX, y+tauY)});
            }
            else
            {
                couleurs.at(id).incr();
            }

            // symétrie
            presence = is_in(couleurs, image_.pixelAbsolute(x+tauX, y+tauY), image_.pixelAbsolute(x,y), id);
            if(not presence)
            {
                couleurs.push_back(color_vois{image_.pixelAbsolute(x+tauX, y+tauY), image_.pixelAbsolute(x,y)});
            }
            else
            {
                couleurs.at(id).incr();
            }

//            }

        }
    }

    return couleurs;
}






//double gaussAC_4D(double moy1, double moy2, double var1, double var2, double AC1, double AC2, double x1, double x2, double y1, double y2)
//{
//    var1 = std::sqrt(var1);
//    var2 = std::sqrt(var2);
//
//    /* inverse de la matrice de covariance pour obtenir la matrice de précision
//     * inverse d'une matrice de la forme
//     * a 0 c 0
//     * 0 b 0 d
//     * c 0 a 0
//     * 0 d 0 b
//     *
//     * résultat de la forme
//     * 1/det * A 0 C 0
//     *         0 B 0 D
//     *         C 0 A 0
//     *         0 D 0 B
//     */
//    double detMC = var1*var1*var2*var2 + AC1*AC1*AC2*AC2 - (AC1*AC1*var2*var2 + AC2*AC2*var1*var1);
//    double A = var1*(var2*var2 - AC2*AC2);
//    double B = var2*(var1*var1 - AC1*AC1);
//    double C = AC1*(AC2*AC2 - var2*var2);
//    double D = AC2*(AC1*AC1 - var1*var1);
//
//    /* produit x^T M x
//     *  x1 x2 y1 y2 * A 0 C 0 * x1
//     *                0 B 0 D   x2
//     *                C 0 A 0   y1
//     *                0 D 0 B   y2
//     *
//     * résultat :
//     * x1(x1*A + y1*C) + x2(x2*B + y2*D) + y1(x1*C + y1*A) + y2(x2*D + y2*B)
//     */
//    double X1 = x1-moy1;
//    double X2 = x2-moy2;
//    double Y1 = y1-moy1;
//    double Y2 = y2-moy2;
//
//    double produit = X1*(X1*A + Y1*C) + X2*(X2*B + Y2*D) + Y1*(X1*C + Y1*A) + Y2*(X2*D + Y2*B);
//
//    /* gaussienne 4D
//     * G = exp(-0.5*(x-moy)^T C^-1 (x-moy))
//     *   = exp(-0.5/det produit)
//     */
//    double normalisation = (2.*M_PI)*(2.*M_PI)*std::sqrt(detMC);
//    double G = std::exp(-(0.5/detMC) * produit) / normalisation;
//
//
//    return G;
//}


double gaussAC_2D(double moy, double var, double AC, double x, double y)
{
    double ecart_type = std::sqrt(var);
    /* matrice covariance :
     * var AC
     * AC  var
     */
    double det = (var*var)-(AC*AC); // déterminant de la matrice de covariance

    double X = x-moy;
    double Y = y-moy;

    /* matrice précision (inverse matrice covariance) :
     * var/det -AC/det
     * -AC/det  var/det
     */
    double var_inv = var/det;
    double AC_inv = -AC/det;

    double produit = X*(X*var_inv + Y*AC_inv) + Y*(Y*var_inv + X*AC_inv);
    double norm = 1.;// 2.*M_PI*std::sqrt(det);

    double G = (1./norm)*std::exp(-0.5*produit);

    return G;
}







double quantified_gaussAC_2D(int x1, int x2, int y1, int y2, int cm_size, double moy1, double moy2, double var1, double var2, double AC1, double AC2)
{
    // x1 ~ Nu, y1 ~ N'u, x2 ~ Nv, y2 ~ N'v
    double X1 = x1/double(cm_size);
    double Y1 = y1/double(cm_size);

    double X2 = x2/double(cm_size);
    double Y2 = y2/double(cm_size);

    double gaussian1 = gaussAC_2D(moy1, var1, AC1, X1, X2);
    double gaussian2 = gaussAC_2D(moy2, var2, AC2, Y1, Y2);

    double Gaussian = gaussian1*gaussian2;//*100;//img_size*img_size;

    return Gaussian;
}





//int quantified_gaussAC(int x1, int x2, int y1, int y2, int cm_size, double moy1, double moy2, double var1, double var2, double AC1, double AC2)
//{
//
//    double X1 = x1/double(cm_size);
//    double X2 = x2/double(cm_size);
//    double Y1 = y1/double(cm_size);
//    double Y2 = y2/double(cm_size);
//
//    double gaussian = gaussAC_4D(moy1, moy2, var1, var2, AC1, AC2, X1, X2, Y1, Y2);
//    int Gaussian = gaussian*10;//img_size*img_size;
//
//    return Gaussian;
//}



std::vector<color_vois> capacityH_vois(ImageGrayd cm, double moy1, double moy2, double var1, double var2, double AC1, double AC2)
{
    std::vector<color_vois> couleurs = {};

    // listage des couleur et compte de leur apparition
    // x1 ~ Nu, y1 ~ N'u, x2 ~ Nv, y2 ~ N'v
    cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P1, int x1, int y1) // parallel_for_all_pixels, for_all_pixels
                      {
                          cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P2, int x2, int y2)
                                            {
                                                int id;
//                                                int qtt = quantified_gaussAC(x1, x2, y1, y2, cm.height(), moy1, moy2, var1, var2, AC1, AC2); // gaussienne 4D quantifiée
                                                double qtt = quantified_gaussAC_2D(x1, x2, y1, y2, cm.height(), moy1, moy2, var1, var2, AC1, AC2);

                                                if(couleurs.empty())
                                                {
                                                    couleurs.push_back(color_vois{P1, P2, qtt});
                                                }
                                                else
                                                {
                                                    bool presence = is_in(couleurs, P1, P2, id);
                                                    if(not presence)
                                                    {
                                                        couleurs.push_back(color_vois{P1, P2, qtt});
                                                    }
                                                    else
                                                    {
                                                        couleurs.at(id).incr(qtt);
                                                    }
                                                }

                                            });
//                          if(x1==0){std::cout<<y1<<" ";}
                      });
    return couleurs;
}








//void autocovariance_iFT(ImageGrayd inputNoise)
//{
//    ImageSpectrald module;
//    ImageSpectrald phase;
//    ImageGrayd autocovariance;
//
//    Fourier::fftForwardModulusAndPhase(inputNoise, module, phase, true); // transformée de fourier du bruit
//
//    IO::save_spectrum(module, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_module.png"); // module
//    IO::save_phase(phase, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_phase.png"); // phase
//
//    module.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y) // carré du module divisé par la surface d'intégration
//                          {
//                              P = (P*P);///(inputNoise.width()*inputNoise.height());
//                          });
//    IO::save_spectrum(module, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_PSD.png"); // PSD
//
//    phase.zero();
//    Fourier::fftInverseModulusAndPhase(module, phase, autocovariance, false); // fourier inverse sur la PSD = AC
//    IO::save(autocovariance, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_AC_Fourier.png");
//
//    std::cout<<autocovariance.pixelAbsolute(1,0)<<std::endl;
//
//
////    Fourier::fftForwardModulusAndPhase(autocovariance, module, phase, false); // transformée de fourier de l'AC
////    IO::save_spectrum(module, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/gabor_PSD_test.png"); // PSD
//}






ImageGrayd histo_2D_AC(ImageGrayd noise1, int tauX, int tauY, int name)
{
    ImageGrayd histo(256, 256, true);
    int max = 0;

    for(int x=0; x<noise1.width()-tauX; x++)
    {
        for(int y=0; y<noise1.height()-tauY; y++)
        {
            double nx = noise1.pixelAbsolute(x, y)*255;
            double ny = noise1.pixelAbsolute(x+tauX, y+tauY)*255;

            if(int(std::round(nx))==127 and int(std::round(ny))==127){max += 1;}

            histo.pixelAbsolute(std::round(nx), std::round(ny)) += 1; // approx à cause de la quantification
        }
    }
//    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_AC_0.png");


    double mX = 0.;
    double mY = 0.;
    double tot = 0.;
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             mX += P*(x/255.); // x/255. : valeur sur 0, 1
//                             mX += P*x;

                             mY += P*(y/255.);
//                             mY += P*y;

                             tot += P;
                         });
    mX = mX/tot;
    mY = mY/tot;
//    std::cout<<"moyenne histo x : "<<mX<<", moyenne histo y : "<<mY<<std::endl;




    double AC = 0.;
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             AC += P*((x/255.)-mX)*((y/255.)-mY); // x/255. : valeur sur 0, 1
//                             AC += P*(x-mX)*(y-mY); // x : valeur sur 0, 255
                         });
    AC = AC/tot;
//    std::cout<<"auto-covariance histo : "<<AC<<std::endl;




    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             P *= 1./double(max);
                         });
    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_AC_" + std::to_string(name) + ".png");
    return histo;
}






void histo_2D_theo_vois(double moy, double var, double AC, int name)
{
    ImageGrayd histo(256, 256, true);

    for(int x=0; x<256; x++)
    {
        for(int y=0; y<256; y++)
        {
            double X = x/256.;
            double Y = y/256.;

            double gaussian = gaussAC_2D(moy, var, AC, X, Y);

            histo.pixelAbsolute(x, y) = gaussian;
        }
    }


    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_theo_vois_" + std::to_string(name) + ".png");
}




#endif //ASTEX_NUMERICAL_TVOIS_H
