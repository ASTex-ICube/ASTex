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

    bool operator<(const color_vois& couleur2)
    {
        if(couleur1_ == couleur2.couleur1_){
            return (couleur2_ < couleur2.couleur2_);
        }
        else{
            return (couleur1_ < couleur2.couleur1_);
        }
//        return (couleur1_ < couleur2.couleur1_);// and (couleur2_ < couleur2.couleur2_);
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
            presence = is_in(couleurs, image_.pixelAbsolute(x,y), image_.pixelAbsolute(x+tauX, y+tauY), id);
            if(not presence)
            {
                couleurs.push_back(color_vois{image_.pixelAbsolute(x,y), image_.pixelAbsolute(x+tauX, y+tauY)});
            }
            else
            {
                couleurs.at(id).incr();
            }

            // symétrie ?
            presence = is_in(couleurs, image_.pixelAbsolute(x+tauX, y+tauY), image_.pixelAbsolute(x,y), id);
            if(not presence)
            {
                couleurs.push_back(color_vois{image_.pixelAbsolute(x+tauX, y+tauY), image_.pixelAbsolute(x,y)});
            }
            else
            {
                couleurs.at(id).incr();
            }

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

    /* matrice précision (inverse matrice covariance) :
     * var/det -AC/det
     * -AC/det  var/det
     */
    double var_inv = var/det;
    double AC_inv = -AC/det;


    double X = x-moy;
    double Y = y-moy;

    double produit = X*(X*var_inv + Y*AC_inv) + Y*(X*AC_inv + Y*var_inv);
    double norm = 2.*M_PI*std::sqrt(det);

    double G = std::exp(-0.5*produit)/norm;

    return G;
}







//double quantified_gaussAC_2D(int x1, int x2, int y1, int y2, int cm_size, double moy1, double moy2, double var1, double var2, double AC1, double AC2)
//{
//    // x1 ~ Nu, y1 ~ N'u, x2 ~ Nv, y2 ~ N'v
//    double X1 = x1/double(cm_size);
//    double Y1 = y1/double(cm_size);
//
//    double X2 = x2/double(cm_size);
//    double Y2 = y2/double(cm_size);
//
//    double gaussian1 = gaussAC_2D(moy1, var1, AC1, X1, X2);
//    double gaussian2 = gaussAC_2D(moy2, var2, AC2, Y1, Y2);
//
//    double Gaussian = gaussian1*gaussian2;//*100;//img_size*img_size;
//
//    return Gaussian;
//}





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
    int cm_size = cm.height();

    // listage des couleur et compte de leur apparition
    // x1 ~ Nu, y1 ~ N'u, x2 ~ Nv, y2 ~ N'v
    cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P1, int x1, int y1) // for_all_pixels // parallel_for_all_pixels
                      {
//                          if(x1==0){std::cout<<y1<<" ";}
//                          std::cout<<x1<<" ";
                          cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P2, int x2, int y2)
                                            {
//                                                if(x2==0){std::cout<<y2<<" ";}
//                                                if(P1 != P2){
                                                    int id;

                                                    // évaluation au centre des pixels
                                                    double X1 = (x1+0.5)/double(cm_size-1); // x1 ~ Nu
                                                    double Y1 = (y1+0.5)/double(cm_size-1); // y1 ~ N'u
                                                    double X2 = (x2+0.5)/double(cm_size-1); // x2 ~ Nv
                                                    double Y2 = (y2+0.5)/double(cm_size-1); // y2 ~ N'v

                                                    double qtt = gaussAC_2D(moy1, var1, AC1, X1, X2) * gaussAC_2D(moy2, var2, AC2, Y1, Y2);

                                                    bool presence = is_in(couleurs, P1, P2, id);
                                                    if(not presence)
                                                    {
                                                        couleurs.push_back(color_vois{P1, P2, qtt});
                                                    }
                                                    else
                                                    {
                                                        couleurs.at(id).incr(qtt);
                                                    }
//                                                }


                                            });

                      });
//    std::cout<<std::endl;
    return couleurs;
}







std::vector<color_vois> capacityH_vois_histo_reel(ImageGrayd cm, ImageGrayd histo_1, ImageGrayd histo_2)
{
    std::vector<color_vois> couleurs = {};

    // listage des couleur et compte de leur apparition
    // x1 ~ Nu, y1 ~ N'u, x2 ~ Nv, y2 ~ N'v
    cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P1, int x1, int y1)
                      {
                          cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P2, int x2, int y2)
                                            {
                                                int id;
                                                // histo_1 : densité couple (Nu,Nv), histo_2 : densité couple (N'u, N'v)
                                                double qtt = histo_1.pixelAbsolute(x1,x2)*histo_2.pixelAbsolute(y1,y2);

                                                bool presence = is_in(couleurs, P1, P2, id);
                                                if(not presence)
                                                {
                                                    couleurs.push_back(color_vois{P1, P2, qtt});
                                                }
                                                else
                                                {
                                                    couleurs.at(id).incr(qtt);
                                                }

                                            });
                      });
    return couleurs;
}



std::vector<color_vois>  capacityH_vois_histo_dist(ImageGrayd cm, ImageGrayd histo_dist_N1, ImageGrayd histo_dist_N2, double moy1, double moy2, double var1, double var2, double AC1, double AC2)
{
    std::vector<color_vois>  couleurs = {};

    // listage des couleur et compte de leur apparition
    // x1 ~ Nu, y1 ~ N'u, x2 ~ Nv, y2 ~ N'v
    cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P1, int x1, int y1)
                      {
                          cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P2, int x2, int y2)
                                            {
                                                int id;
                                                double qtt_eps = histo_dist_N1.pixelAbsolute(x1,x2)*histo_dist_N2.pixelAbsolute(y1,y2);

                                                double X1 = x1/double(cm.height()-1); // x1 ~ Nu
                                                double Y1 = y1/double(cm.width()-1); // y1 ~ N'u
                                                double X2 = x2/double(cm.height()-1); // x2 ~ Nv
                                                double Y2 = y2/double(cm.width()-1); // y2 ~ N'v

                                                double qtt_t = gaussAC_2D(moy1, var1, AC1, X1, X2) * gaussAC_2D(moy2, var2, AC2, Y1, Y2);


                                                bool presence = is_in(couleurs, P1, P2, id);
                                                if(not presence)
                                                {
                                                    couleurs.push_back(color_vois{P1, P2, qtt_t-qtt_eps});
                                                }
                                                else
                                                {
                                                    couleurs.at(id).incr(qtt_t-qtt_eps);
                                                }

                                            });
//                          if(x1==0){std::cout<<y1<<" ";}
                      });
    return couleurs;
}










// ---------------------------------------------------------------------------
ImageGrayd histo_2D_AC(ImageGrayd noise1, int tauX, int tauY, int name)
{
    int histo_size = 256;
    double norme = 0;
    ImageGrayd histo(histo_size, histo_size, true);

    // énumération des paire de valeur de bruits
    for(int x=0; x<noise1.width()-tauX; x++)
    {
        for(int y=0; y<noise1.height()-tauY; y++)
        {
            double nx = noise1.pixelAbsolute(x, y)*(histo_size-1);
            double ny = noise1.pixelAbsolute(x+tauX, y+tauY)*(histo_size-1);

            histo.pixelAbsolute(int(std::round(nx)), int(std::round(ny))) += 1.;
            histo.pixelAbsolute(int(std::round(ny)), int(std::round(nx))) += 1.; // symétrie
        }
    }

    // calcul de l'aire sous la courbe
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             norme += P / (histo_size * histo_size);
                         });


    // estimation des moyennes
    double mX = 0.;
    double mY = 0.;
    double tot = 0.;
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             double X = x/double(histo_size-1);
                             double Y = y/double(histo_size-1);

                             mX += P*X; // x/255. : valeur sur 0, 1
                             mY += P*Y;

                             tot += P;
                         });
    mX = mX/tot;
    mY = mY/tot;
    std::cout<<"moyenne histo x : "<<mX<<", moyenne histo y : "<<mY<<std::endl;




    // estimation de l'auto-covariance
    double AC = 0.;
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             double X = x/double(histo_size-1);
                             double Y = y/double(histo_size-1);

                             AC += P*(X-mX)*(Y-mY); // x/255. : valeur sur 0, 1
//                             AC += P*(x-mX)*(y-mY); // x : valeur sur 0, 255
                         });
    AC = AC/tot;
    std::cout<<"auto-covariance histo : "<<AC<<std::endl;



    // normalisation pour une aire sous la courbe égale à 1
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             P *= 1./double(norme);
                         });
    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_AC_" + std::to_string(name) + ".png");
    return histo;
}




ImageGrayd histo_2D_theo_vois(double moy, double var, double AC, int name)
{
    int histo_size = 256;
    ImageGrayd histo(histo_size, histo_size, true);

    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             double X = x/double(histo_size-1);
                             double Y = y/double(histo_size-1);

                             P = gaussAC_2D(moy, var, AC, X, Y);
                         });


    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_theo_vois_" + std::to_string(name) + ".png");
    return histo;
}



ImageGrayd histo_2D_dist_vois(ImageGrayd histo_reel, double moy, double var, double AC, int name)
{
    int histo_size = 256;
    ImageGrayd histo(histo_size, histo_size, true);

    double max_theo = -1.;
    double max_reel = -1.;
    double max_dist = -1.;

    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             double X = x/double(histo_size-1);
                             double Y = y/double(histo_size-1);

                             double val_theo = gaussAC_2D(moy, var, AC, X, Y) ;
                             double val_reel = histo_reel.pixelAbsolute(x,y);
                             double val_dist = val_theo - val_reel;

                             max_theo = std::max(max_theo,val_theo);
                             max_reel = std::max(max_reel,val_reel);
                             max_dist = std::max(max_dist,std::abs(val_dist));

                             P = val_dist;
                         });

    std::cout<<"bruit "<<std::to_string(name)<<" - val max histo réel : "<<max_reel<<", val max histo théo : "<<max_theo<<", val max distance : "<<max_dist<<std::endl;
    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_dist_vois_" + std::to_string(name) + ".png");
    return histo;
}



#endif //ASTEX_NUMERICAL_TVOIS_H
