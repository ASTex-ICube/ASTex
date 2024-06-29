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


struct color_voisinage_large{
    ImageGrayd::PixelType couleur_1;
    ImageGrayd::PixelType couleur_2;
    std::vector<double> proportion_E;
    std::vector<double> proportion_H_t;
    std::vector<double> proportion_H_r;
    std::vector<double> error_t;
    std::vector<double> error_r;
    std::vector<double> distance_t;
    std::vector<double> distance_r;

    color_voisinage_large(ImageGrayd::PixelType couleur1, ImageGrayd::PixelType couleur2){
        couleur_1 = couleur1;
        couleur_2 = couleur2;
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
                                                    bool presence;

                                                    // évaluation au centre des pixels
                                                    double X1 = (x1+0.5)/double(cm_size-1); // x1 ~ Nu
                                                    double Y1 = (y1+0.5)/double(cm_size-1); // y1 ~ N'u
                                                    double X2 = (x2+0.5)/double(cm_size-1); // x2 ~ Nv
                                                    double Y2 = (y2+0.5)/double(cm_size-1); // y2 ~ N'v

                                                    double qtt = gaussAC_2D(moy1, var1, AC1, X1, X2) * gaussAC_2D(moy2, var2, AC2, Y1, Y2);

//                                                    if(P1 != P2){
                                                        presence = is_in(couleurs, P1, P2, id);
                                                        if(not presence)
                                                        {
                                                            presence = is_in(couleurs, P2, P1, id);
                                                            if (not presence)
                                                            {
                                                                ImageGrayd::PixelType Pi = std::min(P1,P2);
                                                                ImageGrayd::PixelType Pj = std::max(P1,P2);
                                                                couleurs.push_back(color_vois{Pi, Pj, qtt});
                                                            }
                                                            else
                                                            {
                                                                couleurs.at(id).incr(qtt);
                                                            }
                                                        }
                                                        else
                                                        {
                                                            couleurs.at(id).incr(qtt);
                                                        }

//                                                    }
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













#endif //ASTEX_NUMERICAL_TVOIS_H
