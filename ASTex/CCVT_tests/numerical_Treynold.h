//
// Created by grenier on 17/11/23.
//

#ifndef ASTEX_NUMERICAL_TREYNOLD_H
#define ASTEX_NUMERICAL_TREYNOLD_H

#include "ASTex/CCVT_tests/numerical_Tcontent.h"
#include "ASTex/CCVT_tests/numerical_Tvois.h"
#include "ASTex/CCVT_tests/numerical_Tstat.h"










void update_Tcontent(std::vector<color_vois>& couleurs,
                     ImageGrayd::PixelType P1, ImageGrayd::PixelType P2,
                     bool i_isEqual_j)
{
    int id;
    bool presence;

    if(P1 != P2){
        presence = is_in(couleurs, P1, P2, id);
        if(not presence)
        {
            presence = is_in(couleurs, P2, P1, id);
            if(not presence)
            {
                ImageGrayd::PixelType Pi = std::min(P1,P2);
                ImageGrayd::PixelType Pj = std::max(P1,P2);
                couleurs.push_back(color_vois{Pi, Pj});
            }
            else
            {
                couleurs.at(id).incr();
            }
        }
        else
        {
            couleurs.at(id).incr();
        }
    }


    else if(P1 == P2 and i_isEqual_j){
        presence = is_in(couleurs, P1, P2, id);
        if(not presence)
        {
            couleurs.push_back(color_vois{P1, P2});
        }
        else
        {
            couleurs.at(id).incr();
        }
    }


}



std::vector<color_vois> Tcontent_vois_test(ImageGrayd image_, bool i_isEqual_j)
{
    std::vector<color_vois> couleurs = {};

    ImageGrayd::PixelType P2;
    int nb_vois = 0;
    int nb_vois_diff = 0;
    int nb_vois_egl = 0;

//    bool i_isEqual_j = true; // true si on compte les i = j, false si on ne compte pas les i = j

    image_.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                          {
                              if(x+1 < image_.width()){
                                  P2 = image_.pixelAbsolute(x+1, y);
                                  update_Tcontent(couleurs, P, P2, i_isEqual_j);

                                  if(P != P2){nb_vois_diff += 1;}
                                  if(P == P2){nb_vois_egl += 1;}
                                  nb_vois +=1;
                              }

                              if(x-1 >= 0){
                                  P2 = image_.pixelAbsolute(x-1, y);
                                  update_Tcontent(couleurs, P, P2, i_isEqual_j);

                                  if(P != P2){nb_vois_diff += 1;}
                                  if(P == P2){nb_vois_egl += 1;}
                                  nb_vois +=1;
                              }

                              if(y+1 < image_.height()){
                                  P2 = image_.pixelAbsolute(x, y+1);
                                  update_Tcontent(couleurs, P, P2, i_isEqual_j);

                                  if(P != P2){nb_vois_diff += 1;}
                                  if(P == P2){nb_vois_egl += 1;}
                                  nb_vois +=1;
                              }

                              if(y-1 >= 0){
                                  P2 = image_.pixelAbsolute(x, y-1);
                                  update_Tcontent(couleurs, P, P2, i_isEqual_j);

                                  if(P != P2){nb_vois_diff += 1;}
                                  if(P == P2){nb_vois_egl += 1;}
                                  nb_vois +=1;
                              }


                          });
    std::cout<<"nb vois = "<<nb_vois<<", nb vois diff = "<<nb_vois_diff<<", nb vois egl = "<<nb_vois_egl<<std::endl;
    return couleurs;
}















void update_Hcontent_reel(std::vector<color_vois>& couleurs, ImageGrayd::PixelType P1, ImageGrayd::PixelType P2,
                     ImageGrayd histo_, bool i_isEqual_j,
                     int x, int y,  int dx, int dy) {

    int id;
    bool presence;
    double qtt;// = gauss(moy1, moy2, var1, var2, X, Y);
    int cm_size = histo_.height();


    // histo réel
    if(histo_.pixelAbsolute(x, y) < (1./double(cm_size - 1)) or histo_.pixelAbsolute(x + dx, y + dy) < (1./double(cm_size - 1)))
    {
        qtt = 0.;
    }
    else
    {
        qtt = 0.5*(histo_.pixelAbsolute(x, y) + histo_.pixelAbsolute(x + dx, y + dy));
    }


    if(P1 != P2){
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
    }


    else if (P1 == P2 and i_isEqual_j) {

        presence = is_in(couleurs, P1, P2, id);
        if (not presence) {
            couleurs.push_back(color_vois{P1, P2, qtt});
        } else {
            couleurs.at(id).incr(qtt);
        }
    }
}






void update_Hcontent_theo(std::vector<color_vois>& couleurs, ImageGrayd::PixelType P1, ImageGrayd::PixelType P2,
                          int cm_size, bool i_isEqual_j,
                     double moy1, double moy2, double var1, double var2,
                     int x, int y,  int dx, int dy) {

    int id;
    bool presence;
    double X;// = (x+0.5)/double(cm_size-1); // évaluation au centre des pixels
    double Y;// = (y+0.5)/double(cm_size-1);
    double qtt;// = gauss(moy1, moy2, var1, var2, X, Y);


    // histo théorique
    X = (x + 0.5 + double(dx) / 2.) / double(cm_size - 1); // évaluation au centre des bord des pixels
    Y = (y + 0.5 + double(dy) / 2.) / double(cm_size - 1);
    qtt = gauss(moy1, moy2, var1, var2, X, Y);


    if(P1 != P2){
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
    }


    else if (P1 == P2 and i_isEqual_j) {

        presence = is_in(couleurs, P1, P2, id);
        if (not presence) {
            couleurs.push_back(color_vois{P1, P2, qtt});
        } else {
            couleurs.at(id).incr(qtt);
        }
    }
}




std::vector<color_vois> Hcontent_vois_test(ImageGrayd cm_, double moy1, double moy2, double var1, double var2, ImageGrayd histo_, bool i_isEqual_j)
{
    std::vector<color_vois> couleurs = {};

    ImageGrayd::PixelType P2;
    int cm_size = cm_.height();


//    bool i_isEqual_j = true; // true si on compte les i = j, false si on ne compte pas les i = j

    cm_.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                          {
                              if(x+1 < cm_size){
                                  P2 = cm_.pixelAbsolute(x+1, y);
                                  update_Hcontent_theo(couleurs, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, 1, 0);
//                                  update_Hcontent_reel(couleurs, P, P2, histo_, i_isEqual_j, x, y, 1, 0);
                              }

                              if(x-1 >= 0){
                                  P2 = cm_.pixelAbsolute(x-1, y);
                                  update_Hcontent_theo(couleurs, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, -1, 0);
//                                  update_Hcontent_reel(couleurs, P, P2, histo_, i_isEqual_j, x, y, -1, 0);
                              }

                              if(y+1 < cm_size){
                                  P2 = cm_.pixelAbsolute(x, y+1);
                                  update_Hcontent_theo(couleurs, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, 0, 1);
//                                  update_Hcontent_reel(couleurs, P, P2, histo_, i_isEqual_j, x, y, 0, 1);
                              }

                              if(y-1 >= 0){
                                  P2 = cm_.pixelAbsolute(x, y-1);
                                  update_Hcontent_theo(couleurs, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, 0, -1);
//                                  update_Hcontent_reel(couleurs, P, P2, histo_, i_isEqual_j, x, y, 0, -1);
                              }

                          });

    return couleurs;
}









#endif //ASTEX_NUMERICAL_TREYNOLD_H
