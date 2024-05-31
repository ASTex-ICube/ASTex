//
// Created by grenier on 14/05/24.
//

#ifndef ASTEX_MESURE_FULL_H
#define ASTEX_MESURE_FULL_H

#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

#include "ASTex/CCVT_tests/mesure_presence.h"
#include "ASTex/CCVT_tests/mesure_voisinage.h"
#include "ASTex/CCVT_tests/mesure_triplet.h"

#include "ASTex/CCVT_tests/tools.h"

using namespace ASTex;

void mesure_E(ImageRGBu8 image_,
              std::vector<color_info>& presences, std::vector<color_vois>& voisinages, std::vector<color_triple>& triplets,
              bool i_isEqual_j,
              int& nb_vois, int& nb_vois_diff){

    int id;
    bool presence;
    ImageRGBu8::PixelType P2;
    ImageRGBu8::PixelType P3;
    ImageRGBu8::PixelType P4;

    std::vector<color_info> couleurs_pres = {};
    std::vector<color_vois> couleurs_vois = {};
    std::vector<color_triple> couleur_tripl = {};

    int p = 0;
    int d_h = 0;
    int d_v = 0;
    int t = 0;

    image_.for_all_pixels([&] (typename ImageRGBu8::PixelType& P, int x, int y)
                          {
                              // présence
                              color_info col{P};
                              update_Tcontent(couleurs_pres, col);
                              p += 1;

                              // voisinage
//                              if(x+1 < image_.width()){
//                                  P2 = image_.pixelAbsolute(x+1, y);
//                                  update_Tcontent(couleurs_vois, P, P2, i_isEqual_j);
//
//                                  if(P != P2){nb_vois_diff += 1;}
//                                  nb_vois +=1;
//                              }
//
//                              if(x-1 >= 0){
//                                  P2 = image_.pixelAbsolute(x-1, y);
//                                  update_Tcontent(couleurs_vois, P, P2, i_isEqual_j);
//
//                                  if(P != P2){nb_vois_diff += 1;}
//                                  nb_vois +=1;
//                              }
//
//                              if(y+1 < image_.height()){
//                                  P2 = image_.pixelAbsolute(x, y+1);
//                                  update_Tcontent(couleurs_vois, P, P2, i_isEqual_j);
//
//                                  if(P != P2){nb_vois_diff += 1;}
//                                  nb_vois +=1;
//                              }
//
//                              if(y-1 >= 0){
//                                  P2 = image_.pixelAbsolute(x, y-1);
//                                  update_Tcontent(couleurs_vois, P, P2, i_isEqual_j);
//
//                                  if(P != P2){nb_vois_diff += 1;}
//                                  nb_vois +=1;
//                              }



                              if(x+1 < image_.width()){
                                  d_v += 1;
                                  P2 = image_.pixelAbsolute(x+1, y);

                                  color_vois doublet{P, P2};
                                  bool is_doublet = get_doublet(doublet, P, P2);

                                  if(is_doublet){
                                      update_Tcontent(couleurs_vois, doublet);
                                  }
                              }

                              if(y+1 < image_.height()){
                                  d_h += 1;
                                  P2 = image_.pixelAbsolute(x, y+1);

                                  color_vois doublet{P, P2};
                                  bool is_doublet = get_doublet(doublet, P, P2);

                                  if(is_doublet){
                                      update_Tcontent(couleurs_vois, doublet);
                                  }
                              }






                              // triplet
                              if(x+1 < image_.width() and y+1 < image_.height()){
                                  t += 1;
                                  P2 = image_.pixelAbsolute(x+1, y);
                                  P3 = image_.pixelAbsolute(x, y+1);
                                  P4 = image_.pixelAbsolute(x+1, y+1);

                                  color_triple triplet{P2, P3, P4};
                                  bool is_triplet = get_triplet(triplet, P, P2, P3, P4);

                                  if(is_triplet){
                                      update_Tcontent(couleur_tripl, triplet);
                                  }
                              }
                          });

    std::cout<<"pixel : "<<p<<", frontière : "<<d_v+d_h<<", triple : "<<t<<std::endl;

    presences = couleurs_pres;
    voisinages = couleurs_vois;
    triplets = couleur_tripl;
}





void mesure_H(ImageRGBu8 cm,
              std::vector<color_info>& presence_theo, std::vector<color_info>& presence_reel,
              std::vector<color_vois>& voisinage_theo, std::vector<color_vois>& voisinage_reel,
              std::vector<color_triple>& triplet_theo, std::vector<color_triple>& triplet_reel,
              double moy1, double moy2, double var1, double var2,
              ImageGrayd histo, bool i_isEqual_j){

    std::vector<color_info> couleurs_presence_theo = {};
    std::vector<color_info> couleurs_presence_reel = {};

    std::vector<color_vois> couleurs_voisinage_theo = {};
    std::vector<color_vois> couleurs_voisinage_reel = {};

    std::vector<color_triple> couleurs_triple_theo = {};
    std::vector<color_triple> couleurs_triple_reel = {};

    double p = 0;
    double d_h = 0;
    double d_v = 0;
    double t = 0;

    int cm_size = cm.height();
    int id;
    bool presence;
    ImageRGBu8::PixelType P2;
    ImageRGBu8::PixelType P3;
    ImageRGBu8::PixelType P4;

    int tot_pixel = 0;
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y){
        tot_pixel += P;
    });

    double X;
    double Y;

    int qtt_r;
    double qtt_t;

    // listage des couleur et compte de leur apparition
    cm.for_all_pixels([&] (typename ImageRGBu8::PixelType& P, int x, int y)
                      {

                          // presence histo théorique
                          X = (x+0.5);//double(cm_size-1); // évaluation au centre des pixels
                          Y = (y+0.5);//double(cm_size-1);

                          qtt_t = tot_pixel * gauss(moy1, moy2, var1, var2, X, Y, double(cm_size-1));

                          // présence histo réel
                          qtt_r = histo.pixelAbsolute(x,y);

                          p += qtt_t;
//                          std::cout<<"qtt_t : "<<qtt_t<<std::endl;

                          // présence
                          color_info col{P};
                          update_Tcontent(couleurs_presence_reel, col, qtt_r);
                          update_Tcontent(couleurs_presence_theo, col, qtt_t);








                          //voisinages
                          if(x+1 < cm.width()){
                              X = (x+1.0);//double(cm_size-1); // évaluation au centre des pixels
                              Y = (y+0.5);//double(cm_size-1);

                              qtt_t = tot_pixel * gauss(moy1, moy2, var1, var2, X, Y, double(cm_size-1));


                              qtt_r = 0.5*(histo.pixelAbsolute(x,y) + histo.pixelAbsolute(x+1,y));
                              d_v += qtt_r;

                              P2 = cm.pixelAbsolute(x+1, y);

                              color_vois doublet{P, P2};
                              bool is_doublet = get_doublet(doublet, P, P2);

                              if(is_doublet){
                                  update_Tcontent(couleurs_voisinage_reel, doublet, qtt_r);
                                  update_Tcontent(couleurs_voisinage_theo, doublet, qtt_t);
                              }
                          }

                          if(y+1 < cm.height()){
                              X = (x+0.5);//double(cm_size-1); // évaluation au centre des pixels
                              Y = (y+1.0);//double(cm_size-1);

                              qtt_t = tot_pixel * gauss(moy1, moy2, var1, var2, X, Y, double(cm_size-1));

                              qtt_r = 0.5*(histo.pixelAbsolute(x,y) + histo.pixelAbsolute(x,y+1));
                              d_h += qtt_r;

                              P2 = cm.pixelAbsolute(x, y+1);

                              color_vois doublet{P, P2};
                              bool is_doublet = get_doublet(doublet, P, P2);

                              if(is_doublet){
                                  update_Tcontent(couleurs_voisinage_reel, doublet, qtt_r);
                                  update_Tcontent(couleurs_voisinage_theo, doublet, qtt_t);
                              }
                          }
//


//                          // voisinages
//                          if(x+1 < cm_size){
//                              P2 = cm.pixelAbsolute(x+1, y);
//                              update_Hcontent_theo(couleurs_voisinage_theo, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, 1, 0);
//                              update_Hcontent_reel(couleurs_voisinage_reel, P, P2, histo, i_isEqual_j, x, y, 1, 0);
//
////                              if(P != P2){nb_vois_diff += 1;}
////                              nb_vois +=1;
//                          }
//
//                          if(x-1 >= 0){
//                              P2 = cm.pixelAbsolute(x-1, y);
//                              update_Hcontent_theo(couleurs_voisinage_theo, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, -1, 0);
//                              update_Hcontent_reel(couleurs_voisinage_reel, P, P2, histo, i_isEqual_j, x, y, -1, 0);
//
////                              if(P != P2){nb_vois_diff += 1;}
////                              nb_vois +=1;
//                          }
//
//                          if(y+1 < cm_size){
//                              P2 = cm.pixelAbsolute(x, y+1);
//                              update_Hcontent_theo(couleurs_voisinage_theo, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, 0, 1);
//                              update_Hcontent_reel(couleurs_voisinage_reel, P, P2, histo, i_isEqual_j, x, y, 0, 1);
//
////                              if(P != P2){nb_vois_diff += 1;}
////                              nb_vois +=1;
//                          }
//
//                          if(y-1 >= 0){
//                              P2 = cm.pixelAbsolute(x, y-1);
//                              update_Hcontent_theo(couleurs_voisinage_theo, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, 0, -1);
//                              update_Hcontent_reel(couleurs_voisinage_reel, P, P2, histo, i_isEqual_j, x, y, 0, -1);
//
////                              if(P != P2){nb_vois_diff += 1;}
////                              nb_vois +=1;
//                          }



                          // triplet
                          if(x+1 < cm.width() and y+1 < cm.height()){
                              X = (x+1.0);//double(cm_size-1); // évaluation au centre des pixels
                              Y = (y+1.0);//double(cm_size-1);

                              qtt_t = tot_pixel * gauss(moy1, moy2, var1, var2, X, Y, double(cm_size-1));


                              qtt_r = 0.25*(histo.pixelAbsolute(x,y) +
                                           histo.pixelAbsolute(x,y+1) +
                                           histo.pixelAbsolute(x+1,y) +
                                           histo.pixelAbsolute(x+1,y+1));
                              t += qtt_t;

                              P2 = cm.pixelAbsolute(x+1, y);
                              P3 = cm.pixelAbsolute(x, y+1);
                              P4 = cm.pixelAbsolute(x+1, y+1);

                              color_triple triplet{P2, P3, P4};
                              bool is_triplet = get_triplet(triplet, P, P2, P3, P4);

                              if(is_triplet){
                                  update_Tcontent(couleurs_triple_reel, triplet, qtt_r);
                                  update_Tcontent(couleurs_triple_theo, triplet, qtt_t);
                              }
                          }


                      });

    std::cout<<"pixel : "<<p<<", frontière : "<<d_v+d_h<<", triple : "<<t<<std::endl;
//    std::cout<<p<<" "<<p_r<<std::endl;

    presence_theo = couleurs_presence_theo;
    presence_reel = couleurs_presence_reel;

    voisinage_theo = couleurs_voisinage_theo;
    voisinage_reel = couleurs_voisinage_reel;

    triplet_theo = couleurs_triple_theo;
    triplet_reel = couleurs_triple_reel;

}

#endif //ASTEX_MESURE_FULL_H
