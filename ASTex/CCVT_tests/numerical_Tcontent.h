//
// Created by grenier on 11/10/23.
//

#ifndef ASTEX_NUMERICAL_TCONTENT_H
#define ASTEX_NUMERICAL_TCONTENT_H



#include "tools.h"

//-----------------------------------------------------------------------------
struct color_info{
    ImageGrayd::PixelType couleur_;
    int compteur_;
    double compteurD_;

    color_info(ImageGrayd::PixelType couleur){
        couleur_ = couleur;
        compteur_ = 1;
    }

    color_info(ImageGrayd::PixelType couleur, int i){
        couleur_ = couleur;
        compteur_ = i;
    }

    color_info(ImageGrayd::PixelType couleur, double i){
        couleur_ = couleur;
        compteurD_ = i;
    }

    color_info(ImageGrayd::PixelType couleur, double i, int j){
        couleur_ = couleur;
        compteurD_ = i;
        compteur_ = j;
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

    bool operator==(const color_info& couleur2)
    {
        return couleur_ == couleur2.couleur_;
    }

    bool operator<(const color_info& couleur2)
    {
        return couleur_ < couleur2.couleur_;
    }
};



struct color_info_large{
    ImageGrayd::PixelType couleur_;
    std::vector<double> proportion_E;
    std::vector<double> proportion_H_t;
    std::vector<double> proportion_H_r;
    std::vector<double> error_t;
    std::vector<double> error_r;
    std::vector<double> distance_t;
    std::vector<double> distance_r;

    color_info_large(ImageGrayd::PixelType couleur){
        couleur_ = couleur;
    }
};



struct comput_mean{
    ImageGrayd::PixelType couleur_;
    double sum_;
    std::vector<double> values_;

    comput_mean(ImageGrayd::PixelType couleur){
        couleur_ = couleur;
        sum_ = 0.;
    }

    void incr(double value){
        sum_ += value;
        values_.push_back(value);
    }
};



bool is_in(std::vector<color_info>& vecteur, ImageGrayd::PixelType& couleur, int& id)
{
    id = -1;
//    for(ImageRGB8::PixelType col : vecteur)
    for(int i=0; i<vecteur.size(); i++)
    {
        if(vecteur.at(i) == couleur)
        {
            id = i;
            return true;
        }
    }
    return false;
}


std::vector<color_info> Tcontent(ImageGrayd image_)
{
    std::vector<color_info> couleurs = {};


    // listage des couleur et compte de leur apparition
    image_.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                          {
                              int id;
                              bool presence = is_in(couleurs, P, id);

                              if(not presence)
                              {
                                  couleurs.push_back(color_info{P});
                              }
                              else
                              {
                                  couleurs.at(id).incr();
                              }
                          });
    return couleurs;
}









std::vector<color_info> cell_capacity(ImageGrayd cm, double moy1, double moy2, double var1, double var2)
{
    std::vector<color_info> couleurs = {};
    int cm_size = cm.height();

    // listage des couleur et compte de leur apparition
    cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                      {
                          int id;

                          double X = (x+0.5)/double(cm_size-1); // évaluation au centre des pixels
                          double Y = (y+0.5)/double(cm_size-1);
                          double qtt =  gauss(moy1, moy2, var1, var2, X, Y);

                          bool presence = is_in(couleurs, P, id);
                          if(not presence)
                          {
                              couleurs.push_back(color_info{P, qtt});
                          }
                          else
                          {
                              couleurs.at(id).incr(qtt);
                          }
                      });
    return couleurs;
}





std::vector<color_info> cell_capacity_real_histo(ImageGrayd cm, ImageGrayd histo)
{
    std::vector<color_info> couleurs = {};

    // listage des couleur et compte de leur apparition
    cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                      {
                          int id;
                          double qtt = histo.pixelAbsolute(x,y);

                          bool presence = is_in(couleurs, P, id);
                          if(not presence)
                          {
                              couleurs.push_back(color_info{P, qtt});
                          }
                          else
                          {
                              couleurs.at(id).incr(qtt);
                          }
                      });
    return couleurs;
}


std::vector<color_info> cell_dist_real_histo(ImageGrayd cm, ImageGrayd histo_dist, double moy1, double moy2, double var1, double var2)
{
    std::vector<color_info> couleurs = {};

    // listage des couleur et compte de leur apparition
    cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                      {
                          int id;
                          double qtt_eps = histo_dist.pixelAbsolute(x,y);

                          double X = x/double(cm.height()-1);
                          double Y = y/double(cm.width()-1);
                          double qtt_t =  gauss(moy1, moy2, var1, var2, X, Y);

                          bool presence = is_in(couleurs, P, id);
                          if(not presence)
                          {
                              couleurs.push_back(color_info{P, qtt_t-qtt_eps, 1});
                          }
                          else
                          {
                              couleurs.at(id).incr(qtt_t-qtt_eps);
                              couleurs.at(id).incr();
                          }
                      });
    return couleurs;
}





#endif //ASTEX_NUMERICAL_TCONTENT_H
