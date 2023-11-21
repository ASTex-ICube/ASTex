//
// Created by grenier on 17/11/23.
//

#ifndef ASTEX_REYNOLD_H
#define ASTEX_REYNOLD_H

struct Graine{
    double x_;
    double y_;
    double weight_;

    Graine(){
        x_=0.;
        y_=0.;
        weight_ = 1.;
    }

    Graine(double x, double y, double w){
        x_ = x;
        y_ = y;
        weight_ = w;
    }
};



double grain_dist(double couleur1, double couleur2, std::vector<double> H_color, std::vector<Graine> H_seeds){
    int id1 = -1;
    int id2 = -1;
    for(int i=0; i<H_color.size(); i++){
        if(H_color.at(i)==couleur1){id1=i;}
        if(H_color.at(i)==couleur2){id2=i;}
    }

    Graine g1 = H_seeds.at(id1);
    Graine g2 = H_seeds.at(id2);

    return sqrt((g1.x_-g2.x_)*(g1.x_-g2.x_) + (g1.y_-g2.y_)*(g1.y_-g2.y_));
}






std::vector<color_vois> Tcontent_vois_test(ImageGrayd image_)
{
    std::vector<color_vois> couleurs = {};

    int id;
    bool presence;
    ImageGrayd::PixelType P2;

    image_.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                          {
                              if(x+1 < image_.width()){
                                  P2 = image_.pixelAbsolute(x+1, y);

                                  if(P != P2){
                                      presence = is_in(couleurs, P, P2, id);
                                      if(not presence)
                                      {
                                          couleurs.push_back(color_vois{P, P2});
                                      }
                                      else
                                      {
                                          couleurs.at(id).incr();
                                      }
                                  }

                                  if(y+1 < image_.height()){
                                      P2 = image_.pixelAbsolute(x+1, y+1);

                                      if(P != P2){
                                          presence = is_in(couleurs, P, P2, id);
                                          if(not presence)
                                          {
                                              couleurs.push_back(color_vois{P, P2});
                                          }
                                          else
                                          {
                                              couleurs.at(id).incr();
                                          }
                                      }
                                  }

                                  if(y-1 > 0){
                                      P2 = image_.pixelAbsolute(x+1, y-1);

                                      if(P != P2){
                                          presence = is_in(couleurs, P, P2, id);
                                          if(not presence)
                                          {
                                              couleurs.push_back(color_vois{P, P2});
                                          }
                                          else
                                          {
                                              couleurs.at(id).incr();
                                          }
                                      }
                                  }
                              }

                              if(x-1 > 0){
                                  ImageGrayd::PixelType P2 = image_.pixelAbsolute(x-1, y);

                                  if(P != P2){
                                      presence = is_in(couleurs, P, P2, id);
                                      if(not presence)
                                      {
                                          couleurs.push_back(color_vois{P, P2});
                                      }
                                      else
                                      {
                                          couleurs.at(id).incr();
                                      }
                                  }

                                  if(y+1 < image_.height()){
                                      P2 = image_.pixelAbsolute(x-1, y+1);

                                      if(P != P2){
                                          presence = is_in(couleurs, P, P2, id);
                                          if(not presence)
                                          {
                                              couleurs.push_back(color_vois{P, P2});
                                          }
                                          else
                                          {
                                              couleurs.at(id).incr();
                                          }
                                      }
                                  }

                                  if(y-1 > 0){
                                      P2 = image_.pixelAbsolute(x-1, y-1);

                                      if(P != P2){
                                          presence = is_in(couleurs, P, P2, id);
                                          if(not presence)
                                          {
                                              couleurs.push_back(color_vois{P, P2});
                                          }
                                          else
                                          {
                                              couleurs.at(id).incr();
                                          }
                                      }
                                  }
                              }


                              if(y+1 < image_.height()){
                                  ImageGrayd::PixelType P2 = image_.pixelAbsolute(x, y+1);

                                  if(P != P2){
                                      presence = is_in(couleurs, P, P2, id);
                                      if(not presence)
                                      {
                                          couleurs.push_back(color_vois{P, P2});
                                      }
                                      else
                                      {
                                          couleurs.at(id).incr();
                                      }
                                  }
                              }

                              if(y-1 > 0){
                                  ImageGrayd::PixelType P2 = image_.pixelAbsolute(x, y-1);

                                  if(P != P2){
                                      presence = is_in(couleurs, P, P2, id);
                                      if(not presence)
                                      {
                                          couleurs.push_back(color_vois{P, P2});
                                      }
                                      else
                                      {
                                          couleurs.at(id).incr();
                                      }
                                  }
                              }


                          });

    return couleurs;
}




std::vector<color_vois> Hcontent_vois_test(ImageGrayd cm_, double moy1, double moy2, double var1, double var2, ImageGrayd histo_)
{
    std::vector<color_vois> couleurs = {};

    int id;
    bool presence;
    ImageGrayd::PixelType P2;

    cm_.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                          {
                              double X = x/double(cm_.height()-1);
                              double Y = y/double(cm_.width()-1);
                              double qtt = gauss(moy1, moy2, var1, var2, X, Y);
//                              double qtt = histo_.pixelAbsolute(x,y);

                              if(x+1 < cm_.width()){
                                  P2 = cm_.pixelAbsolute(x+1, y);

                                  if(P != P2){
                                      presence = is_in(couleurs, P, P2, id);
                                      if(not presence)
                                      {
                                          couleurs.push_back(color_vois{P, P2, qtt});
                                      }
                                      else
                                      {
                                          couleurs.at(id).incr(qtt);
                                      }
                                  }

                                  if(y+1 < cm_.height()){
                                      P2 = cm_.pixelAbsolute(x+1, y+1);

                                      if(P != P2){
                                          presence = is_in(couleurs, P, P2, id);
                                          if(not presence)
                                          {
                                              couleurs.push_back(color_vois{P, P2, qtt});
                                          }
                                          else
                                          {
                                              couleurs.at(id).incr(qtt);
                                          }
                                      }
                                  }

                                  if(y-1 > 0){
                                      P2 = cm_.pixelAbsolute(x+1, y-1);

                                      if(P != P2){
                                          presence = is_in(couleurs, P, P2, id);
                                          if(not presence)
                                          {
                                              couleurs.push_back(color_vois{P, P2, qtt});
                                          }
                                          else
                                          {
                                              couleurs.at(id).incr(qtt);
                                          }
                                      }
                                  }
                              }

                              if(x-1 > 0){
                                  ImageGrayd::PixelType P2 = cm_.pixelAbsolute(x-1, y);

                                  if(P != P2){
                                      presence = is_in(couleurs, P, P2, id);
                                      if(not presence)
                                      {
                                          couleurs.push_back(color_vois{P, P2, qtt});
                                      }
                                      else
                                      {
                                          couleurs.at(id).incr(qtt);
                                      }
                                  }

                                  if(y+1 < cm_.height()){
                                      P2 = cm_.pixelAbsolute(x-1, y+1);

                                      if(P != P2){
                                          presence = is_in(couleurs, P, P2, id);
                                          if(not presence)
                                          {
                                              couleurs.push_back(color_vois{P, P2, qtt});
                                          }
                                          else
                                          {
                                              couleurs.at(id).incr(qtt);
                                          }
                                      }
                                  }

                                  if(y-1 > 0){
                                      P2 = cm_.pixelAbsolute(x-1, y-1);

                                      if(P != P2){
                                          presence = is_in(couleurs, P, P2, id);
                                          if(not presence)
                                          {
                                              couleurs.push_back(color_vois{P, P2, qtt});
                                          }
                                          else
                                          {
                                              couleurs.at(id).incr(qtt);
                                          }
                                      }
                                  }
                              }


                              if(y+1 < cm_.height()){
                                  P2 = cm_.pixelAbsolute(x, y+1);

                                  if(P != P2){
                                      presence = is_in(couleurs, P, P2, id);
                                      if(not presence)
                                      {
                                          couleurs.push_back(color_vois{P, P2, qtt});
                                      }
                                      else
                                      {
                                          couleurs.at(id).incr(qtt);
                                      }
                                  }
                              }

                              if(y-1 > 0){
                                  P2 = cm_.pixelAbsolute(x, y-1);

                                  if(P != P2){
                                      presence = is_in(couleurs, P, P2, id);
                                      if(not presence)
                                      {
                                          couleurs.push_back(color_vois{P, P2, qtt});
                                      }
                                      else
                                      {
                                          couleurs.at(id).incr(qtt);
                                      }
                                  }
                              }

                          });

    return couleurs;
}








#endif //ASTEX_REYNOLD_H
