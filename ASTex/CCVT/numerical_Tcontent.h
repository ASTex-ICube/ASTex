//
// Created by grenier on 11/10/23.
//

#ifndef ASTEX_NUMERICAL_TCONTENT_H
#define ASTEX_NUMERICAL_TCONTENT_H





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




// ---------------------------------------------------------------------------
double gauss(double moy1, double moy2, double var1, double var2, double x, double y)
{
    // gaussienne en x
    double X = x-moy1;
    double G1 = std::exp(-(X*X)/(2.*var1))/(std::sqrt(2.*M_PI*var1));

    // gaussienne en y
    double Y = y-moy2;
    double G2 = std::exp(-(Y*Y)/(2.*var2))/(std::sqrt(2.*M_PI*var2));

    return G1*G2;
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


// ---------------------------------------------------------------------------
ImageGrayd histo_2D(ImageGrayd noise1, ImageGrayd noise2, int histo_size)
{
//    int histo_size = 256;
    double norme = 0.;
    ImageGrayd histo(histo_size, histo_size, true);

    // énumération des paire de valeur de bruits
    for(int x=0; x<noise1.width(); x++)
    {
        for(int y=0; y<noise1.height(); y++)
        {
            double nx = noise1.pixelAbsolute(x, y)*(histo_size-1);
            double ny = noise2.pixelAbsolute(x, y)*(histo_size-1);

            histo.pixelAbsolute(int(std::round(nx)), int(std::round(ny))) += 1.;
//            histo.pixelAbsolute(int(std::floor(nx)), int(std::floor(ny))) += 1.;
//            histo.pixelAbsolute(int(std::round(nx)), int(std::round(ny))) += 1./(histo_size * histo_size);
        }
    }

    // calcul de l'aire sous la courbe
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             norme += P / (histo_size * histo_size);
                         });


    // normalisation pour une aire sous la courbe égale à 1
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             P *= 1./double(norme);
                         });

    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_reel.png");
    return histo;
}




ImageGrayd histo_2D_theo(double moy1, double moy2, double var1, double var2, int histo_size)
{
//    int histo_size = 256;
    ImageGrayd histo(histo_size, histo_size, true);

    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                      {
                          double X = x/double(histo_size-1);
                          double Y = y/double(histo_size-1);

                          P = gauss(moy1, moy2, var1, var2, X, Y);
                      });


    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_theo.png");
    return histo;
}



ImageGrayd histo_2D_dist(ImageGrayd histo_reel, double moy1, double moy2, double var1, double var2, int histo_size)
{
//    int histo_size = 256;
    ImageGrayd histo(histo_size, histo_size, true);

    double max_theo = -1.;
    double max_reel = -1.;
    double max_dist = -1.;

    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             double X = x/double(histo_size-1);
                             double Y = y/double(histo_size-1);

                             double val_theo = gauss(moy1, moy2, var1, var2, X, Y);
                             double val_reel = histo_reel.pixelAbsolute(x,y);
                             double val_dist = val_theo - val_reel;

                             max_theo = std::max(max_theo,val_theo);
                             max_reel = std::max(max_reel,val_reel);
                             max_dist = std::max(max_dist,std::abs(val_dist));

                             P = val_dist;
                         });

    std::cout<<"val max histo réel : "<<max_reel<<", val max histo théo : "<<max_theo<<", val max distance : "<<max_dist<<std::endl;
    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_dist.png");
    return histo;
}



#endif //ASTEX_NUMERICAL_TCONTENT_H
