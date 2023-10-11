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

                              if(couleurs.empty())
                              {
                                  couleurs.push_back(color_info{P});
                              }
                              else
                              {
                                  bool presence = is_in(couleurs, P, id);
                                  if(not presence)
                                  {
                                      couleurs.push_back(color_info{P});
                                  }
                                  else
                                  {
                                      couleurs.at(id).incr();
                                  }
                              }
                          });
    return couleurs;
}




// ---------------------------------------------------------------------------
double gauss(double moy1, double moy2, double var1, double var2, double x, double y)
{
    double G1 = std::exp(-(x-moy1)*(x-moy1)/(2*var1));///(std::sqrt(2.*M_PI*var1));
    double G2 = std::exp(-(y-moy2)*(y-moy2)/(2*var2));///(std::sqrt(2.*M_PI*var2));

    return G1*G2;
}

double quantified_gauss(int x, int y, int cm_size, int img_size, double moy1, double moy2, double var1, double var2)
{

    double X = x/double(cm_size);
    double Y = y/double(cm_size);

    double gaussian = gauss(moy1, moy2, var1, var2, X, Y);
    int Gaussian = gaussian*1000;//img_size*img_size;

    return gaussian;
}


std::vector<color_info> cell_capacity(ImageGrayd cm, int img_size, double moy1, double moy2, double var1, double var2)
{
    std::vector<color_info> couleurs = {};

    // listage des couleur et compte de leur apparition
    cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                      {
                          int id;
                          double qtt = quantified_gauss(x, y, cm.height(), img_size, moy1, moy2, var1, var2);

                          if(couleurs.empty())
                          {
                              couleurs.push_back(color_info{P, qtt});
                          }
                          else
                          {
                              bool presence = is_in(couleurs, P, id);
                              if(not presence)
                              {
                                  couleurs.push_back(color_info{P, qtt});
                              }
                              else
                              {
                                  couleurs.at(id).incr(qtt);
                              }
                          }
                      });
    return couleurs;
}



// ---------------------------------------------------------------------------
void histo_2D(ImageGrayd noise1, ImageGrayd noise2)
{
    ImageGrayd histo(256, 256, true);
    int max = 0;

    for(int x=0; x<noise1.width(); x++)
    {
        for(int y=0; y<noise1.height(); y++)
        {
            double nx = noise1.pixelAbsolute(x, y)*255;
            double ny = noise2.pixelAbsolute(x, y)*255;

//            nx = std::min(std::max(nx, 0.), 255.);
//            ny = std::min(std::max(ny, 0.), 255.);

            if(int(std::round(nx))==127 and int(std::round(ny))==127){max += 1;}

            histo.pixelAbsolute(int(std::round(nx)), int(std::round(ny))) += 1;
        }
    }


    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             P *= 1./double(max);
                         });

    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_reel.png");
}




void histo_2D_theo(double moy1, double moy2, double var1, double var2)
{
    ImageGrayd histo(256, 256, true);

    for(int x=0; x<256; x++)
    {
        for(int y=0; y<256; y++)
        {
            double X = x/256.;
            double Y = y/256.;

            double gaussian = gauss(moy1, moy2, var1, var2, X, Y);

            histo.pixelAbsolute(x, y) += gaussian;
        }
    }


    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_theo.png");
}




#endif //ASTEX_NUMERICAL_TCONTENT_H
