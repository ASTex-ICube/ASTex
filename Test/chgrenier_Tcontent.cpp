//
// Created by grenier on 30/06/23.
//

#include <ASTex/easy_io.h>
#include <ASTex/image_rgb.h>
using namespace ASTex;


struct color_info{
    ImageRGB8::PixelType couleur_;
    int compteur_;

    color_info(ImageRGB8::PixelType couleur){
        couleur_ = couleur;
        compteur_ = 1;
    }

    void incr()
    {
        compteur_ += 1;
    }

    bool operator==(const color_info& couleur2)
    {
        return couleur_ == couleur2.couleur_;
    }
};



bool is_in(std::vector<color_info>& vecteur, ImageRGB8::PixelType& couleur, int& id)
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




int main()
{
    std::string directory = "/home/grenier/Documents/ASTex_fork/results/T_analysis/";

    int img_size = 512;
    ImageRGB8 image_(img_size, img_size);
    image_.load(directory + "composition_out.png");

    std::vector<color_info> couleurs = {};


    // listage des couleur et compte de leur apparition
    image_.for_all_pixels([&] (typename ImageRGB8::PixelType& P, int x, int y)
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

    // suppression des "bavure" entre les couleurs (non nécessaire avec des donnée synthétique)
    int tot_vol = img_size* img_size;
    for(auto it = couleurs.begin(); it != couleurs.end();)
    {
        if((*it).compteur_ < 0.01*img_size* img_size){
            tot_vol -= (*it).compteur_ ;
            it = couleurs.erase(it);
        }
        else{
            it++;
        }
    }


    // écriture résultats
    for(auto it = couleurs.begin(); it != couleurs.end(); it++)
    {
        std::cout<<(*it).couleur_<<" : "<<(*it).compteur_/float(tot_vol)<<std::endl;
    }


    // analyse des noises
    ImageGrayu8 noise1_(img_size, img_size);
    noise1_.load(directory + "gabor1.png");

    double moy_n1 = 0.;
    double moy_n1_carre = 0.;

    noise1_.for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
                          {
                            double noise = noise1_.pixelAbsolute(x,y)/127. - 1.;
                            moy_n1 += noise;
                            moy_n1_carre += noise*noise;
                          });

    moy_n1 /= img_size*img_size;
    moy_n1_carre /= img_size*img_size;

    double var_n1 = moy_n1_carre - moy_n1*moy_n1;
    std::cout<<"statistiques noise 1 : "<<moy_n1<<", "<<var_n1<<std::endl;


    ImageGrayu8 noise2_(img_size, img_size);
    noise2_.load(directory + "gabor2.png");

    double moy_n2 = 0.;
    double moy_n2_carre = 0.;

    noise2_.for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
                           {
                               double noise = noise2_.pixelAbsolute(x,y)/127. - 1.;
                               moy_n2 += noise;
                               moy_n2_carre += noise*noise;
                           });

    moy_n2 /= img_size*img_size;
    moy_n2_carre /= img_size*img_size;

    double var_n2 = moy_n2_carre - moy_n2*moy_n2;
    std::cout<<"statistiques noise 2 : "<<moy_n2<<", "<<var_n2<<std::endl;


    return 0;
}