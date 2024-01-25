//
// Created by grenier on 17/11/23.
//

#ifndef ASTEX_REYNOLD_H
#define ASTEX_REYNOLD_H

#include "ASTex/CCVT/numerical_Tcontent.h"
#include "ASTex/CCVT/numerical_Tvois.h"
#include "ASTex/CCVT/numerical_Tstat.h"

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









void update_Tcontent(std::vector<color_vois>& couleurs,
                     ImageGrayd::PixelType P1, ImageGrayd::PixelType P2,
                     bool i_isEqual_j)
{
    int id;
    bool presence;

//    if(P1 != P2 or i_isEqual_j){
//        presence = is_in(couleurs, P1, P2, id);
//        if(not presence)
//        {
//            couleurs.push_back(color_vois{P1, P2});
//        }
//        else
//        {
//            couleurs.at(id).incr();
//        }
//    }

//    if(P1 != P2){
//        presence = is_in(couleurs, P1, P2, id);
//        if(not presence)
//        {
//            couleurs.push_back(color_vois{P1, P2});
//        }
//        else
//        {
//            couleurs.at(id).incr();
//        }
//
//        presence = is_in(couleurs, P2, P1, id);
//        if(not presence)
//        {
//            couleurs.push_back(color_vois{P2, P1});
//        }
//        else
//        {
//            couleurs.at(id).incr();
//        }
//    }

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

//    else{
//        std::cout<<"pb E"<<std::endl;
//    }

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















void update_Hcontent(std::vector<color_vois>& couleurs, ImageGrayd::PixelType P1, ImageGrayd::PixelType P2,
                     ImageGrayd histo_, bool i_isEqual_j,
                     double moy1, double moy2, double var1, double var2,
                     int x, int y,  int dx, int dy) {

    int id;
    bool presence;
    double X;// = (x+0.5)/double(cm_size-1); // évaluation au centre des pixels
    double Y;// = (y+0.5)/double(cm_size-1);
    double qtt;// = gauss(moy1, moy2, var1, var2, X, Y);
    int cm_size = histo_.height();


//    if (P1 != P2 or i_isEqual_j) {
//        X = (x + 0.5 + double(dx) / 2.) / double(cm_size - 1); // évaluation au centre des bord des pixels
//        Y = (y + 0.5 + double(dy) / 2.) / double(cm_size - 1);
//
////        qtt = gauss(moy1, moy2, var1, var2, X, Y);
//        qtt = (histo_.pixelAbsolute(x, y) + histo_.pixelAbsolute(x + dx, y + dy));
//
//
//        presence = is_in(couleurs, P1, P2, id);
//        if (not presence) {
//            couleurs.push_back(color_vois{P1, P2, qtt});
//        } else {
//            couleurs.at(id).incr(qtt);
//        }
//
//    }



    X = (x + 0.5 + double(dx) / 2.) / double(cm_size - 1); // évaluation au centre des bord des pixels
    Y = (y + 0.5 + double(dy) / 2.) / double(cm_size - 1);
    qtt = gauss(moy1, moy2, var1, var2, X, Y);
//    qtt = 0.5*(histo_.pixelAbsolute(x, y) + histo_.pixelAbsolute(x + dx, y + dy));


//    if (P1 != P2) {
//
//        presence = is_in(couleurs, P1, P2, id);
//        if (not presence) {
//            couleurs.push_back(color_vois{P1, P2, qtt});
//        } else {
//            couleurs.at(id).incr(qtt);
//        }
//
//        presence = is_in(couleurs, P2, P1, id);
//        if (not presence) {
//            couleurs.push_back(color_vois{P2, P1, qtt});
//        } else {
//            couleurs.at(id).incr(qtt);
//        }
//    }


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

//    else{
//        std::cout<<"pb H"<<std::endl;
//    }
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
                                  update_Hcontent(couleurs, P, P2, histo_, i_isEqual_j, moy1, moy2, var1, var2, x, y, 1, 0);

//                                  if(y+1 < cm_size){
//                                      P2 = cm_.pixelAbsolute(x+1, y+1);
//                                      update_Hcontent(couleurs, P, P2, histo_, i_isEqual_j, moy1, moy2, var1, var2, x, y, 1, 1);
//                                  }
//
//                                  if(y-1 >= 0){
//                                      P2 = cm_.pixelAbsolute(x+1, y-1);
//                                      update_Hcontent(couleurs, P, P2, histo_, i_isEqual_j, moy1, moy2, var1, var2, x, y, 1, -1);
//                                  }
                              }

                              if(x-1 >= 0){
                                  P2 = cm_.pixelAbsolute(x-1, y);
                                  update_Hcontent(couleurs, P, P2, histo_, i_isEqual_j, moy1, moy2, var1, var2, x, y, -1, 0);

//                                  if(y+1 < cm_size){
//                                      P2 = cm_.pixelAbsolute(x-1, y+1);
//                                      update_Hcontent(couleurs, P, P2, histo_, i_isEqual_j, moy1, moy2, var1, var2, x, y, -1, 1);
//                                  }
//
//                                  if(y-1 >= 0){
//                                      P2 = cm_.pixelAbsolute(x-1, y-1);
//                                      update_Hcontent(couleurs, P, P2, histo_, i_isEqual_j, moy1, moy2, var1, var2, x, y, -1, -1);
//                                  }
                              }

                              if(y+1 < cm_size){
                                  P2 = cm_.pixelAbsolute(x, y+1);
                                  update_Hcontent(couleurs, P, P2, histo_, i_isEqual_j, moy1, moy2, var1, var2, x, y, 0, 1);
                              }

                              if(y-1 >= 0){
                                  P2 = cm_.pixelAbsolute(x, y-1);
                                  update_Hcontent(couleurs, P, P2, histo_, i_isEqual_j, moy1, moy2, var1, var2, x, y, 0, -1);
                              }

                          });

    return couleurs;
}







//    // total par couleur1_
//    double last_col = couleur_vois_E.at(0).couleur1_;
//    std::vector<std::array<double, 2>> tot_vois_E{std::array{last_col, 0.}}; // couleur, total
//
//    for(auto it = couleur_vois_E.begin(); it != couleur_vois_E.end(); it++)
//    {
//        if((*it).couleur1_ == last_col){
//            tot_vois_E.at(tot_vois_E.size()-1)[1] += (*it).compteur_;
//        }
//        else{
//            tot_vois_E.push_back(std::array{(*it).couleur1_, double((*it).compteur_)});
//            last_col = (*it).couleur1_;
//        }
//    }
//
////    for(auto it = tot_vois_E.begin(); it != tot_vois_E.end(); it++)
////    {
////        std::cout<<(*it)[0]<<" "<<(*it)[1]<<std::endl;
////    }
//
//    int last_col_id = 0;
//    for(auto it = couleur_vois_E.begin(); it != couleur_vois_E.end(); it++)
//    {
//        if((*it).couleur1_ == tot_vois_E.at(last_col_id)[0]){
//            std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteur_/tot_vois_E.at(last_col_id)[1]<<std::endl;
//        }
//        else{
//            last_col_id += 1;
//            std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") : "<<(*it).compteur_/tot_vois_E.at(last_col_id)[1]<<std::endl;
//        }
//    }



void large_test(int nb_iteration){
    int resolution = 512;
    int img_size = 4096;//4096;//2048; // nombre de pixel dans l'image
    int cm_size = 256;//256;//512;

    float number_of_impulses_per_kernel = 64.;// 64.0;
    unsigned period = 128; // non utilisé

    float K_ = 1.0; //laisser à 1
    float a_ = 0.01; // taille des noyaux (a*a = 1/variance)0.02


    // ---------------------------------------------------------------------------
    std::vector<Graine> H_seeds{Graine(0.15, 0.25, 0.1),
                                Graine(0.81, 0.65, 0.)};
    std::vector<double> H_color{0.2, 0.0};

    ImageGrayd cm_(cm_size, cm_size);
    cm_.parallel_for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                                {
                                    double dist = 100.;
                                    int id_seed = -1;

                                    for(int i=0; i<H_seeds.size(); i++){
                                        double X = double(x)/(cm_size-1);
                                        double Y = double(y)/(cm_size-1);

                                        double new_dist = (H_seeds.at(i).x_ - X)*(H_seeds.at(i).x_ - X) +
                                                          (H_seeds.at(i).y_ - Y)*(H_seeds.at(i).y_ - Y) -
                                                          H_seeds.at(i).weight_;
                                        if(new_dist<dist){
                                            dist = new_dist;
                                            id_seed = i;
                                        }
                                        P = H_color.at(id_seed);
                                    }
                                });

    ImageGrayd histo_N1_N2(cm_size, cm_size);


    for(int i=0; i<nb_iteration; i++){

        unsigned random_offset_ = std::rand();//954248632;
        unsigned seed_ = 1 + std::rand() / ((RAND_MAX + 1u) / (i+1));//4;

        float F_1 = 0.01*(std::rand() % 10)+0.01; // frequ sur [0.01, 0.11]
        float F_2 = 0.01*(std::rand() % 10)+0.01;

        float omega_1 = 3.14*0.01*(std::rand() % 100 + 1);
        float omega_2 = 3.14*0.01*(std::rand() % 100 + 1);

        // ---------------------------------------------------------------------------
        noise noise_1(K_, // isotrope
                      a_,
                      F_1,//0.012, //F_0_,
                      omega_1, //omega_0_,
                      number_of_impulses_per_kernel,
                      period,
                      random_offset_,
                      seed_);
        ImageGrayd image_1 = storing_noise_d(resolution, img_size, noise_1);


        noise noise_2(K_, // anisotrope
                      a_,
                      F_2,//0.02, //F_0_,
                      omega_2, //omega_0_,
                      number_of_impulses_per_kernel,
                      period,
                      random_offset_-4,
                      seed_+12);
        ImageGrayd image_2 = storing_noise_d(resolution, img_size, noise_2);



        // ---------------------------------------------------------------------------
        ImageGrayd res_composition(img_size, img_size);

        res_composition.parallel_for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                                                {
                                                    // valeur sur [0, 1]
                                                    double n1 = image_1.pixelAbsolute(x,y);
                                                    double n2 = image_2.pixelAbsolute(x,y);

                                                    int id_n1 = int(std::round(n1*(cm_size-1)));
                                                    int id_n2 = int(std::round(n2*(cm_size-1)));

                                                    P = cm_.pixelAbsolute(id_n1, id_n2);
                                                });


        // ---------------------------------------------------------------------------
        double mu_1 = moyenne(image_1); // mu
        double mu_2 = moyenne(image_2); // mu'

        ImageGrayd Fourier_1 = covariance_Fourier(image_1);
        ImageGrayd Fourier_2 = covariance_Fourier(image_2);

        double var_1 = Fourier_1.pixelAbsolute(0, 0); // sigma
        double var_2 = Fourier_2.pixelAbsolute(0, 0); // sigma'





//        // ---------------------------------------------------------------------------
//        std::cout<<"mesure présence"<<std::endl;
//
//        // proportion présence dans E
//        std::vector<color_info> couleurs_E = Tcontent(res_composition);
//        std::sort(couleurs_E.begin(), couleurs_E.end());
//
//        double tot_pixel_E = 0.;
//        for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
//        {
//            tot_pixel_E += (*it).compteur_;
//        }
//
//        // capacité présence dans H
//        std::vector<color_info> couleurs_H = cell_capacity(cm_, mu_1, mu_2, var_1, var_2);
//        std::sort(couleurs_H.begin(), couleurs_H.end());
//
//        double tot_pixel_H = 0.;
//        for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
//        {
//            tot_pixel_H += (*it).compteurD_;
//        }
//
//
//        // affichage
//        for(int i=0; i<couleurs_E.size(); i++){
//            double couleur = couleurs_E.at(i).couleur_;
//            double proportion_E = couleurs_E.at(i).compteur_/tot_pixel_E;
//            double proportion_H = couleurs_H.at(i).compteurD_/tot_pixel_H;
//            double err = std::abs(proportion_E - proportion_H)/proportion_E;
//
//            std::cout<<couleur<<" & "<< proportion_E<<" & "<< proportion_H<<" & "<< 100.*err<<" \\\\"<<std::endl;
//        }
//        std::cout<<std::endl;





        // ---------------------------------------------------------------------------
//        std::cout<<"mesure voisinage"<<std::endl;

        // proportion voisinage dans E
        std::vector<color_vois> couleur_vois_E = Tcontent_vois_test(res_composition, true);
        std::sort(couleur_vois_E.begin(), couleur_vois_E.end());

        // total tous voisinage
        double tot_E = 0.;
        for(auto it = couleur_vois_E.begin(); it != couleur_vois_E.end(); it++)
        {
            tot_E += (*it).compteur_;
        }


        // capacité voisinage dans H
        std::vector<color_vois> couleur_vois_H = Hcontent_vois_test(cm_, mu_1, mu_2, var_1, var_2, histo_N1_N2, true);
        std::sort(couleur_vois_H.begin(), couleur_vois_H.end());

        // total tous voisinage
        double tot_H = 0.;
        for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
        {
            tot_H += (*it).compteurD_;
        }



        // affichage
        for(int i=0; i<couleur_vois_E.size(); i++){
            double couleur_1 = couleur_vois_E.at(i).couleur1_;
            double couleur_2 = couleur_vois_E.at(i).couleur2_;

            double proportion_E = couleur_vois_E.at(i).compteur_/tot_E;
            double proportion_H = couleur_vois_H.at(i).compteurD_/tot_H;

            double err = std::abs(proportion_E - proportion_H)/proportion_E;

            std::cout<<"("<<couleur_1<<", "<<couleur_2<<") & "<<proportion_E<<" & "<<proportion_H<<" & "<<100.*err<<" \\\\"<<std::endl;

        }
//        std::cout<<std::endl;



    }
}



#endif //ASTEX_REYNOLD_H
