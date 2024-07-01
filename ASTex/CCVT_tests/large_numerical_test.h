//
// Created by grenier on 09/02/24.
//

#ifndef ASTEX_LARGE_NUMERICAL_TEST_H
#define ASTEX_LARGE_NUMERICAL_TEST_H

#include "ASTex/CCVT_tests/numerical_Tcontent.h"
#include "ASTex/CCVT_tests/numerical_Tvois.h"
#include "ASTex/CCVT_tests/mesure_statistiques.h"
#include "ASTex/CCVT_tests/numerical_Treynold.h"

void mesure_E(ImageGrayd image_, std::vector<color_info>& presences, std::vector<color_vois>& voisinages, bool i_isEqual_j,
              int& nb_vois, int& nb_vois_diff){

    int id;
    bool presence;
    ImageGrayd::PixelType P2;

    std::vector<color_info> couleurs_pres = {};
    std::vector<color_vois> couleurs_vois = {};

    image_.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                          {
                              presence = is_in(couleurs_pres, P, id);

                              if(not presence)
                              {
                                  couleurs_pres.push_back(color_info{P});
                              }
                              else
                              {
                                  couleurs_pres.at(id).incr();
                              }



                              if(x+1 < image_.width()){
                                  P2 = image_.pixelAbsolute(x+1, y);
                                  update_Tcontent(couleurs_vois, P, P2, i_isEqual_j);

                                  if(P != P2){nb_vois_diff += 1;}
                                  nb_vois +=1;
                              }

                              if(x-1 >= 0){
                                  P2 = image_.pixelAbsolute(x-1, y);
                                  update_Tcontent(couleurs_vois, P, P2, i_isEqual_j);

                                  if(P != P2){nb_vois_diff += 1;}
                                  nb_vois +=1;
                              }

                              if(y+1 < image_.height()){
                                  P2 = image_.pixelAbsolute(x, y+1);
                                  update_Tcontent(couleurs_vois, P, P2, i_isEqual_j);

                                  if(P != P2){nb_vois_diff += 1;}
                                  nb_vois +=1;
                              }

                              if(y-1 >= 0){
                                  P2 = image_.pixelAbsolute(x, y-1);
                                  update_Tcontent(couleurs_vois, P, P2, i_isEqual_j);

                                  if(P != P2){nb_vois_diff += 1;}
                                  nb_vois +=1;
                              }
                          });

    presences = couleurs_pres;
    voisinages = couleurs_vois;
}





void mesure_H(ImageGrayd cm, std::vector<color_info>& presence_theo, std::vector<color_info>& presence_reel, std::vector<color_vois>& voisinage_theo, std::vector<color_vois>& voisinage_reel,
              double moy1, double moy2, double var1, double var2, ImageGrayd histo, bool i_isEqual_j){

    std::vector<color_info> couleurs_presence_theo = {};
    std::vector<color_info> couleurs_presence_reel = {};
    std::vector<color_vois> couleurs_voisinage_theo = {};
    std::vector<color_vois> couleurs_voisinage_reel = {};

    int cm_size = cm.height();
    int id;
    bool presence;
    ImageGrayd::PixelType P2;

    double X;
    double Y;

    // listage des couleur et compte de leur apparition
    cm.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                      {

                          // presence histo théorique
                          X = (x+0.5)/double(cm_size-1); // évaluation au centre des pixels
                          Y = (y+0.5)/double(cm_size-1);

                          double qtt_t =  gauss(moy1, moy2, var1, var2, X, Y);

                          presence = is_in(couleurs_presence_theo, P, id);
                          if(not presence)
                          {
                              couleurs_presence_theo.push_back(color_info{P, qtt_t});
                          }
                          else
                          {
                              couleurs_presence_theo.at(id).incr(qtt_t);
                          }


                          // présence histo réel
                          double qtt_r = histo.pixelAbsolute(x,y);
                          presence = is_in(couleurs_presence_reel, P, id);
                          if(not presence)
                          {
                              couleurs_presence_reel.push_back(color_info{P, qtt_r});
                          }
                          else
                          {
                              couleurs_presence_reel.at(id).incr(qtt_r);
                          }



                          // voisinages
                          if(x+1 < cm_size){
                              P2 = cm.pixelAbsolute(x+1, y);
                              update_Hcontent_theo(couleurs_voisinage_theo, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, 1, 0);
                              update_Hcontent_reel(couleurs_voisinage_reel, P, P2, histo, i_isEqual_j, x, y, 1, 0);

//                              if(P != P2){nb_vois_diff += 1;}
//                              nb_vois +=1;
                          }

                          if(x-1 >= 0){
                              P2 = cm.pixelAbsolute(x-1, y);
                              update_Hcontent_theo(couleurs_voisinage_theo, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, -1, 0);
                              update_Hcontent_reel(couleurs_voisinage_reel, P, P2, histo, i_isEqual_j, x, y, -1, 0);

//                              if(P != P2){nb_vois_diff += 1;}
//                              nb_vois +=1;
                          }

                          if(y+1 < cm_size){
                              P2 = cm.pixelAbsolute(x, y+1);
                              update_Hcontent_theo(couleurs_voisinage_theo, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, 0, 1);
                              update_Hcontent_reel(couleurs_voisinage_reel, P, P2, histo, i_isEqual_j, x, y, 0, 1);

//                              if(P != P2){nb_vois_diff += 1;}
//                              nb_vois +=1;
                          }

                          if(y-1 >= 0){
                              P2 = cm.pixelAbsolute(x, y-1);
                              update_Hcontent_theo(couleurs_voisinage_theo, P, P2, cm_size, i_isEqual_j, moy1, moy2, var1, var2, x, y, 0, -1);
                              update_Hcontent_reel(couleurs_voisinage_reel, P, P2, histo, i_isEqual_j, x, y, 0, -1);

//                              if(P != P2){nb_vois_diff += 1;}
//                              nb_vois +=1;
                          }


                      });

    presence_theo = couleurs_presence_theo;
    presence_reel = couleurs_presence_reel;
    voisinage_theo = couleurs_voisinage_theo;
    voisinage_reel = couleurs_voisinage_reel;

}


void large_test(int nb_iteration){
    std::ofstream fichier(TEMPO_PATH+"results/equ_CCVT/test.txt", std::ios::out | std::ios::trunc);


    int resolution = 512;
    int img_size = 4096;//4096;//2048; // nombre de pixel dans l'image
    int cm_size = 256;//256;//512;

    fichier<<"itérations : "<<nb_iteration<<", résolution : "<<resolution<<", img size : "<<img_size<<", cm size : "<<cm_size<<std::endl;

    float number_of_impulses_per_kernel = 64.;// 64.0;
    unsigned period = 128; // non utilisé

    float K_ = 1.0; //laisser à 1
    float a_ = 0.01; // taille des noyaux (a*a = 1/variance)0.02

    // ---------------------------------------------------------------------------
    ImageGrayd cm_(cm_size, cm_size);
    ImageGrayd histo_N1_N2(cm_size, cm_size);
    ImageGrayd image_1(img_size, img_size);
    ImageGrayd image_2(img_size, img_size);
    ImageGrayd res_composition(img_size, img_size);
//    ImageGrayd Fourier_1(img_size, img_size);
//    ImageGrayd Fourier_2(img_size, img_size);


    std::vector<color_info> couleurs_E;
    std::vector<color_info> couleurs_H;
    std::vector<color_info> couleurs_H_r;

    std::vector<color_vois> couleur_vois_E;
    std::vector<color_vois> couleur_vois_H;
    std::vector<color_vois> couleur_vois_H_r;

    // ---------------------------------------------------------------------------
//    std::vector<Graine> H_seeds{Graine(0.15, 0.25, 0.1),
//                                Graine(0.81, 0.65, 0.)};
//    std::vector<double> H_color{0.2, 0.0};

    std::vector<Graine> H_seeds{Graine(0.15, 0.25, 0.1),
                                Graine(0.65, 0.28, 0.),
                                Graine(0.81, 0.65, 0.),
                                Graine(0.52, 0.91 ,0.),
                                Graine(0.12, 0.72, 0.1),
                                Graine(0.45, 0.48, 0.)};
    std::vector<double> H_color{0.2, 0.6, 0.0, 0.4, 0.8, 1.0};


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


    // ---------------------------------------------------------------------------
    std::sort(H_color.begin(), H_color.end());

    std::vector<color_info_large> info_presence;
    for(auto it = H_color.begin(); it != H_color.end(); it++)
    {
        info_presence.push_back(color_info_large(*it));
    }


    std::vector<color_voisinage_large> info_voisinage;
    for(auto it1 = H_color.begin(); it1 != H_color.end(); it1++)
    {
        for(auto it2 = it1; it2 != H_color.end(); it2++)
        {
            info_voisinage.push_back(color_voisinage_large(*it1, *it2));
        }
    }



    // ---------------------------------------------------------------------------
    for(int i=0; i<nb_iteration; i++){
//        std::cout<<nb_iteration-i<<std::endl;
        if(i%(nb_iteration/10)==0){
            std::cout<<i<<std::endl;
        }

        unsigned random_offset_ = std::rand();//954248632;
        unsigned seed_ = 1 + std::rand() / ((RAND_MAX + 1u) / (i+1));//4;

        float F_1 = 0.01*(std::rand() % 10)+0.01; // frequ sur [0.01, 0.11]
        float F_2 = 0.01*(std::rand() % 10)+0.01;

        float omega_1 = 3.14*0.01*(std::rand() % 100 + 1);
        float omega_2 = 3.14*0.01*(std::rand() % 100 + 1);

        // ---------------------------------------------------------------------------
        noise noise_1(K_,
                      a_,
                      F_1,
                      omega_1,
                      number_of_impulses_per_kernel,
                      period,
                      random_offset_,
                      seed_);
        image_1 = storing_noise_d(resolution, img_size, noise_1);


        noise noise_2(K_,
                      a_,
                      F_2,
                      omega_2,
                      number_of_impulses_per_kernel,
                      period,
                      random_offset_-4,
                      seed_+12);
        image_2 = storing_noise_d(resolution, img_size, noise_2);




        // ---------------------------------------------------------------------------
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

        double mu_1_carre = moyenne_carre(image_1);
        double mu_2_carre = moyenne_carre(image_2);

        double var_1 = mu_1_carre - mu_1*mu_1; // sigma
        double var_2 = mu_2_carre - mu_2*mu_2; // sigma'

        histo_N1_N2 = histo_2D(image_1, image_2, cm_size); // histogramme

        int nb_voisinage = 0;
        int nb_voisinage_diff = 0;


        // ---------------------------------------------------------------------------
        mesure_E(res_composition, couleurs_E, couleur_vois_E, true, nb_voisinage, nb_voisinage_diff);
        mesure_H(cm_, couleurs_H, couleurs_H_r, couleur_vois_H, couleur_vois_H_r, mu_1, mu_2, var_1, var_2, histo_N1_N2, false);



        // ---------------------------------------------------------------------------
//        std::cout<<"mesure présence"<<std::endl;



        // proportion présence dans E
        std::sort(couleurs_E.begin(), couleurs_E.end());

        double tot_pixel_E = 0.;
        for(auto it = couleurs_E.begin(); it != couleurs_E.end(); it++)
        {
            tot_pixel_E += (*it).compteur_;
        }

        // capacité présence dans H histo théo
        std::sort(couleurs_H.begin(), couleurs_H.end());

        double tot_pixel_H = 0.;
        for(auto it = couleurs_H.begin(); it != couleurs_H.end(); it++)
        {
            tot_pixel_H += (*it).compteurD_;
        }


        // capacité présence dans H histo réel
        std::sort(couleurs_H_r.begin(), couleurs_H_r.end());

        double tot_pixel_H_r = 0.;
        for(auto it = couleurs_H_r.begin(); it != couleurs_H_r.end(); it++)
        {
            tot_pixel_H_r += (*it).compteurD_;
        }


        // stockage
        for(int i=0; i<info_presence.size(); i++){
            double couleur = info_presence.at(i).couleur_;
            double proportion_E = 0.;
            double proportion_H = 0.;
            double proportion_H_r = 0.;

            int id_E;
            int id_H;
            int id_H_r;


            if(is_in(couleurs_E, couleur, id_E)){
                proportion_E = couleurs_E.at(id_E).compteur_/tot_pixel_E;
            }
            if(is_in(couleurs_H, couleur, id_H)){
                proportion_H = couleurs_H.at(id_H).compteurD_/tot_pixel_H;
            }
            if(is_in(couleurs_H_r, couleur, id_H_r)){
                proportion_H_r = couleurs_H_r.at(id_H_r).compteurD_/tot_pixel_H_r;
            }

            double err_t = 100.*std::abs(proportion_E - proportion_H)/proportion_E;
            double err_r = 100.*std::abs(proportion_E - proportion_H_r)/proportion_E;

            double dist_t = proportion_E - proportion_H;
            double dist_r = proportion_E - proportion_H_r;


            info_presence.at(i).proportion_E.push_back(proportion_E);
            info_presence.at(i).proportion_H_r.push_back(proportion_H_r);
            info_presence.at(i).proportion_H_t.push_back(proportion_H);

            info_presence.at(i).error_r.push_back(err_r);
            info_presence.at(i).error_t.push_back(err_t);

            info_presence.at(i).distance_r.push_back(dist_r);
            info_presence.at(i).distance_t.push_back(dist_t);

//            std::cout<<couleur<<" & "<< proportion_E<<" & "<< proportion_H<<std::endl;//" & "<< 100.*err<<" \\\\"<<std::endl;
        }
//        std::cout<<std::endl;





        // ---------------------------------------------------------------------------
//        std::cout<<"mesure voisinage"<<std::endl;

        // proportion voisinage dans E
//        couleur_vois_E = Tcontent_vois_test(res_composition, true);
        std::sort(couleur_vois_E.begin(), couleur_vois_E.end());

        // total tous voisinage
//        double tot_E = 0.;
//        for(auto it = couleur_vois_E.begin(); it != couleur_vois_E.end(); it++)
//        {
//            tot_E += (*it).compteur_;
//        }

//        for(auto it = couleur_vois_E.begin(); it != couleur_vois_E.end(); it++)
//        {
//            std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") & "<<(*it).compteur_<<std::endl;
//
//        }
//        std::cout<<std::endl;





        // capacité voisinage dans H histo théorique
//        couleur_vois_H = Hcontent_vois_test(cm_, mu_1, mu_2, var_1, var_2, histo_N1_N2, false);
        std::sort(couleur_vois_H.begin(), couleur_vois_H.end());

        // total tous voisinage
        double tot_H = 0.;
        for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
        {
            tot_H += (*it).compteurD_;
        }

//        for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
//        {
//            std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") & "<<(*it).compteurD_<<std::endl;
//        }
//        std::cout<<std::endl;




        // capacité voisinage dans H histo réel
//        couleur_vois_H = Hcontent_vois_test(cm_, mu_1, mu_2, var_1, var_2, histo_N1_N2, false);
        std::sort(couleur_vois_H_r.begin(), couleur_vois_H_r.end());

        // total tous voisinage
        double tot_H_r = 0.;
        for(auto it = couleur_vois_H_r.begin(); it != couleur_vois_H_r.end(); it++)
        {
            tot_H_r += (*it).compteurD_;
        }

//        for(auto it = couleur_vois_H.begin(); it != couleur_vois_H.end(); it++)
//        {
//            std::cout<<"("<<(*it).couleur1_<<", "<<(*it).couleur2_<<") & "<<(*it).compteurD_<<std::endl;
//        }
//        std::cout<<std::endl;



        // affichage
//        for(int i=0; i<couleur_vois_E.size(); i++){
//            double couleur_1 = couleur_vois_E.at(i).couleur1_;
//            double couleur_2 = couleur_vois_E.at(i).couleur2_;
//
//            double proportion_E = couleur_vois_E.at(i).compteur_/tot_E;
//            double proportion_H = couleur_vois_H.at(i).compteurD_/tot_H;
//
//            double err = std::abs(proportion_E - proportion_H)/proportion_E;
//
//            std::cout<<"("<<couleur_1<<", "<<couleur_2<<") & "<<proportion_E<<" & "<<proportion_H<<std::endl;//" & "<<100.*err<<" \\\\"<<std::endl;
//
//        }
//        std::cout<<std::endl;



        // stockage
        for(int i=0; i<info_voisinage.size(); i++){
            double couleur1 = info_voisinage.at(i).couleur_1;
            double couleur2 = info_voisinage.at(i).couleur_2;

            double proportion_E = 0.;
            double proportion_H = 0.;
            double proportion_H_r = 0.;

            int id_E;
            int id_H;
            int id_H_r;

            bool is_in_E = is_in(couleur_vois_E, couleur1, couleur2, id_E);
            bool is_in_H = is_in(couleur_vois_H, couleur1, couleur2, id_H);
            bool is_in_H_r = is_in(couleur_vois_H_r, couleur1, couleur2, id_H_r);

            if(is_in_E){
                if(couleur1 == couleur2){
                    proportion_E = couleur_vois_E.at(id_E).compteur_/double(nb_voisinage);
                }
                else{
                    proportion_E = couleur_vois_E.at(id_E).compteur_/(2.*double(nb_voisinage));
                }
            }

            if(is_in_H){
                if(couleur1 == couleur2){
                    proportion_H = -1.;
                }
                else{
                    proportion_H = (couleur_vois_H.at(id_H).compteurD_/tot_H) * (double(nb_voisinage_diff)/(2.*double(nb_voisinage)));
                }
            }

            if(is_in_H_r){
                if(couleur1 == couleur2){
                    proportion_H_r = -1.;
                }
                else{
                    proportion_H_r = (couleur_vois_H_r.at(id_H_r).compteurD_/tot_H_r) * (double(nb_voisinage_diff)/(2.*double(nb_voisinage)));
                }
            }




            double err_t = 100.*std::abs(proportion_E - proportion_H)/proportion_E;
            double err_r = 100.*std::abs(proportion_E - proportion_H_r)/proportion_E;

            double dist_t = proportion_E - proportion_H;
            double dist_r = proportion_E - proportion_H_r;


            info_voisinage.at(i).proportion_E.push_back(proportion_E);
            info_voisinage.at(i).proportion_H_r.push_back(proportion_H_r);
            info_voisinage.at(i).proportion_H_t.push_back(proportion_H);

            info_voisinage.at(i).error_r.push_back(err_r);
            info_voisinage.at(i).error_t.push_back(err_t);

            info_voisinage.at(i).distance_r.push_back(dist_r);
            info_voisinage.at(i).distance_t.push_back(dist_t);

//            std::cout<<couleur<<" & "<< proportion_E<<" & "<< proportion_H<<std::endl;//" & "<< 100.*err<<" \\\\"<<std::endl;
        }



    }


    fichier<<"couleur & ";
    fichier<<"moy E & med E & min E & max E & Q1 E & Q3 E & ";

    fichier<<"moy Ht & med Ht & min Ht & max Ht & Q1 Ht & Q3 Ht & ";
    fichier<<"moy Hr & med Hr & min Hr & max Hr & Q1 Hr & Q3 Hr & ";

    fichier<<"moy Err t & med Err t & min Err t & max Err t & Q1 Err t & Q3 Err t & ";
    fichier<<"moy Err r & med Err r & min Err r & max Err r & Q1 Err r & Q3 Err r & ";

    fichier<<"moy dist t & med dist t & min dist t & max dist t & Q1 dist t & Q3 dist t & ";
    fichier<<"moy dist r & med dist r & min dist r & max dist r & Q1 dist r & Q3 dist r & ";
    fichier<<std::endl;

    for(auto it = info_presence.begin(); it != info_presence.end(); it++)
    {
        fichier<<(*it).couleur_<<" & ";

        affichage_data((*it).proportion_E, fichier);
        affichage_data((*it).proportion_H_t, fichier);
        affichage_data((*it).proportion_H_r, fichier);
        affichage_data((*it).error_t, fichier);
        affichage_data((*it).error_r, fichier);
        affichage_data((*it).distance_t, fichier);
        affichage_data((*it).distance_r, fichier);

        fichier<<std::endl;


    }

    fichier<<std::endl;
    fichier<<"couleur & ";
    fichier<<"moy E & med E & min E & max E & Q1 E & Q3 E & ";

    fichier<<"moy Ht & med Ht & min Ht & max Ht & Q1 Ht & Q3 Ht & ";
    fichier<<"moy Hr & med Hr & min Hr & max Hr & Q1 Hr & Q3 Hr & ";

    fichier<<"moy Err t & med Err t & min Err t & max Err t & Q1 Err t & Q3 Err t & ";
    fichier<<"moy Err r & med Err r & min Err r & max Err r & Q1 Err r & Q3 Err r & ";

    fichier<<"moy dist t & med dist t & min dist t & max dist t & Q1 dist t & Q3 dist t & ";
    fichier<<"moy dist r & med dist r & min dist r & max dist r & Q1 dist r & Q3 dist r & ";
    fichier<<std::endl;

    for(auto it = info_voisinage.begin(); it != info_voisinage.end(); it++)
    {
        fichier<<(*it).couleur_1<<";"<<(*it).couleur_2<<" & ";

        affichage_data((*it).proportion_E, fichier);
        affichage_data((*it).proportion_H_t, fichier);
        affichage_data((*it).proportion_H_r, fichier);
        affichage_data((*it).error_t, fichier);
        affichage_data((*it).error_r, fichier);
        affichage_data((*it).distance_t, fichier);
        affichage_data((*it).distance_r, fichier);

        fichier<<std::endl;


    }

    fichier.close();
}



#endif //ASTEX_LARGE_NUMERICAL_TEST_H
