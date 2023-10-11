//
// Created by grenier on 22/08/23.
//

#include <ASTex/image_gray.h>

#include "ASTex/Noises/Gabor.h"

#include "ASTex/CCVT/point.h"
#include "ASTex/CCVT/sites.h"
#include "ASTex/CCVT/metric.h"
#include "ASTex/CCVT/optimizer.h"
#include "ASTex/CCVT/tools.h"

using namespace ASTex;

struct color_info{
    double couleur_;
    int compteur_;

    color_info(ImageGrayu8::PixelType couleur){
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


bool is_in(std::vector<color_info>& vecteur, ImageGrayu8::PixelType& couleur, int& id)
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
    // ---------------------------------------------------------------------------
    // création des bruits d'entrée
    // ---------------------------------------------------------------------------
    int resolution = 128; // pour le stockage
    int img_size = 1024; // nombre de pixel dans l'image

    float F_0_ = 0.1;//0.04; // fréquence
    float omega_0_ = 0.; // orientation (seulement dans le cas anisotrope, cf code gabor ligne 160)

    float number_of_impulses_per_kernel = 64.0;
    unsigned period = 128; // non utilisé

    float K_ = 1.0; //laisser à 1
    float a_ = 0.02; // taille des noyaux (a*a = 1/variance)

    unsigned random_offset_ = 954248632;
    unsigned seed_ = 1;


    // ---------------------------------------------------------------------------
    noise noise_x(K_,
                 a_,
                 F_0_,
                 omega_0_,
                 number_of_impulses_per_kernel,
                 period,
                 random_offset_,
                 seed_);

    noise noise_y(K_,
                  a_,
                  F_0_,
                  0.57,
                  number_of_impulses_per_kernel,
                  period,
                  random_offset_,
                  2.4*seed_);

//    float variance_Nx = noise_x.variance()/(3.4*3.4*noise_x.variance());
//    float variance_Ny = noise_y.variance()/(3.4*3.4*noise_y.variance());

    float scale_Nx = 3.4 * std::sqrt(noise_x.variance());
    float scale_Ny = 3.4 * std::sqrt(noise_y.variance());


    // ---------------------------------------------------------------------------
    ImageGrayu8 Check_noise_x = storing_noise(resolution, img_size, noise_x);
    ImageGrayu8 Check_noise_y = storing_noise(resolution, img_size, noise_y);

    Check_noise_x.save("/home/grenier/Documents/ASTex_fork/results/H_ccvt/gabor1.png");
    Check_noise_y.save("/home/grenier/Documents/ASTex_fork/results/H_ccvt/gabor2.png");






    // ---------------------------------------------------------------------------
    // création de la carte H d'entrée
    // ---------------------------------------------------------------------------
    int NUMBER_SITES      = 6; // nombre de graine de cellule
    int TORUS_SIZE        = 256; // taille de l'image (domaine de définition de H, définit sur -1,1, mappé sur 0, 256)
    int NUMBER_POINTS     = TORUS_SIZE;// * TORUS_SIZE * NUMBER_SITES * NUMBER_SITES; // nombre de points pour la contrainte de densité

    typedef MetricEuclidean2 Metric;
    typedef Point2 Point;

    // contrainte de densité
    Point::List points;
    nonconstant_density(points,
                        NUMBER_POINTS,
                        TORUS_SIZE,
                        noise_x.variance()/(scale_Nx*scale_Nx), // variance pour des valeurs de x sur -1, 1
                        noise_y.variance()/(scale_Ny*scale_Ny)); // densité non constante


    // graine voronoi
    int overallCapacity = static_cast<int>(points.size());
    Site<Point>::List sites;
    std::vector<double> colors{0., 42., 212, 170., 127., 85.};

    // donnée initiales
    // initialisation position aléatoire, capacité toutes égales
    for (int i = 0; i < NUMBER_SITES; ++i) {
        double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
        double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
        int capacity = overallCapacity / (NUMBER_SITES - i);
        overallCapacity -= capacity;
        sites.push_back(Site<Point2>(i, capacity, Point2(x, y)));
    }

    // initialisation du ccvt
    Metric metric(Point(TORUS_SIZE, TORUS_SIZE));
    Optimizer<Site<Point>, Point, Metric> optimizer;
    optimizer.initialize(sites, points, metric);

    // optimization
    int iteration = 0;
    bool stable = false;
    do {
//        printf("iteration %d...", ++iteration);
        stable = optimizer.optimize(true);
//        printf("done\n");
    } while (!stable);


    // ---------------------------------------------------------------------------
    const Site<Point>::Vector& Input_CM = optimizer.sites(); // liste des graines décrivant la fonction H

    // ---------------------------------------------------------------------------
    save_res_zone(Input_CM, metric, colors, TORUS_SIZE, "/home/grenier/Documents/ASTex_fork/results/H_ccvt/input_cm.png");




    // ---------------------------------------------------------------------------
    // création de l'exemple d'entrée par une composition (image discrète)
    // ---------------------------------------------------------------------------
    ImageGrayu8 Input_exemple(img_size, img_size);

    Input_exemple.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
           {
               double dist = 12400.;
               int color_id = -1;

               for (unsigned int i = 0; i < Input_CM.size(); ++i)
               {
                   // les noise sont définit sur 0,1
                   float X = x * (float(resolution)/img_size);
                   float Y = y * (float(resolution)/img_size);

                   // valeurs d'intensité sur 0, 1
                   float intensity_Nx = 0.5 + (0.5 * noise_x(X, Y) / scale_Nx);
                   float intensity_Ny = 0.5 + (0.5 * noise_y(X, Y) / scale_Ny);

                   // clamp
//                   intensity_Nx = std::min(std::max(intensity_Nx, 0.f), 1.f);
//                   intensity_Ny = std::min(std::max(intensity_Ny, 0.f), 1.f);

                   // distance
                   double new_dist = metric.distance(Point2(TORUS_SIZE*intensity_Nx,
                                                                TORUS_SIZE*intensity_Ny),
                                                     Input_CM[i].location);
                   if(new_dist < dist){
                       dist = new_dist;
                       color_id = Input_CM[i].id;
                   }
               }

               double color = colors.at(color_id); //255.* color_id / Input_CM.size();

               P = ImageGrayu8::PixelType(color);
           });

    // ---------------------------------------------------------------------------
    Input_exemple.save("/home/grenier/Documents/ASTex_fork/results/H_ccvt/composition_in.png");





    // ---------------------------------------------------------------------------
    // information sur l'image exemple (compte des pixel)
    // ---------------------------------------------------------------------------
    std::vector<color_info> couleurs_info = {};

    // listage des couleur et compte de leur apparition
    Input_exemple.for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
                          {
                              int id;
                              if(couleurs_info.empty())
                              {
                                  couleurs_info.push_back(color_info{P});
                              }
                              else
                              {
                                  bool presence = is_in(couleurs_info, P, id);
                                  if(not presence)
                                  {
                                      couleurs_info.push_back(color_info{P});
                                  }
                                  else
                                  {
                                      couleurs_info.at(id).incr();
                                  }
                              }
                          });

    // écriture résultats
    std::cout<<std::endl;
    for(auto it = couleurs_info.begin(); it != couleurs_info.end(); it++)
    {
        std::cout<<(*it).couleur_<<" : "<<(*it).compteur_<<std::endl;
    }
    std::cout<<std::endl;






//
//    // ---------------------------------------------------------------------------
//    // création de la carte H' de sortie
//    // ---------------------------------------------------------------------------
//    NUMBER_SITES      = couleurs_info.size(); // nombre de graine de cellule
//    NUMBER_POINTS     = img_size * img_size; // nombre de points pour la contrainte de densité
//
//    typedef MetricEuclidean2 Metric;
//    typedef Point2 Point;
//
//    // contrainte de densité
//    Point::List points_out;
//    nonconstant_density(points_out,
//                        NUMBER_POINTS,
//                        TORUS_SIZE,
//                        noise_x.variance()/(scale_Nx*scale_Nx), // variance pour des valeurs de x sur -1, 1
//                        noise_y.variance()/(scale_Ny*scale_Ny)); // densité non constante
//
//
//    // graine voronoi
//    int overallCapacity_out = static_cast<int>(points_out.size());
//    Site<Point>::List sites_out;
//    std::vector<double> couleur_graine;
//
//    // donnée initiales
//    // initialisation position aléatoire, capacité toutes égales
//    for (int i = 0; i < NUMBER_SITES; ++i) {
//        double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//        double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//
//        int capacity = couleurs_info.at(i).compteur_;
//        couleur_graine.push_back(couleurs_info.at(i).couleur_);
//
//        sites_out.push_back(Site<Point2>(i, capacity, Point2(x, y)));
//    }
//
//    // initialisation du ccvt
//    Metric metric_out(Point(TORUS_SIZE, TORUS_SIZE));
//    Optimizer<Site<Point>, Point, Metric> optimizer_out;
//    optimizer_out.initialize(sites_out, points_out, metric_out);
//
//    // optimization
//    iteration = 0;
//    stable = false;
//    do {
////        printf("iteration %d...", ++iteration);
//        stable = optimizer_out.optimize(true);
////        printf("done\n");
//    } while (!stable);
//
//
//    // ---------------------------------------------------------------------------
//    const Site<Point>::Vector& Output_CM = optimizer_out.sites(); // liste des graines décrivant la fonction H
//
//    // ---------------------------------------------------------------------------
//    save_res_zone(Output_CM, metric_out, couleur_graine, TORUS_SIZE, "/home/grenier/Documents/ASTex_fork/results/H_ccvt/output_cm.png");
//
//
//
//
//
//
//
//    // ---------------------------------------------------------------------------
//    // création de l'image de sortie par une composition (image discrète)
//    // ---------------------------------------------------------------------------
//    ImageGrayu8 Output_exemple(img_size, img_size);
//
//    Output_exemple.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
//                  {
//                      double dist = 12400.;
//                      int color_id = -1;
//
//                      for (unsigned int i = 0; i < Output_CM.size(); ++i)
//                      {
//                          // les noise sont définit sur 0,1
//                          float X = x * (float(resolution)/img_size);
//                          float Y = y * (float(resolution)/img_size);
//
//                          // on les force à valeur dans -1, 1
//                          float intensity_Nx =  noise_x(X, Y) / scale_Nx;
//                          float intensity_Ny =  noise_y(X, Y) / scale_Ny;
//
//                          double new_dist = metric.distance(Point2(0.5*(TORUS_SIZE*intensity_Nx + TORUS_SIZE),
//                                                                       0.5*(TORUS_SIZE*intensity_Ny + TORUS_SIZE)),
//                                                            Output_CM[i].location);
//                          if(new_dist < dist){
//                              dist = new_dist;
//                              color_id = Output_CM[i].id;
//                          }
//                      }
//
//                      double color = couleur_graine.at(color_id); //255.* color_id / Output_CM.size();
//
//                      P = ImageGrayu8::PixelType(color);
//                  });
//
//    // ---------------------------------------------------------------------------
//    Output_exemple.save("/home/grenier/Documents/ASTex_fork/results/H_ccvt/composition_out.png");
//
//
//
//
//    // ---------------------------------------------------------------------------
//    // information sur l'image de sortie (compte des pixel)
//    // ---------------------------------------------------------------------------
//    std::vector<color_info> couleurs_info_out = {};
//
//    // listage des couleur et compte de leur apparition
//    Output_exemple.for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
//                                 {
//                                     int id;
//                                     if(couleurs_info_out.empty())
//                                     {
//                                         couleurs_info_out.push_back(color_info{P});
//                                     }
//                                     else
//                                     {
//                                         bool presence = is_in(couleurs_info_out, P, id);
//                                         if(not presence)
//                                         {
//                                             couleurs_info_out.push_back(color_info{P});
//                                         }
//                                         else
//                                         {
//                                             couleurs_info_out.at(id).incr();
//                                         }
//                                     }
//                                 });
//
//    // écriture résultats
//    std::cout<<std::endl;
//    for(auto it = couleurs_info_out.begin(); it != couleurs_info_out.end(); it++)
//    {
//        std::cout<<(*it).couleur_<<" : "<<(*it).compteur_<<std::endl;
//    }
//    std::cout<<std::endl;





    return 0;
}

































