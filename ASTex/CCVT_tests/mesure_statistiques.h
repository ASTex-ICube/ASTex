//
// Created by grenier on 11/10/23.
//

#ifndef ASTEX_MESURE_STATISTIQUES_H
#define ASTEX_MESURE_STATISTIQUES_H

#include "tools.h"


double moyenne(ImageGrayd inputNoise)
{
    double sum = 0.;
    inputNoise.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y) // cm
                              {
                                  sum += P;
                              });
    return sum/(inputNoise.width()*inputNoise.height());
}


//double moyenne(noise inputNoise, int img_size, int resolution)
//{
//    double sum = 0.;
//    double scale_1 = 3.6 * std::sqrt(inputNoise.variance());
//
//    for(int x=0; x<img_size; x++){
//        for(int y=0; y<img_size; y++){
//            float X = x * (float(resolution)/float(img_size));
//            float Y = y * (float(resolution)/float(img_size));
//
//            double intensity = inputNoise(X, Y);// 0.5 + (0.5 * (inputNoise(X, Y) / scale_1));
////            double color = std::min(std::max(intensity, 0.), 1.); // clamping
//
//            sum += intensity;
//        }
//    }
//    return sum/double(img_size*img_size);
//}





double moyenne_carre(ImageGrayd inputNoise)
{
    double sum = 0.;
    inputNoise.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y) // cm
                              {
                                  sum += P*P;
                              });
    return sum/(inputNoise.width()*inputNoise.height());
}




double covariance(ImageGrayd inputNoiseX, int tauX, int tauY)
{
    double mX = 0.; //moyenne(inputNoiseX);
    double mY = 0.; //moyenne(inputNoiseY);

    for(int x=0; x<inputNoiseX.width()-tauX; x++)
    {
        for(int y=0; y<inputNoiseX.height()-tauY; y++)
        {
            mX += inputNoiseX.pixelAbsolute(x,y); // valeurs sur 0, 1
            mY += inputNoiseX.pixelAbsolute(x+tauX,y+tauY);
        }
    }
    mX = mX/((inputNoiseX.width()-tauX)*(inputNoiseX.height()-tauY)); // valeur sur 0, 1
    mY = mY/((inputNoiseX.width()-tauX)*(inputNoiseX.height()-tauY));
//    std::cout<<"moyenne image x : "<<mX<<", moyenne image y : "<<mY<<std::endl;


    double AC = 0.;
    for(int x=0; x<inputNoiseX.width()-tauX; x++)
    {
        for(int y=0; y<inputNoiseX.height()-tauY; y++)
        {
            double nx = inputNoiseX.pixelAbsolute(x,y); // valeurs sur 0, 1
            double ny = inputNoiseX.pixelAbsolute(x+tauX,y+tauY);
//            std::cout<<(nx - mX)<<std::endl;
//            std::cout<<(ny - mY)<<std::endl;

            AC += (nx - mX)*(ny - mY); // valeurs sur 0, 1
        }
    }
//    std::cout<<AC<<std::endl;
    return AC/((inputNoiseX.width()-tauX)*(inputNoiseX.height()-tauY));
}


ImageGrayd covariance_Fourier(ImageGrayd inputNoiseX)
{
    ImageSpectrald module;
    ImageSpectrald phase;
    ImageGrayd result;

    ImageGrayd noise = inputNoiseX;

    // valeurs centrées en 0
    double mean = moyenne(noise);
    noise.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                               {
                                   P -= mean;
                               });

    // transformée de fourier
    Fourier::fftForwardModulusAndPhase(noise, module, phase, false); // bool preserve_energy

    // affichage de test
//    IO::save_spectrum(module, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/module.png");

    // carré du module = PSD
    module.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                          {
                              P = (P*P)/(noise.width()*noise.height());
                          });

    // affichage de test
//    IO::save_spectrum(module, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/psd.png");
//    IO::save_phase(phase, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/phase.png");

    // phase nulle
    phase.initItk(phase.width(), phase.height(), true); // init_to_zero = true


    // transformée inverse = auto-covaraince
    Fourier::fftInverseModulusAndPhase(module, phase, result, false); // bool preserve_energy



//    // affichage de test
//    ImageGrayd result_tmp = result;
//    result_tmp.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
//                          {
//                              P *= 100.;
//                          });
//    IO::save(result_tmp, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/acov.png");

    return result; //result.pixelAbsolute(tauX, tauY);
}






// auto-covarience = auto-corrélation par hypothèse d'ergodicité
double autocorrelation_Fourier(ImageGrayd inputNoiseX, int tauX, int tauY)
{
    ImageGrayd result;
    ImageGrayd noise = inputNoiseX;

    // valeurs centrées en 0
    double mean = moyenne(noise);
    noise.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                          {
                              P -= mean;
                          });

    // auto-corrélation
    result.initItk(noise.width(), noise.height());
    Fourier::autoCorrelation_full_size(noise,result);

    //Normalisation
    result.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                          {
                              P = P/sqrt(noise.width()*noise.height());//sqrt(square_mean-mean+mean);
                          });

    // affichage de test
    IO::save(result, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/acorr.png");

    return result.pixelAbsolute(tauX, tauY);
}

// ---------------------------------------------------------------------------
double moyenne(std::vector<double> data){
    double sum = 0.;
    for(auto it = data.begin(); it != data.end(); it++)
    {
        sum += (*it);
    }
    return sum/data.size();
}


double minimum(std::vector<double> data){
    double mini = 10000.;
    for(auto it = data.begin(); it != data.end(); it++)
    {
        mini = std::min(mini,(*it));
    }
    return mini;
}

double maximum(std::vector<double> data){
    double maxi = -10000.;
    for(auto it = data.begin(); it != data.end(); it++)
    {
        maxi = std::max(maxi,(*it));
    }
    return maxi;
}

double medianne(std::vector<double> data){
    std::sort(data.begin(), data.end());
    return data.at(data.size()/2);
}

double quartil_1(std::vector<double> data){
    std::sort(data.begin(), data.end());
    return data.at(data.size()/4);
}

double quartil_3(std::vector<double> data){
    std::sort(data.begin(), data.end());
    return data.at(3*data.size()/4);
}


void affichage_data(const std::vector<double>& data, std::ofstream& fichier){
    double moy_ = moyenne(data);
    double med_ = medianne((data));
    double min_ = minimum(data);
    double max_ = maximum(data);
    double q1_ = quartil_1(data);
    double q3_ = quartil_3(data);

    fichier<<moy_<<" & "<<med_<<" & "<<min_<<" & "<<max_<<" & "<<q1_<<" & "<<q3_<<" & ";
}



// ---------------------------------------------------------------------------
ImageGrayd dF_dx(ImageGrayd image_)
{
    ImageGrayd Dx(image_.width(), image_.height(), true);

    Dx.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                      {
                          int xph = (x+1)%image_.width();
                          int xmh = (x-1)%image_.width();
                          P = 2.*(image_.pixelAbsolute(xph, y) - image_.pixelAbsolute(xmh, y)) + 0.5;
                      });
    IO::save(Dx, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/test_dx.png");
    return Dx;
}

ImageGrayd dF_dy(ImageGrayd image_)
{
    ImageGrayd Dy(image_.width(), image_.height(), true);

    Dy.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                      {
                          int yph = (y+1)%image_.height();
                          int ymh = (y-1)%image_.height();
                          P = 2.*(image_.pixelAbsolute(x, yph) - image_.pixelAbsolute(x, ymh)) + 0.5;
                      });
    IO::save(Dy, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/test_dy.png");
    return Dy;
}






// ---------------------------------------------------------------------------
ImageGrayd histo_2D(ImageGrayd noise1, ImageGrayd noise2, int histo_size)
{
//    int histo_size = 256;
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
        }
    }
    double norme = getMax(histo);

    // normalisation
//    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
//                         {
//                             P *= 1./norme;  // pour un max à 1
//                         });

    return histo;
}




ImageGrayd histo_2D_theo(double moy1, double moy2, double var1, double var2, int histo_size)
{
//    int histo_size = 256;
    ImageGrayd histo(histo_size, histo_size, true);

    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             double X = x;//double(histo_size-1);
                             double Y = y;//double(histo_size-1);

                             P = gauss(moy1, moy2, var1, var2, X, Y, double(histo_size-1));
                         });
    double norme = getMax(histo);

    // normalisation
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             P *= 1./norme;  // pour un max à 1
                         });

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
                             double X = x;//double(histo_size-1);
                             double Y = y;//double(histo_size-1);

                             double val_theo = gauss(moy1, moy2, var1, var2, X, Y, double(histo_size-1));
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





// ---------------------------------------------------------------------------
ImageGrayd histo_2D_AC(ImageGrayd noise1, int tauX, int tauY, int name)
{
    int histo_size = 256;
    double norme = 0;
    ImageGrayd histo(histo_size, histo_size, true);

    // énumération des paire de valeur de bruits
    for(int x=0; x<noise1.width()-tauX; x++)
    {
        for(int y=0; y<noise1.height()-tauY; y++)
        {
            double nx = noise1.pixelAbsolute(x, y)*(histo_size-1);
            double ny = noise1.pixelAbsolute(x+tauX, y+tauY)*(histo_size-1);

            histo.pixelAbsolute(int(std::round(nx)), int(std::round(ny))) += 1.;
            histo.pixelAbsolute(int(std::round(ny)), int(std::round(nx))) += 1.; // symétrie
        }
    }

    // calcul de l'aire sous la courbe
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             norme += P / (histo_size * histo_size);
                         });


    // estimation des moyennes
    double mX = 0.;
    double mY = 0.;
    double tot = 0.;
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             double X = x/double(histo_size-1);
                             double Y = y/double(histo_size-1);

                             mX += P*X; // x/255. : valeur sur 0, 1
                             mY += P*Y;

                             tot += P;
                         });
    mX = mX/tot;
    mY = mY/tot;
    std::cout<<"moyenne histo x : "<<mX<<", moyenne histo y : "<<mY<<std::endl;




    // estimation de l'auto-covariance
    double AC = 0.;
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             double X = x/double(histo_size-1);
                             double Y = y/double(histo_size-1);

                             AC += P*(X-mX)*(Y-mY); // x/255. : valeur sur 0, 1
//                             AC += P*(x-mX)*(y-mY); // x : valeur sur 0, 255
                         });
    AC = AC/tot;
    std::cout<<"auto-covariance histo : "<<AC<<std::endl;



    // normalisation pour une aire sous la courbe égale à 1
    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             P *= 1./double(norme);
                         });
    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_AC_" + std::to_string(name) + ".png");
    return histo;
}




ImageGrayd histo_2D_theo_vois(double moy, double var, double AC, int name)
{
    int histo_size = 256;
    ImageGrayd histo(histo_size, histo_size, true);

    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             double X = x/double(histo_size-1);
                             double Y = y/double(histo_size-1);

                             P = gaussAC_2D(moy, var, AC, X, Y);
                         });


    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_theo_vois_" + std::to_string(name) + ".png");
    return histo;
}



ImageGrayd histo_2D_dist_vois(ImageGrayd histo_reel, double moy, double var, double AC, int name)
{
    int histo_size = 256;
    ImageGrayd histo(histo_size, histo_size, true);

    double max_theo = -1.;
    double max_reel = -1.;
    double max_dist = -1.;

    histo.for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y)
                         {
                             double X = x/double(histo_size-1);
                             double Y = y/double(histo_size-1);

                             double val_theo = gaussAC_2D(moy, var, AC, X, Y) ;
                             double val_reel = histo_reel.pixelAbsolute(x,y);
                             double val_dist = val_theo - val_reel;

                             max_theo = std::max(max_theo,val_theo);
                             max_reel = std::max(max_reel,val_reel);
                             max_dist = std::max(max_dist,std::abs(val_dist));

                             P = val_dist;
                         });

    std::cout<<"bruit "<<std::to_string(name)<<" - val max histo réel : "<<max_reel<<", val max histo théo : "<<max_theo<<", val max distance : "<<max_dist<<std::endl;
    IO::save(histo, "/home/grenier/Documents/ASTex_fork/results/equ_CCVT/histo_noise_dist_vois_" + std::to_string(name) + ".png");
    return histo;
}


#endif //ASTEX_MESURE_STATISTIQUES_H
