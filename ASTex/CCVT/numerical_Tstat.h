//
// Created by grenier on 11/10/23.
//

#ifndef ASTEX_NUMERICAL_TSTAT_H
#define ASTEX_NUMERICAL_TSTAT_H




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



#endif //ASTEX_NUMERICAL_TSTAT_H
