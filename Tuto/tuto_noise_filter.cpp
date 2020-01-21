#include <cmath>
#include <random>
#include <ASTex/image_gray.h>
#include <ASTex/easy_io.h>

using namespace ASTex;

int n;
double *phases,*frequences,*orientations;

double clamp(double &x,const double &a,const double &b)
{
    if ( x < a)
        x = a;
    if( x > b)
        x = b;
    return x;
}

double noise(const double &x,const double &y)
{
    double sum_cosines = 0.0;
    double p,f,o;
    for (int i = 0; i < n; ++i) {
        p = phases[i];
        f = frequences[i];
        o = orientations[i];
        sum_cosines += std::cos(p + f * std::cos(o) * x + f * std::sin(o) * y );
    }

    /*sum_cosines *= std::sqrt(1.0 / (3.0 * n));
    sum_cosines += 0.5;
    sum_cosines = clamp(sum_cosines,0.0,1.0);*/

    return sum_cosines;
}

int main()
{
    n = 32;
    phases = new double[n];
    frequences = new double[n];
    orientations = new double[n];
    const double Pi = 4.0 * std::atan(1.0);

    std::random_device rd;
    std::mt19937 gen(rd());
    //[0;2Pi[
    std::uniform_real_distribution<double> phases_dis(0.0,2.0 * Pi);
    //[10;20[
    std::uniform_real_distribution<double> frequences_dis(10.0,20.0);
    //[0;Pi[
    std::uniform_real_distribution<double> orientations_dis(0.0,Pi);

    for (int i = 0; i < n; ++i) {
        phases[i] = phases_dis(gen);
        frequences[i] = frequences_dis(gen);
        orientations[i] = orientations_dis(gen);
    }

    ImageGrayd im(512,512);

    im.parallel_for_all_pixels([](double &p,int x,int y)
    {
        double r = noise(double(x),double(y)); // noise value at (x,y)


        p = r;
    });

    IO::save01_in_u8(im,TEMPO_PATH + "noise.png");

    return 0;
}
