#ifndef __MIXMAX_H__
#define __MIXMAX_H__

#include "ASTex/image_rgb.h"
#include "ASTex/image_gray.h"

using namespace ASTex;

class MixMax {
    public: MixMax() {}

    private: double Phi (double x) {
        return 0.5 * (1.0 + erf(x / std::sqrt(2.0)));
    }

    private: double get_mean (ImageGrayd *T) {
        double mean = 0.0;
        double n = 0;

        T->for_all_pixels([&] (double& p) {
            mean += p;
            n += 1;
        });

        mean = mean * (1.0 / n);

        return mean;
    }

    private: ImageGrayd* center (ImageGrayd *T) {
        ImageGrayd *out = new ImageGrayd(T->width(), T->height(), false);

        float mean = get_mean(T);

        out->for_all_pixels([&] (double& p, int x, int y) {
            p = T->pixelAbsolute(x, y) - mean;
        });

        return out;
    }

    public: ImageRGBd* blend (ImageRGBd *T1, ImageRGBd *T2, ImageGrayd *S1, ImageGrayd *S2, ImageGrayd *V, double lambda = 0, double interpolation_field_multiplier = 1.0) {
        ImageRGBd *out = new ImageRGBd(T1->width(), T1->height(), false);

        ImageGrayd* S1c = center(S1);
        ImageGrayd* S2c = center(S2);

        out->parallel_for_all_pixels([&] (itk::RGBPixel<double>& p, int x, int y) {
            itk::RGBPixel<double> t1 = T1->pixelAbsolute(x, y);
            itk::RGBPixel<double> t2 = T2->pixelAbsolute(x, y);
            double s1 = S1c->pixelAbsolute(x, y);
            double s2 = S2c->pixelAbsolute(x, y);
            double v1 = V->pixelAbsolute(x, y);
            double v2 = 1.0 - v1;

            v1 *= interpolation_field_multiplier;
            v2 *= interpolation_field_multiplier;

            double mu = (s2 + v2) - (s1 + v1);
            double sigma = sqrt(2.0 * pow(lambda, 2));

            double w1 = 0;
            if (sigma < 0.001) {
                if (mu < 0) w1 = 1;
            } else {
                w1 = 1.0 - Phi(mu / sigma);
            }
            double w2 = 1.0 - w1;
            
            p = t1 * w1 + t2 * w2;
        });

        delete S1c;
        delete S2c;

        return out;
    }
};


#endif