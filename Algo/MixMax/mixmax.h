#ifndef __MIXMAX_H__
#define __MIXMAX_H__

#include "ASTex/image_rgb.h"
#include "ASTex/image_gray.h"

using namespace ASTex;

namespace MixMax {
    /*
    * PRIVATE FUNCTIONS
    */
    namespace {
        double phi (double x) {
            return 0.5 * (1.0 + erf(x / std::sqrt(2.0)));
        }

        itk::RGBPixel<double> get_mean (ImageRGBd *T, Region& r) {
            itk::RGBPixel<double> mean;
            mean.SetRed(0);
            mean.SetGreen(0);
            mean.SetBlue(0);
            double n = 0;

            T->for_region_pixels(r, [&] (itk::RGBPixel<double>& p) {
                mean += p;
                n += 1;
            });
            mean = mean * (1.0 / n);

            return mean;
        }

        itk::RGBPixel<double> get_mean (ImageRGBd *T) {
            Region r = gen_region(0, 0, T->width(), T->height());
            return get_mean(T, r);
        }

        double get_mean (ImageGrayd *T, Region r) {
            double mean = 0.0;
            double n = 0;

            T->for_region_pixels(r, [&] (double& p) {
                mean += p;
                n += 1;
            });
            mean = mean / n;

            return mean;
        }

        double get_mean (ImageGrayd *T) {
            Region r = gen_region(0, 0, T->width(), T->height());
            return get_mean(T, r);
        }

        double get_variance (ImageGrayd *T, Region r) {
            double mean = get_mean(T, r);
            double variance = 0;
            double n = 0;

            T->for_region_pixels(r, [&] (double& p) {
                variance += (p - mean) * (p - mean);
                n += 1;
            });
            variance = variance * (1.0 / n);

            return variance;
        }

        double get_variance (ImageGrayd *T) {
            Region r = gen_region(0, 0, T->width(), T->height());

            return get_variance(T, r);
        }

        ImageGrayd* center (ImageGrayd *T) {
            ImageGrayd *out = new ImageGrayd(T->width(), T->height(), false);
            float mean = get_mean(T);

            out->for_all_pixels([&] (double& p, int x, int y) {
                p = T->pixelAbsolute(x, y) - mean;
            });

            return out;
        }
    }
    
    /*
    * MIXMAX IMPLEMENTATION
    */
    ImageRGBd* blend (ImageRGBd *T1, ImageRGBd *T2, ImageGrayd *S1, ImageGrayd *S2, ImageGrayd *V, int mip_level = 0, double lambda = 0, double interpolation_field_multiplier = 1.0) {
        int region_size = std::pow(2, mip_level);
        ImageRGBd *out = new ImageRGBd(T1->width() / region_size, T1->height() / region_size, false);
        ImageGrayd* S1c = center(S1);
        ImageGrayd* S2c = center(S2);

        out->parallel_for_all_pixels([&] (itk::RGBPixel<double>& p, int x, int y) {
            Region r = gen_region(x * region_size, y * region_size, region_size, region_size);
            itk::RGBPixel<double> t1_mean = get_mean(T1, r);
            itk::RGBPixel<double> t2_mean = get_mean(T2, r);
            double s1_mean = get_mean(S1c, r);
            double s2_mean = get_mean(S2c, r);
            double s1_variance = get_variance(S1c, r);
            double s2_variance = get_variance(S2c, r);
            double v1_mean = get_mean(V, r);
            double v2_mean = 1.0 - v1_mean;

            v1_mean *= interpolation_field_multiplier;
            v2_mean *= interpolation_field_multiplier;
            double mean = (s2_mean + v2_mean) - (s1_mean + v1_mean);
            double std_dev = sqrt(s1_variance + s2_variance + 2.0 * pow(lambda, 2));
            double w1 = 0;
            if (std_dev < 0.001) {
                if (mean < 0) w1 = 1;
            } else {
                w1 = 1.0 - phi(mean / std_dev);
            }
            double w2 = 1.0 - w1;

            p = t1_mean * w1 + t2_mean * w2;
        });
        delete S1c;
        delete S2c;

        return out;
    }

    ImageRGBd* ground_truth (ImageRGBd *T1, ImageRGBd *T2, ImageGrayd *S1, ImageGrayd *S2, ImageGrayd *V, int mip_level = 0, double lambda = 0, double interpolation_field_multiplier = 1.0) {
        int region_size = std::pow(2, mip_level);
        ImageRGBd *out = new ImageRGBd(T1->width() / region_size, T1->height() / region_size, false);
        ImageRGBd *mixmax = blend(T1, T2, S1, S2, V, 0, lambda, interpolation_field_multiplier);

        out->parallel_for_all_pixels([&] (itk::RGBPixel<double>& p, int x, int y) {
            Region r = gen_region(x * region_size, y * region_size, region_size, region_size);

            p = get_mean(mixmax, r);
        });
        delete mixmax;

        return out;
    }
};


#endif