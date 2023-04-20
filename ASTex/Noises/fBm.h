//
// Created by grenier on 20/04/23.
//

#ifndef ASTEX_FBM_H
#define ASTEX_FBM_H

namespace ASTex
{

    class fBm_noise
    {
//        int x_Resolution_;
//        int y_Resolution_;
        int nb_octave_;

    private:
        double fract(double p)
        {
            return p - std::floor(p);
        }

        Eigen::Vector2d fract(Eigen::Vector2d p)
        {
            double Px = fract(p[0]);
            double Py = fract(p[1]);

            return Eigen::Vector2d{Px, Py};
        }

        Eigen::Vector2d floor(Eigen::Vector2d p)
        {
            double Px = std::floor(p[0]);
            double Py = std::floor(p[1]);

            return Eigen::Vector2d{Px, Py};
        }


    protected:



        double hash(Eigen::Vector2d p)
        {
            Eigen::Vector2d p2 = 50.*fract(p*0.31569742);
            return fract(p2[0]*p2[1]*(p2[0]+p2[1]));
        }

        double noise(Eigen::Vector2d uv)
        {
            Eigen::Vector2d p = floor(uv);
            Eigen::Vector2d w = fract(uv);

            Eigen::Vector2d u = Eigen::Vector2d{w[0]*w[0] * (3.-2.*w[0]), w[1]*w[1] * (3.-2.*w[1])};//w*w*(3.-2.*w);

            double a = hash(p+Eigen::Vector2d{0.,0.});
            double b = hash(p+Eigen::Vector2d{1.,0.});
            double c = hash(p+Eigen::Vector2d{0.,1.});
            double d = hash(p+Eigen::Vector2d{1.,1.});

            double N = 1.*(a + (b-a)*u[0] + (c-a)*u[1] + (a - b - c + d)*u[0]*u[1]);

            return std::clamp(N, 0., 1.);
        }

    public:
        fBm_noise(const int octave):nb_octave_(octave)
        {

        }


        double fBm_pixel(Eigen::Vector2d& uv)
        {
            Eigen::Vector2d UV = uv + Eigen::Vector2d{1.2, 2.4};
            double fbm_noise = 0.;
            double G = 0.5;

            double Ampl = 1.; // amplitude des octaves
            double freq = 24.; // fréquence des octaves

            for(int i=0; i<nb_octave_; i++)
            {
                fbm_noise += Ampl*noise(UV*freq);
                Ampl *= G;
                freq *= 2.;
            }

            return fbm_noise*128.;
        }


        void fBm_image(ImageGrayu8& img_out)
        {
            img_out.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
                                            {
                                                Eigen::Vector2d uv{ double(x) / (img_out.width()), double(y) / (img_out.height()) };
                                                P = fBm_pixel(uv);
                                            });
        }
    };



    fBm_noise compute_fBm(const int nb_octave)
    {
        return fBm_noise(nb_octave);
    }

}

#endif //ASTEX_FBM_H
