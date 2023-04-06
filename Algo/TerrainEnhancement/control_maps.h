//
// Created by grenier on 05/04/23.
//

#ifndef ASTEX_CONTROL_MAPS_H
#define ASTEX_CONTROL_MAPS_H

namespace ASTex{

    class Control_Maps{
//        using IMG = ImageRGBu8;
//        using PIXT = typename IMG::PixelType;
//        using EPIXT = typename IMG::DoublePixelEigen;

        const ImageGrayu8& input_heightfield_;
        const ImageGrayu8& input_gradX_;
        const ImageGrayu8& input_gradY_;

        double p_frequ1_;
        double p_frequ2_;
        double p_ampl1_;
        double p_ampl2_;


    private:
        double fetch(const Eigen::Vector2d& uv, const ImageGrayu8& input)
        {
            // uv mult by map size
            Eigen::Vector2d pix_uv = Eigen::Vector2d{uv[0]*(input.width()-1), uv[1]*(input.height()-1)}; // uv in [0, 1], pix_uv in [0, 255]

            // partie entière et cast en int
            Eigen::Vector2d pix_floor {std::floor(pix_uv[0]), std::floor(pix_uv[1])}; // pix_floor in {0, 255}
            Eigen::Vector2i ipix_floor = pix_floor.cast<int>();

            // partie décimale
            Eigen::Vector2d pix_fract = pix_uv - pix_floor; // pix_fract in [0, 1[

            // accès
            auto input_acces = [&] (int xp, int yp)
            {
                return eigenPixel<double>(input.pixelAbsolute(xp, yp));
            };

            // interpolation bi-linéaire
            double input_interp_1 = (1.-pix_fract[0]) * input_acces(ipix_floor[0], ipix_floor[1])   + pix_fract[0] * input_acces(ipix_floor[0]+1, ipix_floor[1]);
            double input_interp_2 = (1.-pix_fract[0]) * input_acces(ipix_floor[0], ipix_floor[1]+1) + pix_fract[0] * input_acces(ipix_floor[0]+1, ipix_floor[1]+1);

            double input_interp = (1.-pix_fract[1]) * input_interp_1 + pix_fract[1] * input_interp_2;

            return input_interp;
        }





    public:
        Control_Maps(const ImageGrayu8& terrain, const ImageGrayu8& gradX, const ImageGrayu8& gradY):
        input_heightfield_(terrain), input_gradX_(gradX), input_gradY_(gradY)
        {
            p_frequ1_ = 2.;
            p_frequ2_ = 1.;
            p_ampl1_ = 1.;
            p_ampl2_ = 1.;
        }



        double frequency_map(Eigen::Vector2d& uv)
        {
            double h = fetch(uv, input_heightfield_)/255.; // valeurs sur 0,1

            double freq = p_frequ1_*std::pow(h, p_frequ2_) + 1.; // +1 pour ne pas tomber à 0

            return freq; // valeur retournée directement, donc fréquence possible entre 1 et 255 // TODO : à changer ?
        }


        double orientation_map(Eigen::Vector2d& uv)
        {
            double Gx = fetch(uv, input_gradX_)-127; // re-centrage des valeurs dans -127, 127
            double Gy = fetch(uv, input_gradY_)-127;
//            std::cout<<Gx<<std::endl;

            double ori = std::atan2(Gy, Gx);// valeurs dans -pi,pi
//            std::cout<<ori<<std::endl;
//            std::cout<<(ori+M_PI) * (255./(2.*M_PI))<<std::endl;

            return (ori+M_PI) * (255./(2.*M_PI)); // angle sur -pi, pi; mappé sur 0, 255
        }


        double amplitude_map(Eigen::Vector2d& uv)
        {
            double Gx = fetch(uv, input_gradX_);
            double Gy = fetch(uv, input_gradY_);

            Eigen::Vector2d Grad{Gx, Gy};
            double norme = Grad.norm(); // TODO : pb avec la saturation des gradient

            double ampl = p_ampl1_*std::pow(norme,p_ampl2_);

            return ampl*255./std::sqrt(2.); // TODO : quelle plage de valeur ?
        }


        void compute_control(ImageGrayu8& output_frequ, ImageGrayu8& output_or, ImageGrayu8& output_ampl)
        {
            output_frequ.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
                                            {
                                                Eigen::Vector2d uv{ double(x) / (output_frequ.width()), double(y) / (output_frequ.height()) };
                                                P = frequency_map(uv);
                                            });
//
//            output_ampl.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
//                                                 {
//                                                     Eigen::Vector2d uv{ double(x) / (output_ampl.width()), double(y) / (output_ampl.height()) };
//                                                     P = amplitude_map(uv);
//                                                 });

            output_or.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
                                                {
                                                    Eigen::Vector2d uv{ double(x) / (output_or.width()), double(y) / (output_or.height()) };
                                                    P = orientation_map(uv);
                                                });
        }
    };


    Control_Maps create_control_maps(const ImageGrayu8& terrain, const ImageGrayu8& gradX, const ImageGrayu8& gradY)
    {
        return Control_Maps(terrain, gradX, gradY);
    }

}
#endif //ASTEX_CONTROL_MAPS_H
