//
// Created by grenier on 05/04/23.
//

#ifndef ASTEX_CONTROL_MAPS_H
#define ASTEX_CONTROL_MAPS_H

namespace ASTex{

    template<typename ImgGray>
    class Control_Maps{
//        using IMG = ImageRGBu8;
//        using PIXT = typename IMG::PixelType;
//        using EPIXT = typename IMG::DoublePixelEigen;

        const ImgGray& input_heightfield_;
        const ImgGray& input_gradX_;
        const ImgGray& input_gradY_;

        double max_value_;

        double p_frequ1_ = 2.;
        double p_frequ2_ = 1.;
        double p_ampl1_ = 2.;
        double p_ampl2_ = 4.;


    private:
        double Get_max(ImageGrayu16){
            return std::pow(2.,16.)-1.;
        }
        double Get_max(ImageGrayu8){
            return std::pow(2.,8.)-1.;
        }

        double fetch(const Eigen::Vector2d& uv, const ImgGray& input)
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
        Control_Maps(const ImgGray& terrain, const ImgGray& gradX, const ImgGray& gradY):
        input_heightfield_(terrain), input_gradX_(gradX), input_gradY_(gradY)
        {
            max_value_ = Get_max(terrain);
        }

        void Set_Frequency_param(double facteur, double puissance)
        {
            p_frequ1_ = facteur;
            p_frequ2_ = puissance;
        }

        void Set_Amplitude_param(double facteur, double puissance)
        {
            p_ampl1_ = facteur;
            p_ampl2_ = puissance;
        }


        double frequency_map(Eigen::Vector2d& uv) // définit comment la fréquence varie spatialment, pas les valeur min et max
        {
            double h = fetch(uv, input_heightfield_)/max_value_; // valeurs sur 0,1

            double freq = p_frequ1_*std::pow(h, p_frequ2_); // sur 0, p1

            return freq*(max_value_/p_frequ1_); // sortie sur 0, 255
        }


        double orientation_map(Eigen::Vector2d& uv)
        {
            double Gx = fetch(uv, input_gradX_)-max_value_/2.; // re-centrage des valeurs dans -127, 127
            double Gy = fetch(uv, input_gradY_)-max_value_/2.;
//            std::cout<<Gx<<std::endl;

            double ori = std::atan2(Gy, Gx);// valeurs dans -pi,pi
//            std::cout<<ori<<std::endl;
//            std::cout<<(ori+M_PI) * (255./(2.*M_PI))<<std::endl;

            return (ori+M_PI) * (max_value_/(2.*M_PI)); // angle sur -pi, pi; mappé sur 0, 255
        }


        double amplitude_map(Eigen::Vector2d& uv) // définit comment l'amplitude varie spatialment, pas les valeur min et max
        {
            double Gx = (fetch(uv, input_gradX_)-floor(max_value_/2.))/floor(max_value_/2.); // fetch sur 0, 255
            double Gy = (fetch(uv, input_gradY_)-floor(max_value_/2.))/floor(max_value_/2.); // valeurs sur -1, 1

            double norm = std::sqrt(Gx*Gx + Gy*Gy); // sur 0, sqrt(2)

            double ampl = p_ampl1_*std::pow(norm,p_ampl2_); // sur 0, p1*sqrt(2)^p2

            return ampl*(max_value_/(p_ampl1_*std::pow(std::sqrt(2.),p_ampl2_))); // sortie sur 0, 255
        }


        void compute_control(ImgGray& output_frequ, ImgGray& output_or, ImgGray& output_ampl)
        {
            output_frequ.parallel_for_all_pixels([&] (typename ImgGray::PixelType& P, int x, int y)
                                            {
                                                Eigen::Vector2d uv{ double(x) / (output_frequ.width()), double(y) / (output_frequ.height()) };
                                                P = frequency_map(uv);
                                            });

            output_ampl.parallel_for_all_pixels([&] (typename ImgGray::PixelType& P, int x, int y)
                                                 {
                                                     Eigen::Vector2d uv{ double(x) / (output_ampl.width()), double(y) / (output_ampl.height()) };
                                                     P = amplitude_map(uv);
                                                 });

            output_or.parallel_for_all_pixels([&] (typename ImgGray::PixelType& P, int x, int y)
                                                {
                                                    Eigen::Vector2d uv{ double(x) / (output_or.width()), double(y) / (output_or.height()) };
                                                    P = orientation_map(uv);
                                                });
        }
    };


    template<typename ImgGray>
    Control_Maps<ImgGray> create_control_maps(const ImgGray& terrain, const ImgGray& gradX, const ImgGray& gradY)
    {
        return Control_Maps(terrain, gradX, gradY);
    }

}
#endif //ASTEX_CONTROL_MAPS_H
