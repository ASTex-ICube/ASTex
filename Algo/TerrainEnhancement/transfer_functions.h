//
// Created by grenier on 30/03/23.
//

#ifndef ASTEX_TRANSFER_FUNCTIONS_H
#define ASTEX_TRANSFER_FUNCTIONS_H


namespace ASTex
{
//    template<typename IMG>
    class Varying_Profile
    {
//        using PIXT = typename IMG::PixelType;
//        using EPIXT = typename IMG::DoublePixelEigen;

        const ImageRGBu8& vector_noise_; // champ gaussien complexe G
//        const ImageGrayu8& amplitude_input_; // carte d'amplitude a
        const ImageGrayu8& modulation_input_; // modulation du profile N

//        double ampl_max_ = 1.;


    private:
        double smoothstep(double edge0, double edge1, double x)
        {
            double t = std::min(std::max((x - edge0) / (edge1 - edge0), 0.), 1.); // clamp
            return t*t*(3.-2.*t);
        }



        double argument(ImageRGBu8::DoublePixelEigen& vector_value)
        { // pour un uv donné, récupère les donnée du vector noise et retourn l'atan (dans -pi, pi)
            double x = vector_value(0,0)/255. - 0.5;
            double y = vector_value(1,0)/255. - 0.5;
            Eigen::Vector2d XY{x, y};

            // test perturbation radial
            double rho = std::sqrt(x*x+y*y);
            double perturbation = std::sin(48.*rho)/(48.*rho) + std::sin(95.*rho)/(95.*rho);
            Eigen::Matrix2d Perturbation;
            Perturbation << std::cos(perturbation), -std::sin(perturbation), std::sin(perturbation), std::cos(perturbation);
//            XY = Perturbation*XY;

            // phase
            return std::atan2(XY[1],XY[0]);
        }

        double global_profile(double phase)
        { // T : valeur d'entrée dans -pi, pi; valeur de sortie dans 0, 1
//            double profile = std::abs(phase/M_PI);
            double k = 0.01;
            double profile = std::sqrt((phase/M_PI)*(phase/M_PI) + k) + 1. - std::sqrt(1.+k);
            return profile;
        }

        double mask_profile(double phase, double step_param_l, double step_param_u)
        { // M : valeur d'entrée dans -pi, pi (phase) et 0, pi (params); valeur de sortie dans 0, 1
            double profile = smoothstep(step_param_l, step_param_u, std::abs(phase));
            return profile;
        }

    protected:

        ImageRGBu8::DoublePixelEigen fetch(const Eigen::Vector2d& uv, const ImageRGBu8& input)
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
            ImageRGBu8::DoublePixelEigen input_interp_1 = (1.-pix_fract[0]) * input_acces(ipix_floor[0], ipix_floor[1])   + pix_fract[0] * input_acces(ipix_floor[0]+1, ipix_floor[1]);
            ImageRGBu8::DoublePixelEigen input_interp_2 = (1.-pix_fract[0]) * input_acces(ipix_floor[0], ipix_floor[1]+1) + pix_fract[0] * input_acces(ipix_floor[0]+1, ipix_floor[1]+1);

            ImageRGBu8::DoublePixelEigen input_interp = (1.-pix_fract[1]) * input_interp_1 + pix_fract[1] * input_interp_2;

            return input_interp;
        }


        double fetch_map(const Eigen::Vector2d& uv, const ImageGrayu8& map)
        {
            // uv mult by map size
            Eigen::Vector2d pix_uv = Eigen::Vector2d{uv[0]*(map.width()-1), uv[1]*(map.height()-1)}; // uv in [0, 1], pix_uv in [0, 255]

            // partie entière et cast en int
            Eigen::Vector2d pix_floor {std::floor(pix_uv[0]), std::floor(pix_uv[1])}; // pix_floor in {0, 255}
            Eigen::Vector2i ipix_floor = pix_floor.cast<int>();

            // partie décimale
            Eigen::Vector2d pix_fract = pix_uv - pix_floor; // pix_fract in [0, 1[

            // accès
            auto map_acces = [&] (int xp, int yp)
            {
                return eigenPixel<double>(map.pixelAbsolute(xp, yp));
            };

            // interpolation bi-linéaire
            double map_interp_1 = (1.-pix_fract[0]) * map_acces(ipix_floor[0], ipix_floor[1])   + pix_fract[0] * map_acces(ipix_floor[0]+1, ipix_floor[1]);
            double map_interp_2 = (1.-pix_fract[0]) * map_acces(ipix_floor[0], ipix_floor[1]+1) + pix_fract[0] * map_acces(ipix_floor[0]+1, ipix_floor[1]+1);

            double map_interp = (1.-pix_fract[1]) * map_interp_1 + pix_fract[1] * map_interp_2;

            return map_interp;
        }


    public:
        Varying_Profile(const ImageRGBu8& input, const ImageGrayu8& modulation): // constructeur // const ImageGrayu8& amplitude,
        vector_noise_(input), modulation_input_(modulation) // , amplitude_input_(amplitude)
        {

        }


        double create_varying_profile(Eigen::Vector2d uv)
        { // retourne la hauteur des détails pour une coordonnée uv donnée
            // = a * (T + M*N)(uv)
            ImageRGBu8::DoublePixelEigen G_field = fetch(uv, vector_noise_); // champ gaussien
//            double ampl = ampl_max_*fetch_map(uv, amplitude_input_)/255.; // amplitude a, fetch sur 0, 255
            double mod = fetch_map(uv, modulation_input_)/255.; // modulation de profile N, fetch sur 0, 255

            double arg_field = argument(G_field); // champ de phase
//
            double profile_T = global_profile(arg_field); // profile triangle, valeurs dans 0, 1
            double profile_M = mask_profile(arg_field, 0.25*M_PI, M_PI); // mask, valeurs dans 0, 1 // TODO : exposer les param du mask ?

            double blend = (profile_T + mod*profile_M)*0.5; // valeurs sur 0, 1

            return blend *255.; // ampl*
        }


        void details_heighmap(ImageGrayu8& img_out) // récupère la hauteur des détails pour chaque pixel de l'image de sortie
        {
            img_out.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
                                            {
                                                Eigen::Vector2d uv{ double(x) / (img_out.width()), double(y) / (img_out.height()) };
                                                P = create_varying_profile(uv);
                                            });
        }
    };





//    template<typename IMG>
    Varying_Profile create_procedural_details(const ImageRGBu8& input, const ImageGrayu8& modulation) // appel au constructeur
    {
        return Varying_Profile(input, modulation); // input = champ gaussien complexe, controlé en fréquence et en orientation
    }

}

#endif //ASTEX_TRANSFER_FUNCTIONS_H
