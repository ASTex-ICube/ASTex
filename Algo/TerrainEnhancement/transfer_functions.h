//
// Created by grenier on 30/03/23.
//

#ifndef ASTEX_TRANSFER_FUNCTIONS_H
#define ASTEX_TRANSFER_FUNCTIONS_H


namespace ASTex
{
    template<typename IMG>
    class Varying_Profile
    {
        using PIXT = typename IMG::PixelType;
        using EPIXT = typename IMG::DoublePixelEigen;

        const IMG& vector_noise_; // champ gaussien complexe G
        const IMG& amplitude_input_; // carte d'amplitude a
        const IMG& modulation_input_; // modulation du profile N


    private:
        double smoothstep(double edge0, double edge1, double x)
        {
            double t = std::min(std::max((x - edge0) / (edge1 - edge0), 0.), 1.); // clamp
            return t*t*(3.-2.*t);
        }

        double argument(EPIXT& vector_value)
        { // pour un uv donné, récupère les donnée du vector noise et retourn l'atan (dans -pi, pi)
            double x = vector_value(0,0)/255. - 0.5;
            double y = vector_value(1,0)/255. - 0.5;
            return std::atan2(y,x);
        }

        double global_profile(double phase)
        { // T : valeur d'entrée dans -pi, pi; valeur de sortie dans 0, 1
//            if(phase<0.){
//                return -phase/M_PI;
//            }
//            else{
//                return phase/M_PI;
//            }
            return std::abs(phase/M_PI);
        }

        double mask_profile(double phase, double step_param_l, double step_param_u)
        { // M : valeur d'entrée dans -pi, pi (phase) et 0, pi (params); valeur de sortie dans 0, 1
            return smoothstep(step_param_l, step_param_u, std::abs(phase));
        }

    protected:
        inline void clamp_channel(double& c)
        {
            if (std::is_floating_point<typename IMG::DataType>::value)
                c = std::max(0.0,std::min(1.0,c));
            else
                c = std::max(0.0,std::min(255.0,c));
        }

        inline void clamp(Eigen::Vector2d& v)
        {
            clamp_channel(v[0]);
            clamp_channel(v[1]);
        }

        inline void clamp(Eigen::Vector3d& v)
        {
            clamp_channel(v[0]);
            clamp_channel(v[1]);
            clamp_channel(v[2]);
        }

        inline void clamp(Eigen::Vector4d& v)
        {
            clamp_channel(v[0]);
            clamp_channel(v[1]);
            clamp_channel(v[2]);
            clamp_channel(v[3]);
        }

        inline void clamp(double& v)
        {
            clamp_channel(v);
        }

    public:
        Varying_Profile(const IMG& input, const IMG& amplitude, const IMG& modulation): // constructeur
        vector_noise_(input), amplitude_input_(amplitude), modulation_input_(modulation)
        {

        }

        EPIXT fetch(const Eigen::Vector2d& uv, const IMG& input)
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
            EPIXT input_interp_1 = (1.-pix_fract[0]) * input_acces(ipix_floor[0], ipix_floor[1])   + pix_fract[0] * input_acces(ipix_floor[0]+1, ipix_floor[1]);
            EPIXT input_interp_2 = (1.-pix_fract[0]) * input_acces(ipix_floor[0], ipix_floor[1]+1) + pix_fract[0] * input_acces(ipix_floor[0]+1, ipix_floor[1]+1);

            EPIXT input_interp = (1.-pix_fract[1]) * input_interp_1 + pix_fract[1] * input_interp_2;

            return input_interp;
        }


        PIXT create_varying_profile(Eigen::Vector2d uv)
        { // retourne la hauteur des détails pour une coordonnée uv donnée
            // = a * (T + M*N)(uv)
            EPIXT G_field = fetch(uv, vector_noise_); // champ gaussien
            EPIXT ampl = fetch(uv, amplitude_input_); // amplitude a
            EPIXT mod = fetch(uv, modulation_input_); // modulation de profile N

            double arg_field = argument(G_field); // champ de phase

            double profile_T = global_profile(arg_field); // profile triangle
            double profile_M = mask_profile(arg_field, 0.25*M_PI, M_PI); // mask
            double mod_noise = mod(0,0)/255.;

            double blend = (profile_T + mod_noise*profile_M)*0.5;


            return IMG::itkPixel(blend*ampl);
        }


        void details_heighmap(IMG& img_out) // récupère la hauteur des détails pour chaque pixel de l'image de sortie
        {
            img_out.parallel_for_all_pixels([&] (typename IMG::PixelType& P, int x, int y)
                                            {
                                                Eigen::Vector2d uv{ double(x) / (vector_noise_.width()), double(y) / (vector_noise_.height()) };
                                                P = create_varying_profile(uv);
                                            });
        }
    };





    template<typename IMG>
    Varying_Profile<IMG> create_procedural_details(const IMG& input, const IMG& amplitude, const IMG& modulation) // appel au constructeur
    {
        return Varying_Profile<IMG>(input, amplitude, modulation); // input = champ gaussien complexe, controlé en fréquence et en orientation
    }

}

#endif //ASTEX_TRANSFER_FUNCTIONS_H
