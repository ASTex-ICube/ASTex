//
// Created by grenier on 11/04/23.
//

#ifndef ASTEX_TERRAIN_ENHANCEMENT_H
#define ASTEX_TERRAIN_ENHANCEMENT_H

namespace ASTex{

    template<typename ImgGray>
    class Terrain_enhancement{
        const ImgGray& initial_terrain_; // terrain d'entrée
        const ImgGray& generated_details_; // détails générés

        const ImgGray& amplitude_input_; // carte d'amplitude a
        double ampl_max_ = 1.; // amplitude max des détails

        double max_value_;


    private:
        double Get_max(ImageGrayu16){
            return std::pow(2.,16.)-1.;
        }
        double Get_max(ImageGrayu8){
            return std::pow(2.,8.)-1.;
        }

        double fetch_map(const Eigen::Vector2d& uv, const ImgGray& map)
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
        Terrain_enhancement(const ImgGray& terrain, const ImgGray& details, const ImgGray& amplitude):
        initial_terrain_(terrain), generated_details_(details), amplitude_input_(amplitude)
        {
            max_value_ = Get_max(terrain);
        }

        void Set_Amplitude_max (double ampl_max)
        {
            ampl_max_ = ampl_max;
        }

        double terrain_pixel(const Eigen::Vector2d& uv)
        {
            double terrain = fetch_map(uv, initial_terrain_); // fetch sur 0, 255.
            double detail = fetch_map(uv, generated_details_);
            double ampl = ampl_max_*fetch_map(uv, amplitude_input_)/max_value_; // amplitude a, fetch sur 0, 255

            double enhancement = terrain + ampl*detail; // ampl*detail sur 0, ampl_max_

            return enhancement/(1.+ampl_max_);
        }


        void final_terrain_img(ImgGray& img_out)
        {
            img_out.parallel_for_all_pixels([&] (typename ImgGray::PixelType& P, int x, int y)
                                            {
                                                Eigen::Vector2d uv{ double(x) / (img_out.width()), double(y) / (img_out.height()) };
                                                P = terrain_pixel(uv);
                                            });
        }
    };


    template<typename ImgGray>
    Terrain_enhancement<ImgGray> compute_final_terrain(const ImgGray& terrain, const ImgGray& details, const ImgGray& amplitude)
    {
        return Terrain_enhancement(terrain, details, amplitude);
    }

}

#endif //ASTEX_TERRAIN_ENHANCEMENT_H
