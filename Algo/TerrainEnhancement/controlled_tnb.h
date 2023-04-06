//
// Created by grenier on 29/03/23.
//

namespace ASTex
{

/**
 Simple Tiling n blending algo
 Need Gaussian imagea as input
 */
//    template<typename IMG>
    class Tiling_n_Blending
    {
//        using PIXT = typename IMG::PixelType;
//        using EPIXT = typename IMG::DoublePixelEigen;

        const ImageRGBu8& img_input_; // exemple d'entrée
        const ImageGrayu8& frequency_input_; // carte de fréquence f
        const ImageGrayu8& orientation_input_; // carte d'orientation theta

        ImageRGBu8::DoublePixelEigen img_average_;

    protected:

        inline void set_zero(Eigen::Vector2d& v)
        {
            v.setZero();
        }

        inline void set_zero(Eigen::Vector3d& v)
        {
            v.setZero();
        }

        inline void set_zero(Eigen::Vector4d& v)
        {
            v.setZero();
        }

        inline void set_zero(double& v)
        {
            v = 0.0;
        }

        inline void clamp_channel(double& c)
        {
            if (std::is_floating_point<typename ImageRGBu8::DataType>::value)
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


        void TriangleGrid(const Eigen::Vector2d& p_uv, Eigen::Vector3d& Bi, Eigen::Vector2i& vertex1, Eigen::Vector2i& vertex2, Eigen::Vector2i& vertex3)
        {
            const Eigen::Vector2d uv = p_uv * 2.0 * std::sqrt(3.0); // iGridSize

            Eigen::Matrix2d gridToSkewedGrid;
            gridToSkewedGrid << 1.0, -1./sqrt(3.),
                                0.0, 2./sqrt(3.);

            Eigen::Vector2d skewedCoord = gridToSkewedGrid * uv;
            Eigen::Vector2d baseId{ std::floor(skewedCoord[0]), std::floor(skewedCoord[1]) };
            Eigen::Vector3d temp{ skewedCoord[0] - baseId[0], skewedCoord[1] - baseId[1], 0.0 };
            temp[2] = 1.0 - temp[0] - temp[1];

            if (temp[2] > 0.0)
            {
                Bi = Eigen::Vector3d(temp[2], temp[1], temp[0]);
                Eigen::Vector2i ibaseId = baseId.cast<int>();
                vertex1 = ibaseId;
                vertex2 = ibaseId + Eigen::Vector2i(0, 1);
                vertex3 = ibaseId + Eigen::Vector2i(1, 0);
            }
            else
            {
                Bi = Eigen::Vector3d(-temp[2], 1.0 - temp[1], 1.0 - temp[0]);
                Eigen::Vector2i ibaseId = baseId.cast<int>();
                vertex1 = ibaseId + Eigen::Vector2i(1, 1);
                vertex2 = ibaseId + Eigen::Vector2i(1, 0);
                vertex3 = ibaseId + Eigen::Vector2i(0, 1);
            }
        }

        Eigen::Vector2d MakeCenterUV(Eigen::Vector2i& vertex)
        {
            Eigen::Matrix2d invSkewMat;
            invSkewMat << 1.0, 0.5,
                          0.0, sqrt(3.)/2.;

            Eigen::Vector2d center = (invSkewMat * vertex.cast<double>())/(2.0 * std::sqrt(3.0));
            return center;
        }

        Eigen::Matrix2d TransfoMatrix(double theta, double scale)
        {
            Eigen::Matrix2d Rotation;
            Eigen::Matrix2d Scale;

            // theta est sur -pi, pi; mappé sur 0, 255, on veut angle dans -pi,pi
//            std::cout<<theta<<std::endl;
            double angle = (theta*2.*M_PI/255.) - M_PI;
//            std::cout<<angle<<std::endl;

            double cos = std::cos(angle);
            double sin = std::sin(angle);

            Rotation << cos, sin, -sin, cos; // rotation autour de l'axe -z ( (0,0) en haut à gauche avec x vers la droite et y vers le bas)
            Scale << scale, 0., 0., scale;

            return Rotation*Scale;
        }


        //original hash version
        Eigen::Vector2d hash(const Eigen::Vector2i& p)
        {
            Eigen::Matrix2d hashMat;
            hashMat << 127.1, 269.5,
                    311.7, 183.3;

            Eigen::Vector2d q = hashMat * p.cast<double>();
            q[0] = std::sin(q[0]);
            q[1] = std::sin(q[1]);
            q *= 43758.5453;
            return Eigen::Vector2d(q[0] - std::floor(q[0]), q[1] - std::floor(q[1]));

        }

    public:
        Tiling_n_Blending(const ImageRGBu8& in, const ImageGrayu8& freq, const ImageGrayu8& theta): // constructeur : récupère l'image d'entrée et calcul la moyenne
                img_input_(in), frequency_input_(freq), orientation_input_(theta)
        {
            // compute average img value
            ImageRGBu8::DoublePixelEigen sum;
            set_zero(sum);
            img_input_.for_all_pixels([&] ( const ImageRGBu8::PixelType& P)
                                      {
                                          ImageRGBu8::DoublePixelEigen lv = eigenPixel<double>(P);
                                          sum += lv;
                                      });

            img_average_ = sum / double(img_input_.width()*img_input_.height());

        }

        ImageRGBu8::DoublePixelEigen fetch(const Eigen::Vector2d& uv)
        {
            // take fract part of uv mult by image size
            Eigen::Vector2d uvd = Eigen::Vector2d{-0.5 + (uv[0] - std::floor(uv[0]))*img_input_.width(),
                                                  -0.5 + (uv[1] - std::floor(uv[1]))*img_input_.height()};
            Eigen::Vector2d uvfl {std::floor(uvd[0]), std::floor(uvd[1])};
            //  a = coef for linear interpolation
            Eigen::Vector2d a = uvd - uvfl;
            // c = integer texel coord
            Eigen::Vector2i c = uvfl.cast<int>();

            auto acces_repeat = [&] (int xp, int yp)
            {
                int xx = xp < 0 ? img_input_.width() - 1 : ( xp >= img_input_.width() ? 0 : xp );
                int yy = yp < 0 ? img_input_.height() - 1 : ( yp >= img_input_.height() ? 0 : yp );
                return eigenPixel<double>(img_input_.pixelAbsolute(xx, yy));
            };


            return acces_repeat(c[0],c[1]);

            ImageRGBu8::DoublePixelEigen V1 = (1.0 - a[0]) * acces_repeat(c[0],c[1])	+ a[0] * acces_repeat(c[0]+1,c[1]);
            ImageRGBu8::DoublePixelEigen V2 = (1.0 - a[0]) * acces_repeat(c[0],c[1]+1) + a[0] * acces_repeat(c[0]+1,c[1]+1);
            ImageRGBu8::DoublePixelEigen V = (1.0 - a[1]) * V1 + a[1] * V2;

            return V;
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


        ImageRGBu8::PixelType tile_pixel(const Eigen::Vector2d& uv ) // fait le tnb pour un pixel
        {
            // grille
            Eigen::Vector3d B;
            Eigen::Vector2i  vertex1, vertex2, vertex3;
            TriangleGrid(uv, B,	vertex1, vertex2, vertex3);

            // control maps
            double freq = fetch_map(uv, frequency_input_); // fetch sur 0, 255
            double theta = fetch_map(uv, orientation_input_); // fetch sur 0, 255

            Eigen::Matrix2d M_transfo = TransfoMatrix(theta, freq);

            // centers of tiles
            Eigen::Vector2d cen1 = MakeCenterUV(vertex1);
            Eigen::Vector2d cen2 = MakeCenterUV(vertex2);
            Eigen::Vector2d cen3 = MakeCenterUV(vertex3);


            // Assign random offset to each triangle vertex
            Eigen::Vector2d uv1 = M_transfo*(uv-cen1) + cen1 + hash(vertex1);
            Eigen::Vector2d uv2 = M_transfo*(uv-cen2) + cen2 + hash(vertex2);
            Eigen::Vector2d uv3 = M_transfo*(uv-cen3) + cen3 + hash(vertex3);


            ImageRGBu8::DoublePixelEigen t1 = fetch(uv1) - img_average_;
            ImageRGBu8::DoublePixelEigen t2 = fetch(uv2) - img_average_;
            ImageRGBu8::DoublePixelEigen t3 = fetch(uv3) - img_average_;

            auto W = B.normalized();

            ImageRGBu8::DoublePixelEigen P = W[0] * t1 + W[1] * t2 + W[2] * t3 + img_average_ ;

            clamp(P);

            return ImageRGBu8::itkPixel(P);
        }


        void tile_img(ImageRGBu8& img_out) // récupère le tnb pour chaque pixel de l'image de sortie
        {
            img_out.parallel_for_all_pixels([&] (typename ImageRGBu8::PixelType& P, int x, int y)
                                            {
                                                Eigen::Vector2d uv{ double(x) / (img_out.width()), double(y) / (img_out.height()) };
                                                P = tile_pixel(uv);
                                            });
        }

    };




//    template<typename IMG>
    Tiling_n_Blending make_Tiling_n_Blending(const ImageRGBu8& img, const ImageGrayu8& freq, const ImageGrayu8& theta) // appel au constructeur
    {
        return Tiling_n_Blending(img, freq, theta);
    }

}
