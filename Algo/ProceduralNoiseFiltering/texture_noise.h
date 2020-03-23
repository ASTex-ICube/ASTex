#ifndef TEXTURE_NOISE_H
#define TEXTURE_NOISE_H

#include <ASTex/easy_io.h>
#include <ASTex/image_spectral.h>
#include <ASTex/rpn_utils.h>
#include "color_map.h"

namespace ASTex {


template <typename T>
class TextureNoise
{
private:
    ImageGray<T> noise;
public:
    using Vec2 = Eigen::Matrix<T,2,1>;
    using Color = typename Color_map<T>::Color;

    TextureNoise() {}

    ImageGray<T> getNoise() const {
        return noise;
    }

    void setNoise(const ImageGray<T> &im){
        noise = im;
    }

    int width() const {
        return noise.width();
    }

    int height() const {
        return noise.height();
    }

    T get(const int &i, const int &j) const {
        int x = (i % width() + width()) % width();
        int y = (j % height() + height()) % height();

        return noise.pixelAbsolute(x, y);
    }

    T get(const Vec2 &pos) const {
        int x = static_cast<int>(std::floor(pos(0)));
        int y = static_cast<int>(std::floor(pos(1)));
        return get(x,y);
    }

private :
    template<typename T2,typename init, typename func>
    T2 get_over_footprint(const Vec2 &pos, const Vec2 &footprint, const init &initialisation, const func &f) const{
        Vec2 corner = pos - footprint * T(0.5);
        int nb_texel_x = static_cast<int>(std::ceil(footprint(0)));
        int nb_texel_y = static_cast<int>(std::ceil(footprint(1)));
        T2 ret = initialisation();
        for(int i = 0; i < nb_texel_y ; ++i){
            int y = static_cast<int>(std::floor(corner(1))) + i;
            for (int j=0; j < nb_texel_x ; ++j) {
                int x = static_cast<int>(std::floor(corner(0))) + j;
                T value_noise = get(x,y);
                ret += f(value_noise);
            }
        }

        ret /= nb_texel_x * nb_texel_y;

        return ret;

//        Vec2 corner_start = pos - footprint * T(0.5);
//        T2 ret = initialisation();
//        int nb_texel_x = static_cast<int>(std::ceil(footprint(0)));
//        int nb_texel_y = static_cast<int>(std::ceil(footprint(1)));

//        for(int i = 0; i < nb_texel_y ; ++i){
//            //int y = static_cast<int>(std::floor(corner_start(1))) + i;
//            T y  = corner_start(1) + i;
//            T y_floor = std::floor(y);
//            T y_ceil = std::ceil(y);
//            T u = y - y_floor;

//            for (int j=0; j < nb_texel_x; ++j) {
//                //int x = static_cast<int>(std::floor(corner_start(0))) + j;
//                T x = corner_start(0) + j;
//                T x_floor = std::floor(x);
//                T x_ceil = std::ceil(x);
//                T v = x - x_floor;

//                T v1 = (1 - u) * get(x_floor, y_floor) + u * get(x_floor, y_ceil);
//                T v2 = (1 - u) * get(x_ceil, y_floor) + u * get(x_ceil, y_ceil);

//                T value_noise = (1 - v) * v1 + v * v2;
//                ret += f(value_noise);
//            }
//        }

//        ret /= nb_texel_x * nb_texel_y;

//        return ret;
    }

public:
    T get_noise_mean_over_footprint(const Vec2 &pos, const Vec2 &footprint) const{
        return get_over_footprint<T>(pos, footprint, [](){return T(0);}, [](const T&x) {
            return x;
        });
    }

    T get_squared_noise_mean_over_footprint(const Vec2 &pos, const Vec2 &footprint) const{
        return get_over_footprint<T>(pos, footprint, [](){return T(0);}, [](const T&x) {
            return x * x;
        });
    }

    Color get_color_mean_over_footprint(const Vec2 &pos, const Vec2 &footprint, const Color_map<T> &cm) const{
        return get_over_footprint<Color>(pos, footprint, [](){return Color(0,0,0);}, [&](const T&x) {
            return cm.map(x, 0);
        });
    }
};

}
#endif // TEXTURE_NOISE_H
