#ifndef PNF_H
#define PNF_H

#include <Eigen/Eigen>
#include <ASTex/color_map.h>
#include <ASTex/noise.h>

namespace ASTex {

using T = double;
using Vec2 = Eigen::Matrix<T,2,1>;
using Mat22 = Eigen::Matrix<T,2,2>;
using Color = Color_map<T>::Color;

template<typename func>
inline ImageRGB<T> computeIMG(const Vec2 & w_size, const Vec2 &im_size, const func &f){
    Mat22 borns;
    borns << w_size(0), 0, 0, w_size(1);
    Vec2 center = Vec2(borns(0),borns(3)) / T(2);

    ImageRGB<T> im(static_cast<int>(im_size(0)),static_cast<int>(im_size(1)));
    im.parallel_for_all_pixels([&](ImageRGB<T>::PixelType &pix,int i,int j)
    {
        // x = pos(0) between [-center(0),center(0)]
        // y = pos(1) between [-center(1),center(1)]
        Vec2 pos =  borns * Vec2(T(i) / T(im.width()), T(j) / T(im.height())) - center;

        pix = f(pos);
    });

    return im;
}

inline ImageRGB<T> compute_unfiltered_IMG(const Vec2 & w_size,
                                          const Vec2 &im_size,
                                          const Noise<T> &n,
                                          const Color_map<T> &cm)
{
    return computeIMG(w_size, im_size, [&](const Vec2 &pos) {
        return ImageRGB<T>::itkPixel(cm.map(n.basic2D(pos)));
    });
}

inline ImageRGB<T> compute_naive_filter_IMG(const Vec2 & w_size,
                                            const Vec2 &im_size,
                                            const Noise<T> &n,
                                            const Color_map<T> &cm,
                                            const int &nb_sample)
{
    Vec2 footprint(w_size(0) / im_size(0), w_size(1) / im_size(1));
    return computeIMG(w_size, im_size, [&](const Vec2 &pos) {
        T mean = n.get_noise_mean_over_footprint(pos, footprint, nb_sample);
        return ImageRGB<T>::itkPixel(cm.map(mean));
    });
}

inline ImageRGB<T> compute_good_filter_IMG( const Vec2 & w_size,
                                            const Vec2 &im_size,
                                            const Noise<T> &n,
                                            const Color_map<T> &cm,
                                            const int &nb_sample)
{
    Vec2 footprint(w_size(0) / im_size(0), w_size(1) / im_size(1));
    return computeIMG(w_size, im_size, [&](const Vec2 &pos) {
        T mean = n.get_noise_mean_over_footprint(pos, footprint, nb_sample);
        T squared_mean = n.get_squared_noise_mean_over_footprint(pos, footprint, nb_sample);
        T sigma = std::sqrt(squared_mean - mean * mean);
        return cm.map(mean, sigma);
    });
}

}

#endif // PNF_H
