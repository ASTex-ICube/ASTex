#ifndef PNF_H
#define PNF_H

#include <Eigen/Eigen>
#include "color_map.h"
#include "noise.h"
#include "texture_noise.h"

namespace ASTex {

template <typename T>
class Pnf
{
public:
    using Vec2 = Eigen::Matrix<T,2,1>;
    using Mat22 = Eigen::Matrix<T,2,2>;
    using Color = typename Color_map<T>::Color;

private:
    //function noise
    template<typename IMG, typename func>
    static inline IMG computeIMG(const Vec2 & w_size, const Vec2 &im_size, const func &f){
        Mat22 borns;
        borns << w_size(0), 0, 0, w_size(1);
        Vec2 center = Vec2(borns(0),borns(3)) / T(2);

        IMG im(static_cast<int>(im_size(0)),static_cast<int>(im_size(1)));
        im.parallel_for_all_pixels([&](typename IMG::PixelType &pix,int i,int j)
        {
            // x = pos(0) between [-center(0),center(0)]
            // y = pos(1) between [-center(1),center(1)]
            Vec2 pos =  borns * Vec2(T(i) / T(im.width()), T(j) / T(im.height())) - center;

            pix = IMG::itkPixel(f(pos));
        });

        return im;
    }

    //texture noise
    template<typename IMG, typename func>
    static inline IMG computeIMG(const Vec2 &offset, const Vec2 & w_size, const Vec2 &im_size, const func &f){
        Mat22 borns;
        borns << w_size(0), 0, 0, w_size(1);

        IMG im(static_cast<int>(im_size(0)),static_cast<int>(im_size(1)));
        im.parallel_for_all_pixels([&](typename IMG::PixelType &pix,int i,int j)
        {
            Vec2 pos =  borns * Vec2(T(i) / T(im.width()), T(j) / T(im.height())) + offset ;

            pix = IMG::itkPixel(f(pos));
        });

        return im;
    }

public:
    Pnf() {}

    //function noise
    static inline ImageGray<T> compute_noise_unfiltered(const Vec2 & w_size,
                                                        const Vec2 &im_size,
                                                        const Noise<T> &n)
    {
        return computeIMG<ImageGray<T>>(w_size, im_size, [&](const Vec2 &pos) {
            return n.get(pos);
        });
    }

    static inline ImageGray<T> compute_noise_filtered(const Vec2 & w_size,
                                                      const Vec2 &im_size,
                                                      const Noise<T> &n,
                                                      const int &nb_sample)
    {
        Vec2 footprint(w_size(0) / im_size(0), w_size(1) / im_size(1));
        return computeIMG<ImageGray<T>>(w_size, im_size, [&](const Vec2 &pos) {
            return n.get_noise_mean_over_footprint(pos, footprint, nb_sample);
        });
    }
    static inline ImageRGB<T> compute_unfiltered_IMG(const Vec2 & w_size,
                                                     const Vec2 &im_size,
                                                     const Noise<T> &n,
                                                     const Color_map<T> &cm)
    {
        return computeIMG<ImageRGB<T>>(w_size, im_size, [&](const Vec2 &pos) {
            return cm.map(n.get(pos), 0);
        });
    }

    static inline ImageRGB<T> compute_ground_truth_IMG(const Vec2 & w_size,
                                                       const Vec2 &im_size,
                                                       const Noise<T> &n,
                                                       const Color_map<T> &cm,
                                                       const int &nb_sample)
    {
        Vec2 footprint(w_size(0) / im_size(0), w_size(1) / im_size(1));
        return computeIMG<ImageRGB<T>>(w_size, im_size, [&](const Vec2 &pos) {
            return n.get_color_mean_over_footprint(pos, footprint, cm, nb_sample);
        });
    }

    static inline ImageRGB<T> compute_naive_filter_IMG(const Vec2 & w_size,
                                                       const Vec2 &im_size,
                                                       const Noise<T> &n,
                                                       const Color_map<T> &cm,
                                                       const int &nb_sample)
    {
        Vec2 footprint(w_size(0) / im_size(0), w_size(1) / im_size(1));
        return computeIMG<ImageRGB<T>>(w_size, im_size, [&](const Vec2 &pos) {
            T mean = n.get_noise_mean_over_footprint(pos, footprint, nb_sample);
            return cm.map(mean, 0);
        });
    }

    static inline ImageRGB<T> compute_good_filter_IMG(const Vec2 & w_size,
                                                      const Vec2 &im_size,
                                                      const Noise<T> &n,
                                                      const Color_map<T> &cm,
                                                      const int &nb_sample)
    {
        Vec2 footprint(w_size(0) / im_size(0), w_size(1) / im_size(1));
        return computeIMG<ImageRGB<T>>(w_size, im_size, [&](const Vec2 &pos) {
            T mean = n.get_noise_mean_over_footprint(pos, footprint, nb_sample);
            T squared_mean = n.get_squared_noise_mean_over_footprint(pos, footprint, nb_sample);
            T sigma = std::sqrt(squared_mean - mean * mean);
            return cm.map(mean, sigma);
        });
    }

    //function texture noise
    static inline ImageGray<T> compute_noise_unfiltered(const Vec2 & w_size,
                                                        const Vec2 &im_size,
                                                        const TextureNoise<T> &n)
    {
        Vec2 offset(T(n.width()) * T(0.5), T(n.height()) * T(0.5));
        offset -= w_size * T(0.5);

        Vec2 t_size(n.width(), n.height());
        Vec2 footprint(w_size(0) / t_size(0), w_size(1) / t_size(1));

        return computeIMG<ImageGray<T>>(offset, w_size, im_size, [&](Vec2 &pos) {
            Vec2 corner_start = pos - footprint * T(0.5);
            return n.get(corner_start);
        });
    }

    static inline ImageGray<T> compute_noise_filtered(const Vec2 & w_size,
                                                      const Vec2 &im_size,
                                                      const TextureNoise<T> &n)
    {
        Vec2 offset(T(n.width()) * T(0.5), T(n.height()) * T(0.5));
        offset -= w_size * T(0.5);

        Vec2 t_size(n.width(), n.height());
        Vec2 footprint(w_size(0) / t_size(0), w_size(1) / t_size(1));

        return computeIMG<ImageGray<T>>(offset, w_size, im_size, [&](const Vec2 &pos) {
            return n.get_noise_mean_over_footprint(pos, footprint);
        });
    }


    static inline ImageRGB<T> compute_unfiltered_IMG(const Vec2 & w_size,
                                                     const Vec2 &im_size,
                                                     const TextureNoise<T> &n,
                                                     const Color_map<T> &cm)
    {
        Vec2 offset(T(n.width()) * T(0.5), T(n.height()) * T(0.5));
        offset -= w_size * T(0.5);

        Vec2 t_size(n.width(), n.height());
        Vec2 footprint(w_size(0) / t_size(0), w_size(1) / t_size(1));

        return computeIMG<ImageRGB<T>>(offset, w_size, im_size, [&](Vec2 &pos) {
            Vec2 corner_start = pos - footprint * T(0.5);
            return cm.map(n.get(corner_start), 0);
        });
    }

    static inline ImageRGB<T> compute_ground_truth_IMG(const Vec2 & w_size,
                                                       const Vec2 &im_size,
                                                       const TextureNoise<T> &n,
                                                       const Color_map<T> &cm)
    {
        Vec2 offset(T(n.width()) * T(0.5), T(n.height()) * T(0.5));
        offset -= w_size * T(0.5);

        Vec2 t_size(n.width(), n.height());
        Vec2 footprint(w_size(0) / t_size(0), w_size(1) / t_size(1));

        return computeIMG<ImageRGB<T>>(offset, w_size, im_size, [&](const Vec2 &pos) {
            return n.get_color_mean_over_footprint(pos, footprint, cm);
        });
    }

    static inline ImageRGB<T> compute_naive_filter_IMG(const Vec2 & w_size,
                                                       const Vec2 &im_size,
                                                       const TextureNoise<T> &n,
                                                       const Color_map<T> &cm)
    {
        Vec2 offset(T(n.width()) * T(0.5), T(n.height()) * T(0.5));
        offset -= w_size * T(0.5);

        Vec2 t_size(n.width(), n.height());
        Vec2 footprint(w_size(0) / t_size(0), w_size(1) / t_size(1));

        return computeIMG<ImageRGB<T>>(offset, w_size, im_size, [&](const Vec2 &pos) {
            T mean = n.get_noise_mean_over_footprint(pos, footprint);
            return cm.map(mean, 0);
        });
    }

    static inline ImageRGB<T> compute_good_filter_IMG(const Vec2 & w_size,
                                                      const Vec2 &im_size,
                                                      const TextureNoise<T> &n,
                                                      const Color_map<T> &cm)
    {
        Vec2 offset(T(n.width()) * T(0.5), T(n.height()) * T(0.5));
        offset -= w_size * T(0.5);

        Vec2 t_size(n.width(), n.height());
        Vec2 footprint(w_size(0) / t_size(0), w_size(1) / t_size(1));

        return computeIMG<ImageRGB<T>>(offset, w_size, im_size, [&](const Vec2 &pos) {
            T mean = n.get_noise_mean_over_footprint(pos, footprint);
            T squared_mean = n.get_squared_noise_mean_over_footprint(pos, footprint);
            T sigma = std::sqrt(squared_mean - mean * mean);
            return cm.map(mean,sigma);
        });
    }
};

} //end namespace ASTex

#endif // PNF_H
