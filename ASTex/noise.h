#ifndef NOISE_H
#define NOISE_H

#include <ASTex/utils.h>
#include <random>
#include <ASTex/color_map.h>

namespace ASTex {

template <typename T>
class Noise{
private:
    int nb_cosines;
    T *phases,*frequences,*orientations;

public :
    using Vec2 = Eigen::Matrix<T,2,1>;
    using Color = typename Color_map<T>::Color;

    Noise(const int &n, const T &fr_min, const T &fr_max) : nb_cosines(n) {
        phases = new T[nb_cosines];
        frequences = new T[nb_cosines];
        orientations = new T[nb_cosines];

        const T Pi = 4 * std::atan(T(1));

        std::random_device rd;
        std::mt19937 gen(rd());
        //[0;2Pi[
        std::uniform_real_distribution<T> phases_dis(T(0), T(2) * Pi);
        //[10;20[
        std::uniform_real_distribution<T> frequences_dis(fr_min, fr_max);
        //[0;Pi[
        std::uniform_real_distribution<T> orientations_dis(T(0), Pi);

        for (int i = 0; i < n; ++i) {
            phases[i] = phases_dis(gen);
            frequences[i] = frequences_dis(gen);
            orientations[i] = orientations_dis(gen);
        }
    }

    Noise() : Noise(32,T(10),T(20)) {}

    ~Noise(){
        delete[] phases;
        delete[] frequences;
        delete[] orientations;
        phases = frequences = orientations = nullptr;
    }

private:
    template<typename T2,typename init, typename func>
    T2 get_over_footprint(const Vec2 &pos, const Vec2 &footprint,const int & nb_sample,const init &initialisation, const func &f) const{

        T2 ret = initialisation();
        T step(T(1)/T(nb_sample));
        Vec2 corner = pos - footprint * T(0.5);

        for (int i =0; i < nb_sample; ++i) {
            T y = corner(1) + i * step;
            for(int j=0; j < nb_sample; ++j) {
                T x = corner(0) + j * step;

                T value_noise = basic2D(Vec2(x,y)) ;
                ret += f(value_noise);
            }
        }
        ret /= nb_sample * nb_sample;

        return ret;
    }

public:

    T basic2D(const Vec2 &pos) const{
           T sum_cosines(0);
           T ph,fr,o;
           for (int i = 0; i < nb_cosines; ++i) {
               ph = phases[i];
               fr = frequences[i];
               o = orientations[i];
               sum_cosines += std::cos(ph + fr * std::cos(o) * pos(0) + fr * std::sin(o) * pos(1) );
           }

           sum_cosines *= T(1)/ T(6) * std::sqrt(T(2) / T(nb_cosines));
           sum_cosines += T(0.5);
           sum_cosines = clamp_scalar(sum_cosines,T(0),T(1));

           return sum_cosines;
    }

    T get_noise_mean_over_footprint(const Vec2 &pos, const Vec2 &footprint, const int &nb_sample) const{
        return get_over_footprint<T>(pos, footprint, nb_sample, [](){return T(0);}, [](const T&x) {
            return x;
        });
    }

    T get_squared_noise_mean_over_footprint(const Vec2 &pos, const Vec2 &footprint, const int &nb_sample) const{
        return get_over_footprint<T>(pos, footprint, nb_sample, [](){return T(0);}, [](const T&x) {
            return x * x;
        });
    }

    Color get_color_mean_over_footprint(const Vec2 &pos,
                                        const Vec2 &footprint,
                                        const Color_map<T> &cm,
                                        const int &nb_sample) const{
        return get_over_footprint<Color>(pos, footprint, nb_sample, [](){return Color(0,0,0);}, [&](const T&x) {
            return cm.map(x, 0);
        });
    }

};
}

#endif
