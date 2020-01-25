#ifndef NOISE_H
#define NOISE_H

#include <ASTex/utils.h>
#include <random>

namespace ASTex {

template <typename T>
class Noise{
private:
    int nb_cosines;
    T *phases,*frequences,*orientations;

public :
    Noise(const int &n, const T &fr_min, const T &fr_max) : nb_cosines(n) {
        phases = new T[nb_cosines];
        frequences = new T[nb_cosines];
        orientations = new T[nb_cosines];

        const T Pi = 4 * std::atan(T(1));

        std::random_device rd;
        std::mt19937 gen(0);
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

    T basic2D(const T &x, const T &y){
        T sum_cosines(0);
        T p,f,o;
        for (int i = 0; i < nb_cosines; ++i) {
            p = phases[i];
            f = frequences[i];
            o = orientations[i];
            sum_cosines += std::cos(p + f * std::cos(o) * x + f * std::sin(o) * y );
        }

        sum_cosines *= T(1)/ T(6) * std::sqrt(T(2) / T(nb_cosines));
        sum_cosines += T(0.5);
        sum_cosines = clamp_scalar(sum_cosines,T(0),T(1));

        return sum_cosines;
    }
};
}

#endif
