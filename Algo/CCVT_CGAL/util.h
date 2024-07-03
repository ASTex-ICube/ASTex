//
// Created by grenier on 12/10/23.
//

#ifndef CCVT_TEST_CHARLINE_UTIL_H
#define CCVT_TEST_CHARLINE_UTIL_H


#define MINN -1e100
#define MAXN  1e100
#define EPS   1e-12
#define M_PI 3.14159265358979323846

#include <cmath>
#include <vector>

template <class T>
T compute_mean(const std::vector<T>& data)
{
    T mean = 0.0;
//    return std::accumulate(data.begin(), data.end(), mean)/data.size();
    for (unsigned i = 0; i < data.size(); ++i)
        mean += data[i];

    return (mean / data.size());
}

template <class T>
void remove_mean(std::vector<T>& data)
{
    T mean = compute_mean(data);
    for (unsigned i = 0; i < data.size(); ++i)
        data[i] -= mean;
}

template <class T>
double compute_norm(std::vector<T>& data)
{
    double norm2 = 0.0;
    for (unsigned i = 0; i < data.size(); ++i)
        norm2 += data[i]*data[i];
    return std::sqrt(norm2);
}

template <class Point, class OutputIterator>
void generate_regular_polygon(unsigned nb,
                              const double radius,
                              const Point center,
                              OutputIterator out)
{
    if (nb < 3) nb = 3;
    double step = 2*M_PI / nb;
    double start = 0.25*M_PI;
    for (unsigned i = 0; i < nb; ++i)
    {
        double angle = start + i*step;
        Point q(radius*cos(angle) + center.x(),
                radius*sin(angle) + center.y());
        *out++ = q;
    }
}

inline double compute_int01_gauss_t(double mu, double sigma){
    return sqrt(M_PI*sigma*sigma/2.) *
            (
                std::erf(mu/sqrt(2.*sigma*sigma)) -
                std::erf(-(1-mu)/sqrt(2.*sigma*sigma))
            );
}

inline double compute_int01_t_gauss_t(double mu, double sigma){
    return mu*sqrt(M_PI*sigma*sigma/2.) *
           (
                   std::erf(mu/sqrt(2.*sigma*sigma)) -
                   std::erf(-(1-mu)/sqrt(2.*sigma*sigma))
           ) +
           sigma*sigma*
           (
                   exp(-mu*mu/(2.*sigma*sigma))-
                   exp(-(1-mu)*(1-mu)/(2.*sigma*sigma))
           );
}

/* produit de deux fonctions gaussienne définie selon la même variable t
 * exp(-(a*t-mu_1)^2 / (2*sig_1^2)) *  exp(-(b*t-mu_2)^2 / (2*sig_2^2)) = A exp(-(t-mu)^2 / (2*sig^2))
 * product_gaussian_amplitude : calcul de A
 * product_gaussian_mean : calcul de mu
 * product_gaussian_variance : calcul de sig^2
 */

inline double product_gaussian_amplitude(double a, double b, double mu_1, double mu_2, double sig_1, double sig_2){
    return exp(-(b*mu_1 - a*mu_2)*(b*mu_1 - a*mu_2)/(2.*(b*b*sig_1*sig_1 + a*a*sig_2*sig_2)));
}

inline double product_gaussian_mean(double a, double b, double mu_1, double mu_2, double sig_1, double sig_2){
    return (a*sig_2*sig_2*mu_1 + b*sig_1*sig_1*mu_2)/(b*b*sig_1*sig_1 + a*a*sig_2*sig_2);
}

inline double product_gaussian_variance(double a, double b, double mu_1, double mu_2, double sig_1, double sig_2){
    return (sig_1*sig_1 * sig_2*sig_2)/(b*b*sig_1*sig_1 + a*a*sig_2*sig_2);
}

#endif //CCVT_TEST_CHARLINE_UTIL_H
