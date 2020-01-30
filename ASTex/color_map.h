#ifndef COLOR_MAP_H
#define COLOR_MAP_H

#include <ASTex/image_rgb.h>
#include <map>
#include <ASTex/utils.h>

namespace ASTex {


template <typename T>
class Color_map
{
public :
    using Color = Eigen::Matrix<T,3,1>;
    using Vec2 = Eigen::Matrix<T,2,1>;
private:
    std::map<int, Color > palette;
    T step = T(0);

    ImageRGB<T> filtered;
    T sigma_max;
    int width, height;

    T gauss1D(const T &x, const T &mu, const T &sigma) const
    {
        const T Pi = 4 * std::atan(T(1));

        T ret = std::exp(-T(0.5) * std::pow((x-mu)/sigma,2));
        ret /= sigma * std::sqrt(2 * Pi);

        return ret;
    }

    T numeric_integration_gauss1D(const T&a, const T&b, const int &n, const T&mu, const T& sigma){
        T sum(0);
        for(int k= 1; k <= n-1; ++k){
            T v =a + T(k) * (b - a) / T(n);
            sum += gauss1D(v, mu, sigma);
        }
        sum += (gauss1D(a, mu, sigma) + gauss1D(b, mu, sigma)) * T(0.5);
        sum *= (b-a) / T(n);

        return sum;
    }

    Color numeric_integration_col_gauss1D(const T&a, const T&b, const int &n, const T&mu, const T& sigma){
        Color sum(0,0,0);

        for(int k= 1; k <= n-1; ++k){
            T v = a + T(k) * (b - a) / T(n);
            sum += map(v) * gauss1D(v, mu, sigma);
        }
        sum += (map(a) * gauss1D(a, mu, sigma) + map(b) * gauss1D(b, mu, sigma)) * T(0.5);
        sum *= (b-a) / T(n);

        return sum;
    }

public:
    Color_map() {}
    void add_color(const int &pos,const Color &col) {
        palette.insert(std::pair<int, Color>(pos,col));
        step = T(1) / T(palette.rbegin()->first);
    }

    ImageRGB<T> get_filtered() const {
        return filtered;
    }

    void set_filtered(const ImageRGB<T> &i,const T &sm) {
        height = i.height();
        width = i.width();
        filtered = i;
        sigma_max = sm;
    }

    Color map(const T &x) const {
        Color ret;
        if(x < T(0))
            ret = palette.begin()->second;
        else if ( x > T(1))
            ret = palette.rbegin()->second;
        else {
            for(auto it = palette.begin();it!=palette.end();++it){
                T inf = T(it->first) * step;
                T sup = T((++it)->first) * step;
                T t = (x - inf) / (sup -inf);
                it--;
                if (x >= inf && x <= sup) {
                    Color c_inf = it->second;
                    Color c_sup = (++it)->second;
                    it--;
                    ret = (1 - t) * c_inf + t * c_sup;
                    break;
                }
            }
        }
        return ret;
    }

    typename ImageRGB<T>::PixelType map(const T &f, const T &sigma) const {
        int y = f * (height-1);
        int x = sigma * (width-1) / sigma_max;

        return filtered.pixelAbsolute(x,y);
    }

    void filter(const int &w, const int &h,const int &nb_bins, const T &sm){
        width = w;
        height = h;
        sigma_max = sm;

        filtered = ImageRGB<T>(width,height);
        filtered.parallel_for_all_pixels([&](typename ImageRGB<T>::PixelType &p,int i,int j)
        {
            T f = T(j) / T(h);
            T s = T(i) / T(w) * sigma_max;
            Color c(0,0,0);

            if(s != 0)
                c = numeric_integration_col_gauss1D(-3*s+f,3*s+f,nb_bins,f,s);
            else
                c = map(f);

            p = ImageRGB<T>::itkPixel(c);

        });

    }

    std::string to_string() const {
        std::stringstream s;
        auto it = palette.begin();
        s << "(" << it->first << " ";
        s << it->second(0) << " " << it->second(1) << " " << it->second(2);
        it++;
        for (;it!=palette.end();++it) {
            s << ", " << it->first << " ";
            s << it->second(0) << " " << it->second(1) << " " << it->second(2);
        }

        s << ")";

        return s.str();
    }

    void export_palette(const std::string &filename) const {
        std::ofstream fd;
        fd.open(filename);
        fd << "set palette defined ";
        fd << to_string();

        fd.close();
    }

    void export_courbe(const std::string &filename = "data.txt") const {
        std::ofstream fd;
        fd.open(filename);
        for(auto it = palette.begin();it!= palette.end();it++)
        {
            fd << it->second[0] << " " << it->second[1] << " " << it->second[2] << std::endl;
        }
        fd << std::endl << std::endl;
        for (int i =0; i < 100; i++) {
            T x( T(i) / T(100));
            Color c = map(x);
            fd << c[0] << " " << c[1] << " " << c[2] << std::endl;
        }

        fd.close();
    }
};

}

#endif // COLOR_MAP_H
