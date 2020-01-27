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
private:
    std::map<int, Color > palette;
    T step = T(0);

    ImageRGB<T> filtered;
    T sigma_max;
    int width, height;

    T gauss(const T &x, const T &y, const T &mu, const T &sigma) const
    {
        const T Pi = 4 * std::atan(T(1));

        T ret = std::exp(- (std::pow(x-mu,2) + std::pow(y-mu,2)) / (T(2) * std::pow(sigma,2)));
        ret /= sigma * sigma * T(2) * Pi;

        return ret;
    }

    T numeric_integration_gauss(const T&a, const T&b, const int &n, const T&mu, const T& sigma){
        T sum(0);
        for(int i= 1; i <= n-1; ++i)
            for (int j=1; j<=n-1 ;++j)
                sum += gauss(a + T(i) * (b - a) / T(n * n), a + T(j) * (b - a) / T(n * n), mu, sigma);
        sum += (gauss(a, a, mu, sigma) + gauss(b, b, mu, sigma)) / T(2);
        sum *= (b-a) / T(n);

        return sum;
    }

    Color numeric_integration_col_gauss(const T&a, const T&b, const int &n, const T&mu, const T& sigma){
        Color sum(0,0,0);
        //T sum_gauss = numeric_integration_gauss(a,b,n,mu,sigma);
        for(int i= 1; i <= n-1; ++i){
            T x = a + T(i) * (b - a) / T(n);
            for (int j =1;j <= n-1; ++j){
                T y = a + T(j) * (b - a) / T(n);
                sum += map((x+y)/T(2)) * gauss(x, y, mu, sigma);
            }
        }
        sum += (map(a) * gauss(a, a, mu, sigma) + map(b) * gauss(b, b, mu, sigma)) / T(2);
        sum *= (b-a) / T(n);

        return sum;
    }

public:
    Color_map() {}
    void add_color(const int &pos,const Color &col) {
        palette.insert(std::pair<int, Color>(pos,col));
        step = T(1) / T(palette.rbegin()->first);;
    }

    ImageRGB<T> get_filtered() const {
        return filtered;
    }

    Color map(const T &x) const {
        Color ret;
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
        return ret;
    }

    typename ImageRGB<T>::PixelType map(const T &f, const T &sigma) const {
        int y = f * height;
        int x = sigma * width / sigma_max;

        return filtered.pixelAbsolute(x,y);
    }

    void filter(const int &w, const int &h,const int &nb_bins, const T &sm){
        width = w;
        height = h;
        sigma_max = sm;
        Color * tab = new Color[nb_bins];
        for(int i=0;i<nb_bins;++i){
            tab[i] = map(T(i)/T(nb_bins));
        }

        filtered = ImageRGB<T>(width,height);
        filtered.parallel_for_all_pixels([&](typename ImageRGB<T>::PixelType &p,int i,int j)
        {
            T f = T(j) / T(h);
            T sigma = T(i) / T(w) * sigma_max;
            Color c(0,0,0);

            if(sigma != 0){
                //assert(compare_scalar(1.0,numeric_integration_gauss(0,1,nb_bins,f,sigma)));
                c = numeric_integration_col_gauss(0,1,nb_bins,f,sigma);
            }
            else
                c = map(f);

            p = ImageRGB<T>::itkPixel(c);

        });

        delete[] tab;

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
