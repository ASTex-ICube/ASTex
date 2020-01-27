#ifndef COLOR_MAP_H
#define COLOR_MAP_H

#include <ASTex/image_rgb.h>
#include <map>

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

        for(auto it = palette.begin();it!=palette.end();++it){
            T inf = T(it->first) * step;
            T sup = T((++it)->first) * step;
            T t = (x - inf) / (sup -inf);
            it--;
            if (x >= inf && x <= sup) {
                Color c_inf = it->second;
                Color c_sup = (++it)->second;
                //it--;
                return (1 - t) * c_inf + t * c_sup;
            }
        }

        return Color(0,0,0);
    }

    typename ImageRGB<T>::PixelType map(const T &f, const T &sigma) const {
        int y = f * height;
        int x = sigma * width / sigma_max;

        return filtered.pixelAbsolute(x,y);
    }

    void filter(const int &w, const int &h,const int &nb_bins, const T &sm){
        const T Pi = 4 * std::atan(T(1));
        width = w;
        height = h;
        sigma_max = sm;
        Color * tab = new Color[nb_bins];
        for(int i=0;i<nb_bins;++i){
            tab[i] = map(T(i)/T(nb_bins));
        }

        filtered = ImageRGB<T>(width,height);
        filtered.for_all_pixels([&](typename ImageRGB<T>::PixelType &p,int i,int j)
        {
            T f = T(j) / T(h);
            T sigma = T(i) / T(w) * sigma_max;
            Color c(0,0,0);

            if(sigma != 0){
                for(int k = 0; k < nb_bins;k++){
                    T v = T(k) / T(nb_bins);
                    T gauss = std::exp(-std::pow((v-f)/sigma,2)/T(2));
                    gauss /= sigma * std::sqrt(2 * Pi);
                    c += tab[k] * gauss;
                }
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
