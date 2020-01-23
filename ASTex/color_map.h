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
public:
    Color_map() {}
    void add_color(const int &pos,const Color &col) {
        palette.insert(std::pair<int, Color>(pos,col));
    }

    typename ImageRGB<T>::PixelType map(const T &x)
    {
        T pas = T(1) / T(palette.rbegin()->first);

        for(auto it = palette.begin();it!=palette.end();++it)
        {
            T inf = T(it->first) * pas;
            T sup = T((++it)->first) * pas;
            it--;
            if (x >= inf && x <= sup) {
                Color c_inf = it->second;
                Color c_sup = (++it)->second;
                //it--;
                return ImageRGB<T>::itkPixel((1 - x) * c_inf + x * c_sup);
            }
        }

        return ImageRGB<T>::itkPixel(0);
    }

    void export_palette()
    {

    }
};

}

#endif // COLOR_MAP_H
