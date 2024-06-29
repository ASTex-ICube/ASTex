//
// Created by grenier on 12/10/23.
//

#ifndef CCVT_TEST_CHARLINE_PGM_H
#define CCVT_TEST_CHARLINE_PGM_H



#include <vector>
#include <string>
#include <fstream>

class PGMImage
{
private:
    unsigned m_w, m_h;
    double m_max_value;
    std::vector<double> m_value;
    double m_max_original;
    std::vector<double> m_original;

public:
    PGMImage()
    {
        m_w = 0;
        m_h = 0;
        m_max_value = 0.0;
        m_max_original = 0.0;
    }

    PGMImage(unsigned width, unsigned height, bool zero_init = true):
    m_w(width), m_h(height)
    {
        m_max_value = 0.0;
        m_max_original = 0.0;

        double init_value;
        if(zero_init){init_value=0.;}

        m_value.resize(width*height, init_value);
        m_original.resize(width*height, init_value);
    }

    unsigned width () const { return m_w; }

    unsigned height() const { return m_h; }

    bool isNull() const { return m_value.empty(); }

    double max_value() const { return m_max_value; }

    unsigned get_index(unsigned i, unsigned j) const { return (j*width() + i); }

    double pixel(unsigned i, unsigned j) const
    {
        unsigned index = get_index(i, j);
        return m_value[index];
    }

    // set a gaussian density from average and variance
    bool set_data(double moy_x, double moy_y, double variance_x, double variance_y, unsigned size_x, unsigned size_y, double max_val){
        m_original.clear();
        m_h = size_y;// 128; //256;
        m_w = size_x;// 128; //256;
        m_max_original = max_val;// 255.;

        for (unsigned i = 0; i < m_w; ++i)
        {

            double y = double(i)/m_h;
            double Gy = std::exp(-(y-moy_y)*(y-moy_y)/(2.*variance_y));

            for (unsigned j = 0; j  < m_h; j++){
                double x = double(j)/m_w;
                double Gx = std::exp(-(x-moy_x)*(x-moy_x)/(2.*variance_x));

                double value = m_max_original*Gx*Gy;
                m_original.push_back(value);
            }
        }
        tonemap(0.5);
        return true;
    }

    bool load(const std::string& filename)
    {
        m_original.clear();
        std::ifstream input(filename.c_str());

        if(!input.is_open())
        {
            std::cout << "failed to open " << filename << '\n';
            return false;
        }
        else
        {
            std::string header;

            // première ligne = version
            std::getline(input, header);
            if(header.compare("P2") != 0)
            {
                std::cerr << "Version error, should be P2 pgm file" << std::endl;
                return false;
            }

            // deuxième ligne = commentaire
            std::getline(input, header);
//            std::cout << "comment : " << header << std::endl;

            input >> m_w >> m_h;
            input >> m_max_original;
            for (unsigned i = 0; i < m_w*m_h; ++i)
            {
                double value;
                input >> value;
                m_original.push_back(value);
            }

            input.close();
        }

        tonemap(0.5);
        return true;
    }


    void save(const std::string& filename) const
    {
        std::ofstream output(filename.c_str());
        output << "P2" << std::endl;
        output << m_w << " " << m_h << std::endl;
        output << "255" << std::endl;
        for (unsigned i = 0; i < m_w*m_h; ++i)
        {
            int value = static_cast<int>(m_value[i]);
            output << value << std::endl;
        }
        output.close();
    }

    void tonemap(double key)
    {
//        std::cout << " tonemapp with key = "<< key << std::endl;
        m_value = m_original;
        m_max_value = m_max_original;

        double maxv = 0.0;
        double maxo = 0.0;
        for (unsigned i = 0; i < m_value.size(); ++i)
        {
            maxv = std::max(maxv, m_value[i]);
            maxo = std::max(maxo, m_original[i]);
        }
//        std::cout << "MaxV: " << maxv << " ; " << m_max_value << std::endl;
//        std::cout << "MaxO: " << maxo << " ; " << m_max_original << std::endl;
    }
};


#endif //CCVT_TEST_CHARLINE_PGM_H
