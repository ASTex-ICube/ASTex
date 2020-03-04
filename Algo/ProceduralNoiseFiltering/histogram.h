#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <Eigen/Eigen>
#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

namespace ASTex {

template <typename IMG>
class Histogram2{
public:
    using HistoType = Eigen::Matrix<typename IMG::DataType,IMG::NB_CHANNELS,1>;
private :
    HistoType *histo;
    int nb_bins;

    HistoType conversionPixel(const typename ImageGray<typename IMG::DataType>::PixelType & p)
    {
        HistoType ret;
        ret(0) = p;
        return ret;
    }

    HistoType conversionPixel(const typename ImageRGB<typename IMG::DataType>::PixelType & p)
    {
        HistoType ret;
        ret(0) = p.GetRed();
        ret(1) = p.GetGreen();
        ret(2) = p.GetBlue();
        return ret;
    }

public:
	Histogram2() : histo(nullptr),nb_bins(0) {}
	~Histogram2(){delete [] histo;}

    void computeHisto(const IMG &img, const int &n)
    {
        nb_bins = n;
        delete [] histo;
        histo = new HistoType[nb_bins];

        for (int i=0;i<nb_bins;i++) {
            for (int channel =0; channel < IMG::NB_CHANNELS; channel++) {
                histo[i](channel) = 0;
            }
        }

        for (int channel =0;channel < IMG::NB_CHANNELS; channel++) {
            img.parallel_for_all_pixels([&](const typename IMG::PixelType & p){
                HistoType pix = conversionPixel(p);
                int index = static_cast<int>(std::floor(pix(channel) * (nb_bins-1)));
                histo[index](channel)++;
            });
        }

        int nb_pixels = img.width() * img.height();
        for (int i=0;i<nb_bins;i++) {
            for (int channel =0; channel < IMG::NB_CHANNELS; channel++) {
                histo[i](channel) /= nb_pixels;
            }
        }

    }

    void exportHisto(const std::string &filename,const float ymax){
        std::ofstream fd;
        fd.open(filename + ".txt");

        if(fd.is_open()){
            for (int channel =0; channel < IMG::NB_CHANNELS; channel++) {
                for (int i =0; i < nb_bins; i++) {
                    //if(histo[i](channel) != 0)
                        fd << float(i)/float(nb_bins-1) << "\t" << histo[i](channel) << std::endl;
                }

                fd << std::endl << std::endl;
            }

            fd.close();
        }
        else {
            std::cout << "cannot open " << filename << ".txt" << std::endl;
        }

        std::ofstream fd2;
        fd2.open(filename + "_script.txt");

        if(fd2.is_open()){
            fd2 << "set terminal svg size 512,512 enhanced fname \"arial\" butt solid" << std::endl;
            fd2 << "set output '"<< filename << ".svg'" << std::endl;
            fd2 << "set yrange [0:"<< ymax <<"]" << std::endl;
            fd2 << "set size ratio 0.5" << std::endl;

            if(IMG::NB_CHANNELS == 1){
                fd2 << "plot \""<< filename << ".txt\" notitle lt rgb \"blue\" with lines"<< std::endl;
            }

            if(IMG::NB_CHANNELS == 3){
                fd2 << "plot \""<< filename << ".txt\" index 0 notitle lt rgb \"red\" with lines,";
                fd2 << " \"" << filename << ".txt\" index 1 notitle lt rgb \"green\" with lines,";
                fd2 << " \"" << filename << ".txt\" index 2 notitle lt rgb \"blue\" with lines" << std::endl;
            }
            fd2.close();
        }
        else {
            std::cout << "cannot open " << filename << "_script.txt" << std::endl;
        }

    }




};

}

#endif // HISTOGRAM_H
