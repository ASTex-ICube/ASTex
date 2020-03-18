#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include<cmath>
#include<string>
#include <iostream>
#include <fstream>

namespace ASTex {

template <typename IMG>
class Histogram{
public:
    using DataType = typename IMG::DataType;
    using PixelType = typename IMG::PixelType;
private :
    DataType *histo;
    unsigned int nb_bins;

public:
    Histogram() : histo(nullptr),nb_bins(0) {}
    ~Histogram(){delete [] histo;}

    void computeHisto(const IMG &img, const unsigned int &n)
    {
        nb_bins = n;
        delete [] histo;
        histo = new DataType[nb_bins * IMG::NB_CHANNELS];

        for (unsigned int i=0;i<nb_bins;i++) {
            for (unsigned int channel =0; channel < IMG::NB_CHANNELS; channel++) {
                histo[i * IMG::NB_CHANNELS + channel] = DataType(0);
            }
        }

        DataType *pix;
        for (unsigned int channel =0;channel < IMG::NB_CHANNELS; channel++) {
            img.for_all_pixels([&](const PixelType & p){
                PixelType tmp = p;
                pix = reinterpret_cast<DataType*>(&tmp);
                int index = static_cast<int>(std::floor(pix[channel] * (nb_bins-1)));
                histo[index * IMG::NB_CHANNELS + channel]++;
            });
        }

        unsigned int nb_pixels = img.width() * img.height();
        for (unsigned int i=0;i<nb_bins;i++) {
            for (unsigned int channel =0; channel < IMG::NB_CHANNELS; channel++) {
                histo[i * IMG::NB_CHANNELS + channel] /= nb_pixels;
            }
        }

    }

    void exportHisto(const std::string &filename,const float ymax){
        std::ofstream fd;
        fd.open(filename + ".txt");

        if(fd.is_open()){
            for (unsigned int channel =0; channel < IMG::NB_CHANNELS; channel++) {
                for (unsigned int i =0; i < nb_bins; i++) {
                    //if(histo[i](channel) != 0)
                        fd << DataType(i)/DataType(nb_bins-1) << "\t" << histo[i * IMG::NB_CHANNELS + channel] << std::endl;
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
