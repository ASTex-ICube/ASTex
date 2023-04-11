//
// Created by grenier on 03/04/23.
//

#ifndef ASTEX_INPUT_TERRAIN_H
#define ASTEX_INPUT_TERRAIN_H

#include <iostream>
#include <ASTex/image_rgb.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkRGBToLuminanceImageFilter.h"

namespace ASTex{

    class Input_Terrain{
//        using PIXT = typename IMG::PixelType;
//        using EPIXT = typename IMG::DoublePixelEigen;

        const ImageGrayu8& img_input_; // terrain d'entrée
        int input_height_;
        int input_width_;
        int output_height_;
        int output_width_;

    private:
        double fetch_mean(const Eigen::Vector2i& pix_coord) // pix_coord dans l'image de sortie
        {
            double rapport_x = double(input_width_)/output_width_; // combien de pixel d'entrée pour un pixel de sortie
            double rapport_y = double(input_height_)/output_height_;

            double average = 0.;// = eigenPixel<double>(0.);
//            set_zero(average);
            int nb = 0;

            for_indices(rapport_x*pix_coord[0],rapport_x*(pix_coord[0] + 1), rapport_y*pix_coord[1], rapport_y*(pix_coord[1]+1), [&] (int x, int y)
            {
                average += eigenPixel<double>(img_input_.pixelAbsolute(x, y));
                nb ++;
            });
            average /= nb;

            return average;
        }


        double fetch_grad_x(const Eigen::Vector2i& pix_coord) // pix_coord dans l'image de sortie
        {
            int xph = pix_coord[0] >= output_width_-1 ? 0 : 1;
            int xmh = pix_coord[0] <= 0 ? 0 : 1;
            double F_xph = fetch_mean(Eigen::Vector2i{pix_coord[0]+xph, pix_coord[1]});
            double F_xmh = fetch_mean(Eigen::Vector2i{pix_coord[0]-xmh, pix_coord[1]});

            double delta_h = (xph+xmh) * (1./output_width_);
            double diff = (F_xph-F_xmh)/delta_h;

            double saturation = 0.5;
            diff = std::clamp(saturation*diff + 127.5, 0., 255.);
//            std::cout<<diff<<std::endl;

            return diff;///255.;
        }

        double fetch_grad_y(const Eigen::Vector2i& pix_coord) // pix_coord dans l'image de sortie
        {
            int xph = pix_coord[1] >= output_height_-1 ? 0 : 1;
            int xmh = pix_coord[1] <= 0 ? 0 : 1;
            double F_xph = fetch_mean(Eigen::Vector2i{pix_coord[0], pix_coord[1]+xph});
            double F_xmh = fetch_mean(Eigen::Vector2i{pix_coord[0], pix_coord[1]-xmh});

            double delta_h = (xph+xmh) * (1./output_height_);
            double diff = (F_xph-F_xmh)/delta_h;

            double saturation = 0.5;
            diff = std::clamp(saturation*diff + 127.5, 0., 255.);
//            std::cout<<diff<<std::endl;

            return diff;///255.;
        }



    public:
        Input_Terrain(const ImageGrayu8& input, int out_height, int out_width):
        img_input_(input), output_height_(out_height), output_width_(out_width)
        {
            input_height_ = input.height();
            input_width_ = input.width();
        }




        double terrain_pixel(const Eigen::Vector2i& pix_coord) // récupère un pixel du terrain
        {
            double mean = fetch_mean(pix_coord); // fetch
            return mean;
        }

        double gradx_pixel(const Eigen::Vector2i& pix_coord)
        {
            double gradx = fetch_grad_x(pix_coord);
            return gradx;
        }

        double grady_pixel(const Eigen::Vector2i& pix_coord)
        {
            double grady = fetch_grad_y(pix_coord);
            return grady;
        }


        void terrain_img(ImageGrayu8& img_out, ImageGrayu8& out_gradX, ImageGrayu8& out_gradY) // récupère le terrain pour chaque pixel de l'image de sortie
        {
            img_out.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
                                            {
                                                Eigen::Vector2i pix_coord{x,y};//{ double(x) / (img_out.width()), double(y) / (img_out.height()) };
                                                P = terrain_pixel(pix_coord);
                                            });
            out_gradX.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
                                            {
                                                Eigen::Vector2i pix_coord{x,y};//{ double(x) / (img_out.width()), double(y) / (img_out.height()) };
                                                P = gradx_pixel(pix_coord);
                                            });
            out_gradY.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y)
                                              {
                                                  Eigen::Vector2i pix_coord{x,y};//{ double(x) / (img_out.width()), double(y) / (img_out.height()) };
                                                  P = grady_pixel(pix_coord);
                                              });
        }
    };



    Input_Terrain grab_Input_terrain(const ImageGrayu8& input, int out_height, int out_width){
        return Input_Terrain(input, out_height, out_width);
    }




}
#endif //ASTEX_INPUT_TERRAIN_H
