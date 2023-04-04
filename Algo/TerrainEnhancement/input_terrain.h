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
#include "itkDerivativeImageFilter.h"
#include "itkRGBToLuminanceImageFilter.h"

namespace ASTex{

    template<typename IMG>
    class Input_Terrain{
        using PIXT = typename IMG::PixelType;
        using EPIXT = typename IMG::DoublePixelEigen;

        const IMG& img_input_; // terrain d'entrée
        int input_height_;
        int input_width_;
        int output_height_;
        int output_width_;

    protected:
        inline void set_zero(Eigen::Vector2d& v)
        {
            v.setZero();
        }

        inline void set_zero(Eigen::Vector3d& v)
        {
            v.setZero();
        }

        inline void set_zero(Eigen::Vector4d& v)
        {
            v.setZero();
        }

        inline void set_zero(double& v)
        {
            v = 0.0;
        }



    public:
        Input_Terrain(const IMG& input, int out_height, int out_width):
        img_input_(input), output_height_(out_height), output_width_(out_width)
        {
            input_height_ = input.height();
            input_width_ = input.width();
        }


        EPIXT fetch_mean(const Eigen::Vector2i& pix_coord) // pix_coord dans l'image de sortie
        {
            double rapport_x = double(input_width_)/output_width_; // combien de pixel d'entrée pour un pixel de sortie
            double rapport_y = double(input_height_)/output_height_;

            EPIXT average;// = eigenPixel<double>(0.);
            set_zero(average);
            int nb = 0;

            for_indices(rapport_x*pix_coord[0],rapport_x*(pix_coord[0] + 1), rapport_y*pix_coord[1], rapport_y*(pix_coord[1]+1), [&] (int x, int y)
            {
                average += eigenPixel<double>(img_input_.pixelAbsolute(x, y));
                nb ++;
            });
            average /= nb;

            return average;
        }


        PIXT terrain_pixel(const Eigen::Vector2i& pix_coord) // récupère un pixel du terrain
        {
            EPIXT mean = fetch_mean(pix_coord); // fetch
            return IMG::itkPixel(mean);
        }


        void terrain_img(IMG& img_out) // récupère le terrain pour chaque pixel de l'image de sortie
        {
            img_out.parallel_for_all_pixels([&] (typename IMG::PixelType& P, int x, int y)
                                            {
                                                Eigen::Vector2i pix_coord{x,y};//{ double(x) / (img_out.width()), double(y) / (img_out.height()) };
                                                P = terrain_pixel(pix_coord);
                                            });
        }
    };



    template<typename IMG>
    Input_Terrain<IMG> grab_Input_terrain(const IMG& input, int out_height, int out_width){
        return Input_Terrain<IMG>(input, out_height, out_width);
    }



    // Compute gradient images for a given input image
// - image: image to be processed
// - imageGX: horizontal gradients
// - imageGY: vertical gradients
    void compute_Gradient_terrain(ImageRGBu8& image, ImageGrayu8& imageGX, ImageGrayu8& imageGY)
    {
        typedef itk::Image< float, 2 > FloatImageType;
        typedef itk::Image< uint8_t , 2 > UnsignedCharImageType;
        typedef itk::RGBToLuminanceImageFilter< ImageRGBu8::ItkImg, UnsignedCharImageType > ToGrayFilterType;
        typedef itk::DerivativeImageFilter< UnsignedCharImageType, FloatImageType > DerivativeFilterType;
        typedef itk::RescaleIntensityImageFilter< FloatImageType, UnsignedCharImageType > NormalizeFilterType;

        // Convert image to gray
        ToGrayFilterType::Pointer toGrayFilter = ToGrayFilterType::New();
        toGrayFilter->SetInput(image.itk());

        // Create and setup derivative filter along X
        DerivativeFilterType::Pointer derivativeFilterX = DerivativeFilterType::New();
        derivativeFilterX->SetDirection(0);
        derivativeFilterX->SetInput(toGrayFilter->GetOutput());
        NormalizeFilterType::Pointer normalizerX = NormalizeFilterType::New();
        normalizerX->SetOutputMinimum(0);
        normalizerX->SetOutputMaximum(255);
        normalizerX->SetInput(derivativeFilterX->GetOutput());
        normalizerX->Update();
        //imageGX.itk()->CopyInformation(normalizerX->GetOutput());
        imageGX.itk() = ImageGrayu8(normalizerX->GetOutput()).itk();


        // Create and setup derivative filter along Y
        DerivativeFilterType::Pointer derivativeFilterY = DerivativeFilterType::New();
        derivativeFilterY->SetDirection(1);
        derivativeFilterY->SetInput(toGrayFilter->GetOutput());
        NormalizeFilterType::Pointer normalizerY = NormalizeFilterType::New();
        normalizerY->SetOutputMinimum(0);
        normalizerY->SetOutputMaximum(255);
        normalizerY->SetInput(derivativeFilterY->GetOutput());
        normalizerY->Update();
        //imageGY.itk()->CopyInformation(normalizerY->GetOutput());
        imageGY.itk() = ImageGrayu8(normalizerY->GetOutput()).itk();

    }


}
#endif //ASTEX_INPUT_TERRAIN_H
