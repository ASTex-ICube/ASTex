/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/



#include "special_io.h"

#include <cmath>

#include<ASTex/colorspace_filters.h>
#include<ASTex/utils.h>

#include <iostream>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCastImageFilter.h>

#include <itkShiftScaleImageFilter.h>
#include <itkClampImageFilter.h>

#include <itkStatisticsImageFilter.h>
#include <itkSquareImageFilter.h>

namespace ASTex
{

double ASTEX_API getMean (const ImageGrayd& input)
{
	typedef itk::StatisticsImageFilter< ImageGrayd::ItkImg > StatisticsImageFilterType;
    typename StatisticsImageFilterType::Pointer stats = StatisticsImageFilterType::New ();

	stats->SetInput(input.itk());
    stats->Update();
    return stats->GetMean();
}


double ASTEX_API getMax (const ImageGrayd& input)
{
    typedef itk::StatisticsImageFilter< ImageGrayd::ItkImg > StatisticsImageFilterType;
    typename StatisticsImageFilterType::Pointer stats = StatisticsImageFilterType::New ();

    stats->SetInput(input.itk());
    stats->Update();
    return stats->GetMaximum();
}

double ASTEX_API getStDev (const ImageGrayd& input)
{
    typedef itk::StatisticsImageFilter< ImageGrayd::ItkImg > StatisticsImageFilterType;
    typename StatisticsImageFilterType::Pointer stats = StatisticsImageFilterType::New ();

    stats->SetInput(input.itk());
    stats->Update();
    return stats->GetSigma();
}




namespace IO
{

void ASTEX_API monitorStats (const ImageGrayd& input, std::string intro)
{
//    std::cout << intro << std::endl;

	typedef itk::StatisticsImageFilter< ImageGrayd::ItkImg > StatisticsImageFilterType;
    typename StatisticsImageFilterType::Pointer stats = StatisticsImageFilterType::New ();

	stats->SetInput(input.itk());
    stats->Update();
    std::cout << "(" << intro <<") \t Mean: " << stats->GetMean() << "\t Min: " << stats->GetMinimum() << "\t Max: " << stats->GetMaximum() << "\t Sum: " << stats->GetSum()<< "\t Std.: " << stats->GetSigma()<< "\t Var.: " << stats->GetVariance()<< std::endl;

    typedef itk::SquareImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > SquareImageFilterType;
    typename SquareImageFilterType::Pointer squared = SquareImageFilterType::New ();
    squared->SetInput(input.itk());
    stats->SetInput(squared->GetOutput());
    stats->Update();
//    std::cout << "(" << intro <<"^2) \t Mean: " << stats->GetMean() << "\t Min: " << stats->GetMinimum() << "\t Max: " << stats->GetMaximum() << "\t Sum: " << stats->GetSum()<< "\t Std.: " << stats->GetSigma()<< "\t Var.: " << stats->GetVariance()<< std::endl;
    std::cout << "() \t energy: " << stats->GetSum()<< "\t Power: " << stats->GetSum() / ( input.width() * input.height() )<< std::endl;


    // TODO : this section has a strange side effect : try it when monitoring modulusAccu in welch : the mean growth twice faster...
//    typedef itk::MultiplyImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg, ImageGrayd::ItkImg > MultiplyFilterType;
//    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();

//    multiplyFilter->SetInput1(input.itk());
//    multiplyFilter->SetInput2(input.itk());
//    stats->SetInput(multiplyFilter->GetOutput());
//    stats->Update();
//    std::cout << "(squared)\t Mean: " << stats->GetMean() << "\t Std.: " << stats->GetSigma() << "\t Min: " << stats->GetMinimum() << "\t Max: " << stats->GetMaximum() << std::endl;
}

void ASTEX_API load (ImageGrayd& image, const std::string& filename, double min, double max)
{
    // load UCHAR
    typedef itk::ImageFileReader<ImageGrayu8::ItkImg> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(filename);

    // cast to DOUBLE
    typedef itk::CastImageFilter<ImageGrayu8::ItkImg,ImageGrayd::ItkImg> CastType;
    CastType::Pointer cast = CastType::New();
    cast->SetInput(reader->GetOutput());

    // rescale, shift, and convert to DOUBLE
    typedef itk::ShiftScaleImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > RescaleType;
    RescaleType::Pointer rescale = RescaleType::New();
    rescale->SetInput(cast->GetOutput());
    cast->Update();
    ImageGrayd im (cast->GetOutput());
    rescale->SetShift(min * 255 / (max-min));
    rescale->SetScale((max-min)/255.0);
    rescale->Update();
    image.itk() = rescale->GetOutput();
}

void ASTEX_API save (const ImageGrayd& img, const std::string& filename, double min, double max)
{
    typedef itk::ShiftScaleImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > RescaleType;
    RescaleType::Pointer rescale = RescaleType::New();
    rescale->SetInput(img.itk());
    rescale->SetShift(-min);
    rescale->SetScale(255.0 / (max-min));

    typedef itk::ClampImageFilter< ImageGrayd::ItkImg, ImageGrayu8::ItkImg> ClampFilterType;
    ClampFilterType::Pointer clampFilter = ClampFilterType::New();
    clampFilter->SetInput(rescale->GetOutput());
    clampFilter->SetBounds(0,255);

    typedef itk::ImageFileWriter<ImageGrayu8::ItkImg> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(filename);
    writer->SetInput(clampFilter->GetOutput());
    writer->Update();
}

double ASTEX_API load (ImageGrayd& image, const std::string& filename, bool shift_mean_to_zero)
{
    // load UCHAR
	typedef itk::ImageFileReader<ImageGrayu8::ItkImg> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);

    // cast to DOUBLE
	typedef itk::CastImageFilter<ImageGrayu8::ItkImg,ImageGrayd::ItkImg> CastType;
    CastType::Pointer cast = CastType::New();
    cast->SetInput(reader->GetOutput());

    // rescale, shift if required, and convert to DOUBLE
	typedef itk::ShiftScaleImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > RescaleType;
    RescaleType::Pointer rescale = RescaleType::New();
    rescale->SetInput(cast->GetOutput());
    cast->Update();
	ImageGrayd im (cast->GetOutput());
	double mean = getMean(im);
    if (shift_mean_to_zero)
    {
        rescale->SetShift(-mean);
    }
    else
    {
        rescale->SetShift(0.0);
    }
    rescale->SetScale(1.0/255.0);
    rescale->Update();
	image.itk() = rescale->GetOutput();
    return mean/255.0;
}

void ASTEX_API save (const ImageGrayd& img, const std::string& filename, double mean_shift)
{
	typedef itk::ShiftScaleImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > RescaleType;
    RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetInput(img.itk());
    rescale->SetShift(mean_shift);
    rescale->SetScale(255.0);

	typedef itk::ClampImageFilter< ImageGrayd::ItkImg, ImageGrayu8::ItkImg> ClampFilterType;
    ClampFilterType::Pointer clampFilter = ClampFilterType::New();
    clampFilter->SetInput(rescale->GetOutput());
    clampFilter->SetBounds(0,255);

	typedef itk::ImageFileWriter<ImageGrayu8::ItkImg> WriterType;
    WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filename);
    writer->SetInput(clampFilter->GetOutput());
    writer->Update();
}


void ASTEX_API load_spectrum (ImageSpectrald& image, const std::string& filename)
{
    load(image,filename, false);
    assert (image.width() == image.height());
    const int T = image.width();

	image.for_all_pixels([&] (double& P)
	{
		P = P*P*T;
	});
//    image.setCenter(T/2,T/2); // done in ImageSpectrald constructor
}

void ASTEX_API save_spectrum (const ImageSpectrald& image, const std::string& filename)
{
    assert (image.width() == image.height());
    const int T = image.width();
	ImageSpectrald result(T);

	for_all_pixels(image,result,[&](double P, double& Q)
	{
		Q = std::sqrt(P/double(T));
	});

    save(result,filename, false);
}


void ASTEX_API load_phase (ImageSpectrald& image, const std::string& filename)
{
    load(image,filename, -M_PI, M_PI);
    assert (image.width() == image.height());
    image.setCenter(image.width()/2,image.height()/2);
}

void ASTEX_API save_phase (const ImageSpectrald& image, const std::string& filename)
{
    assert (image.width() == image.height());
    save(image,filename, -M_PI, M_PI);
}

void ASTEX_API load_coordinates (ImageGrayd& image, const std::string& filename)
{
    load(image,filename, -1, 2);
}

void ASTEX_API save_coordinates (const ImageGrayd& image, const std::string& filename)
{
    save(image,filename, -1, 2);
}

void ASTEX_API load_PCA_coordinates (ImageGrayd& image, const std::string& filename)
{
    load(image,filename, -10.0/7.0, 10.0/7.0);
}

void ASTEX_API save_PCA_coordinates (const ImageGrayd& image, const std::string& filename)
{
    save(image,filename, -10.0/7.0, 10.0/7.0);
}



void ASTEX_API  load_RGB_2_gray(ImageGrayd& output, const std::string& filename)
{
	ImageRGBd input_color;
	input_color.load(filename);
	output.initItk(input_color.width(), input_color.height(),true);

	typedef ImageRGBd::ItkImg IMG_DBL;

	// first transform [0,255] double -> [0,1] double
	ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::Pointer filter0 =
			ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::New();
	filter0->SetInput(input_color.itk());
	filter0->Update();


	input_color.itk()=filter0->GetOutput() ;

	output.initItk(input_color.width(),input_color.height(),true);

	for (int i = 0; i < input_color.width(); ++i)
		for (int j = 0; j < input_color.height() ; ++j)
		{
			if(input_color.pixelAbsolute(i,j)[0]  == input_color.pixelAbsolute(i,j)[0])
			{
				output.pixelAbsolute(i,j)= input_color.pixelAbsolute(i,j)[0] ;
			}
			else
			{
				output.pixelAbsolute(i,j)=0;
			}
		}
}

void ASTEX_API  load_RGB_2_luminance(ImageGrayd& output, const std::string& filename)
{
	typedef ImageRGBd::ItkImg IMG_DBL;

	ImageRGBd input_color;
	input_color.load(filename);

	const int im_w = input_color.width();
	const int im_h = input_color.height();
	output.initItk(im_w, im_h,true);


	// first transform [0,255] double -> [0,1] double
	ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::Pointer filter0 =
			ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::New();
	filter0->SetInput(input_color.itk());

	// RGB double -> XYZ float
	ColorSpace::FilterRGBtoXYZ<IMG_DBL,IMG_DBL>::Pointer filter1 =
			ColorSpace::FilterRGBtoXYZ<IMG_DBL,IMG_DBL>::New();
	filter1->SetInput(filter0->GetOutput());
	filter1->Update();

	ImageRGBd input_XYZ;
	input_XYZ.itk()= filter1->GetOutput() ;

	for (int i = 0; i < im_w; ++i)
		for (int j = 0; j < im_h ; ++j)
		{
			if(input_XYZ.pixelAbsolute(i,j)[1]  == input_XYZ.pixelAbsolute(i,j)[1])
			{
				output.pixelAbsolute(i,j)= input_XYZ.pixelAbsolute(i,j)[1] ;
			}
			else
			{
				output.pixelAbsolute(i,j)=0;
			}
		}
}

void ASTEX_API  load_RGB_2_lightness(ImageGrayd& output, const std::string& filename)
{
	typedef ImageRGBd::ItkImg IMG_DBL;

	ImageRGBd input_color;
	input_color.load(filename);

	const int im_w = input_color.width();
	const int im_h = input_color.height();
	output.initItk(im_w, im_h,true);


	// first transform [0,255] double -> [0,1] double
	ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::Pointer filter0 =
			ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::New();
	filter0->SetInput(input_color.itk());

	// RGB double -> XYZ float
	ColorSpace::FilterRGBtoXYZ<IMG_DBL,IMG_DBL>::Pointer filter1 =
			ColorSpace::FilterRGBtoXYZ<IMG_DBL,IMG_DBL>::New();
	filter1->SetInput(filter0->GetOutput());

	// XYZ float -> LUV float
	ColorSpace::FilterXYZtoLUV<IMG_DBL,IMG_DBL>::Pointer filter2 =
			ColorSpace::FilterXYZtoLUV<IMG_DBL,IMG_DBL>::New();
	filter2->SetInput(filter1->GetOutput());
	filter2->Update();

	ImageRGBd input_LUV;
	input_LUV.itk()= filter2->GetOutput() ;

	for (int i = 0; i < im_w; ++i)
		for (int j = 0; j < im_h ; ++j)
		{
			if(input_LUV.pixelAbsolute(i,j)[0]  == input_LUV.pixelAbsolute(i,j)[0])
			{
				output.pixelAbsolute(i,j)= input_LUV.pixelAbsolute(i,j)[0] / 100.0 ;
			}
			else
			{
				output.pixelAbsolute(i,j)=0;
			}
		}
}



void ASTEX_API  save(const ImageRGBd& input, const std::string&  name, double n,  double mean_shift_0,  double mean_shift_1,  double mean_shift_2)
{
	ImageRGBu8 output;
	output.initItk(input.width(),input.height(),true);

	for_all_pixels(input,output,[&](const ImageRGBd::PixelType& P, ImageRGBu8::PixelType& Q)
	{
		Q[0] = clamp_scalar((P[0]+mean_shift_0)*255/n,0.0,255.0);
		Q[1] = clamp_scalar((P[1]+mean_shift_1)*255/n,0.0,255.0);
		Q[2] = clamp_scalar((P[2]+mean_shift_2)*255/n,0.0,255.0);

	});
	output.save(name);
}


} // namespace
} // namespace
