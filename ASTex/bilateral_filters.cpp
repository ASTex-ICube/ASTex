
#include <fstream>
#include <cmath>
#include <algorithm>

#include "bilateral_filters.h"

//includes pour la convert d'espace de couleur
#include "colorspace_filters.h"
#include "itkCastImageFilter.h"

#include "utils.h"
#include "fourier.h"


namespace ASTex
{



void ASTEX_API bilateral_filter(const ImageGrayd& input, ImageGrayd& output, int filtre_size, double id, double cd)
{
	input.for_all_pixels([&] (int i, int j)
	{
		double sumWeight = 0;
		double sum_g = 0;

		for(int y = -filtre_size/2; y <= filtre_size/2; y++)
			for(int x = -filtre_size/2; x <= filtre_size/2; x++)
			{
				if(i+x >= 0 && j+y>=0 && i+x < input.width()&& j+y < input.height() )
				{
					double currWeight;
					// define bilateral filter kernel weights
					float imageDist = (float)((x)*(x) + (y)*(y)) ;
					float colorDist =   (input.pixelAbsolute(i,j) - input.pixelAbsolute(i+x,j+y)) * (input.pixelAbsolute(i,j) - input.pixelAbsolute(i+x,j+y)) ;
					currWeight = std::exp(-(imageDist/(id*id*2)))*std::exp(-(colorDist/(cd*cd*2)) );
					sumWeight += currWeight;
					sum_g += currWeight*input.pixelAbsolute(i+x,j+y);
				}
			}

		sum_g /= sumWeight;
		output.pixelAbsolute(i,j) = (sum_g);
	});
}


void ASTEX_API bilateral_filter(const ImageRGBd& input, ImageRGBd& output, int filtre_size, double id, double cd)
{
	input.for_all_pixels([&] (int i, int j)
	{

		double sumWeight = 0;
		double sum_r = 0;
		double sum_g = 0;
		double sum_b = 0;

		for(int y = -filtre_size/2; y <= filtre_size/2; y++)
			for(int x = -filtre_size/2; x <= filtre_size/2; x++)
			{
				if(i+x >= 0 && j+y>=0 && i+x < input.width()&& j+y < input.height() )
				{

					double currWeight;

					// define bilateral filter kernel weights
					float imageDist = (float)((x)*(x) + (y)*(y)) ;

					float colorDist =   (input.pixelAbsolute(i,j)[0] - input.pixelAbsolute(i+x,j+y)[0]) * (input.pixelAbsolute(i,j)[0] - input.pixelAbsolute(i+x,j+y)[0]) +
							(input.pixelAbsolute(i,j)[1] - input.pixelAbsolute(i+x,j+y)[1]) * (input.pixelAbsolute(i,j)[1] - input.pixelAbsolute(i+x,j+y)[1]) +
							(input.pixelAbsolute(i,j)[2] - input.pixelAbsolute(i+x,j+y)[2]) * (input.pixelAbsolute(i,j)[2] - input.pixelAbsolute(i+x,j+y)[2]) ;

					currWeight = std::exp(-(imageDist/(id*id*2)))*std::exp(-(colorDist/(cd*cd*2)) );
					sumWeight += currWeight;

					sum_r += currWeight*input.pixelAbsolute(i+x,j+y)[0];
					sum_g += currWeight*input.pixelAbsolute(i+x,j+y)[1];
					sum_b += currWeight*input.pixelAbsolute(i+x,j+y)[2];

				}
			}

		sum_r /= sumWeight;
		sum_g /= sumWeight;
		sum_b /= sumWeight;

		output.pixelAbsolute(i,j)[0] = ((sum_r));
		output.pixelAbsolute(i,j)[1] = ((sum_g));
		output.pixelAbsolute(i,j)[2] = ((sum_b));

	});
}


namespace Fourier
{


void ASTEX_API frequency_bilateral_filter(const ImageGrayd& input, ImageGrayd& output, int filtre_size, LocalSpectrum& lsp , double id, double cd)
{
	int left_step = -(filtre_size-1) / 2;
	int right_step = (filtre_size) / 2;

	//lsp.compute_local_guidance_map(filtre_size);
	if (lsp.get_used_filter_size()!= filtre_size)
		std::cerr << "Warning in frequency_bilateral_filter lsp not fully computed"<< std::endl;

	input.for_all_pixels([&] (int i, int j)
	{
		float sumWeight = 0;
		output.pixelAbsolute(i,j)= 0;
		double sum_g = 0;
		for(int y = left_step; y <= right_step; y++)
		{
			for(int x = left_step; x <= right_step; x++)
			{
				if(i+x >= 0 && j+y >=0 && i+x < input.width() && j+y < input.height())
				{
					double currWeight;

					//define bilateral filter kernel weights
					float imageDist = (float)((x)*(x) + (y)*(y)) ;

					float freqDist = lsp.local_guidance(i,j).pixelAbsolute(x-left_step,y-left_step);
					currWeight = std::exp(-(imageDist/(id*id*2)));
					currWeight *= std::exp(-(freqDist/(cd*cd*2)));
					sumWeight += currWeight;
					sum_g += currWeight*input.pixelAbsolute(i+x,j+y);
				}
			}
		}
		sum_g /= sumWeight;
		output.pixelAbsolute(i,j) = sum_g;
	});
}




void ASTEX_API frequency_bilateral_filter(const ImageRGBd& input, ImageRGBd& output, int filtre_size, LocalSpectrum& lsp , double id, double cd)
{
	int left_step = -(filtre_size-1) / 2;
	int right_step = (filtre_size) / 2;

	// CHECK LSP
	if (lsp.get_used_filter_size()!= filtre_size)
		std::cerr << "Warning in frequency_bilateral_filter lsp not fully computed"<< std::endl;

	input.for_all_pixels([&] (int i, int j)
	{
		double sumWeight = 0;
		RGBd sum_col(0,0,0);

		for(int y = left_step; y <= right_step; y++)
		{
			for(int x = left_step; x <= right_step; x++)
			{
				if(i+x >= 0 && j+y >=0 && i+x < input.width() && j+y < input.height())
				{

					double currWeight;

					//define bilateral filter kernel weights
					float imageDist = (float)((x)*(x) + (y)*(y)) ;
					float freqDist = lsp.local_guidance(i,j).pixelAbsolute(x-left_step,y-left_step);

					currWeight = std::exp(-(imageDist/(id*id*2)));
					currWeight *= std::exp(-(freqDist/(cd*cd*2)));

					sumWeight += currWeight;
					sum_col += input.pixelAbsolute(i+x,j+y) * currWeight;
				}
			}
		}
		sum_col /=sumWeight;
		output.pixelAbsolute(i,j) = sum_col;
	});

}




void ASTEX_API frequency_joint_bilateral_filter(const ImageGrayd& input, ImageGrayd& output,const std::vector<ImageGrayd>& guidance_maps, int filtre_size, double id, double cd)
{
	int left_step = -(filtre_size-1) / 2;
	int right_step = (filtre_size) / 2;

// // FREQ BF
	input.for_all_pixels([&] (int i, int j)
	{
		float sumWeight = 0;
		output.pixelAbsolute(i,j)= 0;
		double sum_g = 0;

		for(int y = left_step; y <= right_step; y++)
		{
			for(int x = left_step; x <= right_step; x++)
			{
				if(i+x >= 0 && j+y >=0 && i+x < input.width() && j+y < input.height())
				{
					double currWeight;

					//define bilateral filter kernel weights
					float imageDist = (float)((x)*(x) + (y)*(y)) ;

					//Calcul dist
					float sum = 0;
					for (uint32_t w = 0; w < guidance_maps.size(); ++w)
					{
						sum += (guidance_maps[w].pixelAbsolute(i+x,j+y) - guidance_maps[w].pixelAbsolute(i,j)) * (guidance_maps[w].pixelAbsolute(i+x,j+y) - guidance_maps[w].pixelAbsolute(i,j)) ;
					}

					float freqDist = std::sqrt(sum) ;
					 currWeight = std::exp(-(imageDist/(id*id*2)));
					currWeight *= std::exp(-(freqDist/(cd*cd*2)));
					sumWeight += currWeight;
					sum_g += currWeight*input.pixelAbsolute(i+x,j+y);
				}
			}
		}
	   sum_g /= sumWeight;
	  output.pixelAbsolute(i,j) = sum_g;
	});
 }


} //namespace Fourier
} //namespace ASTex
