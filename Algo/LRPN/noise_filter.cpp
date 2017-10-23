
#include <ASTex/special_io.h>
#include <ASTex/fourier.h>
#include <ASTex/bilateral_filters.h>
#include <ASTex/utils.h>
#include <ASTex/colorspace_filters.h>
#include <stdio.h>
#include <time.h>

namespace ASTex
{


/**
 * @brief RPnoise
 * @param inputfile
 * @param path_out (must finish with / and exist)
 * @param size
 * @param size_fft
 * @param sig_freq
 */
void RPnoise(const std::string& inputfile, const std::string& path_out,int size, int size_fft, float sig_freq, int nb_step_filt)
{
	using IMG_DBL = ImageRGBd::ItkImg ;
	std::string name_file = IO::remove_ext(IO::remove_path(inputfile));

	srand(time(NULL));

	// LOAD INPUT
	ImageRGBd input_color;
	input_color.load(inputfile);
	// first transform [0,255] double -> [0,1] double
	ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::Pointer filter0 =
			ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::New();
	filter0->SetInput(input_color.itk());
	filter0->Update();
	input_color.itk() = filter0->GetOutput();

	ImageGrayd input_grey;
	IO::load_RGB_2_luminance(input_grey, inputfile);

	int color_width = input_color.width();
	int color_height = input_color.height();

	int FFT_SIZE = size_fft;
	float value_z = sig_freq;
	int filter_size = size;

	double sigma_dist = filter_size;
	double sigma_freq = value_z;


	Fourier::LocalSpectrum lsp (input_grey,FFT_SIZE,'p');
	lsp.welch_post_loading((FFT_SIZE*2)-1, 1);

	int step =0;
//	std::string gpath = path_out+name_file+"_Filtersize_"+std::to_string(size)+"_FFTsize_"+std::to_string(FFT_SIZE)+"_SigFreq_"+std::to_string(value_z)+"_STEP_";
	std::string gpath = path_out+name_file+"_filtered_";

	// Local compute function
	auto compute = [&] (ImageRGBd& in_color, ImageRGBd& out_color) -> void
	{
		Fourier::frequency_bilateral_filter(in_color, out_color, filter_size, lsp , sigma_dist, sigma_freq);
		if (step == nb_step_filt)
			IO::save(out_color, gpath+"S.png",1,0,0,0);
		ImageRGBd RES(color_width,color_height,true);

		// SOUSTRACT
		RES.for_all_pixels([&] (ImageRGBd::PixelType& R, int i, int j)
		{
			R = input_color.pixelAbsolute(i,j) - out_color.pixelAbsolute(i,j);
		});

		ImageRGBd RES_final(color_width,color_height,false);
		// ADD 0.5
		RES_final.for_all_pixels([&RES] (ImageRGBd::PixelType& R, int i, int j)
		{
			R[0] = clamp_scalar(RES.pixelAbsolute(i,j)[0]+0.5,0,1);
			R[1] = clamp_scalar(RES.pixelAbsolute(i,j)[1]+0.5,0,1);
			R[2] = clamp_scalar(RES.pixelAbsolute(i,j)[2]+0.5,0,1);
		});

		if (step == nb_step_filt)
			IO::save(RES_final, gpath+"N.png",1,0,0,0);
	};


	ImageRGBd test_color(input_grey.width(),input_grey.height(),true);
	ImageRGBd loop_c(input_grey.width(),input_grey.height(),true);

	lsp.compute_local_guidance_map(filter_size);

	compute(input_color, test_color);
	step++;
	//for (int i = 1; i < 10; ++i)
	while (step <= nb_step_filt)
	{
		compute(test_color, loop_c);
		step++;
		compute(loop_c,test_color);
		step++;
	}
}

}
