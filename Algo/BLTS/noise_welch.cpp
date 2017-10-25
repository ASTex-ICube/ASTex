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




#include "io.h"
#include "fourier.h"
#include "itkJPEGImageIO.h"
#include "mode_seeking.h"
#include "output_geoffrey.h"
#include "image_treatment.h"
#include "utilities.h"
#include "tests.h"
#include "colorspace_filters.h"


using namespace ASTex;



/**
 * @brief RPnoise random phase noise
 * @param inputfile a grayscale input example
 * @param outputfile is a patchwork image containing input + spectrum + output (pure random phase)
 * @return error code
 */
int RPnoise(std::string inputfile)
{
	std::string delimiter = ".";
    std::string name_file = inputfile.substr(0, inputfile.find(delimiter));

	 srand48(time(NULL));

	std::cout << "TRAITEMENT DE  " <<inputfile<< std::endl;
	
	std::string path_in ="/home/guingo/Images/Test_methods/sample/"; 
	std::string path_out ="/home/guingo/Images/Test_methods/result_noise/";

		// LOAD INPUT
    ImageGrayd input;
	load2lightness(inputfile, input);

// ***************************************
	// std::string path_in_b ="/home/guingo/Images/Test_methods/bench_article/"+name_file+"/";
		std::string path_in_b ="/home/guingo/Images/Test_methods/bench_article/"; 
	// system(std::string("mkdir -p "+path_in_b).c_str());

	// ImageGrayd test;
	// test.initItk(input.width(),input.height(),true);

	// ImageGrayd loop;
	// loop.initItk(input.width(),input.height(),true);


	// ImageRGBd test_color;
	// test_color.initItk(input.width(),input.height());

	// ImageRGBd loop_c;
	// loop_c.initItk(input.width(),input.height());

	int welchsize = 63, welchstep = 4, fft_size=32;
			// int welchsize = 15, welchstep = 4, fft_size=8;


	ImageGrayd test_before;
	test_before.initItk(input.width(),input.height(),true);
	ImageGrayd test_after;
	test_after.initItk(input.width(),input.height(),true);
	ImageGrayd test_constructor;
	test_constructor.initItk(input.width(),input.height(),true);

	ImageGrayd test_dist;
	test_dist.initItk(input.width(),input.height(),true);

	ImageGrayd test_dist_mosa;
	test_dist_mosa.initItk(input.width(),input.height(),true);


		


	local_spectrum lsp;
	new (&lsp) local_spectrum(input,fft_size,'g');
	lsp.mosaic(test_before);
	lsp.welch_post_loading(welchsize, welchstep);
	lsp.mosaic(test_after);

	// local_spectrum lsp_w;
	// new (&lsp_w) local_spectrum (input,fft_size,welchsize,welchstep);
	// lsp_w.mosaic(test_constructor);
	
	double v_max_before = percentile (test_before,0.8);
	double v_max_after = percentile (test_after,0.8);
	// double v_max_constructor = percentile (test_constructor,0.80);

	print_image(test_before,path_in_b,"MOSA_WELCH_BEFORE" ,v_max_before);
	print_image(test_after,path_in_b,"MOSA_WELCH_AFTER" ,v_max_after);
	// print_image(test_constructor,path_in_b,"MOSA_WELCH_CONSTRUCTOR" ,v_max_constructor);


	const int dx = input.width()/4;
	const int dy = input.height()/4;
	

						for (int pix_x = dx; pix_x < input.width(); pix_x+=dx){
							for (int pix_y = dy; pix_y < input.height(); pix_y+=dy){
//										lsp.distance_map_to_pixel(test_dist, pix_x, pix_y);
				
								double v_max = percentile (test_dist,0.8);
									std::string pix_loc= "_pix_"+std::to_string(pix_x)+"_"+std::to_string(pix_y);
								print_image(test_dist,path_in_b,name_file+"_wsize_"+std::to_string(lsp.get_windowsize())+ pix_loc+"_VMAX_"+std::to_string(v_max) ,v_max);
								

							}
						}

	lsp.distance_mosaic (test_dist_mosa,32);
	double v_max = percentile (test_dist_mosa,0.99);
	std::cout << " vmax : " << v_max << std::endl; 
	print_image(test_dist_mosa,path_in_b,"Dist_Mosa" ,v_max);


}

int main( int argc, char ** argv )
{
    if( argc < 2 )
    { 
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " Source" << std::endl;
    return EXIT_FAILURE;
    }

    std::string filename_source = argv[1];



	return RPnoise(filename_source);
	// return RPnoise("/home/guingo/Images/Test_methods/sample/bright_black.png","/home/guingo/Images/Test_methods/result_noise/result_bright_black.png");

	// return RPnoise("/home/guingo/Images/Test_methods/sample/lave2.png","/home/guingo/Images/Test_methods/result_noise/result_lave2.png");

//		return testPhase(TEMPO_PATH+"refnoise.png","/tmp/result.png");
}
