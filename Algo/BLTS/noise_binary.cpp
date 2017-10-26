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
int binarize_one(std::string inputfile, std::string outputfile, double ratio)
{
	std::string delimiter = ".";
    std::string name_file = inputfile.substr(0, inputfile.find(delimiter));

	 srand48(time(NULL));

	std::cout << "TRAITEMENT DE  " <<inputfile<< std::endl;
	
	// std::string path_in ="/home/guingo/Images/Test_methods/sample/"; 
	// std::string path_out ="/home/guingo/Images/Test_methods/result_noise/";


	std::string path_in_b ="/home/guingo/Images/Test_methods/bench_article/"; 

	ImageGrayd input;
	load2lightness(inputfile,input);
	
	ImageGrayd output;
	output.initItk(input.width(),input.height(),true);	
	binarize(input, output, ratio);

	// std::vector<ImageGrayd&> input;
	// input.resize( nb_image );	

	// for (int i = 0; i < nb_image; ++i){
	// 	load2gray(inputfile+to_string(i),input));
	// }


	print_image(output, path_in_b,outputfile, 1 );


}


int main( int argc, char ** argv )
{
    if( argc < 5 )
    { 
	    std::cerr << "Usage: " << std::endl;
	    std::cerr << argv[0] << " <input grey 1> <output bin 1>  <input grey 2> <output bin 2> [ <input grey> <output bin> ]*" << std::endl;

	    return EXIT_FAILURE;
    }

    const int nb_maps = (argc - 1) / 2;
    
    std::vector<ImageGrayd> input;
    input.resize(nb_maps);

    std::vector<ImageGrayd> output;
    output.resize(nb_maps);

    std::vector<std::string> output_names;

    std::string delimiter = ".";
    for (int i=0; i<nb_maps; ++i)
    {
    	input[i].load(argv[2*i+1]);
	    std::string name_file = std::string(argv[2*i+2]).substr(0, std::string(argv[2*i+2]).find(delimiter));
    	output_names.push_back(name_file);
    } 

	binarize_distance_maps(input,output);

    for (int i=0; i<nb_maps; ++i)
    {
        print_image(output[i], "./", output_names[i], 255.0);
    }
    return 0;
}
