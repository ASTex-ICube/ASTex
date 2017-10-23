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
int Synthesis( std::string inputfile,const ImageGrayd N_input,const ImageGrayd S_input,const std::vector<ImageGrayd> crops, const std::vector<ImageGrayd> guidance)
{
	std::string delimiter = ".p";
    std::string name_file = inputfile.substr(0, inputfile.find(delimiter));

	srand48(time(NULL));

	std::cout << "TRAITEMENT DE  " <<inputfile<< std::endl;

	std::string path_in_b ="/home/guingo/Images/Test_methods/bench_article/Synthesis_LRP"; 
	system(std::string("mkdir -p "+path_in_b).c_str());
	path_in_b ="/home/guingo/Images/Test_methods/bench_article/Synthesis_LRP/"+name_file+"/";
	system(std::string("mkdir -p "+path_in_b).c_str());

	// const int im_w = size;
	// const int im_h = size;

	// int m_left_bnd = -(im_w-1) / 2;
 //   	int m_right_bnd = (im_h) / 2;


 //   	std::vector<std::vector< double>> W;
	// W.resize( im_w , std::vector< double>( im_h ) );

	std::vector<int> nb_pix_cluster;
	nb_pix_cluster.resize(guidance.size());

	// std::vector<ImageGrayd> modulus;
 //    modulus.resize(guidance.size());

 //    std::vector<ImageGrayd> phase;
 //    phase.resize(guidance.size());

    ImageGrayd noise_localised;
    noise_localised.initItk(N_input.width(),N_input.height(),true);
        ImageGrayd input_synthesis;
    input_synthesis.initItk(N_input.width(),N_input.height(),true);

    std::vector<ImageGrayd> noise_synthesis;
    noise_synthesis.resize(guidance.size());

    for (int k=0; k<guidance.size(); ++k){

		// LOAD INPUT
		ImageGrayd input = crops[k];

		// COLLECTION OF IMAGES
		image_collector collec;
		collec.add(input, 0.0); // add input to collection

		// COMPUTE FFT
		ImageGrayd modulus;
		ImageGrayd phase;
		fftForwardModulusAndPhase(input,modulus,phase);	
		welch(input, modulus, 1);
		collec.add(modulus); // add modulus to collection
	//	monitorStats(modulus, "modulus"); // look stats

		double prop [8] = {0,0.01,0.02, 0.05, 0.1, 0.2, 0.5, 1};
	// 	for (int i=0; i<8; ++i)
	// 	{
	// 		// RANDOM PHASE
	// 		// ImageGrayd phase;
	// 		// phase.initItk(modulus.width(),modulus.height(),true);
	
	// //		randomPhase(phase,mask_smallest_values(modulus,1-prop[i]));
	// 		randomPhase(phase,mask_largest_values(modulus,prop[i]));
	// //		monitorStats(phase, "phase"); // look stats
	// 		collec.add(phase,0.0); // add result to collection

	// 		// INVERSE FFT
	// 		ImageGrayd result;
	// 		fftInverseModulusAndPhase(modulus, phase, result);
	// 		monitorStats(result, "Output (space)"); // look stats
	// 		collec.add(result,0.0230799); // add result to collection
	// 	}
		for (int i=7; i>=0; --i)
		{
			// RANDOM PHASE
			fftForwardModulusAndPhase(input,modulus,phase);	

			// ImageGrayd phase;
			// phase.initItk(modulus.width(),modulus.height(),true);
			randomPhase(phase,mask_smallest_values(modulus,1-prop[i]));
	//		randomPhase(phase,mask_largest_values(modulus,prop[i]));
	//		monitorStats(phase, "phase"); // look stats
			collec.add(phase,0.0); // add result to collection

			// INVERSE FFT
			ImageGrayd result;
			fftInverseModulusAndPhase(modulus, phase, result);
	//		monitorStats(result, "Output (space)"); // look stats
			collec.add(result,0); // add result to collection
		}
		// SAVE
		save(collec.collect(),path_in_b+std::to_string(k)+".png", 0.0230799);
    }

}

int main( int argc, char ** argv )
{
    if( argc < 9 )
    { 
	    std::cerr << "Usage: " << std::endl;
	    std::cerr << argv[0] << "<input_name> <input_N><input_S> <size1> <X position patch1><Y position patch1> <size2><X position patch2><Y position patch2> <guidance grey 1>  <guidance grey 2> [ <guidance grey> ]*" << std::endl;

	    return EXIT_FAILURE;
    }

    const int nb_maps = (argc - 10) ;


    std::string filename_source = argv[1];
    
    ImageGrayd N_input; 	
    load(N_input,argv[2]);

    ImageGrayd S_input; 	
    load(S_input,argv[3]);

    std::vector<ImageGrayd> crops;
    crops.resize(2);


    for (int i=0; i<2; ++i){ 
    	int size = atoi(argv[(i*3)+4]);
    	int im_padding_left = (size-1)/2;
    	int im_padding_right = size/2;
    	
    	int x = atoi(argv[(i*3)+5]);
    	int y = atoi(argv[(i*3)+6]);

    	crops[i].initItk(size,size,true);

	  	crop_image(N_input, crops[i], x-im_padding_left, x+im_padding_right, y-im_padding_left, y+im_padding_right);

	} 

    std::vector<ImageGrayd> guidance;
    guidance.resize(nb_maps);

	  for (int i=0; i<nb_maps; ++i)
    {
        load(guidance[i],argv[i+10],0,2);
	} 
	Synthesis(filename_source,N_input,S_input,crops,guidance);

	return 0;
	// return RPnoise("/home/guingo/Images/Test_methods/sample/bright_black.png","/home/guingo/Images/Test_methods/result_noise/result_bright_black.png");

	// return RPnoise("/home/guingo/Images/Test_methods/sample/lave2.png","/home/guingo/Images/Test_methods/result_noise/result_lave2.png");

//		return testPhase(TEMPO_PATH+"refnoise.png","/tmp/result.png");
}
