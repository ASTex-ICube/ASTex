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


//**********************************



	//EGALISATION HISTOGRAMME 
	// ImageRGBd input_egalised;
	// input_egalised.initItk(input.width(),input.height());
	// egalisation_histogramme(input_color,input_egalised,inputfile);



	// Gaussian de l'input 
		// ImageRGBd Input_blured;

		//CREATION DU BLUR
		// Gaussian_filter(input_color,Input_blured,size_filter,input_color.width()/4);
		// Input_blured.save(path_blured+"Gaussian_"+inputfile);

		//LOADING DE L'IMAGE
	  	// std::string path_blured ="/home/guingo/Images/Test_methods/result_noise/blured/result_strutextract/";
	  	// Input_blured.load(path_blured+inputfile);
		
	  	// Soustract_img_plus_average(input_color,Input_blured,inputfile);
	  	// Soustract_img_plus(input_color,Input_blured,127,inputfile);
	
	


	// COLLECTION OF IMAGES
	// image_collector collec;
	// collec.add(input, mean); // add input to collection

	// COMPUTE FFT
	ImageGrayd modulus;
	ImageGrayd phase;

	welch(input, modulus, 8);

	// fftForwardModulusAndPhase(input, modulus, phase);

	//Plus grande taille image multiple de 2 
	const int W = input.width();
	const int H = input.height();

    const int M = std::min(W,H);
    int T = 1;
    while ((2*T) < M)
    {
        T*=2;
    }

	// collec.add(modulus); // add modulus to collection

	// monitorStats(modulus, "modulus"); // look stats


  		 
	//BIDOUILLE

	//incrementation successive de la taille du patch / du seuil 

	// double distance = compute_square_size_v2(modulus,0.05);
	// 
	// std::string savefile = "input_egalised_"+inputfile;	
		
	// input_egalised.save (path_out+savefile);



	float distance=1;
	
	// PERIODO TXT
	std::map<double, int> map;
	// std::map<double, double> values_map = print_periodogram_CIRCLE(modulus,inputfile);
	std::map<double, double> values_map = print_periodogram_reverse_size(modulus,inputfile,MIN(input.width(),input.height()) );
	std::vector<double> tab_mode;
	
	// //MODE SEEKING
	for (int i = 1; i < 10000; ++i){
		tab_mode = compute_mode_kmean(values_map, 4);
		
		for(int j=1;j <tab_mode.size(); j++)
		map[tab_mode[j]]++;
	}

	// PRINT MAP
	std::set<double> seeds;
	std::map<double, int>::iterator it = map.begin();
	
	for (it=map.begin(); it!=map.end(); ++it){
		if(it->second < 10)
			map.erase(it);
		else{
	    std::cout << it->first << " => " << it->second << '\n';
			seeds.insert(it->first);
		}

	}
	std::cout << " MAX OCCURE POST RECTIF" << distance << std ::endl;

	distance = compute_mode_meanshift(values_map,seeds,5);


//  TROVUER DISTANCE VERSION 1 ****************************
	// std::vector<double> index;
	// std::vector<double> values;



	// //Split des valeurs
	// for (std::pair<double,double> entry : values_map) {
 //    	index.push_back(entry.first);
 //    	values.push_back(entry.second);
 //    }

 //    float max_combo = 0;

 //    for (it=map.begin(); it!=map.end(); ++it){
    
 //    //FIND LE PLUS PROCHE
 //    int k = 0;
 //    while(index [k]< it->first )
 //    	k++;

 //    int closest ;
    
 //    //Check entre le suivant et l'actuel, garde la plus petite distance
 //    // std::cout << " CLOSEST " << index [k] << " ou " << index [k-1] << " de " << it->first << std ::endl;
 //    ( std::fabs(index[k] - it->first) < std::fabs(index[k-1] - it->first) )? closest = k : closest = k-1;
    	
 //    // std::cout << " CLOSEST " << index [closest] << "bla " << values[closest]<<  std ::endl;

 //    //Ponderation empirique de départ
 //    // if(index [closest] * values [closest] *it->second  > max_combo){
 //    // 	max_combo = index [closest] * values [closest] *it->second ;
 //    // 	distance = it->first;
 //    // }

 //    //Pour chaque seed resultat, regarder les values 
 //    int filtre_size = 3; 
 //    float actual_combo = 0;

 //    for(int i = -filtre_size/2; i <= filtre_size/2; i++){
	//     if(i+closest >= 0 && i+closest < index.size() ){
	// 		actual_combo += values[ i+closest];    	
	//     }
	
 //    }
 //    if (actual_combo > max_combo){
 //     	distance = index[closest];//it->first;
 //     	max_combo = actual_combo;
 //    }
    

	
	// }

    		    // std::cout << entry.first << " => " << entry.second << '\n';

//***********************************************************************
	
  	double periode = 1/distance;
  	
  	std:: cout <<"DISTANCE "<< distance<<	" patch size "<< T*periode <<" "<<T*periode <<" plus grand 2^ "<< T<<  std::endl; 
	
	int patch_size = std::round(T*periode);
		
    	// Drawing result
	  draw_result(modulus,inputfile,distance,patch_size);

  	// QUILTING 
		// std:: cout <<	"patch size * "<< i <<  std::endl; 
  		for (int i = 1; i < 10; ++i){

  		quilting(patch_size*i,inputfile);
  		}

  		quilting(patch_size    ,inputfile);
  		quilting(patch_size*1.5,inputfile);
  		quilting(patch_size*2  ,inputfile);
  		quilting(patch_size*3.5,inputfile);
  		quilting(patch_size*5  ,inputfile);

  	// OUTPUT PERIODO
  	 std::cout << "entreé" << std::endl;

	gnuplot_output(inputfile, (float)(distance));

 	std::cout << "sortie" << std::endl<< std::endl;



	// // RANDOM PHASE
	// phase.initItk(modulus.width(),modulus.height());
	// randomPhase(phase);
	// monitorStats(phase, "phase"); // look stats

	// // INVERSE FFT
	// ImageGrayd result;
	// fftInverseModulusAndPhase(modulus, phase, result);
	// collec.add(result,mean); // add result to collection

	// monitorStats(result, "Output (space)"); // look stats
	
	// // SAVE
	// save(collec.collect(),path_out+inputfile, 0.0);

	// save(result,"/home/guingo/Images/Test_methods/result_noise/result.png", 0.0);
	// modulus.save("/home/guingo/Images/Test_methods/result_noise/FFT.png");
	// std::cout << "TRAITEMENT DE  " <<inputfile<< "SUCCESS" <<std::endl<< std::endl;


	return EXIT_SUCCESS;
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
