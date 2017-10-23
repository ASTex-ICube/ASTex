
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
int Filter( std::string inputfile,const int size, int size_fft, float sig_freq, const std::vector<ImageGrayd> guidance)
{
    std::string delimiter = ".p";
    std::string name_file = inputfile.substr(0, inputfile.find(delimiter));

	 srand48(time(NULL));

	std::cout << "TRAITEMENT DE  " <<inputfile<< std::endl;
	
	std::string path_in ="/home/guingo/Images/Test_methods/sample/"; 
	std::string path_out ="/home/guingo/Images/Test_methods/result_noise/";

	// LOAD INPUT
    ImageGrayd input;
	// load2luminance(inputfile, input);
//	load2lightness(inputfile, input);
    load2luminance(inputfile, input);

// ***************************************
	std::vector<double> v;
	v.push_back(0.01);
	v.push_back(0.02);
	v.push_back(0.05);
    int FFT_SIZE = size_fft;


	// for (int i = 0; i < 2; ++i)
	// {
	// 	for (int j = 0; j < v.size(); ++j)
	// 	{
	// 		 code 
  		// float value_z = v[j]* std::pow(10,i);
//				float value_z = 0.25;
                float value_z = sig_freq;


    std::string path_in_b ="/home/guingo/Images/Test_methods/bench_article/Joint_Filtered/"+name_file+"/size_"+std::to_string(FFT_SIZE)+"/nb_cluster"+std::to_string(guidance.size())+"/"+std::to_string(value_z)+"/";
	system(std::string("mkdir -p "+path_in_b).c_str());

    std::string path_in_histo ="/home/guingo/Images/Test_methods/bench_article/Joint_Filtered/"+name_file+"/size_"+std::to_string(FFT_SIZE)+"/nb_cluster"+std::to_string(guidance.size())+"/"+std::to_string(value_z)+"/histo/";
	system(std::string("mkdir -p "+path_in_histo+"S/").c_str());
	system(std::string("mkdir -p "+path_in_histo+"N/").c_str());


	ImageGrayd test;
	test.initItk(input.width(),input.height(),true);
	
	test.initItk(input.width(),input.height(),true);

    float filter_size = size;
    float sigma_dist = size/4;
	float sigma_freq = value_z;
    frequency_joint_bilateral_filter(input, test,guidance ,filter_size ,FFT_SIZE, name_file ,sigma_dist,sigma_freq);

		std::cout<<"FILTERING DONE"<<std::endl;
	// 	}
	// }
		

	// print_image(test,"/home/guingo/Images/Test_methods/bench_article/","TestFilter" );


	// for (int window_type = 0; window_type < 3; ++window_type){
	// 	std::cout<<"local_spectrum window_TYPE: "<<window_type<<std::endl;

	// 	for (int window_size = 8; window_size < 64; window_size*=2){
			
	// 		std::cout<<"local_spectrum window_size: "<<window_size<<std::endl;
	// 		switch(window_type){
	// 			case 0 : 	    	//GAUSSIAN
	// 				new (&lsp) local_spectrum(input,window_size,'g');
	//         	break;
	//         	case 1 : 	    	//HAMMING
	// 				new (&lsp) local_spectrum(input,window_size,'h');
	//         	break;
	//         	case 2 : 	    	//PORTE
	// 				new (&lsp) local_spectrum(input,window_size,'p');
	//         	break;
	// 		}
				
	// 		int welchsize = window_size, welchstep = window_size/4;
	// 		// local_spectrum lsp (input,16,welchsize,welchstep);
	// 		std::cout<<"done compute"<<std::endl;

	// 		for(int switch_welch=0; switch_welch <2 ; switch_welch ++){				
				
	// 			std::string welch_type_s = "FFT";
	// 			 if(switch_welch == 1){
	// 			 		lsp.welch_post_loading(welchsize, welchstep);
	// 			 		welch_type_s = "WELCH";
	// 			 }

	// 			const int dx = input.width()/4;
	// 			const int dy = input.height()/4;
	// 			for (int dist_type = 0; dist_type < 2; ++dist_type){
	// 						std::cout<<"local_spectrum dist_type: "<<dist_type<<std::endl;
	// 						std::string dist_type_s;

	// 				for (int pix_x = dx; pix_x < input.width(); pix_x+=dx){
	// 					for (int pix_y = dy; pix_y < input.height(); pix_y+=dy)
	// 					{
	// 						// std::cout<<"distance_map_to_pixel"<<std::endl;
	// 						switch(dist_type){
	// 							case 0 :
	// 								lsp.distance_map_to_pixel(test, pix_x, pix_y);
	// 				        		dist_type_s = "Distsimple";
	// 				        	break;
	// 				        	case 1 :	    	
	// 								lsp.distance_map_to_pixel_voisinage(test, pix_x, pix_y);
	// 				        		dist_type_s = "Distvoisinage";
	// 				        	break;
	// 				        }

	// 						double v_max = percentile (test,0.80);
	// 						// std::string pix_loc= "_pix_"+(pix_x<10?"0":"")+(pix_x<100?"0":"")+std::to_string(pix_x)+"_"+(pix_y<10?"0":"")+(pix_y<100?"0":"")+std::to_string(pix_y);
	// 						std::string pix_loc= "_pix_"+std::to_string(pix_x)+"_"+std::to_string(pix_y);
	// 						print_image(test,path_in_b,name_file+"_"+welch_type_s+"_"+dist_type_s+"_wtype_"+std::string(1,lsp.get_windowtype())+"_wsize_"+std::to_string(lsp.get_windowsize())+ pix_loc ,v_max);

	// 					}
	// 				}
	// 				std::cout<<"local_spectrum MOSAIC DIST: "<<std::endl;
	// 				lsp.distance_mosaic (test,16);
	// 				double v_max = percentile (test,0.99);
	// 				print_image(test,path_in_b,name_file+"_DistanceMosaic"+"_"+welch_type_s+"_"+dist_type_s+"_wtype_"+std::string(1,lsp.get_windowtype())+"_wsize_"+std::to_string(lsp.get_windowsize()) ,v_max);
	// 				std::cout<<"DONE"<<std::endl;
	// 			}
	// 				std::cout<<"local_spectrum MOSAIC SPECTRUM: "<<std::endl;
	// 				lsp.mosaic(test);
	// 				double v_max = percentile (test,0.80);
	// 				print_image(test,path_in_b,name_file+"_Mosaic"+"_"+welch_type_s+"_wsize_"+std::to_string(lsp.get_windowsize()) ,v_max);
	// 				std::cout<<"DONE"<<std::endl;
	// 		}
	// 	}
	// }





	// print_image(test,name_file+"_Mosaic_"+"_wsize_"+std::to_string(lsp.get_windowsize())+"_welchsize_"+std::to_string(welchsize)+"_wstep_"+std::to_string(welchstep) ,v_max);



	//loop integré pour concerver les ftt 
	// frequency_bilateral_filter(input,test,31,8,0.05);

	// bilateral_filter(input,test,8,8,0.05);
// for (int i = 0; i < 500; ++i)
// {
	// bilateral_filter(test,loop,8,8,0.05);
	// bilateral_filter(loop,test,8,8,0.05);
// }

//     color_bilateral_filter(input_color,test_color,8,8,0.05*3);
// for (int i = 0; i < 500; ++i)
// {
// 	color_bilateral_filter(test_color,loop_c,8,8,0.05*3);
// 	color_bilateral_filter(loop_c,test_color,8,8,0.05*3);
// }


	return EXIT_SUCCESS;
}

int main( int argc, char ** argv )
{
    if( argc < 7 )
    { 
	    std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << "<input grey 1> <FILERSIZE> <FFTSIZE> <SIGMA FREQ> <guidance grey 1>  <guidance grey 2> [ <guidance grey> ]*" << std::endl;

	    return EXIT_FAILURE;
    }

    const int nb_maps = (argc - 5) ;


    std::string filename_source = argv[1];
    int size = atoi(argv[2]);
    int size_fft = atoi(argv[3]);
    float sig_freq = atof(argv[4]);


    std::vector<ImageGrayd> guidance;
    guidance.resize(nb_maps);

    for (int i=0; i<nb_maps; ++i)
    {
        load(guidance[i],argv[i+5]);
		std::cout<<guidance[i].pixelAbsolute(50,50)<<std::endl;

    } 

    Filter(filename_source,size,size_fft,sig_freq,guidance);

	return 0;
	// return RPnoise("/home/guingo/Images/Test_methods/sample/bright_black.png","/home/guingo/Images/Test_methods/result_noise/result_bright_black.png");

	// return RPnoise("/home/guingo/Images/Test_methods/sample/lave2.png","/home/guingo/Images/Test_methods/result_noise/result_lave2.png");

//		return testPhase(TEMPO_PATH+"refnoise.png","/tmp/result.png");
}
