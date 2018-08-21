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
//    std::string path_out ="/home/guingo/Images/Test_methods/bench_article/standard_bilateral_filter/";
    std::string path_out ="/home/guingo/Images/Test_methods/bench_article/GREY_INPUT/";

    system(std::string("mkdir -p "+path_out).c_str());
//    path_out ="/home/guingo/Images/Test_methods/bench_article/standard_bilateral_filter/"+name_file+"/";


    system(std::string("mkdir -p "+path_out).c_str());

	// LOAD INPUT
		
    ImageGrayd input;

//	ImageGrayd modulus;
//	ImageGrayd phase;
    ImageRGBd input_color;
    input_color.load(path_in+inputfile);

    typedef ImageRGBd::ItkImg IMG_DBL;
    typedef ImageRGBf::ItkImg IMG_FLT;
    typedef ImageRGBu8::ItkImg IMG_U8;

    // first transform [0,255] double -> [0,1] double
    ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::Pointer filter0 =
            ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::New();
    filter0->SetInput(input_color.itk());
            filter0->Update();

    input_color.itk()=filter0->GetOutput() ;

    load2luminance(inputfile, input);
    print_image(input,path_out,name_file+"_LUMINANCE" );

//	fftForwardModulusAndPhase(input,modulus, phase);

//	print_image(modulus,path_out,"fft_noise2",1);
	
// ***************************************
// std::string path_in_b ="/home/guingo/Images/Test_methods/";
// system(std::string("mkdir -p "+path_out).c_str());

	// for (int i = 0; i < input.width(); ++i)
	// {
	// 	for (int j = 0; j < input.height(); ++j)
	// 	{
	// 		input.pixelAbsolute(i,j) = 0.5+ ((i/double(input.width()))*(input.pixelAbsolute(i,j)-0.5));
	// 	}
	// }
	// print_image(input,path_in_b,"Var_intensite",1);

	// ImageGrayd loop;
	// loop.initItk(input.width(),input.height(),true);


//    ImageRGBd test_color;
//    test_color.initItk(input_color.width(),input_color.height(),true);

//    ImageRGBd loop_c;
//    loop_c.initItk(input_color.width(),input_color.height(),true);

// 	local_spectrum* lsp;
// lsp = new local_spectrum(input,8,'g');

	// ImageGrayd test;
	// 	test.initItk(input.width(),input.height(),true);


// 	ImageGrayd mosa;


// 	for (int window_type = 0; window_type < 3; ++window_type){
// 		std::cout<<"local_spectrum window_TYPE: "<<window_type<<std::endl;

// 		for (int window_size = 8; window_size < 64; window_size*=2){
// 			test.initItk(input.width(),input.height(),true);

// 			std::cout<<"local_spectrum window_size: "<<window_size<<std::endl;
// 			switch(window_type){
// 				case 0 : 	    	//GAUSSIAN
// 					delete lsp;
// 					lsp = new local_spectrum(input,window_size,'g');
// 	        	break;
// 	        	case 1 : 	    	//HAMMING
// 	        	delete lsp;
// 					lsp = new local_spectrum(input,window_size,'h');
// 	        	break;
// 	        	case 2 : 	    	//PORTE
// 					delete lsp;
// 					lsp = new local_spectrum(input,window_size,'p');
// 	        	break;
// 			}
				
// 			int welchsize = (window_size*2)-1, welchstep = 4;
// 			// local_spectrum lsp (input,16,welchsize,welchstep);
// 			std::cout<<"done compute"<<std::endl;

// 			for(int switch_welch=0; switch_welch <2 ; switch_welch ++){				
				
// 				std::string welch_type_s = "FFT";
// 				 if(switch_welch == 1){
// 				 		lsp->welch_post_loading(welchsize, welchstep);
// 				 		welch_type_s = "WELCH";
// 				 }

// 				const int dx = input.width()/4;
// 				const int dy = input.height()/4;
// 				for (int dist_type = 0; dist_type < 2; ++dist_type){
// 							std::cout<<"local_spectrum dist_type: "<<dist_type<<std::endl;
// 							std::string dist_type_s;
// 					if(window_size==32 && dist_type ==1){}
// 					else{			
// 						for (int pix_x = dx; pix_x < input.width(); pix_x+=dx){
// 							for (int pix_y = dy; pix_y < input.height(); pix_y+=dy){
// 								// std::cout<<"distance_map_to_pixel"<<std::endl;
// 								switch(dist_type){
// 									case 0 :
// 										lsp->distance_map_to_pixel(test, pix_x, pix_y);
// 						        		dist_type_s = "Distsimple";
// 						        	break;
// 						        	case 1 :	    	
// 										lsp->distance_map_to_pixel_voisinage(test, pix_x, pix_y);
// 						        		dist_type_s = "Distvoisinage";
// 						        	break;
// 						        }

// 								double v_max = percentile (test,0.80);
// 								// std::string pix_loc= "_pix_"+(pix_x<10?"0":"")+(pix_x<100?"0":"")+std::to_string(pix_x)+"_"+(pix_y<10?"0":"")+(pix_y<100?"0":"")+std::to_string(pix_y);
// 								std::string pix_loc= "_pix_"+std::to_string(pix_x)+"_"+std::to_string(pix_y);
// 								print_image(test,path_in_b,name_file+"_"+welch_type_s+"_"+dist_type_s+"_wtype_"+std::string(1,lsp->get_windowtype())+"_wsize_"+std::to_string(lsp->get_windowsize())+ pix_loc ,v_max);
								

// 							}
// 						}
// 						if(dist_type != 1){	
// 							std::cout<<"local_spectrum MOSAIC DIST: "<<std::endl;
// 							mosa.initItk(input.width(),input.height(),true);
// 							lsp->distance_mosaic (mosa,16);
						
// 							double v_max;
// 							if(welchstep==1){
// 								v_max= percentile (mosa,0.80);

// 							}else{
// 								v_max= percentile (mosa,0.80);						
// 							}
// 							print_image(mosa,path_in_b,name_file+"_DistanceMosaic"+"_"+welch_type_s+"_"+dist_type_s+"_wtype_"+std::string(1,lsp->get_windowtype())+"_wsize_"+std::to_string(lsp->get_windowsize()) ,v_max);

// 							std::cout<<"DONE"<<std::endl;
// 						}
// 					}
// 				}
// 						mosa.initItk(input.width(),input.height(),true);
// 					std::cout<<"local_spectrum MOSAIC SPECTRUM: "<<std::endl;
// 					lsp->mosaic(mosa);
// 					double v_max = percentile (mosa,0.80);
// 					print_image(mosa,path_in_b,name_file+"_Mosaic"+"_"+welch_type_s+"_wsize_"+std::to_string(lsp->get_windowsize()) ,v_max);

// 					std::cout<<"DONE"<<std::endl;
// 			}
// 		}
// 	}





	// print_image(test,name_file+"_Mosaic_"+"_wsize_"+std::to_string(lsp.get_windowsize())+"_welchsize_"+std::to_string(welchsize)+"_wstep_"+std::to_string(welchstep) ,v_max);



	//loop integré pour concerver les ftt 
	// frequency_bilateral_filter(input,test,31,8,0.05);

//    bilateral_filter(input,test,8,8,0.05);
//    print_image(test,path_in_b,name_file+"bilateral_filter"+std::to_string(0) ,1);

// for (int i = 1; i < 9; i+=2)
// {
//    bilateral_filter(test,loop,8,8,0.05);
//        print_image(loop,path_in_b,name_file+"bilateral_filter"+std::to_string(i) ,1);

//    bilateral_filter(loop,test,8,8,0.05);
//        print_image(test,path_in_b,name_file+"bilateral_filter"+std::to_string(i+1) ,1);

// }

//     color_bilateral_filter(input_color,test_color,8,8,0.05*3);
//    print_image(test_color,path_out,name_file+"bilateral_filter"+std::to_string(0)+"_S" ,1);


//    // SOUSTRACT
//        ImageRGBd RES;
//        ImageRGBd RES_final;
//        RES.initItk(test_color.width(),test_color.height(),true);
//        RES_final.initItk(test_color.width(),test_color.height(),true);
//        for (int i = 0; i < test_color.width(); ++i)
//            for (int j = 0; j < test_color.height(); ++j){
//                RES.pixelAbsolute(i,j)[0] = input_color.pixelAbsolute(i,j)[0]- test_color.pixelAbsolute(i,j)[0];
//                RES.pixelAbsolute(i,j)[1] = input_color.pixelAbsolute(i,j)[1]- test_color.pixelAbsolute(i,j)[1];
//                RES.pixelAbsolute(i,j)[2] = input_color.pixelAbsolute(i,j)[2]- test_color.pixelAbsolute(i,j)[2];
//            }
//    //ADD 127 TO substract
//        for (int i = 0; i < test_color.width(); ++i)
//            for (int j = 0; j < test_color.height(); ++j){
//                RES_final.pixelAbsolute(i,j)[0] =clamp( RES.pixelAbsolute(i,j)[0]+0.5,0,1);
//                RES_final.pixelAbsolute(i,j)[1] =clamp( RES.pixelAbsolute(i,j)[1]+0.5,0,1);
//                RES_final.pixelAbsolute(i,j)[2] =clamp( RES.pixelAbsolute(i,j)[2]+0.5,0,1);
//            }
//        print_image(RES_final,path_out,name_file+"bilateral_filter"+std::to_string(0)+"_N" );



// for (int i = 0; i < 5; ++i)
// {
//    color_bilateral_filter(test_color,loop_c,8,8,0.05*3);
//            print_image(loop_c,path_out,name_file+"bilateral_filter"+std::to_string(i)+"_S" ,1);

//                    RES.initItk(test_color.width(),test_color.height(),true);
//                    RES_final.initItk(test_color.width(),test_color.height(),true);
//                    for (int i = 0; i < test_color.width(); ++i)
//                        for (int j = 0; j < test_color.height(); ++j){
//                            RES.pixelAbsolute(i,j)[0] = input_color.pixelAbsolute(i,j)[0]- loop_c.pixelAbsolute(i,j)[0];
//                            RES.pixelAbsolute(i,j)[1] = input_color.pixelAbsolute(i,j)[1]- loop_c.pixelAbsolute(i,j)[1];
//                            RES.pixelAbsolute(i,j)[2] = input_color.pixelAbsolute(i,j)[2]- loop_c.pixelAbsolute(i,j)[2];
//                        }
//                //ADD 127 TO substract
//                    for (int i = 0; i < test_color.width(); ++i)
//                        for (int j = 0; j < test_color.height(); ++j){
//                            RES_final.pixelAbsolute(i,j)[0] =clamp( RES.pixelAbsolute(i,j)[0]+0.5,0,1);
//                            RES_final.pixelAbsolute(i,j)[1] =clamp( RES.pixelAbsolute(i,j)[1]+0.5,0,1);
//                            RES_final.pixelAbsolute(i,j)[2] =clamp( RES.pixelAbsolute(i,j)[2]+0.5,0,1);
//                        }
//                    print_image(RES_final,path_out,name_file+"bilateral_filter"+std::to_string(i)+"_N" );





//    color_bilateral_filter(loop_c,test_color,8,8,0.05*3);
//            print_image(test_color,path_out,name_file+"bilateral_filter"+std::to_string(i+1)+"_S" ,1);


//                    RES.initItk(test_color.width(),test_color.height(),true);
//                    RES_final.initItk(test_color.width(),test_color.height(),true);
//                    for (int i = 0; i < test_color.width(); ++i)
//                        for (int j = 0; j < test_color.height(); ++j){
//                            RES.pixelAbsolute(i,j)[0] = input_color.pixelAbsolute(i,j)[0]- test_color.pixelAbsolute(i,j)[0];
//                            RES.pixelAbsolute(i,j)[1] = input_color.pixelAbsolute(i,j)[1]- test_color.pixelAbsolute(i,j)[1];
//                            RES.pixelAbsolute(i,j)[2] = input_color.pixelAbsolute(i,j)[2]- test_color.pixelAbsolute(i,j)[2];
//                        }
//                //ADD 127 TO substract
//                    for (int i = 0; i < test_color.width(); ++i)
//                        for (int j = 0; j < test_color.height(); ++j){
//                            RES_final.pixelAbsolute(i,j)[0] =clamp( RES.pixelAbsolute(i,j)[0]+0.5,0,1);
//                            RES_final.pixelAbsolute(i,j)[1] =clamp( RES.pixelAbsolute(i,j)[1]+0.5,0,1);
//                            RES_final.pixelAbsolute(i,j)[2] =clamp( RES.pixelAbsolute(i,j)[2]+0.5,0,1);
//                        }
//                    print_image(RES_final,path_out,name_file+"bilateral_filter"+std::to_string(i+1)+"_N" );

// }


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
