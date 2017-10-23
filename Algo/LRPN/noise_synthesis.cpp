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

int spectrum_accumulate (ImageGrayd& accu, const ImageGrayd& s)
{
	for (int fx=0; fx < s.width(); ++fx)
	{
		for (int fy=0; fy < s.height(); ++fy)
		{
			accu.pixelAbsolute(fx,fy) += s.pixelAbsolute(fx,fy);
		}
	}
	return 0;
}
int spectrum_rescale (ImageGrayd& s, double ratio)
{
	for (int fx=0; fx < s.width(); ++fx)
	{
		for (int fy=0; fy < s.height(); ++fy)
		{
			s.pixelAbsolute(fx,fy) *= ratio;
		}
	}
	return 0;
}

/**
 * @brief RPnoise random phase noise
 * @param inputfile a grayscale input example
 * @param outputfile is a patchwork image containing input + spectrum + output (pure random phase)
 * @return error code
 */
int Synthesis( std::string inputfile,std::string synthese_name,const int A,const int B, std::vector<ImageGrayd> crops,  std::vector<ImageGrayd> guidance)
{


    // LOAD INPUT
    ImageGrayd input;
//    load2gray(inputfile, input);


    std::string delimiter = ".p";
    std::string name_file = inputfile.substr(0, inputfile.find(delimiter));

	srand48(time(NULL));

	std::cout << "TRAITEMENT DE  " <<inputfile<< std::endl;

    std::string path_in_b ="/home/guingo/Images/Test_methods/bench_article/Synthesis/"+name_file+"/"+synthese_name+"/"+std::to_string(crops.size())+"/";
	system(std::string("mkdir -p "+path_in_b).c_str());
	
//    const int im_w = crops[0].width();
//    const int im_h = crops[0].height();

//    const int im_w = 32;
//    const int im_h = 32;

//    int m_left_bnd = -(im_w-1) / 2;
//    int m_right_bnd = (im_h) / 2;


//    std::vector<std::vector< double>> W;
//    W.resize( im_w , std::vector< double>( im_h ) );

//    m_left_bnd = -(im_w-1) / 2;
//    m_right_bnd = im_h / 2;

//    for(int x = m_left_bnd; x < m_right_bnd ; x++){
//        for(int y = m_left_bnd; y < m_right_bnd; y++){
////            W[x - m_left_bnd][y - m_left_bnd] = g_non_norm(x, y, m_windowsize/3);

//            std::cout << x - m_left_bnd << " "<< y -m_left_bnd << std::endl;
//            for (int var = 0; var < 2; ++var) {
//                crops[var].pixelAbsolute(x - m_left_bnd,y - m_left_bnd) *=  g_non_norm(x, y, im_w/2);
//            }
//        }
//    }

    std::cout << " ok tree file "<< std::endl;


//	std::vector<int> nb_pix_cluster;
//	nb_pix_cluster.resize(guidance.size());

	std::vector<ImageGrayd> spectrums;
    spectrums.resize(guidance.size());

    ImageGrayd noise_localised;
    noise_localised.initItk(guidance[0].width(),guidance[0].height(),true);

    ImageGrayd input_synthesis;
    input_synthesis.initItk(guidance[0].width(),guidance[0].height(),true);

    std::vector<ImageGrayd> noise_synthesis;
    noise_synthesis.resize(guidance.size());

    std::vector<ImageGrayd> noise_synthesis_localised;
    noise_synthesis_localised.resize(guidance.size());

    for (int i=0; i<guidance.size(); ++i){
//    	nb_pix_cluster[i]=0;
//    	spectrums[i].initItk(crops[i].width(),crops[i].height(),true);
//        spectrums[i].initItk(im_w,im_w,true);

        noise_synthesis[i].initItk(guidance[0].width(),guidance[0].height(),true);
        noise_synthesis_localised[i].initItk(guidance[0].width(),guidance[0].height(),true);

    }


//  Calcul des spectres par la moyenne des FFT des crop dans chaque cluster binarisé

//     for (int i = 0; i < N_input.width(); ++i)
//     {
//        for (int j = 0; j < N_input.height(); ++j)
//        {
//            int max_indice = 0;
//            float max_val = 0;

//            for (int c = 0; c < guidance.size(); ++c)
//            {
//                if (max_val < guidance[c].pixelAbsolute(i,j))
//                {
//                    max_indice=c;
//                    max_val = guidance[c].pixelAbsolute(i,j);
//                }
//            }

////            if(max_val>0.90){

//                const int ic = std::max(std::min(i,N_input.width()- m_right_bnd-1),-m_left_bnd); // the crop is not centered on i when too close to image boundaries
//                const int jc = std::max(std::min(j,N_input.height()- m_right_bnd-1),-m_left_bnd); // the crop is not centered on j when too close to image boundaries

//                ImageGrayd temp_crop;
//                temp_crop.initItk(im_w,im_h,true);
				
//                ImageGrayd sp;
//                ImageGrayd phase;

//                crop_image(N_input,temp_crop,ic+m_left_bnd,ic+m_right_bnd,jc+m_left_bnd,jc+m_right_bnd);
//                fftForwardModulusAndPhase(temp_crop, sp, phase);

//                spectrum_accumulate(spectrums[max_indice],sp);
//                nb_pix_cluster[max_indice]++;
////            }
//        }
//     }
    //ICI on load les spectres
        spectrums[0].load("/home/guingo/Seafile/2016_BiLayerTexturing_figures/spectrum_"+std::to_string(A) +".png");
        spectrums[1].load("/home/guingo/Seafile/2016_BiLayerTexturing_figures/spectrum_"+std::to_string(B) +".png");
        spectrum_rescale (spectrums[0],1/255.f);
        spectrum_rescale (spectrums[1],1/255.f);


        //

        monitorStats (guidance[0], "guidance 0 pst scale ");
        monitorStats (guidance[1], "guidance 1 pst scale ");

//        load_spectrum(spectrums[0],"/home/guingo/Seafile/2016_BiLayerTexturing_figures/spectrum_"+std::to_string(A) +".png");
//        load_spectrum(spectrums[1],"/home/guingo/Seafile/2016_BiLayerTexturing_figures/spectrum_"+std::to_string(B) +".png");

        monitorStats (spectrums[0], "spectrums 0");
        float mean = getMean(spectrums[0]);
        float fac = mean/0.5;

//        spectrum_rescale (spectrums[0],1/fac);
        monitorStats (spectrums[0], "spectrums 0 pst scale ");

         monitorStats (spectrums[1], "spectrums 1");
         mean = getMean(spectrums[1]);
         fac = mean/0.5;

//         spectrum_rescale (spectrums[1],1/fac);

         monitorStats (spectrums[1], "spectrums 1 pst scale ");


    std::cout << " ok load spec "<< std::endl;

    for (int i=0; i<guidance.size(); ++i)
    {
//        std::cout << " NB PIX "<< nb_pix_cluster[i] << std::endl;
        // if (nb_pix_cluster[i]>0)
//		ImageGrayd phase;

		// ICI, spectre calculer sur toute l'image
//         fftForwardModulusAndPhase(guidance[i],spectrums[i], phase);
    	
    	//ICI spectre calculer par des crops 
//        welch(crops[i],spectrums[i],1);
//         fftForwardModulusAndPhase(crops[i],spectrums[i], phase);

        //ICI spectre calculer par des crops
//        welch(crops[i],spectrums[i],1);
//

//    	print_image(spectrums[i],path_in_b,"SPECTRUM "+std::to_string(i) ,1);
        std::cout << " bodeore mosa "<< std::endl;

    	RPnoise_mosaic(spectrums[i],noise_synthesis[i],4);
        std::cout << " ok mosa "<< std::endl;


        for (int x = 0; x < guidance[0].width(); ++x)
        {
            for (int y = 0; y < guidance[0].height(); ++y)
            {
                    noise_synthesis[i].pixelAbsolute(x,y) += 0.5;

                }
            }

        print_image(noise_synthesis[i],path_in_b,"RESULT"+std::to_string(i) ,1);

	}

    for (int i = 0; i < guidance[0].width(); ++i)
	{
        for (int j = 0; j < guidance[0].height(); ++j)
		{
			for (int c = 0; c < guidance.size(); ++c)
			{
                //conservation des bord
//                if(guidance[c].pixelAbsolute(i,j)< 0.6 && guidance[c].pixelAbsolute(i,j)> 0.4){
//                           noise_localised.pixelAbsolute(i,j) = input.pixelAbsolute(i,j)-0.15f;
//                }else{
                    noise_localised.pixelAbsolute(i,j) += noise_synthesis[c].pixelAbsolute(i,j)*(double)guidance[c].pixelAbsolute(i,j);
//                }
            }
            noise_synthesis_localised[0].pixelAbsolute(i,j) = noise_synthesis[0].pixelAbsolute(i,j)*(double)guidance[0].pixelAbsolute(i,j);
            noise_synthesis_localised[1].pixelAbsolute(i,j) = noise_synthesis[1].pixelAbsolute(i,j)*(double)guidance[1].pixelAbsolute(i,j);

		}
	}

    
    print_image(noise_localised,path_in_b,"RESULT_localised ",1);

    print_image(noise_synthesis_localised[0],path_in_b,"mask_colorised_0",1);
    print_image(noise_synthesis_localised[1],path_in_b,"mask_colorised_1",1);


        // AJOUT N localised +S

//	for (int i = 0; i < N_input.width(); ++i)
//	{
//		for (int j = 0; j < N_input.height(); ++j)
//		{
//            input_synthesis.pixelAbsolute(i,j) += noise_localised.pixelAbsolute(i,j)+S_input.pixelAbsolute(i,j)-0.5;//-0.5f
//		}
//	}
////       double v_max = percentile (input_synthesis,0.80);

//        print_image(input_synthesis,path_in_b,"RESULT_SYNTHESIS " ,1);

	// RESULT
   

	// float filter_size = FFT_SIZE;
	// float sigma_dist = FFT_SIZE/4;
	// float sigma_freq = value_z;
	// frequency_joint_bilateral_filter(input, test,guidance ,filter_size , name_file ,sigma_dist,sigma_freq);

	// 	std::cout<<"FILTERING DONE"<<std::endl;


	// compute local spectra 








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
        std::cerr << argv[0] << "<input_name> <synthese name> <input_N> <input_S> <guidance grey 1>  <guidance grey 2> [ <guidance grey> ]*" << std::endl;

	    return EXIT_FAILURE;
    }

    const int nb_maps = (argc - 5) ;


    std::string filename_source = argv[1];
    std::string synthese_name = argv[2];

    int A = atoi(argv[3]);
//   load(N_input,argv[3]);

    int B = atoi(argv[4]);
//   load(S_input,argv[4]);

    std::vector<ImageGrayd> crops;
    crops.resize(nb_maps);


    std::string delimiter = ".p";
    std::string name_file = filename_source.substr(0, filename_source.find(delimiter));

//    std::string path_in_b ="/home/guingo/Seafile/Geoffrey/N_et_S_COLOR/"+name_file+"/size_4/0.030000/";
        std::string path_in_b ="/home/guingo/Seafile/Geoffrey/N_et_S_COLOR/"+name_file+"/size_4/0.008000/";


//        for (int i = 0; i < nb_maps; ++i) {

//            load2luminance(name_file+"_nclusters_"+std::to_string(nb_maps)+"_noise_cl00"+std::to_string(i)+".inpaint.png", crops[i],path_in_b);
//        }

    std::vector<ImageGrayd> guidance;
    guidance.resize(nb_maps);


//    if(synthese_name.compare("distance")==0){
        for (int i=0; i<nb_maps; ++i)
           load(guidance[i],argv[i+5]);
//    }else{
//        if(synthese_name.compare("coord")==0){
//            std::cout <<synthese_name<<" ok coord"<< std::endl;

//            if(nb_maps==2){
//                for (int i=0; i<nb_maps; ++i)
//                   load(guidance[i],argv[i+5],-1,3);
//            }
//            if(nb_maps==3){
//                for (int i=0; i<nb_maps; ++i)
//                   load(guidance[i],argv[i+5],-1,3);
//            }

//        }else{
//            if(nb_maps==2){
//                for (int i=0; i<nb_maps; ++i)
//                   load(guidance[i],argv[i+5],0,2);
//            }
//            if(nb_maps==3){
//                for (int i=0; i<nb_maps; ++i)
//                   load(guidance[i],argv[i+5],0,2);
//            }
//        }
//    }

//      load_coordinates(guidance[i],argv[i+10]);


    Synthesis(filename_source,synthese_name,A,B,crops,guidance);

	return 0;
	// return RPnoise("/home/guingo/Images/Test_methods/sample/bright_black.png","/home/guingo/Images/Test_methods/result_noise/result_bright_black.png");

	// return RPnoise("/home/guingo/Images/Test_methods/sample/lave2.png","/home/guingo/Images/Test_methods/result_noise/result_lave2.png");

//		return testPhase(TEMPO_PATH+"refnoise.png","/tmp/result.png");
}

