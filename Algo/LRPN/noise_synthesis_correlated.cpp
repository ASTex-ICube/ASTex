
#include <ASTex/special_io.h>
#include <ASTex/fourier.h>
#include <ASTex/local_spectrum.h>
#include <ASTex/utils.h>
#include <ASTex/colorspace_filters.h>
#include <ASTex/mask.h>

using namespace ASTex;

/**
 * @brief RPnoise random phase noise
 * @param inputfile a grayscale input example
 * @param outputfile is a patchwork image containing input + spectrum + output (pure random phase)
 * @return error code
 */
void Synthesis(const std::string& out_path, const std::string& inputfile, const std::string& synthese_name,
			  const ImageGrayd& N_input, const ImageGrayd S_input,
			  std::vector<ImageSpectrald>& spectrum, const std::vector<ImageGrayd>& guidance)
{


    std::string filename_plus_ext = basename(inputfile.c_str());
//    std::cout << "filename+ext "<<filename_plus_ext << std::endl;

    std::string delimiter = ".p";
    std::string name_file = filename_plus_ext.substr(0, filename_plus_ext.find(delimiter));

	create_directory(out_path);
	std::string path_in_b = out_path+name_file+"/"+synthese_name+"/"+
			std::to_string(guidance.size())+"/"+
			std::to_string(spectrum[0].width())+"/";
	create_directory(path_in_b);

    std::cout << "Creation file tree ok" <<name_file<< std::endl;

	// LOAD INPUT
	ImageGrayd input;
	IO::load_RGB_2_gray(input,inputfile);

	ImageGrayd noise_localised;
	noise_localised.initItk(N_input.width(),N_input.height(),true);

	ImageGrayd input_synthesis;
	input_synthesis.initItk(N_input.width(),N_input.height(),true);

	std::vector<ImageSpectrald> noise_synthesis;
	noise_synthesis.resize(guidance.size());

	for (uint i=0; i<guidance.size(); ++i)
	{
		noise_synthesis[i].initItk(N_input.width(),N_input.height(),true);
	}


	std::cout << "Print spectrum et synthese noise" << std::endl;
	srand (time(NULL));
	for (uint i=0; i<guidance.size(); ++i)
	{

		IO::save(spectrum[i], path_in_b+"SPECTRUM"+std::to_string(i)+".png");
		Fourier::RPnoise_mosaic(spectrum[i],noise_synthesis[i],4);
		IO::save(noise_synthesis[i], path_in_b+"RESULT"+std::to_string(i)+".png");
	}

	std::cout << "Synthese noise localisé" << std::endl;
	for (int i = 0; i < N_input.width(); ++i)
	{
		for (int j = 0; j < N_input.height(); ++j)
		{
			for (uint c = 0; c < guidance.size(); ++c)
			{

				noise_localised.pixelAbsolute(i,j) += noise_synthesis[c].pixelAbsolute(i,j)*guidance[c].pixelAbsolute(i,j);
			}
		}
	}

	IO::save(noise_localised,path_in_b+"RESULT_localised.png");


	std::cout << "Synthese Finale " << std::endl;
	for (int i = 0; i < N_input.width(); ++i)
	{
		for (int j = 0; j < N_input.height(); ++j)
		{
			input_synthesis.pixelAbsolute(i,j) += noise_localised.pixelAbsolute(i,j)+S_input.pixelAbsolute(i,j);
		}
	}

	IO::save(input_synthesis, path_in_b+"RESULT_SYNTHESIS.png");

}


int main( int argc, char ** argv )
{
	if( argc < 8 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <out_path> <input_filename> <synthese name> <input_N> <input_S> <guidance grey 1>  <guidance grey 2> [ <guidance grey> ]*" << std::endl;

		return EXIT_FAILURE;
	}

	const int nb_maps = (argc - 6) ;

	std::string out_path = argv[1];
	std::string filename_source = argv[2];
	std::string synthese_name = argv[3];

	ImageGrayd input;
	/*double mean = */
	IO::load(input,filename_source, true);


	ImageGrayd N_input;
	IO::load(N_input,argv[4],true);

	ImageGrayd S_input;
	IO::load(S_input,argv[5]);

	std::vector<ImageSpectrald> spectrum;
	spectrum.resize(nb_maps);


	std::vector<ImageGrayd> guidance;
	guidance.resize(nb_maps);


	if(synthese_name.compare("distance")==0)
	{
		for (int i=0; i<nb_maps; ++i)
			IO::load(guidance[i],argv[i+6]);
	}
	else
	{
		if(synthese_name.compare("coord")==0)
		{
			std::cout <<synthese_name<<" ok coord"<< std::endl;

			if(nb_maps==2)
			{
				for (int i=0; i<nb_maps; ++i)
					IO::load(guidance[i],argv[i+6],-1,3);
			}
			if(nb_maps==3)
			{
				for (int i=0; i<nb_maps; ++i)
					IO::load(guidance[i],argv[i+6],-1,3);
			}

		}
		else
		{
			if(nb_maps==2)
			{
				for (int i=0; i<nb_maps; ++i)
					IO::load(guidance[i],argv[i+6],0,2);
			}
			if(nb_maps==3)
			{
				for (int i=0; i<nb_maps; ++i)
					IO::load(guidance[i],argv[i+6],0,2);
			}
		}
	}

    std::vector<int> nb_pix_mask;
    nb_pix_mask.resize(nb_maps);

    std::vector<double> prop_selection_masque;
    prop_selection_masque.resize(nb_maps);


    for (uint k = 0; k < guidance.size(); ++k){
        prop_selection_masque[k] = 0;
        nb_pix_mask[k] = 0;
    }

    for (int i = 0; i < N_input.width(); ++i)
        for (int j = 0; j < N_input.height(); ++j)
        {
           double max = guidance[0].pixelAbsolute(i,j);
            int indice = 0;
                for (uint k = 1; k < guidance.size(); ++k){
                    if (max< guidance[k].pixelAbsolute(i,j))
                    {
                        max = guidance[k].pixelAbsolute(i,j);
                        indice = k;
                    }
                }
            nb_pix_mask[indice]++;
        }
    int nb_pix_max = N_input.width()*N_input.height();
    for (uint k = 0; k < guidance.size(); ++k){
        prop_selection_masque[k] = double(nb_pix_mask[k])/double(nb_pix_max);
        std::cout <<"NBPIX "<<  nb_pix_mask[k] << " PROP MASK " << k << " " << double(nb_pix_mask[k])/double(nb_pix_max)<< std::endl;
    }


	// Charger les spectres par auto-covariance
	for (int i = 0; i < nb_maps; ++i)
	{
        MaskLargestValue<double> mask(guidance[i],prop_selection_masque[i]);
		spectrum[i].initItk(64,64);
        Fourier::spectrum_by_autocorrelation_small_size(N_input,spectrum[i],mask,0.7,64);
	}


	Synthesis(out_path,filename_source,synthese_name,N_input,S_input,spectrum,guidance);

	return 0;
}

