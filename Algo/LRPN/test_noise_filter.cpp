#include <iostream>
#include <ASTex/internal.h>
#include "noise_filter.h"

int main( int argc, char ** argv )
{
	if( argc < 1 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " Source " << std::endl;
		return EXIT_FAILURE;
	}
	std::string filename_source = argv[1];
	std::string base_dir = ASTex::IO::remove_ext(filename_source) + "_noise_filter/";

	int filt_size = 16;
	int size_fft = 8;
	float sig_freq = 0.03;
	int nb_step_filt = 5;

	ASTex::RPnoise(filename_source, base_dir, filt_size, size_fft, sig_freq, nb_step_filt);


	return 0;
}
