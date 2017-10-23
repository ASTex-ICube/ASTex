
#include <iostream>
#include <ASTex/internal.h>
#include "noise_synthesis_correlated_color.h"

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

	int nb_clusters = 2;

	ASTex::Synthesis_corel(filename_source, base_dir, nb_clusters);

	return EXIT_SUCCESS;
}

