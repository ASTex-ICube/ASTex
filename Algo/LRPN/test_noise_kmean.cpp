
#include <iostream>
#include <ASTex/internal.h>

#include "noise_kmean.h"

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
	int size_fft = 8;

	ASTex::KMnoise(filename_source, base_dir, nb_clusters, size_fft);

        return EXIT_SUCCESS;
}
