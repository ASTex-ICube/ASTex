
#include <iostream>
#include <ASTex/internal.h>

#include "biscalenoisepatchexg.h"

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

        int nb_clusters = 5;

        ASTex::biscalenoisepatchexg(filename_source, base_dir, nb_clusters, 15, 3, 5, 3, 2);

        return EXIT_SUCCESS;
}


