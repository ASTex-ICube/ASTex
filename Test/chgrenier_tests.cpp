//
// Created by grenier on 28/03/23.
//
#include <ASTex/image_gray.h>

#include "Algo/CCVT_CGAL/ccvt.h"
#include "Algo/CCVT_CGAL/types.h"

using namespace ASTex;


int main()
{
    std::string working_directory = "/home/grenier/Documents/ASTex_fork/results/CCVT_CGAL/";
    CCVT main_ccvt;
    main_ccvt.load_image(working_directory + "density_theo.pgm");

    return EXIT_SUCCESS;
}


