//
// Created by grenier on 28/03/23.
//
#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include "ASTex/Noises/fBm.h"

using namespace ASTex;



int main()
{

    std::string img_to_gen_name = "/home/grenier/Documents/ASTex_fork/results/test_result";


    ImageGrayu8 img_out{1024, 1024, false};

    auto fBm = compute_fBm(12);
    fBm.fBm_image(img_out);

    img_out.save(img_to_gen_name + "_fBm.png");

    return EXIT_SUCCESS;
}


