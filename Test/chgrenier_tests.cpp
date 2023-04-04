//
// Created by grenier on 28/03/23.
//
#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>
#include <ASTex/TilingBlending/controlled_tnb.h>
#include <Algo/TerrainEnhancement/transfer_functions.h>
#include <Algo/TerrainEnhancement/input_terrain.h>

using namespace ASTex;



int main(int argc, char** argv)
{
    if (argc<2)
    {
        std::cout << argv[0] << " img_terrain out_name " << std::endl;
        return 1;
    }

    using IMGT = ImageRGBu8;
    using IMGG = ImageGrayu8;
    // input terrain
    IMGT img_terrain;
    img_terrain.load(std::string(argv[1]));


    // output image
    int w = 256;
    int h = 256;
    IMGT img_out{w, h, false};
    IMGG gradX_out{w, h, false};
    IMGG gradY_out{w, h, false};


    // test récupération terrain et des gradient
    auto terrain = grab_Input_terrain(img_terrain, img_out.height(), img_out.width());
    terrain.terrain_img(img_out);
    compute_Gradient_terrain(img_out, gradX_out, gradY_out);

    img_out.save(std::string(argv[2])+".png");
    gradX_out.save(std::string(argv[2])+"_gradX.png");
    gradY_out.save(std::string(argv[2])+"_gradY.png");


    return EXIT_SUCCESS;
}


