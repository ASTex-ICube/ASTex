//
// Created by grenier on 28/03/23.
//
#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>
#include <ASTex/TilingBlending/controlled_tnb.h>
#include <Algo/TerrainEnhancement/transfer_functions.h>

using namespace ASTex;



int main(int argc, char** argv)
{
    if (argc<8)
    {
        std::cout << argv[0] << " img_example ctrl_freq ctrl_or ctrl_ampl ctrl_mod img_to_gen width height" << std::endl;
        return 1;
    }

    using IMGT = ImageRGBu8;
    // input example
    IMGT img_ex;
    img_ex.load(std::string(argv[1]));

    // control maps
    IMGT control_freq;
    control_freq.load(std::string(argv[2]));

    IMGT control_or;
    control_or.load(std::string(argv[3]));

    IMGT control_ampl;
    control_ampl.load(std::string(argv[4]));

    IMGT control_mod;
    control_mod.load(std::string(argv[5]));

    // output image
    int w = atoi(argv[7]);
    int h = atoi(argv[8]);
    IMGT img_out_tnb{w, h, false};
    IMGT img_out{w, h, false};


    // tiling and blending
    auto tnb = make_Tiling_n_Blending(img_ex, control_freq, control_or);
    tnb.tile_img(img_out_tnb);
//    img_out_tnb.save(std::string(argv[6]));


    // test transfer function
    auto tr_func = create_procedural_details(img_out_tnb, control_ampl, control_mod);
    tr_func.details_heighmap(img_out);
    img_out.save(std::string(argv[6]));


    return EXIT_SUCCESS;
}


