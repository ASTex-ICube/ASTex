//
// Created by grenier on 30/03/23.
//
#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>

#include <Algo/TerrainEnhancement/controlled_tnb.h>
#include <Algo/TerrainEnhancement/transfer_functions.h>
#include <Algo/TerrainEnhancement/input_terrain.h>

using namespace ASTex;



int main(int argc, char** argv)
{
//    if (argc<8)
//    {
//        std::cout << argv[0] << " img_example ctrl_freq ctrl_or ctrl_ampl ctrl_mod img_to_gen width height" << std::endl;
//        return 1;
//    }


    std::string tnb_example_name = "/home/grenier/Documents/ASTex_fork/results/bi_chanel_noise.png";
    std::string terrain_name = "/home/grenier/Documents/ASTex_fork/results/terrain_HR.png";

    std::string ctrl_fr_name = "/home/grenier/Documents/ASTex_fork/results/grad_x.png";
    std::string ctrl_or_name = "/home/grenier/Documents/ASTex_fork/results/grad_y.png";
    std::string ctrl_ampl_name = "/home/grenier/Documents/ASTex_fork/results/noise_test.png";
    std::string ctrl_modulation_name = "/home/grenier/Documents/ASTex_fork/results/BF_HF.png";

    std::string img_to_gen_name = "/home/grenier/Documents/ASTex_fork/results/terrain_result";

    using IMGT = ImageRGBu8;
    using IMGG = ImageGrayu8;

    // input example
    IMGT img_ex;
    img_ex.load(tnb_example_name); // TODO : noise ?

    // control maps
    IMGT control_freq;
    control_freq.load(ctrl_fr_name);

    IMGT control_or;
    control_or.load(ctrl_or_name);

    IMGT control_ampl;
    control_ampl.load(ctrl_ampl_name);

    IMGT control_mod;
    control_mod.load(ctrl_modulation_name); // TODO : noise

    // input terrain
    IMGT img_terrain;
    img_terrain.load(terrain_name);


    // output image
    int w_coarse = 256;
    int h_coarse = 256;
    int w_fine = 1024;
    int h_fine = 1024;

    IMGT img_out_terrain{w_coarse, h_coarse, false};
    IMGG gradX_out{w_coarse, h_coarse, false};
    IMGG gradY_out{w_coarse, h_coarse, false};
    IMGT img_out_tnb{w_fine, h_fine, false};
    IMGT img_out_details{w_fine, h_fine, false};


    // terrain information
    auto terrain = grab_Input_terrain(img_terrain, h_coarse, w_coarse);
    terrain.terrain_img(img_out_terrain);
    compute_Gradient_terrain(img_out_terrain, gradX_out, gradY_out);

    img_out_terrain.save(img_to_gen_name+"_terrain.png");
    gradX_out.save(img_to_gen_name+"_gradX.png");
    gradY_out.save(img_to_gen_name+"_gradY.png");


    // tiling and blending
    auto tnb = make_Tiling_n_Blending(img_ex, control_freq, control_or);
    tnb.tile_img(img_out_tnb);
    img_out_tnb.save(img_to_gen_name+"_tnb.png");


    // transfer function
    auto tr_func = create_procedural_details(img_out_tnb, control_ampl, control_mod);
    tr_func.details_heighmap(img_out_details);
    img_out_details.save(img_to_gen_name+"_details.png");


    return EXIT_SUCCESS;
}
