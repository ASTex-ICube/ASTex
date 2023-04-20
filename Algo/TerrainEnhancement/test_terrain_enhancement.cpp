//
// Created by grenier on 28/03/23.
//
#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include "ASTex/Noises/fBm.h"

#include <Algo/TerrainEnhancement/input_terrain.h>
#include <Algo/TerrainEnhancement/control_maps.h>
#include <Algo/TerrainEnhancement/controlled_tnb.h>
#include <Algo/TerrainEnhancement/transfer_functions.h>
#include <Algo/TerrainEnhancement/terrain_enhancement.h>

using namespace ASTex;



int main()
{
//    using IMGT = ImageRGBu8;
//    using IMGG = ImageGrayu8;

    /* principe :
     * en entrée j'ai un terrain,
     * je récupère le champ de hauteur et les gradients
     * en sortie j'ai les carte de fréquence, orientation et amplitude
     * que j'utilise pour générer des détails à ajouter au terrain
     */

    std::string tnb_example_name = "/home/grenier/Documents/ASTex_fork/results/bi_chanel_noise.png";
    std::string terrain_name = "/home/grenier/Documents/ASTex_fork/results/terrain_ubisoft_u8.png";

    std::string img_to_gen_name = "/home/grenier/Documents/ASTex_fork/results/test_result";

//    std::string ctrl_fr_name = "/home/grenier/Documents/ASTex_fork/results/un.png";
//    std::string ctrl_or_name = "/home/grenier/Documents/ASTex_fork/results/grad_y.png";
//    std::string ctrl_ampl_name = "/home/grenier/Documents/ASTex_fork/results/un.png";
//    std::string ctrl_modulation_name = "/home/grenier/Documents/ASTex_fork/results/BF_HF.png";

    // génération d'un gradient de test
//    ImageGrayu8 img_grad{256, 256, false};
//    img_grad.parallel_for_all_pixels([&] (typename ImageGrayu8::PixelType& P, int x, int y){P = ImageGrayu8::PixelType(std::clamp(x+y, 0, 255));});
//    img_grad.save(img_to_gen_name+"_sample.png");


    // input terrain
    ImageGrayu8 img_terrain; // TODO : u8 ou u16 ?
    img_terrain.load(terrain_name);
//    img_terrain.save(img_to_gen_name+"_load_terrain.png");

    // input example
    ImageRGBu8 img_ex;
    img_ex.load(tnb_example_name);


    // control maps
//    ImageGrayu8 control_freq;
//    control_freq.load(ctrl_fr_name);

//    ImageGrayu8 control_or;
//    control_or.load(ctrl_or_name);

//    ImageGrayu8 control_ampl;
//    control_ampl.load(ctrl_ampl_name);

//    ImageGrayu8 control_mod;
//    control_mod.load(ctrl_modulation_name);





    // output image
    int w_coarse = img_terrain.width()/64.;
    int h_coarse = img_terrain.height()/64.;
    int w_fine = img_terrain.width()*1.;
    int h_fine = img_terrain.height()*1.;

    ImageGrayu8 img_out_terrain{w_coarse, h_coarse, false};
    ImageGrayu8 gradX_out{w_coarse, h_coarse, false};
    ImageGrayu8 gradY_out{w_coarse, h_coarse, false};

    ImageGrayu8 control_freq{w_coarse, h_coarse, false};
    ImageGrayu8 control_or{w_coarse, h_coarse, false};
    ImageGrayu8 control_ampl{w_coarse, h_coarse, false};

    ImageGrayu8 control_mod{w_fine, h_fine, false};
    ImageRGBu8 img_out_tnb{w_fine, h_fine, false};
    ImageGrayu8 img_out_details{w_fine, h_fine, false};
    ImageGrayu8 img_out_final{w_fine, h_fine, false};




    // terrain information
    auto terrain = grab_Input_terrain(img_terrain, h_coarse, w_coarse);
    terrain.terrain_img(img_out_terrain, gradX_out, gradY_out);

    img_out_terrain.save(img_to_gen_name+"_terrain.png");
    gradX_out.save(img_to_gen_name+"_gradX.png");
    gradY_out.save(img_to_gen_name+"_gradY.png");


    // cartes de controle
    auto control_maps = create_control_maps(img_out_terrain, gradX_out, gradY_out);
    control_maps.Set_Frequency_param(2., 2.);
    control_maps.Set_Amplitude_param(2., 6.);
    control_maps.compute_control(control_freq, control_or, control_ampl);

    control_freq.save(img_to_gen_name+"_frequ.png");
    control_ampl.save(img_to_gen_name+"_ampl.png");
    control_or.save(img_to_gen_name+"_or.png");



    // tiling and blending
    auto tnb = make_Tiling_n_Blending(img_ex, control_freq, control_or);
    tnb.Set_Frequency_max(2.);
    tnb.Set_Lattice_resolution(2.); // attention : triangle petit = meilleur orientation mais pb avec les basse fréquences
    tnb.tile_img(img_out_tnb);

    img_out_tnb.save(img_to_gen_name+"_tnb.png");


    // generate fBm for profile modulation
    auto fBm = compute_fBm(6);
    fBm.fBm_image(control_mod);

    control_mod.save(img_to_gen_name + "_fBm.png");


    // transfer function
    auto tr_func = create_procedural_details(img_out_tnb, control_mod);
    tr_func.details_heighmap(img_out_details);

    img_out_details.save(img_to_gen_name+"_details.png");


    // terrain amplifié
    auto final_terrain = compute_final_terrain(img_terrain, img_out_details, control_ampl);
    final_terrain.Set_Amplitude_max(.1);
    final_terrain.final_terrain_img(img_out_final);

    img_out_final.save(img_to_gen_name+"_final.png");

    /* détail d'implem :
     * j'ai finalement mis l'utilisation de l'amplitude max au moment de l'application des détails
     * la carte de contrôle est calculée avant mais ne donne que les variation spatiale de l'amplitude (valeurs normées)
     * définir l'amplitude max à la fin permet de garder la séparation des étpaes dans des images séparées (y'a peut être mieux comme méthode mais c'ets la seule que j'ai vu)
     */



    return EXIT_SUCCESS;
}


