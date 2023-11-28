#include "ASTex/image_rgb.h"
#include "ASTex/image_gray.h"
#include "Algo/MixMax/mixmax.h"
#include "ASTex/easy_io.h"

using namespace ASTex;

int main() {
    ImageRGBd T1, T2;
    ImageGrayd S1, S2, V;
    IO::loadu8_in_01(T1, TEMPO_PATH+"PBRMaterials/dirt_color.png");
    IO::loadu8_in_01(T2, TEMPO_PATH+"PBRMaterials/paving_color.png");
    IO::loadu8_in_01(S1, TEMPO_PATH+"PBRMaterials/dirt_height.png");
    IO::loadu8_in_01(S2, TEMPO_PATH+"PBRMaterials/paving_height.png");
    IO::loadu8_in_01(V, TEMPO_PATH+"Mask1.png");

    ImageRGBd out_texture;
    for(int i = 0; i <= 8; i++) {
        MixMax::blend(out_texture, T1, T2, S1, S2, V, i, 0.01, 0.4, false);
        IO::save01_in_u8(out_texture, TEMPO_PATH+"test_mixmax_" + std::to_string(i) + ".png");
        MixMax::ground_truth(out_texture, T1, T2, S1, S2, V, i, 0.01, 0.4, false);
        IO::save01_in_u8(out_texture, TEMPO_PATH+"test_gound_truth_" + std::to_string(i) + ".png");
    }

    return 0;
}
