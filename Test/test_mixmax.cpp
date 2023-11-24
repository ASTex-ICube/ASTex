#include "ASTex/image_rgb.h"
#include "ASTex/image_gray.h"
#include "Algo/MixMax/mixmax.h"
#include "ASTex/easy_io.h"

using namespace ASTex;

int main() {
    MixMax mixmax;

    ImageRGBd T1, T2;
    ImageGrayd S1, S2, V;

    IO::loadu8_in_01(T1, TEMPO_PATH+"/PBRMaterials/grass_color.png");
    IO::loadu8_in_01(T2, TEMPO_PATH+"/PBRMaterials/rock_color.png");
    IO::loadu8_in_01(S1, TEMPO_PATH+"/PBRMaterials/grass_height.png");
    IO::loadu8_in_01(S2, TEMPO_PATH+"/PBRMaterials/rock_height.png");
    IO::loadu8_in_01(V, TEMPO_PATH+"/Mask2.png");

    std::cout << S1.pixelAbsolute(0,0) << std::endl;
    std::cout << S1.pixelAbsolute(484,658) << std::endl;
    std::cout << S1.pixelAbsolute(918,33) << std::endl;

    ImageRGBd *OUT = mixmax.blend(&T1, &T2, &S1, &S2, &V, 0.01, 0.4);

    IO::save01_in_u8(*OUT, TEMPO_PATH+"test_mixmax.png");

    return 0;
}