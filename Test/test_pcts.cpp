#include <stdlib.h>
#include "ASTex/PCTS/pcts.h"
#include "ASTex/easy_io.h"

using namespace ASTex;

int main(int argc, char **argv)
{
    std::cout << "PCTS started" << std::endl;

    ImageRGBd input, guidanceSynthesis, guidanceInput;
    ImageRGBd inputLabel, inputSkeleton, inputSkeletonInv, synthesisSkeleton, synthesisSkeletonInv;
    ImageRGBd synthesis, synthesisMask, stencil;

    std::string input_directory = std::string(argv[1]);

    std::string inputFilename = input_directory + "/" + "input.png";
    std::string inputLabelFilename = input_directory + "/" + "input_label.png";
    std::string inputSkeletonFilename = input_directory + "/" + "input_skeleton.png";
    std::string inputSkeletonInvFilename = input_directory + "/" + "input_skeleton_inv.png";
    std::string stencilFilename = input_directory + "/" + "stencil.png";

    std::string synthesisFilename = input_directory + "/" + "synthesis.png";
    std::string synthesisSkeletonFilename = input_directory + "/" + "synthesis_skeleton.png";
    std::string synthesisSkeletonInvFilename = input_directory + "/" + "synthesis_skeleton_inv.png";
    std::string synthesisMaskFilename = input_directory + "/" + "synthesis_mask.png";

    //set this image and macro PCTS_DEBUG_DIRECTORY in pcts.h
    IO::loadu8_in_01(input, inputFilename);
    IO::loadu8_in_01(inputLabel, inputLabelFilename);
    IO::loadu8_in_01(inputSkeleton, inputSkeletonFilename);
    IO::loadu8_in_01(inputSkeletonInv, inputSkeletonInvFilename);
    IO::loadu8_in_01(stencil, stencilFilename);
    guidanceInput.initItk(inputSkeleton.width(), inputSkeleton.height());
    guidanceInput.for_all_pixels([&](ImageRGBd::PixelType &pix, int x, int y)
    {
        double col[3];
        col[0] = pow(inputSkeleton.pixelAbsolute(x, y)[0],0.25);
        col[1] = pow(inputSkeletonInv.pixelAbsolute(x, y)[0],0.25);
        col[2] = 0.0;
        pix = ImageRGBd::PixelType(col);
    });
    IO::loadu8_in_01(synthesisSkeleton, synthesisSkeletonFilename);
    IO::loadu8_in_01(synthesisSkeletonInv, synthesisSkeletonInvFilename);
    guidanceSynthesis.initItk(synthesisSkeleton.width(), synthesisSkeleton.height());
    guidanceSynthesis.for_all_pixels([&](ImageRGBd::PixelType &pix, int x, int y)
    {
        double col[3];
        col[0] = pow(synthesisSkeleton.pixelAbsolute(x, y)[0],0.25);
        col[1] = pow(synthesisSkeletonInv.pixelAbsolute(x, y)[0],0.25);
        col[2] = 0.0;
        pix = ImageRGBd::PixelType(col);
    });
    IO::loadu8_in_01(synthesis, synthesisFilename);
    IO::loadu8_in_01(synthesisMask, synthesisMaskFilename);

    ASTex::Pcts<ImageRGBd> pcts;
    pcts.setTexture(input);
    pcts.setWidth(800);
    pcts.setHeight(800);

    unsigned    seed=0,
                blockSize = 8,
                sizeLog = 5,
                samplesNNM = 30,
                nbRefinementsNNM = 2,
                radiusNNM = 60;

     double     labelW = 0.9,
                guidanceW = 0.8,
                stencilW = 1.0;

    srand(seed);
    pcts.setLvl0BlockSize(blockSize);
    pcts.setMinimumSizeLog(sizeLog);
    pcts.setNbSamplesNNM(samplesNNM);
    pcts.setNbRefinementsNNM(nbRefinementsNNM);
    pcts.setRadiusScaleNNM(radiusNNM);

    pcts.setLabel(inputLabel, labelW);
//    pcts.setGuidance(guidanceSynthesis, guidanceInput, guidanceW, 1.0);
//    pcts.setSynthesis(synthesis, synthesisMask);
    pcts.setStencil(stencil, stencilW); // weight between 0 and 1

    IO::save01_in_u8(pcts.generate(), input_directory + "/out_pcts_seed" + std::to_string(seed)
                     + "_bs" + std::to_string(blockSize)
                     + "_sl" + std::to_string(sizeLog)
                     + "_samples" + std::to_string(samplesNNM)
                     + "_ref" + std::to_string(nbRefinementsNNM)
                     + "_labelW" + std::to_string(labelW)
                     + "_guidanceW" + std::to_string(guidanceW)
                     + "_stencilW" + std::to_string(stencilW)
                     + ".png");

    std::cout << "PCTS ended" << std::endl;
    return 0;
}
