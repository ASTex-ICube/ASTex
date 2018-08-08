#include <stdlib.h>
#include "ASTex/easy_io.h"
#include "ASTex/utils.h"
#include "Stamping/stamper.h"
#include "histogram.h"
#include "rpn_utils.h"
#include "mipmap.h"
#include "texton_io.h"
#include "imageviewer.h"
#include "ContentExchange/patchProcessor.h"
#include "ContentExchange/atlas.h"

// v define your own
#define MY_PATH std::string("/home/nlutz/img/")

using namespace ASTex;

int test_contentExchange(int argc, char **argv)
{
    if(argc < 3)
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " <in_texture> <in_patchMap>" << std::endl;
        return EXIT_FAILURE;
    }

    ImageRGBd im_in, im_patches;

    if(!IO::loadu8_in_01(im_in, std::string(MY_PATH)+argv[1]))
        return 1;
    if(!IO::loadu8_in_01(im_patches, std::string(MY_PATH)+argv[2]))
        return 1;

    std::string out_dir=std::string(MY_PATH)+"contentMipmaps/";
    create_directory(out_dir);

    ContentExchange::PatchProcessor<ImageRGBd> pProcessor(im_in);
    //pProcessor.setFilteringMode(ANISOTROPIC);
    pProcessor.setFilteringMode(NO_FILTER);
    //add some patches here
    pProcessor.debug_setPatchFromImageRGB(im_patches);
    //patchProcessor.initializePatchesRegularGrid(64);
    pProcessor.initializeContents();
    //add some contents there
    pProcessor.debug_setRandomContents(5);

    ImageRGBd im_out;
    Mipmap<ImageRGBd> output(pProcessor.generate(512, 512));
    IO::save01_in_u8(output.texture(), MY_PATH + "cntexch_output.png");

    return 0;
}

int main( int argc, char **argv )
{
    return test_contentExchange(argc, argv);
}

