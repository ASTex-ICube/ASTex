
#include "io.h"
#include "fourier.h"
#include "itkJPEGImageIO.h"
#include "mode_seeking.h"
#include "output_geoffrey.h"
#include "image_treatment.h"
#include "utilities.h"
#include "tests.h"
#include "colorspace_filters.h"


using namespace ASTex;





int main( int argc, char ** argv )
{
    if( argc < 5 )
    { 
	    std::cerr << "Usage: " << std::endl;
	    std::cerr << argv[0] << " <output> <input> <mask 1> <colorized mask 1>  [ <mask> <colorized mask> ]*" << std::endl;

	    return EXIT_FAILURE;
    }

    const int nb_maps = (argc - 2) / 2;
    
    ImageRGBd im;
    im.load(argv[2]);

    const int im_w = im.width();
    const int im_h = im.height();

    ImageRGBd mosaic;
    mosaic.initItk(im_w*(nb_maps+1),im_h*2,true);

    fulfill_crop_image(im,mosaic,0,0);

    for (int i = 0; i < nb_maps; ++i)
    {
	    im.load(argv[2*i+3]);
    	fulfill_crop_image(im,mosaic,im_w*(i+1),0);
	    im.load(argv[2*i+4]);
    	fulfill_crop_image(im,mosaic,im_w*(i+1),im_h);
    }

  	print_image(mosaic, "./", argv[1], 255.0);

    return 0;
}
