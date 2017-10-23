
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



int colorize(const ImageGrayd& mask, ImageRGBd& image)
{
	for (int i = 0; i < image.width(); ++i)
	{
		for (int j = 0; j < image.height(); ++j)
		{
			image.pixelAbsolute(i,j)[0] *=  mask.pixelAbsolute(i,j);
			image.pixelAbsolute(i,j)[1] *=  mask.pixelAbsolute(i,j);
			image.pixelAbsolute(i,j)[2] *=  mask.pixelAbsolute(i,j);
		}
	}
}


int main( int argc, char ** argv )
{
    if( argc != 4 )
    { 
	    std::cerr << "Usage: " << std::endl;
	    std::cerr << argv[0] << " <input mask> <input color>  <output colorized mask> " << std::endl;

	    return EXIT_FAILURE;
    }

    ImageGrayd mask;
    load(mask,argv[1],0,2);
    ImageRGBd image;
    image.load(argv[2]);

	colorize(mask,image);

	std::string delimiter = ".p";
	std::string name_file = std::string(argv[3]).substr(0, std::string(argv[3]).find(delimiter));
    print_image(image, "./", name_file,255);
    return 0;
}
