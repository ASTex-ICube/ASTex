/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/




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
