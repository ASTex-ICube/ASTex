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
