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



#include "PatchProcessor.h"
#include <stdlib.h>
#include <ASTex/slic.h>
#include <ASTex/easy_io.h>


using namespace ASTex;
using namespace ASTex::ContentExchg;



ASTex::ImageRGBu8 contentExchangeTest( const ASTex::ImageRGBu8 &img )
{
    int fragmentMinSize        = 150;
    int fragmentMaxSize        = 500;
    int fragmentColorThreshold = 40;
    int requiredPatchNumber    = 16;

    // Fragments generation.

    std::cout << "* Creating fragments..." << std::flush;
    ContentExchg::FragmentProcessor fProc( img );
    int nFrag = fProc.createFragments( fragmentMaxSize, fragmentColorThreshold );
    std::cout << " " << nFrag << " fragments created!" << std::endl;

    std::cout << "* Cleaning smallest fragments..." << std::flush;
    nFrag = fProc.cleanupSmallestFragments( fragmentMinSize );
    std::cout << " " << nFrag << " fragments created!" << std::endl;

    std::cout << "* Updating fragment attributes..." << std::flush;
    fProc.updateFragmentsAttributes();
    std::cout << " Done!" << std::endl;

#if 0

    ASTex::ImageRGBu8::PixelType *idColors = new ASTex::ImageRGBu8::PixelType [ fProc.fragmentCount() ];
    for( int i=0; i<fProc.fragmentCount(); ++i )
    {
        idColors[i].SetRed  ( 64+(rand()%192) );
        idColors[i].SetGreen( 64+(rand()%192) );
        idColors[i].SetBlue ( 64+(rand()%192) );
    }

    ASTex::ImageRGBu8 *imgFrag = new ASTex::ImageRGBu8( img.width(), img.height() );
    for( int i=0; i<img.height(); ++i )
        for( int j=0; j<img.width(); ++j )
            imgFrag->pixelAbsolute(j,i) = idColors[ fProc.fragmentIdAt(j,i) ];

	imgFrag->save("/tmp/testImg/frag.png");
#endif


    // Patches generation.

    std::cout << "* Creating patches..." << std::flush;
    ContentExchg::PatchProcessor pProc( fProc );
    pProc.createPatches( requiredPatchNumber );
    std::cout << " " << pProc.patchCount() << " patches created!" << std::endl;

    std::cout << "* Computing patch boundaries..." << std::flush;
    pProc.computePatchBoundaries();
    std::cout << " Done!" << std::endl;

#if 0
    ASTex::ImageRGBu8 *imgPatch = new ASTex::ImageRGBu8( img.width(), img.height() );
    for( auto &p : pProc.patches() )
        for( auto &frag : p.fragments )
            for( auto &pixel : fProc.fragmentById(frag).pixels )
                imgPatch->pixelAbsolute(pixel) = idColors[p.id];

    for( auto &p : pProc.patchById(0).boundary )
        imgPatch->pixelAbsolute(p).SetRed( 0 );

    delete [] idColors;
	imgPatch->save("/tmp/testImg/patch.png");
#endif


    // Look for alternative contents.

    std::vector<double> rotations;
    rotations.push_back( 0.0 );
    rotations.push_back( 0.5*M_PI );
    rotations.push_back( M_PI );
    rotations.push_back( -0.5*M_PI );

    std::vector<double> scales;
    scales.push_back( 1.0000 );
    //scales.push_back( 1.0125 );
    //scales.push_back( 1.0250 );
    //scales.push_back( 1.0375 );
    //scales.push_back( 1.0500 );

    std::cout << "* Looking for alternative contents..." << std::flush;
	pProc.findAlternativeContents( /*img,*/ rotations, scales, 10, 128 );
	std::cout << "Analyse Done!" << std::endl;


    // Perform synthesis test.

	std::cout << "Synthetising (4000x4000)" << std::endl;

    ASTex::ImageRGBAf patchMap;
    pProc.getPatchMap( patchMap, 0.005f * std::max(img.width(),img.height()) );

	ASTex::ImageRGBu8 synthesized( 4000, 4000 );

	synthesized.for_all_pixels( [&] (ASTex::ImageRGBu8::PixelType &pixel, int x, int y)
    {
        int sourceX = x % img.width();
        int sourceY = y % img.height();
        auto &pm = patchMap.pixelAbsolute( sourceX, sourceY );

        auto &patch = pProc.patchById( (int) pm.GetRed() );

        int centroidX = x - (int) pm.GetGreen();
        int centroidY = y - (int) pm.GetBlue ();

		std::srand( centroidY * img.width() + centroidX );
		auto &content = patch.contents[ std::rand() % patch.contents.size() ];

        Eigen::Matrix2d transform;
        transform(0,0) =  std::cos(content.angle) * content.scale;
        transform(1,0) =  std::sin(content.angle) * content.scale;
        transform(0,1) = -transform(1,0);
        transform(1,1) =  transform(0,0);

        Eigen::Vector2d transformedPixel = (transform * Eigen::Vector2d(pm.GetGreen(),pm.GetBlue()) ) + Eigen::Vector2d(content.offset[0]+0.5,content.offset[1]+0.5);
        itk::Index<2> shiftedPixel;
        shiftedPixel[0] = ((int) transformedPixel[0] + 4*img.width ()) % img.width ();
        shiftedPixel[1] = ((int) transformedPixel[1] + 4*img.height()) % img.height();

        float alpha = pm.GetAlpha();

        for( int i=0; i<3; ++i )
            pixel[i] = alpha * img.pixelAbsolute( shiftedPixel )[i] + (1.0-alpha) * img.pixelAbsolute(sourceX,sourceY)[i];
    });

	return synthesized;
}



int main( int argc, char **argv )
{
	if( argc != 2 )
        return EXIT_FAILURE;

    ASTex::ImageRGBu8 img;
	bool ok = img.load( argv[1] );

	if (ok)
	{
		ASTex::ImageRGBu8 out = contentExchangeTest( img );
		std::string out_name = ASTex::IO::remove_ext(argv[1])+"content_exchange.png";
		out.save(out_name);
		std::cout << out_name << " generated"<< std::endl;
	}

    return EXIT_SUCCESS;
}
