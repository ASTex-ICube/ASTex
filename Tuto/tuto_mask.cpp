#include <iostream>

#include <ASTex/image_rgb.h>
#include <ASTex/mask.h>
#include <ASTex/easy_io.h>
#include <ASTex/region_traversor.h>

using namespace ASTex;


int main()
{
	ImageRGBu8 image;
	bool ok = image.load(TEMPO_PATH+"simpleRGB.png");
	if (!ok)
		return EXIT_FAILURE;

	//
	//	a MaskBool is an image of bool (store in uint8)
	//  with operator ()(x,y) -> bool (true if diff of 0)
	//
	MaskBool mb(image.width(), image.height());
	// fill randomly with true at 10%
	mb.random(0.1);

	// logical operators &= and |= apply on MaskBool
	// parameter is something like f(x,y)->bool (Mask or lambda)

	// example of &= |= operator with lambda-mask
	mb &= [&](int i, int /*j*/) { return i<image.width()/4;};

	(mb |= [&](int i, int) { return i> 3*image.width()/4;}) |= [&](int, int j) { return j> 3*image.width()/4;};

	//

	ImageGrayu8 mb_im(image.width(), image.height());

	// any kind of mask can be exported in an image uc(0/255)
	mb.export_binary_image(mb_im);
	mb_im.save(TEMPO_PATH+"mask.png");


	// Usage of mask is pincipaly to customize traversal
	image.for_all_pixels([&] (ImageRGBu8::PixelType& pix)
	{
		pix=RGBu8(0,0,0);
	},
	mb // mask
	);

	image.save(TEMPO_PATH+"simpleRGB_masked1.png");

	//
	// Example of usage of a MaskAboveThreshold
	//
	// create image
	ImageGrayd gd(image.width(), image.height());
	gd.for_all_pixels([&] (ImageGrayd::PixelType& v, int i, int j)
	{
		int dx = i - image.width()/2;
		int dy = j - image.height()/2;
		v =std::min(1.0,std::sqrt(double(dx*dx+dy*dy))/(image.width()/2));
	});

	// create mask from image
	MaskAboveThreshold<double> mat(gd,0.5);
	mat.export_binary_image(mb_im);
	mb_im.save(TEMPO_PATH+"mask_thre.png");


	return EXIT_SUCCESS;
}



