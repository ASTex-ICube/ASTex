#include <stdlib.h>

#include "ASTex/Stamping/stamper.h"
#include "ASTex/histogram.h"
#include "ASTex/easy_io.h"
#include "ASTex/utils.h"
#include "ASTex/Stamping/stamper.h"
#include "ASTex/histogram.h"
#include "ASTex/rpn_utils.h"
#include "ASTex/texton_io.h"

/**
 * @brief this is a template for making texton noise.
 * Change the variable MY_PATH to the patch you have your images.
 * You can safely change the following lines:
 * sampler.setNbPoints(300);
 *      ^ which specifies the number of times texton noise hits the output texture
 * stamp.setInterpolationRule(Stamping::StampDiscrete<ImageRGBd>::BILINEAR);
 *      ^ which specifies the interpolation rule (we shoot between pixels)
 * tamponneur.setPeriodicity(false);
 *      ^ which specifies whether the output image is allowed to be periodic or not
 * @return 0
 */
int main(int argc, char **argv)
{
	if( argc < 5 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << "<input image> <output filename> <output width> <output height>" << std::endl;

		return EXIT_FAILURE;
	}

	std::string filename_source=argv[1];
	std::string name_file = IO::remove_path(filename_source);
	std::string name_noext = IO::remove_ext(name_file);

	std::string filename_output=argv[2];

	unsigned width = std::atoi(argv[3]);
	unsigned height = std::atoi(argv[4]);

	//Loading im_in

	ImageRGBd im_in, im_out, im_sample, im_texton;
	IO::loadu8_in_01(im_in, filename_source);

	//these import functions turn a texton file into a vizualizable image
	if(!import_texton(im_texton, filename_source))
	{
		import_texton_from_png(im_texton, filename_source);
	}

	HistogramRGBd histo_texton(im_texton);
	ImageRGBd::PixelType mean;
	for(int i=0; i<3; ++i)
		mean[i] = histo_texton.mean(i);
	//transformation image -> texton
	im_texton.for_all_pixels([&] (ImageRGBd::PixelType &pix) {
		pix -= mean;
		pix = pix * (1.0/std::sqrt(im_texton.width()*im_texton.height()));
	});

	//Testing

	Stamping::SamplerUniform sampler; //you can make it a PoissonGrid but it won't be as good
	sampler.setNbPoints(400); //< more than 30 with a PoissonSampler, more than, maybe, 100 with a uniform sampler
	Stamping::StampDiscrete<ImageRGBd> stamp(im_texton);
	stamp.setInterpolationRule(Stamping::StampDiscrete<ImageRGBd>::BILINEAR); //< you can change that too

	Stamping::StamperTexton<ImageRGBd> stamper(&sampler, &stamp);

	stamper.setPeriodicity(false);
	stamper.setUseMargins(true);

	im_out = stamper.generate(width, height);

	//transformation texton -> image
	im_out.for_all_pixels([&] (ImageRGBd::PixelType &pix) {
		//pix = pix * std::sqrt(im_texton.width()*im_texton.height()); //mistake: no need to re-normalize on top
		pix += mean;
	});

	im_sample.initItk(width, height, true);


	std::vector<Eigen::Vector2f> pointArray = sampler.generate();
	for(std::vector<Eigen::Vector2f>::iterator it=pointArray.begin(); it!=pointArray.end(); ++it)
	{
		int i = im_sample.width() * (*it)[0]; //i & j: single point coordinates in im_out
		int j = im_sample.height() * (*it)[1];

		im_sample.pixelAbsolute(i, j)=1.0;
	}

	//colored_RPN(im_in, im_rpn, RGB_SPACE, NORMAL, 0, 0, false, true, false, 1.0);

	auto clamp = [&] (ImageRGBd::PixelType &pix) { for(int i=0; i<3; ++i) pix[i] = pix[i] > 1.0 ? 1.0 : (pix[i] < 0.0 ? 0.0 : pix[i]); };
	im_out.for_all_pixels(clamp);

	//Saving sample
	//Saving im_out
	//IO::save01_in_u8(im_sample, some_output);
	IO::save01_in_u8(im_out, filename_output);
//	IO::save01_in_u8(im_texton, std::string(MY_PATH) + name_noext + "_texton.png");

	return 0;
}
