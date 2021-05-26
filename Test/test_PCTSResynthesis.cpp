#include "ASTex/PCTS/pcts.h"

using namespace ASTex;

//Authors: Jean-Michel Dischler (dischler@unistra.fr), Nicolas Lutz (nicolas.nlutz@gmail.com)
//This file is not to be comitted to the main branch until it undergoes a complete rework,
//As it is hard to use for anyone but whoever used it in the state that it is.
//Additionally, every single map derived from the rigidity or the wrapped rigidity should be computed by this program,
//which include: mask.png, init_Binary_warped_specific_DT.png, init_Binary_warped_specific_DT_neg.png,
//Binary_warped_specific_DT.png, Binary_warped_specific_DT_neg.png.
//Finally, the program should take 5 arguments: <input> <wrapped input> <rigidity/stencil> <wrapped rigidity> <argument file>.
//The argument file should be read to provide the PCTS class arguments such as seed, block size, number of samples, etc.
//Comments made by Nicolas Lutz

int main(int argc, char **argv)
{
	std::cout << "PCTS started" << std::endl;

	ImageRGBd image, guid, seg;
	ImageRGBd mask, Ipos, Ineg, I2pos, I2neg;
	ImageRGBd synth, binary, stencil;

	char fname[256];

	std::string input_directory = std::string(argv[1]); //the working directory.
	char *name = argv[2];

	sprintf(fname, "%s/%s.png", input_directory.c_str(), name);
	IO::loadu8_in_01(image, fname); //image represents the input, called input.png (at workingDirectory/input.png).
	sprintf(fname, "%s/mask.png", input_directory.c_str());
	IO::loadu8_in_01(mask, fname);	//this represents the mask called mask.png.
									//I don't remember what it was for but it's supposed to be a grayscale distance map
									//to the center of the rigid regions, but with increased luminosity..????
									//For the last results I simply inverted the rigidity and had no problem.
	sprintf(fname, "%s/init_Binary_warped_specific_DT.png", input_directory.c_str());
	IO::loadu8_in_01(Ipos, fname);	//This is supposed to be a distance map to the center of the rigid regions.
	sprintf(fname, "%s/init_Binary_warped_specific_DT_neg.png", input_directory.c_str());
	IO::loadu8_in_01(Ineg, fname);  //This is supposed to be a distance map to the center of the non-rigid regions.
	sprintf(fname, "%s/rigidity_RGB.png", input_directory.c_str());
	IO::loadu8_in_01(stencil, fname);	//This is supposed to be the rigidity map, or an eroded version of it.
										//The stencil itself determines regions that the algorithm
										//is strictly forbidden to copy from the input,
										//So erode if the rigidity is taking ambiguous parts of the texture (sort of half rigid).
	seg.initItk(Ipos.width(), Ipos.height());
	seg.for_all_pixels([&](ImageRGBd::PixelType &pix, int x, int y)
	{
		double col[3];
		col[0] = pow(Ipos.pixelAbsolute(x, y)[0],0.25);
		col[1] = pow(Ineg.pixelAbsolute(x, y)[0],0.25);
		col[2] = 0.0;
		pix = ImageRGBd::PixelType(col);
	});
	sprintf(fname, "%s/Binary_warped_specific_DT.png", input_directory.c_str());
	IO::loadu8_in_01(I2pos, fname); //Same as init* but for the wrapped rigidity.
	sprintf(fname, "%s/Binary_warped_specific_DT_neg.png", input_directory.c_str());
	IO::loadu8_in_01(I2neg, fname); //Same as init* but for the wrapped rigidity.
	guid.initItk(I2pos.width(), I2pos.height());
	guid.for_all_pixels([&](ImageRGBd::PixelType &pix, int x, int y)
	{
		double col[3];
		col[0] = pow(I2pos.pixelAbsolute(x, y)[0],0.25);
		col[1] = pow(I2neg.pixelAbsolute(x, y)[0],0.25);
		col[2] = 0.0;
		pix = ImageRGBd::PixelType(col);
	});
	sprintf(fname, "%s/%s_C_10.png", input_directory.c_str(), name);
	IO::loadu8_in_01(synth, fname); //This is simply the wrapped texture to modify.
	sprintf(fname, "%s/%s_C_10_bin.png", input_directory.c_str(), name);
	IO::loadu8_in_01(binary, fname); //This is the wrapped rigidity that tells us which regions to modify.
									 //You can erode/dilate/modify this map offline if the result is not satisfying.

	ASTex::Pcts<ImageRGBd> pcts;	//You can use any image type. The mse used for distance computations was made
									//for ANY ImageType under 64 bits/channel.
									//Type mse() and F2 on it to discover the tricks that your doctor does not want you to know.

	//pcts.setWidth(800); //those arguments are used when creating an output, but we're just modifying one.
	//pcts.setHeight(800);

	//For the following arguments, the function in which they are used is self-explanatory and/or commented in pcts.h.
	unsigned    seed=13,
				bs = 20,
				sl = 6,
				samples = 15,
				ref = 1,
				radius = 60;

	 double     labelW = 0.5,
				guidanceW = 0.5,
				stencilW = 1.0;

	for(bs=14; bs<20; bs+=2)	//I made this loop to generate several results at once
								//but it would be better to do this sort of things outside of the program
		for(seed=15; seed<19;++seed)
		{
			srand(seed);
			pcts.setLvl0BlockSize(bs);
			pcts.setMinimumSizeLog(sl);
			pcts.setNbSamplesNNM(samples);
			pcts.setNbRefinementsNNM(ref);
			pcts.setRadiusScaleNNM(radius);
			pcts.setTextureLabel(mask, labelW);
			pcts.setGuidance(guid, seg, guidanceW, 1.0);
			pcts.setSynthesis(synth, binary);
			pcts.setStencil(stencil, stencilW); // weight between 0 and 1
			IO::save01_in_u8(pcts.generate(), input_directory + "/out_pcts_seed" + std::to_string(seed)
							 + "_bs" + std::to_string(bs)
							 + "_sl" + std::to_string(sl)
							 + "_samples" + std::to_string(samples)
							 + "_ref" + std::to_string(ref)
							 + "_labelW" + std::to_string(labelW)
							 + "_guidanceW" + std::to_string(guidanceW)
							 + "_stencilW" + std::to_string(stencilW)
							 + ".png");

			std::cout << "PCTS ended" << std::endl;
		}

	return 0;
}
