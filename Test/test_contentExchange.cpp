#include <stdlib.h>
#include "ASTex/ContentExchange/atlas.h"

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> <out_directory> [save rendering pack]" << std::endl;
		return EXIT_FAILURE;
	}

	ImageRGBu8 im_in;
	im_in.load(std::string(argv[1]));

	std::string out_dir = argv[2];
	std::string name_file = IO::remove_path(argv[1]);
	std::string name_noext = IO::remove_ext(name_file);
	create_directory(out_dir);

    ContentExchange::PatchProcessor<ImageRGBu8> pProcessor(im_in);
    pProcessor.setFilteringMode(ISOTROPIC);
	pProcessor.setNbContentsPerPatch(5);
	pProcessor.setSeed(6);
//	pProcessor.patches_initRandom(32);
//	pProcessor.contents_initDefault();
//	pProcessor.contents_initRandom();
	pProcessor.fullProcess_oldMethod();

	create_directory(out_dir);
	if(argc>3 && std::atoi(argv[3])!=0)
	{
		std::string renderingPackDir = out_dir + "/" + name_noext + "_renderingPack";
		create_directory(renderingPackDir);
		pProcessor.saveRenderingPack(renderingPackDir);
	}
	pProcessor.setOutputSize(2*im_in.width(), 2*im_in.height());
	pProcessor.generate().texture().save(out_dir + "/" + name_noext + "_" + std::to_string(time(0)) + ".png");

	return 0;
}
