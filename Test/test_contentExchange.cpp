#include <stdlib.h>
#include "ASTex/ContentExchange/atlas.h"
#include "ASTex/ContentExchange/benchmarker.h"

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> <out_directory>" << std::endl;
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
	pProcessor.fullProcess_oldMethod();

	std::string renderingDirectory = out_dir + "/" + name_noext + "_" + std::to_string(std::time(0)) + "/";
	create_directory(renderingDirectory);
	pProcessor.saveRenderingPack(renderingDirectory);

	return 0;
}
