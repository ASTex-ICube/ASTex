#include <iostream>
#include <ASTex/image_gray.h>

using namespace ASTex;


int main()
{
	ImageGrayu8 image;

	bool ok = image.load(TEMPO_PATH+"simpleGray.png");
	if (!ok)
		return 1;

	std::cout << " Taille =" << image.width() << "/" << image.height() << std::endl;

	uint32_t W = image.width()/4;
	uint32_t H = image.height()/4;

	for(uint32_t j=0; j<H  ; ++j)
	{
		for(uint32_t i=0; i< W ; ++i)
		{
			image.pixelAbsolute(i,j) = 128;
		}
	}

	image.save(TEMPO_PATH+"simpleG2.png");

  return EXIT_SUCCESS;
}



