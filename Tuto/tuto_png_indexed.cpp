#include <iostream>

#include <ASTex/image_indexed.h>

using namespace ASTex;



int main()
{
	ImageIndexedu8 img;
	bool ok = img.loadIndexedPNG(TEMPO_PATH+"small4.png");
	if (!ok)
		return 1;

	std::vector<ImageIndexedu8::PaletteColorType>& pal = img.palette();

	for (size_t i= 0; i<pal.size(); ++i)
	{
		std::cout << pal[i] << std::endl;
	}



	uint32_t HW = std::min(img.height(),img.width());

	for(uint32_t j=0; j<HW; ++j)
	{
		img.pixelAbsolute(j,j)=3;
	}

	for (ImageIndexedu8::IteratorIndexed it = img.beginIteratorIndexed(); !it.IsAtEnd(); ++it)
	{
		std::cout << "Pixel["<<it.GetIndex()[0]<< ","<< it.GetIndex()[1]<<"] >>> index = "<<int(it.Get()) << " >>> color = "<<img.color(it.Get())<< std::endl;
	}

	img.save(TEMPO_PATH+"indexed2.png");

	img.saveIndexedPNG(TEMPO_PATH+"ind4.png");

	ImageRGBu8 img2 = img.createRGB();
	img2.save(TEMPO_PATH+"indexedRGB.png");

	return EXIT_SUCCESS;

}

