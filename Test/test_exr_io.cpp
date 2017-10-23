#include <iostream>
#include <iomanip>

#include <ASTex/image_rgb.h>

#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>
#include <ASTex/exr_io.h>



using namespace ASTex;


void test_rgba()
{
	ImageRGBAf img1;
	img1.initItk(4,4);
	img1.for_all_pixels([&](ImageRGBAf::PixelType& p)
	{
		p[0] = std::numeric_limits<float>::min();
		p[1] = std::numeric_limits<float>::max();
		p[2] = std::numeric_limits<float>::lowest();
		p[3] = std::numeric_limits<float>::epsilon();
	});

	IO::EXR::save(img1, TEMPO_PATH+"img_exr_test1.exr");
	ImageRGBAf img2;
	IO::EXR::load(img2,TEMPO_PATH+"img_exr_test1.exr");

	std::cout << "== RGBAf == " << std::endl;
	img1.for_all_pixels([&](ImageRGBAf::PixelType& p, int x, int y)
	{
		std::cout << std::scientific <<p << " => " <<img2.pixelAbsolute(x,y)<< std::endl;
	});
}

void test_gray()
{
	ImageGrayf img1;
	img1.initItk(4,4);
	img1.for_all_pixels([&](ImageGrayf::PixelType& p)
	{
		p = std::numeric_limits<float>::min();
	});

	IO::EXR::save(img1, TEMPO_PATH+"img_exr_test2.exr");
	ImageGrayf img2;
	IO::EXR::load(img2,TEMPO_PATH+"img_exr_test2.exr");

	std::cout << "== GRAYf == " << std::endl;
	img1.for_all_pixels([&](ImageGrayf::PixelType& p, int x, int y)
	{
		std::cout << std::scientific <<p << " => " <<img2.pixelAbsolute(x,y)<< std::endl;
	});

}


void test_grayd()
{
	ImageGrayd img1;
	img1.initItk(4,4);
	img1.for_all_pixels([&](ImageGrayd::PixelType& p)
	{
		p = std::numeric_limits<double>::min();
	});

	IO::EXR::save(img1, TEMPO_PATH+"img_exr_test3.exr");
	ImageGrayd img2;
	IO::EXR::load(img2,TEMPO_PATH+"img_exr_test3.exr");

	std::cout << "== GRAYd == " << std::endl;
	img1.for_all_pixels([&](ImageGrayd::PixelType& p, int x, int y)
	{
		std::cout << std::scientific <<p << " == " <<img2.pixelAbsolute(x,y)<< std::endl;
	});

}


void test_rgbd()
{
	ImageRGBd img1;
	img1.initItk(4,4);
	img1.for_all_pixels([&](ImageRGBd::PixelType& p)
	{
		p[0] = std::numeric_limits<double>::min();
		p[1] = std::numeric_limits<double>::max();
		p[2] = std::numeric_limits<double>::lowest();
	});

	IO::EXR::save(img1, TEMPO_PATH+"img_exr_test4.exr");
	ImageRGBd img2;
	IO::EXR::load(img2,TEMPO_PATH+"img_exr_test4.exr");

	std::cout << "== RGBd == " << std::endl;
	img1.for_all_pixels([&](ImageRGBd::PixelType& p, int x, int y)
	{
		std::cout << std::scientific <<p << " => " <<img2.pixelAbsolute(x,y)<< std::endl;
	});

}



int main()
{
	test_rgba();
	test_gray();

	test_grayd();
	test_rgbd();

	ImageRGBf img;
	IO::loadu8_in_01(img, TEMPO_PATH+"simpleRGB.png");
	IO::EXR::save(img, TEMPO_PATH+"simpleRGB.exr");

  return EXIT_SUCCESS;
}

