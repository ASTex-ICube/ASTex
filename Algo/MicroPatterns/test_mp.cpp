#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include <ASTex/easy_io.h>

#include <Algo/MicroPatterns/micropatterns.h>
#include <Algo/MicroPatterns/tiling.h>
#include <Algo/MicroPatterns/blurred_colorMap.h>
#include <Algo/MicroPatterns/mipmaped_noise.h>

#include<iomanip>

using namespace ASTex;



template<typename IMGD>
void mipmap_sqp2_image(const IMGD& img, std::vector<IMGD>& mipmap)
{
	using TD = typename IMGD::PixelType;
	using TED = typename IMGD::DoublePixelEigen;
	int nbl = std::log2(img.width())+1;
	mipmap.clear();
	mipmap.reserve(nbl);
	int w = img.width();
	const IMGD* mipmap_up = &img;
	while(w>1)
	{
		w /= 2;
		mipmap.emplace_back();
		mipmap.back().initItk(w, w);
		mipmap.back().parallel_for_all_pixels([&](TD& p, int x, int y)
		{
			int x2 = x * 2;
			int y2 = y * 2;

			TED v = eigenPixel<double>(mipmap_up->pixelAbsolute(x2, y2));
			++x2;
			v += eigenPixel<double>(mipmap_up->pixelAbsolute(x2, y2));
			++y2;
			v += eigenPixel<double>(mipmap_up->pixelAbsolute(x2, y2));
			--x2;
			v += eigenPixel<double>(mipmap_up->pixelAbsolute(x2, y2));
			v/=4.0;
			p = IMGD::itkPixel(v);
		});
		mipmap_up = &(mipmap.back());
	}
}

void repete(const PackedInputNoises& noises, double level, const Eigen::Vector2d& uv, Eigen::Vector2d& mean, Eigen::Vector2d& variance)
{
	Eigen::Vector4d n = noises.fetch(uv, level);
	mean = Eigen::Vector2d(n[0], n[1]);
	variance = Eigen::Vector2d(n[2], n[3]);
}


#define COLOR

#ifdef  COLOR
using T_IMG = ImageRGBu8;
using T_IMG_D = ImageRGBd;
#else
using T_IMG = ImageGray8;
using T_IMG_D = ImageGrayd;
#endif

int main_all(int argc, char** argv, const std::vector<std::array < std::string, 4>>& names)
{
	bool bs = false;
	bs = std::string{ argv[argc-1] } == std::string{ "save" };

	std::string dir {argv[1]};
	int conf = std::atoi(argv[2]);

	T_IMG_D cm;
	IO::loadu8_in_01(cm, dir+"colormap/"+names[conf][3] + ".png");
	auto start_chrono = std::chrono::system_clock::now();
	BlurredColorMaps<T_IMG_D> bcm{ &cm }; //{0.375, 0.25, 0.0625}, 7};

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "Color maps blurred in " << elapsed_seconds.count() << " s." << std::endl;
	if (bs)
		bcm.save(names[conf][0]);

	ImageGrayd nA;
	IO::loadu8_in_01(nA, dir+names[conf][1]+".png");
	ImageGrayd nB;
	IO::loadu8_in_01(nB, dir+names[conf][2]+".png");
	start_chrono = std::chrono::system_clock::now();
	PackedInputNoises pin(nA, nB);
	elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "noise var mipmaped in " << elapsed_seconds.count() << " s." << std::endl;
	if (bs)
		pin.save(names[conf][0]);
	int sz = std::atoi(argv[3]);;
	double scale = std::atof(argv[4]);;

	T_IMG_D img_patterns;
	Tiling_n_Blendinng tnb;

	int i=0; //for the name
	while(sz>0)
	{
		start_chrono = std::chrono::system_clock::now();
		std::cout << "compute "<<i<<" => sz="<<sz<<" => sc="<<scale<< std::endl;
		patterns(img_patterns, sz, pin, bcm, scale, tnb);
		std::cout << i << " computed "<<std::endl;
		std::string fn = names[conf][0] + "_.png";
		fn[fn.length() - 5] = (i <= 9) ? '0' + i : 'A' + i - 10;
		elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
		std::cout << fn << " of " << sz<<" x "<< sz << " pixels generated in " << elapsed_seconds.count() << " s." << std::endl;
		IO::save01_in_u8(img_patterns, fn);
		
		std::cout << "file exported "<< std::endl;
		sz /= 2;
		scale *= 2.0;
		i++;
	}
	
	return EXIT_SUCCESS;
}

int main_mipmap(int argc, char** argv, const std::vector<std::array < std::string, 4>>& names)
{
	bool bs = false;
	bs = std::string{ argv[argc - 1] } == std::string{ "save" };

	double scale = std::atof(argv[4]);
	int sz = std::atoi(argv[3]);

	int conf = std::atoi(argv[2]);
	std::string dir{argv[1]};
	T_IMG_D cm;
	IO::loadu8_in_01(cm, dir + "colormap/" + names[conf][3] + ".png");
	auto start_chrono = std::chrono::system_clock::now();
	BlurredColorMaps<T_IMG_D> bcm{ &cm };
	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "Color maps blurred in " << elapsed_seconds.count() << " s." << std::endl;
	if (bs)
		bcm.save(names[conf][0]);

	ImageGrayd nA;
	IO::loadu8_in_01(nA, dir + names[conf][1] + ".png");
	ImageGrayd nB;
	IO::loadu8_in_01(nB, dir + names[conf][2] + ".png");
	start_chrono = std::chrono::system_clock::now();
	PackedInputNoises pin(nA, nB);
	elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "noise var mipmaped in " << elapsed_seconds.count() << " s." << std::endl;
	if (bs)
		pin.save(names[conf][0]);

	T_IMG_D img_patterns;

	start_chrono = std::chrono::system_clock::now();
	Tiling_n_Blendinng tnb;
	patterns(img_patterns, sz, pin, bcm, scale, tnb);
	std::string fn = names[conf][0] + "_mm_0.png";
	elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << fn << " of " << sz << " x " << sz << " pixels generated in " << elapsed_seconds.count() << " s." << std::endl;
	std::vector<T_IMG_D> mm_im;
	mipmap_sqp2_image(img_patterns, mm_im);
	IO::save01_in_u8(img_patterns, fn);
	int i=1;
	for(const auto& im : mm_im)
	{
		std::string fn = names[conf][0] + "_mm_X.png";
		fn[fn.length() - 5] = (i <= 9) ? '0' + i : 'A' + i - 10;
							IO::save01_in_u8(im, fn);
		++i;
	}
	std::cout << "file exported " << std::endl;

	return EXIT_SUCCESS;
}

int main_stat(int argc, char** argv, const std::vector<std::array < std::string, 4>>& names)
{
	double scale = std::atof(argv[4]);
	int sz = std::atoi(argv[3]);

	int conf = std::atoi(argv[2]);
	std::string dir{argv[1]};
	T_IMG_D cm;
	IO::loadu8_in_01(cm, dir + "colormap/" + names[conf][3] + ".png");
	BlurredColorMaps<T_IMG_D> bcm{ &cm };

	ImageGrayd nA;
	IO::loadu8_in_01(nA, dir + names[conf][1] + ".png");
	ImageGrayd nB;
	IO::loadu8_in_01(nB, dir + names[conf][2] + ".png");
	PackedInputNoises pin(nA, nB);


	T_IMG_D img_patterns;

	Tiling_n_Blendinng tnb;
	double mult = (argc>4) ? std::atof(argv[4]) : 256.0;
	patterns_stat(img_patterns, sz, pin, bcm, scale, tnb ,mult); 
		
	return EXIT_SUCCESS;
}


int main(int argc, char** argv)
{
	std::vector<std::array < std::string, 4>> names = { {"stone", "noise_pierre_1", "noise_pierre_2", "cm11"},
		{"bark", "noise_ecorce_1", "noise_ecorce_2", "cm10"},
		{"camouflage", "noise_1024_2", "noise_1024_4", "cm1"},
		{"green", "noise_1024_4", "noise_1024_2", "cm2"},
		{"blue", "noise_256_2", "noise_256_4", "cm5"},
		{"lava", "noise_1024_2", "noise_1024_4", "cm8"},
		{"water", "n1", "n2", "eau_cm"},
		{"hexa", "n1", "n2", "hexa_cm"},
		{"phasor_sand", "noise_sin_3", "noise_cos_3", "colormap_phasor_sand"},
		{"phasor_sin", "noise_sin_1", "noise_cos_1", "colormap_phasor_sin"},
		{"phasor_square", "noise_sin_1", "noise_cos_1", "colormap_phasor_square"},
		{"grey_cm", "n1", "n2", "grey_cm_1"} };

	if (argc == 1)
	{
		std::cout << "all dir config size scale_initial" << std::endl;
		std::cout << "mipmap dir config size scale_initial" << std::endl;
		std::cout << "stat dir config size scale" <<std::endl;

		std::cout << "config: ";
		for (int i = 0; i < int(names.size()); ++i)
		
			std::cout << i << ":" << names[i][0] << " / ";
		std::cout << std::endl;
		return EXIT_SUCCESS;
	}


	if (std::string(argv[1])=="mipmap")
		main_mipmap(argc-1, argv+1,names);
	if (std::string(argv[1])=="all")
		main_all(argc-1, argv+1,names);
	if (std::string(argv[1])=="stat")
		main_stat(argc-1, argv+1,names);
}

//["stone", "noise_pierre_1", "noise_pierre_2", "cm11"],
//["bark", "noise_ecorce_1", "noise_ecorce_2", "cm10"],
//["camouflage", "noise_1024_2", "noise_1024_4", "cm1"],
//["green", "noise_1024_4", "noise_1024_2", "cm2"],
//["blue", "noise_256_2", "noise_256_4", "cm5"],
//["lava", "noise_1024_2", "noise_1024_4", "cm8"],
//["water", "n1", "n2", "eau_cm"],
//["hexa", "n1", "n2", "hexa_cm"],
//["phasor_sand", "noise_sin_3", "noise_cos_3", "colormap_phasor_sand"],
//["phasor_sin", "noise_sin_1", "noise_cos_1", "colormap_phasor_sin"],
//["phasor_square", "noise_sin_1", "noise_cos_1", "colormap_phasor_square"],