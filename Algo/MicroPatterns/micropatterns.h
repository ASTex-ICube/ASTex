#ifndef _MP_PATTERNS_
#define _MP_PATTERNS_


#include <Algo/MicroPatterns/blurred_colorMap.h>
#include <Algo/MicroPatterns/mp_inputs.h>

namespace ASTex
{

// 
//  void tile(const PackedInputNoises& noises, double level, const Eigen::Vector2d& uv, Eigen::Vector2d& mean, Eigen::Vector2d& variance)
//  noises : MIPMAP de N1,N2,ecart N1/moyenne ,ecart N2/moyenne 
//  

template <typename IMGD, typename FUNC>
void patterns(IMGD& pat, int w, const MicroPatternsInput& noises, const BlurredColorMaps<IMGD>& bcm, double scale, FUNC tile)
{
	using TD = typename IMGD::PixelType;

	pat.initItk(w,w);
	double level = std::min(double(noises.depth()-1), std::log2(scale));
	pat.parallel_for_all_pixels([&](typename IMGD::PixelType& p, int x, int y)
		{
			Eigen::Vector3d color;
			Eigen::Vector2d uv_cm;
			Eigen::Vector2d sigm2;
			Eigen::Vector2d uv{ (double(x)+0.5) / (noises.width()), (double(y)+0.5) / (noises.width()) };
			Eigen::Vector2d uvs = uv * scale;
			tile(noises, level, uvs, uv_cm, sigm2);

			Eigen::Vector2d sigma{ sqrt(sigm2[0]),sqrt(sigm2[1]) };

			Eigen::Vector2d var = sigma * 256.0;
			Eigen::Vector2d flod = Eigen::Vector2d(std::log2(1.0 + var[0]), std::log2(1.0 + var[1]));

			p = bcm.fetch(uv_cm, std::min(bcm.table_width(),int(std::round(flod[0]))), std::min(bcm.table_width(),int(std::round(flod[1]))));
	//		p = bcm.fetch(uv_cm, flod[0], flod[1]);
		});
}



template <typename IMGD, typename FUNC>
void patterns_stat(IMGD& pat, int w, const MicroPatternsInput& noises, const BlurredColorMaps<IMGD>& bcm, double scale, FUNC tile, double mult)
{
	using TD = typename IMGD::PixelType;

	std::vector < double> lx;
	lx.resize(bcm.table_width() + 4, 0.0);
	std::vector < double> ly;
	ly.resize(bcm.table_width() + 4, 0.0); 

	pat.initItk(w, w);
	double level = std::min(double(noises.depth()), std::log2(scale));
	std::cout << "level = "<< level<< std::endl;
	std::cout << "level max = "<< noises.depth()<< std::endl;

	pat.for_all_pixels([&](typename IMGD::PixelType& p, int x, int y)
		{
			Eigen::Vector3d color;
			Eigen::Vector2d uv_cm;
			Eigen::Vector2d sigm2;
			Eigen::Vector2d uv{ (double(x)+0.5) / (noises.width()), (double(y)+0.5) / (noises.width()) };
			Eigen::Vector2d uvs = uv * scale;
			tile(noises, level, uvs, uv_cm, sigm2);

			Eigen::Vector2d sigma{ sqrt(sigm2[0]),sqrt(sigm2[1]) };
			Eigen::Vector2d var = sigma * mult;//256.0 ;
			Eigen::Vector2d flod = Eigen::Vector2d(std::log2(1.0 + var[0]), std::log2(1.0 + var[1]));

			p = bcm.fetch(uv_cm, int(std::round(flod[0])), int(std::round(flod[1])));
			//p = bcm.fetch(uv_cm, flod[0], flod[1]);

			double a = flod[0]-std::floor(flod[0]);
			lx[int(std::round(flod[0]))] += 1.0;
			ly[int(std::round(flod[1]))] += 1.0;	
		});

	for (int i = 0; i < 8; ++i)
	{
		double xl = 0.001 * int(double(lx[i]) / (0.001 * w * w));
		double yl = 0.001 * int(double(ly[i]) / (0.001 * w * w));
		std::cout << "L[" << i << "] = " << std::setprecision(3) << xl << " / " << yl << std::endl;
	}
}

} //namespace ASTex

#endif