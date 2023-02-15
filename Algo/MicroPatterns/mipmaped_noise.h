#ifndef _MP_MIPMAPED_NOISES_
#define _MP_MIPMAPED_NOISES_


#include <ASTex/image_gray.h>
#include <Algo/MicroPatterns/mp_inputs.h>
#include <Algo/MicroPatterns/vec_ops.h>


namespace ASTex
{

class PackedInputNoises: public MicroPatternsInput
{
	using IMG = ImageRGBAd;
	std::vector<IMG> mipmap;
	using PIXT = IMG::PixelType;
	Eigen::Vector2d	average_value_;

public:
	PackedInputNoises(const ImageGrayd& noiseA, const ImageGrayd& noiseB)
	{
		int nbl = std::log2(noiseA.width())+1;
		mipmap.resize(nbl);
		int lp2 = 1;
		int w = noiseA.width();

		mipmap[0].initItk(w, w);
		mipmap[0].parallel_for_all_pixels([&](PIXT& p, int x, int y)
		{
			double eA = noiseA.pixelAbsolute(x, y);
			double eB = noiseB.pixelAbsolute(x, y);
			p = ImageRGBAd::itkPixel(eA, eB, 0.0, 0.0);
		});
		
		for (int ni = 1; ni < nbl; ++ni)
		{
			lp2 *= 2;
			w /= 2;
			mipmap[ni].initItk(w, w);
			mipmap[ni].	parallel_for_all_pixels([&](PIXT& p, int x, int y)
				{
					int x2 = x * 2;
					int y2 = y * 2;	
					const auto& v0 = mipmap[ni - 1].pixelAbsolute(x2, y2);
					double eA = v0[0];
					double eB = v0[1];
					++x2;
					const auto& v1 = mipmap[ni - 1].pixelAbsolute(x2, y2);
					eA += v1[0];
					eB += v1[1];
					++y2;
					const auto& v3 = mipmap[ni - 1].pixelAbsolute(x2, y2);
					eA += v3[0];
					eB += v3[1];
					--x2;
					const auto& v2 = mipmap[ni - 1].pixelAbsolute(x2, y2);
					eA += v2[0];
					eB += v2[1];
					eA /= 4.0;
					eB /= 4.0;

					int xx = x * lp2;
					int yy = y * lp2;

					double varA = 0.0;
					double varB = 0.0;
					for (int j = 0; j < lp2; ++j)
					{
						for (int i = 0; i < lp2; ++i)
						{
							int xxx = xx + i;
							int yyy = yy + j;
							double v = noiseA.pixelAbsolute(xxx, yyy) - eA;
							varA += v * v;
							v = noiseB.pixelAbsolute(xxx, yyy) - eB;
							varB += v * v;
						}
					}
					double nb = double(lp2 * lp2);
					p = ImageRGBAd::itkPixel(eA, eB, varA / nb, varB / nb);
				});
			
		}
	
		const auto& p = mipmap.back().pixelAbsolute(0, 0);
		average_value_ = Eigen::Vector2d{p[0], p[1]};

		for (int ni = 0; ni < nbl; ++ni)
		{
			mipmap[ni].parallel_for_all_pixels([&](PIXT& p)
				{
					p[0] -= average_value_[0];
					p[1] -= average_value_[1];
				});
		}

	}

	int width() const override
	{
		return mipmap[0].width();
	}

	int depth() const override
	{
		return mipmap.size();
	}


	void save(const std::string& base)
	{
		for (int i = 0; i < mipmap.size(); ++i)
		{
			const auto& mi = mipmap[i];
			ImageGrayd imA;
			imA.initItk(mi.width(), mi.height());
			ImageGrayd imB;
			imB.initItk(mi.width(), mi.height());
			ImageGrayd imVA;
			imVA.initItk(mi.width(), mi.height());
			ImageGrayd imVB;
			imVB.initItk(mi.width(), mi.height());

			mi.parallel_for_all_pixels([&](const PIXT& p, int x, int y)
				{
					imA.pixelAbsolute(x, y) =  p[0] + average_value_[0];
					imB.pixelAbsolute(x, y) =  p[1] + average_value_[1];
					imVA.pixelAbsolute(x, y) = p[2] * 16.0;
					imVB.pixelAbsolute(x, y) = p[3] * 16.0;
				});

			std::string  fn = base + std::string("noiseA_X.png");
			fn[fn.length() - 5] = (i<=9) ? '0' + i : 'A' + i - 10;
			IO::save01_in_u8(imA, fn);

			fn = base + std::string("noiseB_X.png");
			fn[fn.length() - 5] = (i<=9) ? '0' + i : 'A' + i - 10;
			IO::save01_in_u8(imB, fn);

			fn = base + std::string("varA_X.png");
			fn[fn.length() - 5] = (i<=9) ? '0' + i : 'A' + i - 10;
			IO::save01_in_u8(imVA, fn);

			fn = base + std::string("varB_X.png");
			fn[fn.length() - 5] = (i<=9) ? '0' + i : 'A' + i - 10;
			IO::save01_in_u8(imVB, fn);

		}
	}


	Eigen::Vector4d fetch(const Eigen::Vector2d& uv, double level) const override
	{
		level = std::max(0.0,std::min(level, double(mipmap.size())));
		int l = int(std::floor(level));
		
		double w = double(mipmap[l].width());
		Eigen::Vector2d uvd =  w * vec_fract(uv) - Eigen::Vector2d{ 0.5,0.5 } ;
		Eigen::Vector2d uvfl = vec_floor(uvd);//{std::floor(uvd[0]), std::floor(uvd[1])};
		Eigen::Vector2d a = uvd - uvfl;
		Eigen::Vector2i c = uvfl.cast<int>();

		auto acces_repeat = [&] (int lp, int xp, int yp) ->Eigen::Vector4d
		{
			int w = mipmap[lp].width();
			int xx = xp < 0 ? w - 1 : ( xp >= w ? 0 : xp );
			int yy = yp < 0 ? w - 1 : ( yp >= w ? 0 : yp );
			return eigenPixel<double>(mipmap[lp].pixelAbsolute(xx, yy));
		};


		Eigen::Vector4d V1 = (1.0 - a[0]) * acces_repeat(l,c[0],c[1])	+ a[0] * acces_repeat(l,c[0]+1,c[1]);
		Eigen::Vector4d V2 = (1.0 - a[0]) * acces_repeat(l,c[0],c[1]+1) + a[0] * acces_repeat(l,c[0]+1,c[1]+1);
		Eigen::Vector4d V = (1.0 - a[1]) * V1 + a[1] * V2;

		int le = int(std::ceil(level));
		if (le == l)
			return  V;

		l = le; //for usage in access_repeat

		w = double(mipmap[le].width());
		uvd = Eigen::Vector2d {-0.5, -0.5} + w * vec_fract(uv);//Eigen::Vector2d { -0.5+w*(uv[0] - floor(uv[0])),-0.5+w*(uv[1] - floor(uv[1])) };
		uvfl = vec_floor(uvd); //Eigen::Vector2d { std::floor(uvd[0]), std::floor(uvd[1]) };
		a = uvd - uvfl;
		c = Eigen::Vector2i{int(uvfl[0]),int(uvfl[1])};

	
		// bilinear interpolation
		V1 = (1.0 - a[0]) * acces_repeat(l,c[0],c[1])	+ a[0] * acces_repeat(l,c[0]+1,c[1]);
		V2 = (1.0 - a[0]) * acces_repeat(l,c[0],c[1]+1)	+ a[0] * acces_repeat(l,c[0]+1,c[1]+1);
		Eigen::Vector4d VV = (1.0 - a[1]) * V1 + a[1] * V2;

		double b = level - std::floor(level);

		return (1.0 - b) * V + b * VV;
	}


	Eigen::Vector2d fetch_average() const override
	{
		return average_value_;
	}
};

} //namspace ASTex
#endif