#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>
#include <random>


using namespace ASTex;


// inline double convert_pix(uint8xp>=_)
// {
// 	return doubxp>=l) / 255.0;
// }

// inline itk::RGBPixel<double> convert_pix(const itk::RGBPixel<uint8_txp>=>)
// {
// 	double d[3] = { doubxp>=l[0]) / 255.0,doubxp>=l[1]) / 255.0,doubxp>=l[2]) / 255.0 };
// 	return itk::RGBPixel<double>(d);
// }



template <typename IMGD, int N>
class BlurringColorMaps
{
	using TD = typename IMGD::PixelType;
	using TED = typename IMGD::DoublePixelEigen;

	const IMGD* color_map_;
	IMGD buffer1_;
	IMGD buffer2_;

	IMGD* inbuf_;
	IMGD* outbuf_;
	Region reg_;

	//    std::array<double,BR+1> convol_;
	std::vector<double> convol_;
	int BR;

	std::array<std::array<IMGD, N>, N> res_;

	int width_;

	void blur1D()
	{
		if constexpr (std::is_same_v<TD, double>)
			outbuf_->parallel_for_region_pixels(reg_, [&](TD& p, int x, int y)
				{
					auto pe = convol_[0] * double(inbuf_->pixelAbsolute(x, y));
					for (int i = 1; i <= BR; ++i)
						pe += convol_[i] * (double(inbuf_->pixelAbsolute(x - i, y)) + double(inbuf_->pixelAbsolute(x + i, y)));
					p = IMGD::itkPixel(pe);
				});
		else
			outbuf_->parallel_for_region_pixels(reg_, [&](TD& p, int x, int y)
				{
					auto pe = (convol_[0] * eigenPixel<double>(inbuf_->pixelAbsolute(x, y))).eval();
					for (int i = 1; i <= BR; ++i)
						pe += convol_[i] * (eigenPixel<double>(inbuf_->pixelAbsolute(x - i, y))
							+ eigenPixel<double>(inbuf_->pixelAbsolute(x + i, y)));
					p = IMGD::itkPixel(pe);
				});
	}

public:
	BlurringColorMaps(const IMGD* colormap, const std::vector<double>& convol) :
		color_map_(colormap), convol_(convol)
	{
		BR = convol.size() - 1;
		width_ = color_map_->width();
		buffer1_.initItk(width_ + 2 * BR, width_);
		buffer2_.initItk(width_ + 2 * BR, width_);
		reg_ = gen_region(BR, 0, width_, width_);
		inbuf_ = &buffer1_;
		outbuf_ = &buffer2_;
		for (int j = 0; j < N; ++j)
			for (int i = 0; i < N; ++i)
				res_[j][i].initItk(color_map_->width(), color_map_->width());

	}

	void compute_column0()
	{
		res_[0][0].for_all_pixels([&](TD& p, int x, int y)
			{
				p = color_map_->pixelAbsolute(x, y);
			});

		for (int x = 0; x < width_; ++x)
		{
			auto bound = color_map_->pixelAbsolute(x, 0);
			for (int k = 0; k < BR; ++k)
			{
				inbuf_->pixelAbsolute(k, x) = bound;
				outbuf_->pixelAbsolute(k, x) = bound;
			}
		}
		inbuf_->for_region_pixels(reg_, [&](TD& p, int x, int y)
			{
				p = color_map_->pixelAbsolute(y, x - BR);
			});
		for (int x = 0; x < width_; ++x)
		{
			auto bound = color_map_->pixelAbsolute(x, width_ - 1);
			for (int k = 0; k < BR; ++k)
			{
				inbuf_->pixelAbsolute(width_ + BR + k, x) = color_map_->pixelAbsolute(x, width_ - 1);
				outbuf_->pixelAbsolute(width_ + BR + k, x) = color_map_->pixelAbsolute(x, width_ - 1);
			}
		}
//		IO::save01_in_u8(*inbuf_, "buff_in.png");

		int ip = 0;
		for (int p = 1; p < N; ++p )
		{
			int nb = 1 << (2 * (p - 1));
			while (ip < nb)
			{
				blur1D();
				std::swap(inbuf_, outbuf_);
				++ip;
			}
			// std::string fn{ "buff_in_X.png" };
			// fn[fn.length()-5] = 'xp>=0;
			// IO::save01_in_u8(*inbuf_, fn);

			//                Region r =gen_region(0,0,1,width_);
			
			res_[p][0].for_all_pixels([&](TD& p, int x, int y)
				{
					p = inbuf_->pixelAbsolute(y + BR, x);
				});
		}
	}

	void compute_rows()
	{
		for (int i = 0; i < N; ++i)
		{
			for (int y = 0; y < width_; ++y)
			{
				for (int k = 0; k < BR; ++k)
				{
					inbuf_->pixelAbsolute(k, y) = res_[i][0].pixelAbsolute(0, y);
					outbuf_->pixelAbsolute(k, y) = res_[i][0].pixelAbsolute(0, y);
				}
			}

			inbuf_->for_region_pixels(reg_, [&](TD& p, int x, int y)
				{
					p = res_[i][0].pixelAbsolute(x - BR, y);
				});

			for (int y = 0; y < width_; ++y)
			{
				for (int k = 0; k < BR; ++k)
				{
					inbuf_->pixelAbsolute(width_ + BR + k, y) = res_[i][0].pixelAbsolute(width_ - 1, y);
					outbuf_->pixelAbsolute(width_ + BR + k, y) = res_[i][0].pixelAbsolute(width_ - 1, y);
				}
			}
			int ip = 0;
			for (int j = 1; j < N; ++j)
			{
				int nb = 1 << (2 * (j - 1));
				while (ip < nb)
				{
					blur1D();
					std::swap(inbuf_, outbuf_);
					++ip;
				}
				//                   Region r =gen_region(0,gy,width_,1);
				res_[i][j].for_all_pixels([&](TD& p, int x, int y)
					{
						p = inbuf_->pixelAbsolute(x + BR, y);
					});
			}
		}

	}

	double max_nbc() const
	{
		return N;
	}

	void save()
	{
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				std::string fn("cm_blur_X_Y.png");
				fn[fn.length() - 7] = '0' + i;
				fn[fn.length() - 5] = '0' + j;
				IO::save01_in_u8(res_[i][j], fn);
			}
		}
	}


	TD fetch(const Eigen::Vector2d& uv, int i, int j) const
	{
		i = std::min(N - 1, i);
		j = std::min(N - 1, j);
		const auto& cm = res_[j][i];

		double w = double(width_);
		Eigen::Vector2d uvd{-0.5 + w * uv[0] , -0.5 + w * uv[1] };
		Eigen::Vector2d uvfl {std::floor(uvd[0]), std::floor(uvd[1])};
		Eigen::Vector2d a = uvd - uvfl;
		Eigen::Vector2i c{int(uvfl[0]),int(uvfl[1])};

		auto acces_repeat = [&] (int xp, int yp)
		{
			int xx = xp < 0 ? 0 : ( xp >= width_ ? width_-1 : xp );
			int yy = yp < 0 ? 0 : ( yp >= width_ ? width_-1 : yp );
			return eigenPixel<double>(res_[j][i].pixelAbsolute(xx, yy));
		};

		

		TED V1 = (1.0 - a[0]) * acces_repeat(c[0],c[1])	+ a[0] * acces_repeat(c[0]+1,c[1]);
		TED V2 = (1.0 - a[0]) * acces_repeat(c[0],c[1]+1) + a[0] * acces_repeat(c[0]+1,c[1]+1);
		TED V = (1.0 - a[1]) * V1 + a[1] * V2;

	//	std::cout <<j<<","<<i<< " ->" << c[0]<<" , "<<c[1]<< " -> "<< V.transpose()<<std::endl;
		
		
		return IMGD::itkPixel(V);
	}

};


class PackedInputNoises
{
	using IMG = ImageRGBAd;
	using PIXT = IMG::PixelType;
	std::vector<IMG> mipmap;
public:
	PackedInputNoises(const ImageGrayd& noiseA, const ImageGrayd& noiseB)
	{
		int nbl = std::log2(noiseA.width())+1;
		mipmap.resize(nbl);
		int lp2 = 1;
		int w = noiseA.width();
		for (int ni = 0; ni < nbl; ++ni)
		{
			mipmap[ni].initItk(w, w);
			mipmap[ni].parallel_for_all_pixels([&](PIXT& p, int x, int y)
				{

					double eA;
					double eB;
					if (ni == 0)
					{
						eA = noiseA.pixelAbsolute(x, y);
						eB = noiseB.pixelAbsolute(x, y);
					}
					else
					{
						int x2 = x * 2;
						int y2 = y * 2;

						const auto& v0 = mipmap[ni - 1].pixelAbsolute(x2, y2);
						eA = v0[0];
						eB = v0[1];
						const auto& v1 = mipmap[ni - 1].pixelAbsolute(x2 + 1, y2);
						eA += v1[0];
						eB += v1[1];
						const auto& v2 = mipmap[ni - 1].pixelAbsolute(x2, y2 + 1);
						eA += v2[0];
						eB += v2[1];
						const auto& v3 = mipmap[ni - 1].pixelAbsolute(x2 + 1, y2 + 1);
						eA += v3[0];
						eB += v3[1];
						eA /= 4.0;
						eB /= 4.0;
					}

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
			lp2 *= 2;
			w /= 2;
		}
	}

	int width() const
	{
		return mipmap[0].width();
	}

	void save()
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
					imA.pixelAbsolute(x, y) =  p[0];
					imB.pixelAbsolute(x, y) =  p[1];
					imVA.pixelAbsolute(x, y) = p[2] * 64.0;
					imVB.pixelAbsolute(x, y) = p[3] * 64.0;
				});

			std::string fn("noiseA_X.png");
			fn[fn.length() - 5] = (i<=9) ? '0' + i : 'A' + i - 10;
			IO::save01_in_u8(imA, fn);

			fn = std::string("noiseB_X.png");
			fn[fn.length() - 5] = (i<=9) ? '0' + i : 'A' + i - 10;
			IO::save01_in_u8(imB, fn);

			fn = std::string("varA_X.png");
			fn[fn.length() - 5] = (i<=9) ? '0' + i : 'A' + i - 10;
			IO::save01_in_u8(imVA, fn);

			fn = std::string("varB_X.png");
			fn[fn.length() - 5] = (i<=9) ? '0' + i : 'A' + i - 10;
			IO::save01_in_u8(imVB, fn);

		}
	}

	Eigen::Vector4d fetch_nearest(const Eigen::Vector2d& uv, double level) const
	{
		level = std::max(level, 0.0);
		level = std::min(level, double(mipmap.size()-1));

		int l = int(std::round(level));
		double w = double(mipmap[l].width());
		Eigen::Vector2d uvd{-0.5+w*(uv[0] - floor(uv[0])),-0.5+w*(uv[1] - floor(uv[1])) };
		int x = std::round(uvd[0]);
		int y = std::round(uvd[0]);

		return eigenPixel<double>(mipmap[l].pixelAbsolute(x, y));
	}

	Eigen::Vector4d fetch(const Eigen::Vector2d& uv, double level) const
	{
		level = std::max(level, 0.0);
		level = std::min(level, double(mipmap.size()));

		int l = int(std::floor(level));
		double w = double(mipmap[l].width());
		Eigen::Vector2d uvd{-0.5 + w * (uv[0] - floor(uv[0])), -0.5 + w * (uv[1] - floor(uv[1]))};
		Eigen::Vector2d uvfl {std::floor(uvd[0]), std::floor(uvd[1])};
		Eigen::Vector2d a = uvd - uvfl;
		Eigen::Vector2i c{int(uvfl[0]),int(uvfl[1])};

		auto acces_repeat = [this] (int lp, int xp, int yp) ->Eigen::Vector4d
		{
			int xx = xp < 0 ? mipmap[lp].width() - 1 : ( xp >= mipmap[lp].width() ? 0 : xp );
			int yy = yp < 0 ? mipmap[lp].width() - 1 : ( yp >= mipmap[lp].width() ? 0 : yp );
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
		uvd = Eigen::Vector2d { -0.5+w*(uv[0] - floor(uv[0])),-0.5+w*(uv[1] - floor(uv[1])) };
		uvfl = Eigen::Vector2d { std::floor(uvd[0]), std::floor(uvd[1]) };
		a = uvd - uvfl;
		c = Eigen::Vector2i{int(uvfl[0]),int(uvfl[1])};

		
		// bilinear interpolation
		V1 = (1.0 - a[0]) * acces_repeat(l,c[0],c[1])	+ a[0] * acces_repeat(l,c[0]+1,c[1]);
		V2 = (1.0 - a[0]) * acces_repeat(l,c[0],c[1]+1)	+ a[0] * acces_repeat(l,c[0]+1,c[1]+1);
		Eigen::Vector4d VV = (1.0 - a[1]) * V1 + a[1] * V2;


		double b = level - std::floor(level);

		return (1.0 - b) * V + b * VV;
	}


	Eigen::Vector4d fetch_average() const
	{
		const auto& p = mipmap.back().pixelAbsolute(0, 0);
		return Eigen::Vector4d	{p[0], p[1], p[2], p[3]};
	}
};


void TriangleGrid(const Eigen::Vector2d& p_uv, Eigen::Vector3d& Bi, Eigen::Vector2d& vertex1, Eigen::Vector2d& vertex2, Eigen::Vector2d& vertex3)
{
	const Eigen::Vector2d uv = p_uv * 2.0 * std::sqrt(3.0);
	Eigen::Matrix2d gridToSkewedGrid;
	gridToSkewedGrid << 1.0, 0.0, -0.57735027, 1.15470054;
	Eigen::Vector2d skewedCoord = gridToSkewedGrid * uv;
	Eigen::Vector2d baseId{ std::floor(skewedCoord[0]), std::floor(skewedCoord[1]) };
	Eigen::Vector3d temp{ skewedCoord[0] - baseId[0], skewedCoord[1] - baseId[1], 0.0 };
	temp[2] = 1.0 - temp[0] - temp[1];

	if (temp[2] > 0.0)
	{
		Bi = Eigen::Vector3d(temp[2], temp[1], temp[0]);
		vertex1 = baseId;
		vertex2 = baseId + Eigen::Vector2d(0.0, 1.0);
		vertex3 = baseId + Eigen::Vector2d(1.0, 0.0);
	}
	else
	{
		Bi = Eigen::Vector3d(-temp[2], 1.0 - temp[1], 1.0 - temp[0]);
		vertex1 = baseId + Eigen::Vector2d(1.0, 1.0);
		vertex2 = baseId + Eigen::Vector2d(1.0, 0.0);
		vertex3 = baseId + Eigen::Vector2d(0.0, 1.0);
	}
}

//std::random_device rd;  // Will be used to obtain a seed for the random number engine
//std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
//std::uniform_real_distribution<> dis(0.0, 1.0);

Eigen::Vector2d hash(const Eigen::Vector2d& p)
{
//	return Eigen::Vector2d(dis(gen), dis(gen));
	Eigen::Matrix2d hashMat;
	hashMat << 127.1, 311.7, 269.5, 183.3;
	Eigen::Vector2d q = hashMat * p ;
	q[0] = std::sin(q[0]);
	q[1] = std::sin(q[1]);
	q *= 43758.5453;
	return Eigen::Vector2d(q[0] - std::floor(q[0]), q[1] - std::floor(q[1]));
}

void Tile_n_blend(const PackedInputNoises& noises, double level, const Eigen::Vector2d& uv, Eigen::Vector2d& mean, Eigen::Vector2d& variance)
{
	Eigen::Vector3d B;
	Eigen::Vector2d  vertex1, vertex2, vertex3;
	TriangleGrid(uv, B,//b1, b2, b3,
		vertex1, vertex2, vertex3);

	// Assign random offset to each triangle vertex
	Eigen::Vector2d uv1 = uv + hash(vertex1);
	Eigen::Vector2d uv2 = uv + hash(vertex2);
	Eigen::Vector2d uv3 = uv + hash(vertex3);
	
	Eigen::Vector4d n1 = noises.fetch(uv1, level);
	Eigen::Vector4d n2 = noises.fetch(uv2, level);
	Eigen::Vector4d n3 = noises.fetch(uv3, level);

	Eigen::Vector4d nu = noises.fetch_average();

	B.normalize();
	Eigen::Matrix<double, 2, 3> M;// = mat3x2(n1.xy - nu, n2.xy - nu, n3.xy - nu);
	M << n1[0] - nu[0], n2[0] - nu[0], n3[0] - nu[0],
		 n1[1] - nu[1], n2[1] - nu[1], n3[1] - nu[1];
	mean = M * B + Eigen::Vector2d(nu[0], nu[1]);
	mean[0] = std::max(0.0,std::min(1.0,mean[0]));
	mean[1] = std::max(0.0,std::min(1.0,mean[1]));

	B = B.cwiseProduct(B);
	Eigen::Matrix<double, 2, 3> S; // = mat3x2(n1.zw, n2.zw, n3.zw) 
	S << n1[2], n2[2], n3[2],
		 n1[3], n2[3], n3[3];
	variance = S * B;
}

void repete(const PackedInputNoises& noises, double level, const Eigen::Vector2d& uv, Eigen::Vector2d& mean, Eigen::Vector2d& variance)
{
	Eigen::Vector4d n = noises.fetch(uv, level);
	mean = Eigen::Vector2d(n[0], n[1]);
	variance = Eigen::Vector2d(n[3], n[4]);
}



template <typename IMGD, int N, typename FUNC>
void patterns(IMGD& pat, int w, const PackedInputNoises& noises, const BlurringColorMaps<IMGD, N >& bcm, double scale, FUNC tile)
{
	using TD = typename IMGD::PixelType;

	Eigen::Vector2d uv_cm;
	Eigen::Vector2d sigm2;

	pat.initItk(w,w);
	double level = std::log2(scale);
	pat.for_all_pixels([&](typename IMGD::PixelType& p, int x, int y)
		{
			Eigen::Vector2d uv{ double(x) / (noises.width()), double(y) / (noises.width()) };
			Eigen::Vector2d uvs = uv*scale;
			//Tile_n_blend(noises, level, uvs, uv_cm, sigm2);
			tile(noises, level, uvs, uv_cm, sigm2);

// if constexpr (std::is_same_v<TD, double>)
// 	p = IMGD::itkPixel(uv_cm[0]);
// else
// 	p = IMGD::itkPixel(uv_cm[0],uv_cm[0],uv_cm[0]);
// return;
			Eigen::Vector2d sigma{ sqrt(sigm2[0]),sqrt(sigm2[1]) };
			Eigen::Vector2d var{ sigma[0] * 256.0, sigma[1] * 256.0 };
			var[0] = std::max(1.0, var[0]);
			var[1] = std::max(1.0, var[1]);
			double pnbcm = std::pow(2.0, bcm.max_nbc()) - 0.1;
			var[0] = std::min(pnbcm, var[0]);
			var[1] = std::min(pnbcm, var[1]);


			Eigen::Vector2d flod = Eigen::Vector2d(std::log2(var[0]), std::log2(var[1]));
			// flod[0] = std::max(0.0, flod[0]);
			// flod[1] = std::max(0.0, flod[1]);
			// double pnbcm =  bcm.max_nbc() - 1.0;
			// var[0] = std::min(pnbcm, flod[0]);
			// var[1] = std::min(pnbcm, flod[1]);



			//Eigen::Vector2i ilod = ivec2(floor(flod));
			//Eigen::Vector2d t = fract(flod);
			//Eigen::Vector2i ilodr = ivec2(round(flod));

			int ilodrx = int(std::round(flod[0]));
			int ilodry = int(std::round(flod[1]));
			p = bcm.fetch(uv_cm, ilodrx, ilodry);
			//p = IMGD::itkPixel(flod[0]/8.0,flod[1]/8.0,0.0);
		});

}



#define COLOR

#ifdef  COLOR
using T_IMG = ImageRGBu8;
using T_IMG_D = ImageRGBd;
#else
using T_IMG = ImageGray8;
using T_IMG_D = ImageGrayd;
#endif

int main(int argc, char** argv)
{
	T_IMG_D cm;
	IO::loadu8_in_01(cm, std::string(argv[1]));
	BlurringColorMaps<T_IMG_D, 7> bcm{ &cm,{0.375,0.25,0.0625} };
	bcm.compute_column0();
	bcm.compute_rows();
	std::cout << "Color maps blurred" << std::endl;
	//bcm.save();

	ImageGrayd nA;
	IO::loadu8_in_01(nA, std::string(argv[2]));
	ImageGrayd nB;
	IO::loadu8_in_01(nB, std::string(argv[3]));

	PackedInputNoises pin(nA, nB);
	std::cout << "noise var mipmaped" << std::endl;
	//pin.save();
	int SZGEN =4096;
	double scale = 1.0;
	T_IMG_D img_patterns;
	for(int i=0;i<8;i++)
	{
		patterns(img_patterns, SZGEN, pin, bcm, scale,Tile_n_blend);	
		std::string fn("patterns0.png");
		fn[fn.length() - 5] = (i<=9) ? '0' + i : 'A' + i - 10;
		IO::save01_in_u8(img_patterns, fn);
		std::cout << fn<<" generated" << std::endl;
		SZGEN/=2;
		scale*=2.0;
	}

	// patterns(img_patterns, SZGEN, pin, bcm, 1.0,Tile_n_blend);
	// std::cout << "micro patterns generated" << std::endl;
	// IO::save01_in_u8(img_patterns, "patterns0.png");

	// patterns(img_patterns, SZGEN, pin, bcm, 2.0,Tile_n_blend);
	// std::cout << "micro patterns generated" << std::endl;
	// IO::save01_in_u8(img_patterns, "patterns1.png");
	// patterns(img_patterns, SZGEN, pin, bcm, 4.0,Tile_n_blend);
	// std::cout << "micro patterns generated" << std::endl;
	// IO::save01_in_u8(img_patterns, "patterns2.png");
	// patterns(img_patterns, SZGEN, pin, bcm, 8.0,Tile_n_blend);
	// std::cout << "micro patterns generated" << std::endl;
	// IO::save01_in_u8(img_patterns, "patterns3.png");
	// patterns(img_patterns, SZGEN, pin, bcm, 16.0,Tile_n_blend);
	// std::cout << "micro patterns generated" << std::endl;
	// IO::save01_in_u8(img_patterns, "patterns4.png");
	// patterns(img_patterns, SZGEN, pin, bcm, 32.0,Tile_n_blend);
	// std::cout << "micro patterns generated" << std::endl;
	// IO::save01_in_u8(img_patterns, "patterns5.png");
	// patterns(img_patterns, SZGEN, pin, bcm, 64.0,Tile_n_blend);
	// std::cout << "micro patterns generated" << std::endl;
	// IO::save01_in_u8(img_patterns, "patterns6.png");



	return EXIT_SUCCESS;
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