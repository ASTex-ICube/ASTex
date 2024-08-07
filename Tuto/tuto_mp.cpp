#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>

#include <random>
#include<iomanip>


using namespace ASTex;

template <typename VEC, typename FUNC>
VEC applyFuncCompo(const VEC& u , const FUNC& f)
{ 
	VEC v;
	for (int i = 0; i < v.size(); ++i)
		v[i] = f(u[i]);
	return v;
}


template <typename VEC>
VEC vec_floor(const VEC& u)
{
	return applyFuncCompo(u, [](double d) { return std::floor(d); });
}

template <typename VEC>
VEC vec_fract(const VEC& u)
{
	return applyFuncCompo(u, [](double d) { return d - std::floor(d); });
}


template <typename VEC>
VEC vec_clamp01(const VEC& u)
{
	return applyFuncCompo(u, [](double d) { return std::max(0.0,std::min(1.0,d)); });
}


class MicroPatternsInput
{
	/**
	* @param uv coordinate to fetch
	* @param level level of mipmap in which to fetch
	* @return Eigen::Vector4d(noiseA-average(NoiseA),noiseB-average(noiseB),var(noiseA),var(noiseB)
	*/
	virtual Eigen::Vector4d fetch(const Eigen::Vector2d& uv, double level) const = 0;

	/**
	* @return average value (noiseA,noiseB)
	*/
	virtual Eigen::Vector2d fetch_average() const = 0;

	/**
	* @return width of based squared textures of mipmap
	*/
	virtual int width() const = 0;

	/**
	* @return return number of level in the mipmap
	*/
	virtual int depth() const = 0;
};



template <typename IMGD>
class BlurredColorMaps
{
	using TD = typename IMGD::PixelType;
	using TED = typename IMGD::DoublePixelEigen;

	const IMGD* color_map_;
	IMGD buffer1_;
	IMGD buffer2_;

	IMGD* inbuf_;
	IMGD* outbuf_;
	Region reg_;

	std::vector<double> convol_;
	int BR;

	std::vector<IMGD> blurred_;

	int nb_blurred_;

	int width_;

	void blur1D()
	{
		outbuf_->parallel_for_region_pixels(reg_, [&](TD& p, int x, int y)
		{
			TED pe = (convol_[0] * eigenPixel<double>(inbuf_->pixelAbsolute(x, y)));
			for (int i = 1; i <= BR; ++i)
				pe += convol_[i] * (eigenPixel<double>(inbuf_->pixelAbsolute(x - i, y))
					+ eigenPixel<double>(inbuf_->pixelAbsolute(x + i, y)));
			p = IMGD::itkPixel(pe);
		});
	}
	
	void average1D_X(const IMGD& img_in, IMGD& img_out)
	{
		TED z;
		z.setZero();
		std::vector<TED> aver{ std::size_t(width_), z };

		img_in.for_all_pixels([&](const TD& p, int x, int y)
			{
				aver[y] += eigenPixel<double>(p);
			});
		for (auto& a : aver)
			a /= width_;
		img_out.for_all_pixels([&](TD& p, int x, int y)
			{
				p = IMGD::itkPixel(aver[y]);
			});
	}


	void average1D_Y(const IMGD& img_in, IMGD& img_out)
	{
		TED z;
		z.setZero();
		std::vector<TED> aver{ std::size_t(width_), z };

		img_in.for_all_pixels([&](const TD& p, int x, int y)
			{
				aver[x] += eigenPixel<double>(p);
			});
		for (auto& a : aver)
			a /= width_;
		img_out.for_all_pixels([&](TD& p, int x, int y)
			{
				p = IMGD::itkPixel(aver[x]);
			});
	}

	void compute_column0()
	{
		int N_1 = nb_blurred_ - 1;
		blurred_[0].for_all_pixels([&](TD& p, int x, int y)
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

		//0 1 4 16 64 256 1024
		int ip = 0;
		for (int j = 1; j < N_1; ++j)
		{
			//int nb = 1 << (2 * (j - 1));
			int nb = 1 << (2 * (j - 1)); //WARNINNNNG
			while (ip < nb)
			{
				blur1D();
				std::swap(inbuf_, outbuf_);
				++ip;
			}

			blurred_[j * nb_blurred_].for_all_pixels([&](TD& p, int x, int y)
				{
					p = inbuf_->pixelAbsolute(y + BR, x);
				});
		}
	}

	void compute_rows()
	{
		int N_1 = nb_blurred_ - 1;
		for (int i = 0; i < N_1; ++i)
		{
			auto& firstColBlurred = blurred_[i * nb_blurred_];
			for (int y = 0; y < width_; ++y)
			{
				for (int k = 0; k < BR; ++k)
				{
					inbuf_->pixelAbsolute(k, y) = firstColBlurred.pixelAbsolute(0, y);
					outbuf_->pixelAbsolute(k, y) = firstColBlurred.pixelAbsolute(0, y);
				}
			}

			inbuf_->for_region_pixels(reg_, [&](TD& p, int x, int y)
				{
					p = firstColBlurred.pixelAbsolute(x - BR, y);
				});

			for (int y = 0; y < width_; ++y)
			{
				for (int k = 0; k < BR; ++k)
				{
					inbuf_->pixelAbsolute(width_ + BR + k, y) = firstColBlurred.pixelAbsolute(width_ - 1, y);
					outbuf_->pixelAbsolute(width_ + BR + k, y) = firstColBlurred.pixelAbsolute(width_ - 1, y);
				}
			}
			int ip = 0;
			for (int j = 1; j < N_1; ++j)
			{
				//int nb = 1 << (2 * (j - 1));
				int nb = 1 << (2 * (j - 1)); //WARNINNNNG
				while (ip < nb)
				{
					blur1D();
					std::swap(inbuf_, outbuf_);
					++ip;
				}
				//                   Region r =gen_region(0,gy,width_,1);
				blurred_[i * nb_blurred_ + j].for_all_pixels([&](TD& p, int x, int y)
					{
						p = inbuf_->pixelAbsolute(x + BR, y);
					});
			}
		}

		for (int i = 0; i < N_1; ++i)
			average1D_X(blurred_[i * nb_blurred_ + N_1 - 1], blurred_[i * nb_blurred_ + N_1]);

		for (int i = 0; i < nb_blurred_; ++i)
			average1D_Y(blurred_[(N_1 - 1) * nb_blurred_ + i], blurred_[N_1 * nb_blurred_ + i]);

	}


public:
	/**
	* @param colormap colormap double image gray or RGB or RGBA
	* @param convol the convolution mask, default {0.375 0.25 0.0625}
	* @param nb number of blur img computed in the 2 directions, last is average (defaault 7)
	*/
	BlurredColorMaps(const IMGD* colormap,
	 const std::vector<double>& convol = { 0.375,0.25,0.0625 },
	 int nb=8) :
		color_map_(colormap), convol_(convol),nb_blurred_(nb)
	{
		BR = convol.size() - 1;
		width_ = color_map_->width();
		buffer1_.initItk(width_ + 2 * BR, width_);
		buffer2_.initItk(width_ + 2 * BR, width_);
		reg_ = gen_region(BR, 0, width_, width_);
		inbuf_ = &buffer1_;
		outbuf_ = &buffer2_;
		blurred_.resize(nb_blurred_ * nb_blurred_);
		for(auto& bl: blurred_)
			bl.initItk(color_map_->width(), color_map_->width());
		compute_column0();
		compute_rows();
	}


	int table_width() const
	{
		return nb_blurred_;
	}

	void save(const std::string& base)
	{
		for (int i = 0; i < nb_blurred_; ++i)
		{
			for (int j = 0; j < nb_blurred_; ++j)
			{
				std::string fn = base + std::string("_cm_blur_X_Y.png");
				fn[fn.length() - 7] = '0' + i;
				fn[fn.length() - 5] = '0' + j;
				IO::save01_in_u8(blurred_[i * nb_blurred_ + j], fn);
			}
		}
	}


	TD fetch(const Eigen::Vector2d& uv, int i, int j) const
	{
		i = std::min(nb_blurred_ - 1, i);
		j = std::min(nb_blurred_ - 1, j);

		double w = double(width_);
		Eigen::Vector2d uvd =  w * uv - Eigen::Vector2d{ 0.5,0.5 } ;
		Eigen::Vector2d uvfl {std::floor(uvd[0]), std::floor(uvd[1])};
		Eigen::Vector2d a = uvd - uvfl;
		Eigen::Vector2i c{int(uvfl[0]),int(uvfl[1])};

		auto acces_edge = [&] (int xp, int yp)
		{
			int xx = xp < 0 ? 0 : ( xp >= width_ ? width_-1 : xp );
			int yy = yp < 0 ? 0 : ( yp >= width_ ? width_-1 : yp );
			return eigenPixel<double>(blurred_[j * nb_blurred_ + i].pixelAbsolute(xx, yy));
		};

		TED V1 = (1.0 - a[0]) * acces_edge(c[0],c[1])	+ a[0] * acces_edge(c[0]+1,c[1]);
		TED V2 = (1.0 - a[0]) * acces_edge(c[0],c[1]+1) + a[0] * acces_edge(c[0]+1,c[1]+1);
		TED V = (1.0 - a[1]) * V1 + a[1] * V2;

		return IMGD::itkPixel(V);
	}

	TD fetch(const Eigen::Vector2d& uv, double ui, double uj) const
	{

		auto acces_edge = [&](int xp, int yp, int ip, int jp)
		{
			int xx = xp < 0 ? 0 : (xp >= width_ ? width_ - 1 : xp);
			int yy = yp < 0 ? 0 : (yp >= width_ ? width_ - 1 : yp);
			return eigenPixel<double>(blurred_[jp*nb_blurred_+ip].pixelAbsolute(xx, yy));
		};

		double w = double(width_);
		Eigen::Vector2d uvd{ -0.5 + w * uv[0] , -0.5 + w * uv[1] };
		Eigen::Vector2d uvfl{ std::floor(uvd[0]), std::floor(uvd[1]) };
		Eigen::Vector2d a = uvd - uvfl;
		Eigen::Vector2i c{ int(uvfl[0]),int(uvfl[1]) };

		auto fetch_cm = [&](int xp, int yp)
		{
			int xx = xp < 0 ? 0 : (xp >= nb_blurred_ ? nb_blurred_ - 1 : xp);
			int yy = yp < 0 ? 0 : (yp >= nb_blurred_ ? nb_blurred_ - 1 : yp);

			TED V1 = (1.0 - a[0]) * acces_edge(c[0], c[1], xx, yy) + a[0] * acces_edge(c[0] + 1, c[1], xx, yy);
			TED V2 = (1.0 - a[0]) * acces_edge(c[0], c[1] + 1, xx, yy) + a[0] * acces_edge(c[0] + 1, c[1] + 1, xx, yy);
			TED V = (1.0 - a[1]) * V1 + a[1] * V2;

			return V;
		};

		double uifl = floor(ui);
		double ai = ui - uifl;
		int i = int(uifl);

		double ujfl = floor(uj);
		double aj = uj - ujfl;
		int j = int(ujfl);

		TED V1 = (1.0 - ai) * fetch_cm(i, j) + ai * fetch_cm(i + 1, j);
		TED V2 = (1.0 - ai) * fetch_cm(i, j + 1) + ai * fetch_cm(i + 1, j + 1);
		TED V = (1.0 - aj) * V1 + aj * V2;

		return IMGD::itkPixel(V);
	}
};


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

class Tiling_n_Blendinng
{

	void TriangleGrid(const Eigen::Vector2d& p_uv, Eigen::Vector3d& Bi, Eigen::Vector2i& vertex1, Eigen::Vector2i& vertex2, Eigen::Vector2i& vertex3)
	{
		const Eigen::Vector2d uv = p_uv * 2.0 * std::sqrt(3.0);

		Eigen::Matrix2d gridToSkewedGrid;
		gridToSkewedGrid << 1.0, -0.57735027,
							0.0, 01.15470054;

		Eigen::Vector2d skewedCoord = gridToSkewedGrid * uv;
		Eigen::Vector2d baseId{ std::floor(skewedCoord[0]), std::floor(skewedCoord[1]) };
		Eigen::Vector3d temp{ skewedCoord[0] - baseId[0], skewedCoord[1] - baseId[1], 0.0 };
		temp[2] = 1.0 - temp[0] - temp[1];


		if (temp[2] > 0.0)
		{																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																
			Bi = Eigen::Vector3d(temp[2], temp[1], temp[0]);
			Eigen::Vector2i ibaseId = baseId.cast<int>();
			vertex1 = ibaseId;
			vertex2 = ibaseId + Eigen::Vector2i(0, 1);
			vertex3 = ibaseId + Eigen::Vector2i(1, 0);
		}
		else
		{
			Bi = Eigen::Vector3d(-temp[2], 1.0 - temp[1], 1.0 - temp[0]);
			Eigen::Vector2i ibaseId = baseId.cast<int>();
			vertex1 = ibaseId + Eigen::Vector2i(1, 1);
			vertex2 = ibaseId + Eigen::Vector2i(1, 0);
			vertex3 = ibaseId + Eigen::Vector2i(0, 1);
		}
	}


	//original hash version
	Eigen::Vector2d hash(const Eigen::Vector2i& p)
	{
		Eigen::Matrix2d hashMat;
		hashMat << 127.1, 269.5,
				311.7, 183.3;

		Eigen::Vector2d q = hashMat * p.cast<double>();
		q[0] = std::sin(q[0]);
		q[1] = std::sin(q[1]);
		q *= 43758.5453;
		return Eigen::Vector2d(q[0] - std::floor(q[0]), q[1] - std::floor(q[1]));
	}

	public:
	Tiling_n_Blendinng()
	{}

	void operator()(const PackedInputNoises& noises, double level, const Eigen::Vector2d& uv, Eigen::Vector2d& mean, Eigen::Vector2d& variance)
	{
		Eigen::Vector3d B;
		Eigen::Vector2i  vertex1, vertex2, vertex3;
		TriangleGrid(uv, B,	vertex1, vertex2, vertex3);

		// Assign random offset to each triangle vertex
		Eigen::Vector2d uv1 = uv + hash(vertex1);
		Eigen::Vector2d uv2 = uv + hash(vertex2);
		Eigen::Vector2d uv3 = uv + hash(vertex3);

		Eigen::Vector4d n1 = noises.fetch(uv1, level);
		Eigen::Vector4d n2 = noises.fetch(uv2, level);
		Eigen::Vector4d n3 = noises.fetch(uv3, level);

		Eigen::Vector2d nu = noises.fetch_average();

		B.normalize();
		Eigen::Matrix<double, 2, 3> M;
		M << n1[0], n2[0], n3[0],
			 n1[1], n2[1], n3[1];
		mean = M * B +nu;
		mean[0] = std::max(0.0,std::min(1.0,mean[0]));
		mean[1] = std::max(0.0,std::min(1.0,mean[1]));

		B = B.cwiseProduct(B);
		Eigen::Matrix<double, 2, 3> S; 
		S << n1[2], n2[2], n3[2],
			 n1[3], n2[3], n3[3];
		variance = S * B;
	}
};



// 
//  void tile(const PackedInputNoises& noises, double level, const Eigen::Vector2d& uv, Eigen::Vector2d& mean, Eigen::Vector2d& variance)
//  noises : MIPMAP de N1,N2,ecart N1/moyenne ,ecart N2/moyenne 
//  

template <typename IMGD, typename FUNC>
void patterns(IMGD& pat, int w, const PackedInputNoises& noises, const BlurredColorMaps<IMGD>& bcm, double scale, FUNC tile)
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
void patterns_stat(IMGD& pat, int w, const PackedInputNoises& noises, const BlurredColorMaps<IMGD>& bcm, double scale, FUNC tile, double mult)
{
	using TD = typename IMGD::PixelType;

	std::vector < double> lx;
	lx.resize(bcm.table_width() + 4, 0.0);
	std::vector < double> ly;
	ly.resize(bcm.table_width() + 4, 0.0); 
	//{ bcm.table_width() + 4, 0.0 };

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

int main_n(int argc, char** argv, const std::vector<std::array < std::string, 4>>& names)
{
	bool bs = false;
	bs = std::string{ argv[argc-1] } == std::string{ "save" };

	int conf = std::atoi(argv[1]);
	std::string dir = "../../data/";
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
	int sz = std::atoi(argv[2]);;
	double scale = std::atof(argv[3]);;

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

int main_one_mipmap(int argc, char** argv, const std::vector<std::array < std::string, 4>>& names)
{
	bool bs = false;
	bs = std::string{ argv[argc - 1] } == std::string{ "save" };

	double scale = std::atof(argv[3]);
	int sz = std::atoi(argv[2]);

	int conf = std::atoi(argv[1]);
	std::string dir = "../../data/";
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
	double scale = std::atof(argv[3]);
	int sz = std::atoi(argv[2]);

	int conf = std::atoi(argv[1]);
	std::string dir = "../../data/";
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
		std::cout << "all config size scale_initial" << std::endl;
		std::cout << "mipmap config size scale_initial" << std::endl;
		std::cout << "stat config size scale" <<std::endl;

		std::cout << "config: ";
		for (int i = 0; i < names.size(); ++i)
		
			std::cout << i << ":" << names[i][0] << " / ";
		std::cout << std::endl;
		return EXIT_SUCCESS;
	}


	if (std::string(argv[1])=="mipmap")
		main_one_mipmap(argc-1, argv+1,names);
	if (std::string(argv[1])=="all")
		main_n(argc-1, argv+1,names);
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