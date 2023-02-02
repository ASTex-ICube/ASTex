#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>

#include <random>
#include<iomanip>


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

template <typename VEC, typename FUNC>
VEC applyFuncCompo(const VEC& u , const FUNC& f)
{ 
	VEC v;
	for (int i = 0; i < v.size(); ++i)
		v[i] = f(u[i]);
	return v;
}


template <typename VEC>
VEC vecfloor(const VEC& u)
{
	return applyFuncCompo(u, [](double d) { return std::floor(d); });
}

template <typename VEC>
VEC vecfract(const VEC& u)
{
	return applyFuncCompo(u, [](double d) { return d - std::floor(d); });
}


class PreHashed
{
	int width_;
	int shift_;
	int nb_;
	std::vector < Eigen::Vector2d > rands_;

public:
	void init(uint32_t seed)
	{
		std::mt19937 gen{seed};
		std::uniform_real_distribution<> dis(0.0, 1.0);
		rands_.clear();
		for (int i = 0; i < nb_; ++i)
			rands_.emplace_back(dis(gen), dis(gen));
	std::cout << nb_ << " random couples gnerated"<<std::endl;
	}

	PreHashed(int w, uint32_t seed) : width_(6*w+1),shift_(2*w)
	{
		nb_ = width_ * width_;
		rands_.reserve(nb_);
		init(seed);
	}
	
	const Eigen::Vector2d& get(const Eigen::Vector2i& p) const
	{
		int i = p[0] + shift_;
		int j = p[1] + shift_;
		return rands_[j * width_ + i];
	}
};

template <typename IMGD>
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

//	std::array<std::array<IMGD, N>, N> res_;
	std::vector<IMGD> blurred_;

	int N;

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

public:
	BlurringColorMaps(const IMGD* colormap, const std::vector<double>& convol, int nb) :
		color_map_(colormap), convol_(convol),N(nb+1)
	{
		BR = convol.size() - 1;
		width_ = color_map_->width();
		buffer1_.initItk(width_ + 2 * BR, width_);
		buffer2_.initItk(width_ + 2 * BR, width_);
		reg_ = gen_region(BR, 0, width_, width_);
		inbuf_ = &buffer1_;
		outbuf_ = &buffer2_;
		blurred_.resize(N * N);
		for(auto& bl: blurred_)
			bl.initItk(color_map_->width(), color_map_->width());
	}

	void compute_column0()
	{
		int N_1 = N - 1;
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
		for (int j = 1; j < N_1; ++j )
		{
			int nb = 1 << (2 * (j - 1));
			while (ip < nb)
			{
				blur1D();
				std::swap(inbuf_, outbuf_);
				++ip;
			}
			
			blurred_[j*N].for_all_pixels([&](TD& p, int x, int y)
				{
					p = inbuf_->pixelAbsolute(y + BR, x);
				});
		}
	}

	void compute_rows()
	{
		int N_1 = N - 1;
		for (int i = 0; i < N_1; ++i)
		{
			auto& firstColBlurred = blurred_[i * N];
			for (int y = 0; y < width_; ++y)
			{
				for (int k = 0; k < BR; ++k)
				{
					inbuf_->pixelAbsolute(k, y) =  firstColBlurred.pixelAbsolute(0, y);
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
				int nb = 1 << (2 * (j - 1));
				while (ip < nb)
				{
					blur1D();
					std::swap(inbuf_, outbuf_);
					++ip;
				}
				//                   Region r =gen_region(0,gy,width_,1);
				blurred_[i * N + j].for_all_pixels([&](TD& p, int x, int y)
					{
						p = inbuf_->pixelAbsolute(x + BR, y);
					});
			}
		}

		for (int i = 0; i < N_1; ++i)
			average1D_X(blurred_[i * N + N_1 - 1], blurred_[i * N + N_1]);

		for (int i = 0; i < N; ++i)
			average1D_Y(blurred_[(N_1-1) * N + i], blurred_[N_1 * N + i]);

	}

	int table_width() const
	{
		return N;
	}

	void save(const std::string& base)
	{
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				std::string fn = base + std::string("_cm_blur_X_Y.png");
				fn[fn.length() - 7] = '0' + i;
				fn[fn.length() - 5] = '0' + j;
				IO::save01_in_u8(blurred_[i * N + j], fn);
			}
		}
	}


	TD fetch(const Eigen::Vector2d& uv, int i, int j) const
	{
		i = std::min(N - 1, i);
		j = std::min(N - 1, j);

		double w = double(width_);
		Eigen::Vector2d uvd{-0.5 + w * uv[0] , -0.5 + w * uv[1] };
		Eigen::Vector2d uvfl {std::floor(uvd[0]), std::floor(uvd[1])};
		Eigen::Vector2d a = uvd - uvfl;
		Eigen::Vector2i c{int(uvfl[0]),int(uvfl[1])};

		auto acces_edge = [&] (int xp, int yp)
		{
			int xx = xp < 0 ? 0 : ( xp >= width_ ? width_-1 : xp );
			int yy = yp < 0 ? 0 : ( yp >= width_ ? width_-1 : yp );
			return eigenPixel<double>(blurred_[j * N + i].pixelAbsolute(xx, yy));
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
			return eigenPixel<double>(blurred_[jp*N+ip].pixelAbsolute(xx, yy));
		};

		double w = double(width_);
		Eigen::Vector2d uvd{ -0.5 + w * uv[0] , -0.5 + w * uv[1] };
		Eigen::Vector2d uvfl{ std::floor(uvd[0]), std::floor(uvd[1]) };
		Eigen::Vector2d a = uvd - uvfl;
		Eigen::Vector2i c{ int(uvfl[0]),int(uvfl[1]) };

		auto fetch_cm = [&](int xp, int yp)
		{
			int xx = xp < 0 ? 0 : (xp >= N ? N - 1 : xp);
			int yy = yp < 0 ? 0 : (yp >= N ? N - 1 : yp);

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
					imA.pixelAbsolute(x, y) =  p[0];
					imB.pixelAbsolute(x, y) =  p[1];
					imVA.pixelAbsolute(x, y) = p[2] * 64.0;
					imVB.pixelAbsolute(x, y) = p[3] * 64.0;
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


	Eigen::Vector4d fetch(const Eigen::Vector2d& uv, double level) const
	{
		level = std::max(0.0,std::min(level, double(mipmap.size())));
		int l = int(std::floor(level));
		Eigen::Vector2d uvm{uv[0] - floor(uv[0]), uv[1] - floor(uv[1])};
		double w = double(mipmap[l].width());
		Eigen::Vector2d uvd = Eigen::Vector2d{ -0.5,-0.5 } + w * uvm;
		Eigen::Vector2d uvfl {std::floor(uvd[0]), std::floor(uvd[1])};
		Eigen::Vector2d a = uvd - uvfl;
		Eigen::Vector2i c = uvfl.cast<int>();

		auto acces_repeat = [&] (int lp, int xp, int yp) ->Eigen::Vector4d
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

class Tiling_n_Blendinng
{

	PreHashed prH_;

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
//		std::cout << p_uv.transpose()  <<" => " << vertex1.transpose() << " & " << vertex2.transpose() << " & " << vertex3.transpose() << std::endl;
	}


	//original hash version
	Eigen::Vector2d hash(const Eigen::Vector2i& p)
	{
		Eigen::Matrix2d hashMat;
		hashMat << 127.1, 269.5,
				311.7, 183.3;

		//hashMat << 127.1, 311.7, 269.5, 183.3;
		Eigen::Vector2d q = hashMat * p.cast<double>();
		q[0] = std::sin(q[0]);
		q[1] = std::sin(q[1]);
		q *= 43758.5453;
		return Eigen::Vector2d(q[0] - std::floor(q[0]), q[1] - std::floor(q[1]));
	}

	public:
	Tiling_n_Blendinng(double scale_max, bool default_seed): prH_(int(std::ceil(scale_max)),default_seed)
	{
	}

	void operator()(const PackedInputNoises& noises, double level, const Eigen::Vector2d& uv, Eigen::Vector2d& mean, Eigen::Vector2d& variance)
	{
		Eigen::Vector3d B;
		Eigen::Vector2i  vertex1, vertex2, vertex3;
		TriangleGrid(uv, B,	vertex1, vertex2, vertex3);

		// Assign random offset to each triangle vertex
		//Eigen::Vector2d uv1 = uv + hash(vertex1);
		//Eigen::Vector2d uv2 = uv + hash(vertex2);
		//Eigen::Vector2d uv3 = uv + hash(vertex3);
		Eigen::Vector2d uv1 = uv + prH_.get(vertex1);
		Eigen::Vector2d uv2 = uv + prH_.get(vertex2);
		Eigen::Vector2d uv3 = uv + prH_.get(vertex3);


		Eigen::Vector4d n1 = noises.fetch(uv1, level);
		Eigen::Vector4d n2 = noises.fetch(uv2, level);
		Eigen::Vector4d n3 = noises.fetch(uv3, level);

		Eigen::Vector4d nu = noises.fetch_average();

		B.normalize();
		Eigen::Matrix<double, 2, 3> M;
		M << n1[0] - nu[0], n2[0] - nu[0], n3[0] - nu[0],
			 n1[1] - nu[1], n2[1] - nu[1], n3[1] - nu[1];
		mean = M * B + Eigen::Vector2d(nu[0], nu[1]);
		mean[0] = std::max(0.0,std::min(1.0,mean[0]));
		mean[1] = std::max(0.0,std::min(1.0,mean[1]));

		B = B.cwiseProduct(B);
		Eigen::Matrix<double, 2, 3> S; 
		S << n1[2], n2[2], n3[2],
			 n1[3], n2[3], n3[3];
		variance = S * B;
	}
};



void repete(const PackedInputNoises& noises, double level, const Eigen::Vector2d& uv, Eigen::Vector2d& mean, Eigen::Vector2d& variance)
{
	Eigen::Vector4d n = noises.fetch(uv, level);
	mean = Eigen::Vector2d(n[0], n[1]);
	variance = Eigen::Vector2d(n[2], n[3]);
}


// 
//  void tile(const PackedInputNoises& noises, double level, const Eigen::Vector2d& uv, Eigen::Vector2d& mean, Eigen::Vector2d& variance)
//  noises : MIPMAP de N1,N2,ecart N1/moyenne ,ecart N2/moyenne 
//  



template <typename IMGD, typename FUNC>
void patterns(IMGD& pat, int w, const PackedInputNoises& noises, const BlurringColorMaps<IMGD>& bcm, double scale, FUNC tile)
{
	using TD = typename IMGD::PixelType;

	//std::array < long, N> lx;
	//std::array < long, N> ly;
	//for (int i = 0; i < N; ++i)
	//{
	//	lx[i] = 0;
	//	ly[i] = 0;
	//}


	pat.initItk(w,w);
	double level = std::log2(scale);
	pat.parallel_for_all_pixels([&](typename IMGD::PixelType& p, int x, int y)
		{
			Eigen::Vector3d color;
			Eigen::Vector2d uv_cm;
			Eigen::Vector2d sigm2;
			Eigen::Vector2d uv{ double(x) / (noises.width()), double(y) / (noises.width()) };
			Eigen::Vector2d uvs = uv * scale;
			//Tile_n_blend(noises, level, uvs, uv_cm, sigm2);
			tile(noises, level, uvs, uv_cm, sigm2);

			Eigen::Vector2d sigma{ sqrt(sigm2[0]),sqrt(sigm2[1]) };
			//Eigen::Vector2d sigma{ sigm2[0],sigm2[1] };

			Eigen::Vector2d var = sigma * 256.0 * 2;
			//var[0] = std::max(1.0, var[0]);
			//var[1] = std::max(1.0, var[1]);
			//double pnbcm = std::pow(2.0, double(bcm.table_width())) - 0.1;
			//var[0] = std::min(pnbcm, var[0]);
			//var[1] = std::min(pnbcm, var[1]);


			Eigen::Vector2d flod = Eigen::Vector2d(std::log2(1.0 + var[0]), std::log2(1.0 + var[1]));

	//		p = bcm.fetch(uv_cm, int(std::round(flod[0])), int(std::round(flod[1])));

			p = bcm.fetch(uv_cm, flod[0], flod[1]);

	//		p = IMGD::itkPixel(uv_cm[0],uv_cm[0],uv_cm[0]);

			//lx[ilodrx]++;
			//ly[ilodry]++;	
		});

	//for (int i = 0; i < lx.size(); ++i)
	//{
	//	double xl = 0.001 * int(double(lx[i]) / (0.001 * w * w));
	//	double yl = 0.001 * int(double(ly[i]) / (0.001 * w * w));
	//	std::cout << "L[" << i << "] = " << std::setprecision(3) << xl << " / " << yl << std::endl;
	//}
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
			v += eigenPixel<double>(mipmap_up->pixelAbsolute(x2 + 1, y2));
			v += eigenPixel<double>(mipmap_up->pixelAbsolute(x2, y2 + 1));
			v += eigenPixel<double>(mipmap_up->pixelAbsolute(x2 + 1, y2 + 1));
			v/=4.0;
			p = IMGD::itkPixel(v);
		});
		mipmap_up = &(mipmap.back());
	}
}



#define COLOR

#ifdef  COLOR
using T_IMG = ImageRGBu8;
using T_IMG_D = ImageRGBd;
#else
using T_IMG = ImageGray8;
using T_IMG_D = ImageGrayd;
#endif

int main_n(int argc, char** argv)
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
		{"hexa_cm", "n1", "n2", "grey_cm_1"}};

	if (argc == 1)
	{
		for (int i = 0; i < names.size(); ++i)
			std::cout << i << ":" << names[i][0] << " / ";
		std::cout << std::endl;
		return EXIT_SUCCESS;
	}

	bool bs = false;
	bs = std::string{ argv[argc-1] } == std::string{ "save" };

	int conf = std::atoi(argv[1]);
	std::string dir = "../../data/";
	T_IMG_D cm;
	IO::loadu8_in_01(cm, dir+"colormap/"+names[conf][3] + ".png");
	auto start_chrono = std::chrono::system_clock::now();
	BlurringColorMaps<T_IMG_D> bcm{ &cm,{0.375,0.25,0.0625}, 7 };
	bcm.compute_column0();
	bcm.compute_rows();
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

	Tiling_n_Blendinng tnb((scale * std::log2(sz) * sz) / pin.width(), 753418);
	int i=0; //for the name
	while(sz>0)
	{
		start_chrono = std::chrono::system_clock::now();

		patterns(img_patterns, sz, pin, bcm, scale, tnb);
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

int main_one(int argc, char** argv)
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
		for (int i = 0; i < names.size(); ++i)
			std::cout << i << ":" << names[i][0] << " / ";
		std::cout << std::endl;
		return EXIT_SUCCESS;
	}

	bool bs = false;
	bs = std::string{ argv[argc - 1] } == std::string{ "save" };

	double scale = std::atof(argv[3]);
	int sz = std::atoi(argv[2]);

	int conf = std::atoi(argv[1]);
	std::string dir = "../../data/";
	T_IMG_D cm;
	IO::loadu8_in_01(cm, dir + "colormap/" + names[conf][3] + ".png");
	auto start_chrono = std::chrono::system_clock::now();
	BlurringColorMaps<T_IMG_D> bcm{ &cm,{0.375,0.25,0.0625}, (scale<2.0)?4:7};
	bcm.compute_column0();
	bcm.compute_rows();
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
	Tiling_n_Blendinng tnb((scale * std::log2(sz) * sz) / pin.width(), 753418);
	// Tiling_n_Blendinng tnb((scale * sz) / pin.width(),753418);
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

int main(int argc, char** argv)
{
	main_n(argc, argv);
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