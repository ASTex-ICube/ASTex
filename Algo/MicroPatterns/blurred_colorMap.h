#ifndef _MP_BLURRED_COLOR_MAP_
#define _MP_BLURRED_COLOR_MAP_


#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>


namespace ASTex
{

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

}
#endif