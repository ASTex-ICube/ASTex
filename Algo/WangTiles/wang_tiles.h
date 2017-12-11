/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/

#ifndef _WANG_TILES_H_
#define _WANG_TILES_H_

#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>
#include <ASTex/easy_io.h>

#include "min_path_cut.h"


namespace ASTex
{


template <typename IMG>
class WangTiles
{
	std::vector<IMG> tiles_;

	int32_t nbcolors_;

	uint16_t nord(uint32_t v)
	{
		return v/(nbcolors_*nbcolors_*nbcolors_);
	}

	uint16_t sud(uint32_t v)
	{
		return (v/(nbcolors_*nbcolors_))%nbcolors_;
	}

	uint16_t east(uint32_t v)
	{
		return (v/nbcolors_)%nbcolors_;
	}

	uint16_t west(uint32_t v)
	{
		return v%nbcolors_;
	}


	uint16_t rand_west(uint32_t e)
	{
		uint32_t r;
		do
		{
			r = (rand()%(nbcolors_*nbcolors_*nbcolors_)) * nbcolors_ + east(e);
		} while (r == e);
		return r;
	}

	uint16_t rand_nord(uint32_t s)
	{
		uint32_t r;
		do
		{
			r = sud(s)*nbcolors_*nbcolors_*nbcolors_ + rand()%(nbcolors_*nbcolors_*nbcolors_);
		} while (r == s);
		return r;
	}

	uint16_t rand_nw(uint32_t s, uint32_t e)
	{
		uint32_t r;
		do
		{
			r = sud(s)*nbcolors_*nbcolors_*nbcolors_ + (rand()%(nbcolors_*nbcolors_))*nbcolors_ + east(e);
		} while ((r == s)||(r==e));
		return r;
	}

public:

	WangTiles(const std::vector<IMG>& tiles, int32_t nb):
		tiles_(tiles), nbcolors_(nb)
	{}

	WangTiles(std::vector<IMG>&& tiles, int32_t nb):
		tiles_(tiles), nbcolors_(nb)
	{}

	IMG compose(int32_t nbw, int32_t nbh)
	{
		std::vector<uint16_t> compo_idx(nbw*nbh,0);
		auto idx = [&compo_idx, nbw] (int32_t i, int32_t j) -> uint16_t& { return compo_idx[i+ j*nbw];};

		if (nbcolors_ > 1)
		{
			idx(0,0) = rand()%(nbcolors_*nbcolors_*nbcolors_*nbcolors_);
			for(int32_t i=1;i<nbw;++i)
				idx(i,0) = rand_west(idx(i-1,0));

			for(int32_t j=1;j<nbh;++j)
			{
				idx(0,j) = rand_nord(idx(0,j-1));
				for(int32_t i=1;i<nbw;++i)
					idx(i,j) = rand_nw(idx(i,j-1), idx(i-1,j));
			}
		}

		int32_t w = tiles_[0].width();
		int32_t h = tiles_[0].height();

		IMG compo(w*nbw, h*nbh);
		for(int32_t j=0;j<nbh;++j)
			for(int32_t i=0;i<nbw;++i)
			{
				int32_t k = idx(i,j);
				compo.copy_pixels(gen_index(i*w,j*h), tiles_[k], gen_region(0,0,w,h));
			}
		return compo;
	}

	IMG all_tiles()
	{
		std::vector<typename IMG::PixelType> colors;
		colors.push_back(IMG::itkPixelNorm(1,0,0));
		colors.push_back(IMG::itkPixelNorm(0,1,0));
		colors.push_back(IMG::itkPixelNorm(0,0,1));
		colors.push_back(IMG::itkPixelNorm(1,1,0));
		colors.push_back(IMG::itkPixelNorm(1,0,1));
		colors.push_back(IMG::itkPixelNorm(0,1,1));
		colors.push_back(IMG::itkPixelNorm(0.5,1,0));
		colors.push_back(IMG::itkPixelNorm(0.5,0,1));
		colors.push_back(IMG::itkPixelNorm(0,0.5,1));
		colors.push_back(IMG::itkPixelNorm(1,0.5,0));
		colors.push_back(IMG::itkPixelNorm(1,0,0.5));
		colors.push_back(IMG::itkPixelNorm(0,1,0.5));
		colors.push_back(IMG::itkPixelNorm(1,1,0.5));
		colors.push_back(IMG::itkPixelNorm(1,0.5,1));
		colors.push_back(IMG::itkPixelNorm(0.5,1,1));

		int32_t w = tiles_[0].width();
		int32_t h = tiles_[0].height();

		int32_t NB = nbcolors_*nbcolors_;

		IMG compo((w+2)*NB, (h+2)*NB);
		int32_t k =0;
		for(int32_t j=0;j<NB;++j)
			for(int32_t i=0;i<NB;++i)
			{
				int32_t c0 = j/nbcolors_;
				int32_t c1 = j%nbcolors_;
				int32_t c2 = i/nbcolors_;
				int32_t c3 = i%nbcolors_;
				compo.copy_pixels(gen_index(i*(w+2), j*(h+2)), tiles_[k], gen_region(0,0,w,h));
				for (int32_t ll=0; ll<4; ++ll)
				for (int32_t l=w/3+ll; l<(2*w)/3-ll; ++l)
				{
					compo.pixelAbsolute(i*(w+2)+l, j*(h+2)+ll) = colors[c0];
					compo.pixelAbsolute(i*(w+2)+l, j*(h+2)+h-1-ll) = colors[c1];
					compo.pixelAbsolute(i*(w+2)+ll, j*(h+2)+l) = colors[c2];
					compo.pixelAbsolute(i*(w+2)+w-1-ll, j*(h+2)+l) = colors[c3];
				}
				k++;
			}
		return compo;
	}
};





template<typename IMG>
class WangTilesGenerator
{

public:

	/**
	 * @brief WangTilesGenerator
	 * @param img input image
	 * @param tw tile width
	 * @param nbcol nb of color per edge (tile algo)
	 */
	WangTilesGenerator(const IMG& img, int32_t tw, int32_t nb) :
		input_img_(img), nbcolors_(nb), tw_(tw + tw / 5), ov_(tw / 5)
	{
		error_func_ = ssd_error_pixel<IMG>;
	}

	WangTiles<IMG> create()
	{
		rotate45();
		return WangTiles<IMG>(create_tiles(),nbcolors_);
	}

	static WangTiles<IMG> create(const IMG& img, int32_t tw, int32_t nb)
	{
		WangTilesGenerator wtg(img,tw,nb);
		wtg.rotate45();
		return WangTiles<IMG>(wtg.create_tiles(),nb);
	}

	template <typename ERROR_PIX>
	void set_error_func(const ERROR_PIX& ef)
	{
		error_func_ = ef;
	}
	ImageGrayu8 choosen_img_;
private:

	static const int32_t NB_RAND= 2000;
	static const int32_t NB_BEST= 10;

	using PIX = typename IMG::PixelType;
	using T = typename IMG::DataType;

	const IMG& input_img_;

	/// nb color edge
	int32_t nbcolors_;

	int32_t tw_;

	int32_t ov_;

	/// rotated orginal image
	IMG rot_img_;

	/// buffer for tile creation (min-cut)
	IMG buff_img_;

	std::function<double(const PIX&, const PIX&)> error_func_;

	std::vector<Index> random_pos;

	static inline int32_t random_int(int32_t min, int32_t max)
	{
		return  min + std::rand()%(max-min);
	}

	Index invRotPix(const Index& ind)
	{
		int32_t yb = input_img_.width()-1;
		int32_t x = (ind[0]-ind[1]+1+yb)/2;
		int32_t y = (ind[1]+ind[0]+1-yb)/2;

		return gen_index(x,y);
	}


	void rotate45()
	{
		using PIX = typename IMG::PixelType;
		using DPIX = typename IMG::DoublePixelEigen;
		using T = typename IMG::DataType;
		int32_t s = input_img_.width()+input_img_.height();
		rot_img_.initItk(s,s,true);

		int32_t yb = input_img_.width()-1;
		input_img_.for_all_pixels([&] (const PIX& P, int32_t x, int32_t y)
		{
			int32_t xx = x+y;
			int32_t yy = yb+y-x;
			rot_img_.pixelAbsolute(xx,yy) = P;
		});

		for (int32_t j = 1; j < input_img_.height(); ++j)
		{
			for (int32_t i = 1; i < input_img_.width(); ++i)
			{
				int32_t x = i + j - 1;
				int32_t y = yb + j - i;

				DPIX n = rot_img_.pixelEigenAbsolute(x,y-1);
				DPIX s = rot_img_.pixelEigenAbsolute(x,y+1);
				DPIX w = rot_img_.pixelEigenAbsolute(x-1,y);
				DPIX e = rot_img_.pixelEigenAbsolute(x+1,y);
				rot_img_.pixelEigenAbsolute(x, y) = (n + s + w + e) / 4;
			}
		}
	}

	IMG invRot45(const IMG& im, int32_t w)
	{
		IMG res;
		using PIX = typename IMG::PixelType;
		res.initItk(w,w,true);

		int32_t xb = im.width()/2 + 1 - w; // ? +1 check
		int32_t yb = im.width()/2;
		res.for_all_pixels([&] (PIX& P, int32_t x, int32_t y)
		{
			int32_t xx = xb+ x + y;
			int32_t yy = yb + y - x;
			P = im.pixelAbsolute(xx, yy);
		});
		return res;
	}


	void random_index_patches45(int32_t nb)
	{
		Size sz = input_img_.size();

		int32_t wm = sz[0]-tw_;
		int32_t hm = sz[1]-tw_;

		std::vector<Index>& vect = random_pos;
		vect.clear();
		vect.reserve(nb);


		for (int32_t i=0;i<nb;++i)
		{
			int32_t x = random_int(tw_-tw_/2,wm+tw_/2);
			int32_t y = random_int(0,hm);
			vect.push_back(gen_index(x+y, sz[0]+y-x));
		}
	}

	void update_choosen(Index ind)
	{
		int32_t RADIUS = 7;
		int32_t cx = ind[0] / 8;
		int32_t cy = ind[1] / 8;
		int32_t x = std::max(cx - RADIUS, 0);
		int32_t y = std::max(cy - RADIUS, 0);
		int32_t w = std::min(choosen_img_.width() - x, RADIUS) + RADIUS;
		int32_t h = std::min(choosen_img_.height() - y, RADIUS) + RADIUS;

		Region reg = gen_region(x,y,w,h);
		choosen_img_.for_region_pixels(reg, [&] (uint8_t& p, int32_t i, int32_t j)
		{
			p = std::max(p, uint8_t(std::max(RADIUS*RADIUS -((i-cx)*(i-cx)+(j-cy)*(j-cy)),0)));
		});

	}

	inline uint8_t choosen_val(Index ind)
	{
		return choosen_img_.pixelAbsolute(ind[0]/8, ind[1]/8);
	}

	double best_matching_corner(const Index& A,  Index& Ax)
	{
		Region rAnw = gen_region(A[0],A[1],ov_,ov_);
		Region rAse = gen_region(A[0]+tw_-ov_,A[1]+tw_-ov_,ov_,ov_);

		std::vector<double> errors(random_pos.size());
		for (std::size_t j=0; j<random_pos.size(); ++j)
		{
			const Index& X = random_pos[j];
			Region rXnw = gen_region(X[0],X[1],ov_,ov_);
			Region rXse = gen_region(X[0]+tw_-ov_,X[1]+tw_-ov_,ov_,ov_);
			errors[j] = computeErrorOverlap(rot_img_, rAnw, rot_img_, rXse, error_func_) + computeErrorOverlap(rot_img_, rAse,rot_img_, rXnw, error_func_);
		}

		// find the NB_BEST min an place them first
		std::vector<uint32_t> best(NB_BEST);
		for(uint32_t j=0; j<NB_BEST; ++j)
		{
			uint32_t m_i = j;
			for(uint32_t i=j+1; i<errors.size(); ++i)
				if ( errors[i] < errors[m_i] )
					m_i = i;
			best[j] = m_i;
			std::swap(errors[m_i],errors[j]);
		}

		// place further first
		int32_t imin = best[0];
		for(uint32_t j=1; j<NB_BEST; ++j)
		{
			if (choosen_val(random_pos[best[j]])<  choosen_val(random_pos[imin]))
				imin = best[j];
		}

		update_choosen(random_pos[imin]);

		// and create positions
		Ax = random_pos[imin];
		random_pos[imin] = random_pos.back();
		random_pos.pop_back();

		return errors[best[0]];
	}


	double best_matching(const std::vector<Index>& A, std::vector<Index>& B)
	{
		Size ho = gen_size(tw_,ov_);
		Size vo = gen_size(ov_,tw_);

		std::vector<Region> reg_buff(nbcolors_ * 4);
		Region* rAn = reg_buff.data();
		Region* rAs = rAn + nbcolors_;
		Region* rAw = rAs + nbcolors_;
		Region* rAe = rAw + nbcolors_;

		for (int32_t i=0; i<nbcolors_; ++i)
		{
			const Index& iA = A[i];
			rAn[i] = {iA,ho};
			rAs[i] = {gen_index(iA[0],iA[1]+tw_-ov_),ho};
			rAw[i] = {iA,vo};
			rAe[i] = {gen_index(iA[0]+tw_-ov_,iA[1]),vo};
		}

		std::vector<double> errors(random_pos.size(),0.0);
		for (std::size_t j=0; j<random_pos.size(); ++j)
		{
			const Index& iB = random_pos[j];
			Region rBn = {iB,ho};
			Region rBs = {gen_index(iB[0],iB[1]+tw_-ov_),ho};
			Region rBw = {iB,vo};
			Region rBe = {gen_index(iB[0]+tw_-ov_,iB[1]),vo};

			for (int32_t i=0; i<nbcolors_; ++i)
			{
				errors[j] += computeErrorOverlap(rot_img_, rAn[i], rot_img_, rBs, error_func_);
				errors[j] += computeErrorOverlap(rot_img_, rAs[i],rot_img_, rBn, error_func_);
				errors[j] += computeErrorOverlap(rot_img_, rAw[i],rot_img_, rBe, error_func_);
				errors[j] += computeErrorOverlap(rot_img_, rAe[i],rot_img_, rBw, error_func_);
			}
		}

		// find the NB_BEST*nbcolors_ min an place them first
		const uint32_t nb_best = NB_BEST*nbcolors_;
		std::vector<uint32_t> best(nb_best);
		for(uint32_t j=0; j<nb_best; ++j)
		{
			uint32_t m_i = j;
			for(uint32_t i=j+1; i<errors.size(); ++i)
				if ( errors[i] < errors[m_i] )
					m_i = i;
			best[j] = m_i;
			std::swap(errors[m_i],errors[j]);
		}


		for(int32_t i=0; i<nbcolors_; ++i)
		{
			int32_t jmin = i;
			for(uint32_t j=i; j<nb_best; ++j)
			{
				if (choosen_val(random_pos[best[j]]) < choosen_val(random_pos[best[jmin]]))
					jmin = j;
			}
			std::swap(best[i],best[jmin]);
			update_choosen(random_pos[best[i]]);
		}

		double err = 0.0;
		// and create positions
		for (int32_t i=0; i<nbcolors_; ++i)
		{
			B[i] = random_pos[best[i]];
			err += errors[i];
		}

		return err;
	}


	IMG create_tile(const Index& nw, const Index& ne, const Index& se, const Index& sw)
	{
		IMG img(2*tw_-ov_,2*tw_-ov_,true);
		MinCutBuffer<IMG,1> mcbh(rot_img_,rot_img_,tw_,ov_);
		int32_t ti_sz = tw_-ov_;

		img.copy_pixels(gen_index(0,0),     rot_img_, gen_region(nw[0],nw[1],         ti_sz, ti_sz));
		img.copy_pixels(gen_index(tw_,0),   rot_img_, gen_region(ne[0]+ov_,ne[1],     ti_sz, ti_sz));
		img.copy_pixels(gen_index(0,tw_),   rot_img_, gen_region(sw[0],sw[1]+ov_,     ti_sz, ti_sz));
		img.copy_pixels(gen_index(tw_,tw_), rot_img_, gen_region(se[0]+ov_,se[1]+ov_, ti_sz, ti_sz));

		mcbh.fusion(gen_index(nw[0]+ti_sz,nw[1]),ne,buff_img_,gen_index(0,0));
		mcbh.fusion(gen_index(sw[0]+ti_sz,sw[1]),se,buff_img_,gen_index(ov_,0));

		img.copy_pixels(gen_index(ti_sz,0),   buff_img_, gen_region(0,0,ov_,ti_sz));
		img.copy_pixels(gen_index(ti_sz,tw_), buff_img_, gen_region(ov_,ov_,ov_,ti_sz));


		buff_img_.copy_pixels(gen_index(0,tw_),     rot_img_,  gen_region(nw[0],nw[1]+ti_sz,ti_sz,ov_));
		buff_img_.copy_pixels(gen_index(ti_sz,tw_), buff_img_, gen_region(0,ti_sz,ov_,ov_));
		buff_img_.copy_pixels(gen_index(tw_,tw_),   rot_img_,  gen_region(ne[0]+ov_,ne[1]+ti_sz,ti_sz,ov_));

		buff_img_.copy_pixels(gen_index(0,tw_+ov_),     rot_img_,  gen_region(sw[0],sw[1],ti_sz,ov_));
		buff_img_.copy_pixels(gen_index(ti_sz,tw_+ov_), buff_img_, gen_region(ov_,0,ov_,ov_));
		buff_img_.copy_pixels(gen_index(tw_,tw_+ov_),   rot_img_,  gen_region(se[0]+ov_,se[1],ti_sz,ov_));

		MinCutBuffer<IMG,0> mcbv(buff_img_,buff_img_,2*tw_-ov_,ov_);
		mcbv.fusion(gen_index(0,tw_),gen_index(0,tw_+ov_), img, gen_index(0,ti_sz));

		return invRot45(img,(img.width()-ov_)/2);
	}


	double choose_tiles_pos(std::vector<Index>& A, std::vector<Index>& B)
	{
		A[0] = random_pos.back();
		random_pos.pop_back();
		update_choosen(A[0]);

		for(int32_t j=1; j<nbcolors_; ++j)
			best_matching_corner(A[0],A[j]);

		random_index_patches45(NB_RAND);

		// find the best B
		double err = best_matching(A,B);

		return err;
	}


	std::vector<IMG> create_tiles()
	{
		srand(time(nullptr));
		random_index_patches45(NB_RAND);

		choosen_img_.initItk(rot_img_.width()/8, rot_img_.height()/8, true);

		std::vector<Index> A(nbcolors_);
		std::vector<Index> B(nbcolors_);

		choose_tiles_pos(A,B);

		buff_img_.initItk(2*tw_,2*tw_,true);

		std::vector<IMG> tiles;
		tiles.reserve(nbcolors_*nbcolors_*nbcolors_*nbcolors_);

		for(Index nw : A)
			for(Index se : A)
				for(Index ne : B)
					for(Index sw : B)
						tiles.push_back(create_tile(nw,ne,se,sw));


		std::cout << "Tiles:"<<tiles[0].size() << " * "<< tiles.size() <<std::endl;
		return tiles;

	}



};

}

#endif
