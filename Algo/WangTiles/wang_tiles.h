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


template <typename IMG, int NBC>
class WangTiles
{
	std::vector<IMG> tiles_;

	uint16_t nord(uint32_t v)
	{
		return v/(NBC*NBC*NBC);
	}

	uint16_t sud(uint32_t v)
	{
		return (v/(NBC*NBC))%NBC;
	}

	uint16_t east(uint32_t v)
	{
		return (v/NBC)%NBC;
	}

	uint16_t west(uint32_t v)
	{
		return v%NBC;
	}


	uint16_t rand_west(uint32_t e)
	{
		uint32_t r;
		do
		{
			r = (rand()%(NBC*NBC*NBC)) * NBC + east(e);
		} while (r == e);
		return r;
	}

	uint16_t rand_nord(uint32_t s)
	{
		uint32_t r;
		do
		{
			r = sud(s)*NBC*NBC*NBC + rand()%(NBC*NBC*NBC);
		} while (r == s);
		return r;
	}

	uint16_t rand_nw(uint32_t s, uint32_t e)
	{
		uint32_t r;
		do
		{
			r = sud(s)*NBC*NBC*NBC + (rand()%(NBC*NBC))*NBC + east(e);
		} while ((r == s)||(r==e));
		return r;
	}

public:

	WangTiles(const std::vector<IMG>&& tiles):
		tiles_(tiles)
	{}

	WangTiles(std::vector<IMG>&& tiles):
		tiles_(tiles)
	{}

	IMG compose(int nbw, int nbh)
	{
		std::vector<uint16_t> compo_idx(nbw*nbh,0);
		auto idx = [&compo_idx, nbw] (int i, int j) -> uint16_t& { return compo_idx[i+ j*nbw];};

		if (NBC > 1)
		{
			idx(0,0) = rand()%(NBC*NBC*NBC*NBC);
			for(int i=1;i<nbw;++i)
				idx(i,0) = rand_west(idx(i-1,0));

			for(int j=1;j<nbh;++j)
			{
				idx(0,j) = rand_nord(idx(0,j-1));
				for(int i=1;i<nbw;++i)
					idx(i,j) = rand_nw(idx(i,j-1), idx(i-1,j));
			}
		}

		int32_t w = tiles_[0].width();
		int32_t h = tiles_[0].height();

		IMG compo(w*nbw, h*nbh);
		for(int j=0;j<nbh;++j)
			for(int i=0;i<nbw;++i)
			{
				int32_t k = idx(i,j);
//				if  (i==0)
//					std::cout << std::endl;
//				std::cout << k << "  ";
				compo.copy_pixels(gen_index(i*w,j*h), tiles_[k], gen_region(0,0,w,h));
			}
		return compo;
	}

	IMG all_tiles()
	{
		int32_t w = tiles_[0].width();
		int32_t h = tiles_[0].height();

		const int NB = NBC*NBC;

		IMG compo((w+2)*NB, (h+2)*NB);
		int k =0;
		for(int j=0;j<NB;++j)
			for(int i=0;i<NB;++i)
				compo.copy_pixels(gen_index(i*(w+2), j*(h+2)), tiles_[k++], gen_region(0,0,w,h));
		return compo;
	}
};





template<typename IMG, int NBC>
class WangTilesGenerator
{

public:

	/**
	 * @brief WangTilesGenerator
	 * @param img input image
	 * @param tw tile width
	 * @param nbcol nb of color per edge (tile algo)
	 */
	WangTilesGenerator(const IMG& img, int tw):
	input_img_(img),tw_(tw+tw/5),ov_(tw/5)
	{
		error_func_ = ssd_error_pixel<IMG>;
	}

	WangTiles<IMG,NBC> create()
	{
		rotate45();
		return WangTiles<IMG,NBC>(create_tiles());
	}

	static WangTiles<IMG,NBC> create(const IMG& img, int tw)
	{
		WangTilesGenerator wtg(img,tw);
		wtg.rotate45();
		return WangTiles<IMG,NBC>(wtg.create_tiles());
	}

	template <typename ERROR_PIX>
	void set_error_func(const ERROR_PIX& ef)
	{
		error_func_ = ef;
	}

private:

	static const int NB_RAND= 2000;
	static const int NB_BEST= 5;

	using PIX = typename IMG::PixelType;
	using T = typename IMG::DataType;

	const IMG& input_img_;

	/// rotated orginal image
	IMG rot_img_;

	/// buffer for tile creation (min-cut)
	IMG buff_img_;

	std::function<double(const PIX&, const PIX&)> error_func_;

	int tw_;

	int ov_;

	std::vector<Index> random_pos;

	static inline int random_int(int min, int max)
	{
		return  min + std::rand()%(max-min);
	}

	Index invRotPix(const Index& ind)
	{
		int yb = input_img_.width()-1;
		int x = (ind[0]-ind[1]+1+yb)/2;
		int y = (ind[1]+ind[0]+1-yb)/2;

		return gen_index(x,y);
	}


	void rotate45()
	{
		using PIX = typename IMG::PixelType;
		using DPIX = typename IMG::DoublePixelEigen;
		using T = typename IMG::DataType;
		int s = input_img_.width()+input_img_.height();
		rot_img_.initItk(s,s,true);

		rot_img_.for_all_pixels([&] (PIX& P)
		{
			P=PIX(T(0));
		});


		int yb = input_img_.width()-1;
		input_img_.for_all_pixels([&] (const PIX& P, int x, int y)
		{
			int xx = x+y;
			int yy = yb+y-x;
			rot_img_.pixelAbsolute(xx,yy) = P;
		});

		for (int j = 1; j < input_img_.height(); ++j)
		{
			for (int i = 1; i < input_img_.width(); ++i)
			{
				int x = i + j - 1;
				int y = yb + j - i;

				DPIX n = rot_img_.pixelEigenAbsolute(x,y-1);
				DPIX s = rot_img_.pixelEigenAbsolute(x,y+1);
				DPIX w = rot_img_.pixelEigenAbsolute(x-1,y);
				DPIX e = rot_img_.pixelEigenAbsolute(x+1,y);
				rot_img_.pixelEigenAbsolute(x, y) = (n + s + w + e) / 4;
			}
		}
	}

	IMG invRot45(const IMG& im, int w)
	{
		IMG res;
		using PIX = typename IMG::PixelType;
		res.initItk(w,w,true);

		int xb = im.width()/2 + 1 - w; // ? +1 check
		int yb = im.width()/2;
		res.for_all_pixels([&] (PIX& P, int x, int y)
		{
			int xx = xb+ x + y;
			int yy = yb + y - x;
			P = im.pixelAbsolute(xx, yy);
		});
		return res;
	}


	void random_index_patches45(int nb)
	{
		Size sz = input_img_.size();

		int wm = sz[0]-tw_;
		int hm = sz[1]-tw_;

		std::vector<Index>& vect = random_pos;
		vect.clear();
		vect.reserve(nb);


		for (int i=0;i<nb;++i)
		{
			int x = random_int(tw_-tw_/2,wm+tw_/2);
			int y = random_int(0,hm);
			vect.push_back(gen_index(x+y, sz[0]+y-x));
		}
	}


	double best_matching_corner(const Index& A,  Index& Ax)
	{
		Region rAnw = gen_region(A[0],A[1],ov_,ov_);
		Region rAse = gen_region(A[0]+tw_-ov_,A[1]+tw_-ov_,ov_,ov_);

//		auto pix_err = [] (const PIX& P,const PIX& Q) {return (IMG::template normalized<double>(Q)-IMG::template normalized<double>(P)).squaredNorm();};

		std::vector<double> errors(random_pos.size(),0.0);
		for (std::size_t j=0; j<random_pos.size(); ++j)
		{
			const Index& X = random_pos[j];
			Region rXnw = gen_region(X[0],X[1],ov_,ov_);
			Region rXse = gen_region(X[0]+tw_-ov_,X[1]+tw_-ov_,ov_,ov_);

						errors[j] += computeErrorOverlap(rot_img_, rAnw, rot_img_, rXse, error_func_);
						errors[j] += computeErrorOverlap(rot_img_, rAse,rot_img_, rXnw, error_func_);
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

		// place a random choice first
			std::swap(best[0],best[std::rand()%best.size()]);

		// and create positions
		Ax = random_pos[best[0]];
		random_pos[best[0]] = random_pos.back();
		random_pos.pop_back();

		return errors[best[0]];
	}


	double best_matching(const std::array<Index,NBC>& A, std::array<Index,NBC>& B)
	{
		Size ho = gen_size(tw_,ov_);
		Size vo = gen_size(ov_,tw_);

		std::array<Region,NBC> rAn;
		std::array<Region,NBC> rAs;
		std::array<Region,NBC> rAw;
		std::array<Region,NBC> rAe;
		for (int i=0; i<NBC; ++i)
		{
			const Index& iA = A[i];
			rAn[i] = {iA,ho};
			rAs[i] = {gen_index(iA[0],iA[1]+tw_-ov_),ho};
			rAw[i] = {iA,vo};
			rAe[i] = {gen_index(iA[0]+tw_-ov_,iA[1]),vo};
		}

//		auto pix_err = [] (const PIX& P,const PIX& Q) {return (IMG::template normalized<double>(Q)-IMG::template normalized<double>(P)).squaredNorm();};

		std::vector<double> errors(random_pos.size(),0.0);
		for (std::size_t j=0; j<random_pos.size(); ++j)
		{
			const Index& iB = random_pos[j];
			Region rBn = {iB,ho};
			Region rBs = {gen_index(iB[0],iB[1]+tw_-ov_),ho};
			Region rBw = {iB,vo};
			Region rBe = {gen_index(iB[0]+tw_-ov_,iB[1]),vo};

			for (int i=0; i<NBC; ++i)
			{
								errors[j] += computeErrorOverlap(rot_img_, rAn[i], rot_img_, rBs, error_func_);
								errors[j] += computeErrorOverlap(rot_img_, rAs[i],rot_img_, rBn, error_func_);
								errors[j] += computeErrorOverlap(rot_img_, rAw[i],rot_img_, rBe, error_func_);
								errors[j] += computeErrorOverlap(rot_img_, rAe[i],rot_img_, rBw, error_func_);
			}
		}

		// find the NB_BEST*NBC min an place them first
		const uint32_t nb_best = NB_BEST*NBC;
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


		// place a random choice first
		for(int i=0; i<NBC; ++i)
			std::swap(best[i],best[i+std::rand()%(best.size()-i)]);

		double err = 0.0;
		// and create positions
		for (uint32_t i=0; i<NBC; ++i)
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
		int ti_sz = tw_-ov_;

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


	double choose_tiles_pos(std::array<Index,NBC>& A, std::array<Index,NBC>& B)
	{
		A[0] = random_pos.back();
		random_pos.pop_back();
		for(uint32_t j=1; j<NBC; ++j)
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

		std::array<Index,NBC> A,B;
		choose_tiles_pos(A,B);

		buff_img_.initItk(2*tw_,2*tw_,true);

		std::vector<IMG> tiles;
		tiles.reserve(NBC*NBC*NBC*NBC);

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
