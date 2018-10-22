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



#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>
#include "min_path_cut.h"


#include "imageviewer.h"

namespace ASTex
{

//template <typename IMG, typename EF>
//double computeErrorOverlap(const IMG& imgA, const Region& rA, const IMG& imgB, const Region& rB, const EF& error_func)
//{
//	using ITER = typename IMG::ConstIterator;

//	double total=0.0;

//	ITER it2 = imgB.beginConstIterator(rB);
//	for (ITER it = imgA.beginConstIterator(rA); !it.IsAtEnd(); ++it, ++it2)
//	{
//		total += error_func_(it.Get(), it2.Get());
//	}

//	return total;
//}

template <typename IMG, typename EF>
auto computeErrorOverlap(const IMG& imgA, const Region& rA, const IMG& imgB, const Region& rB, const EF& error_func)
-> typename std::enable_if<(function_traits<EF>::arity==2), double>::type
{
	assert_msg(rA.GetSize() == rB.GetSize(),"computeErrorOverlap: regions must have same size");

	// shift between rA & rB
	int dx = rB.GetIndex()[0] - rA.GetIndex()[0];
	int dy = rB.GetIndex()[1] - rA.GetIndex()[1];

	// one error sum for each thread
	std::vector<double> totals(nb_launched_threads(),0.0);

	//
	imgA.parallel_for_region_pixels(rA, [&] (const typename IMG::PixelType& P,int x, int y, uint16_t t)
	{
		const auto& Q = imgB.pixelAbsolute(x+dx,y+dy);
		totals[t] += error_func(P,Q);
	});

	double total=0.0;
	for(double t: totals)
		total+= t;

	return total;
}

template <typename IMG, typename EF>
auto computeErrorOverlap(const IMG& imgA, const Region& rA, const IMG& imgB, const Region& rB, const EF& error_func)
-> typename std::enable_if<(function_traits<EF>::arity==4), double>::type
{
	assert_msg(rA.GetSize() == rB.GetSize(),"computeErrorOverlap: regions must have same size");

	int dx = rB.GetIndex()[0] - rA.GetIndex()[0];
	int dy = rB.GetIndex()[1] - rA.GetIndex()[1];

	std::vector<double> totals(nb_launched_threads(),0.0);

	imgA.parallel_for_region_pixels(rA, [&] (int x, int y, uint16_t t)
	{
		totals[t] += error_func(imgA,gen_index(x,y),imgB,gen_index(x+dx,y+dy));
	});
	double total=0.0;
	for(double t: totals)
		total+= t;

	return total;
}


template<typename IMG, int NBC>
class WangTilesAlgo
{
	using PIX = typename IMG::PixelType;
	using T = typename IMG::DataType;

	const IMG& input_img_;

	IMG rot_img_;

	std::function<double(const IMG&, const Index&, const IMG&, const Index&)> error_func_;

	int tw_;

	int ov_;

	std::vector<Index> random_pos;



	static inline int random_int(int min, int max)
	{
//		double d = double(std::rand())/RAND_MAX;
//		return  min + int(d*(max-min));
		return  min + std::rand()%(max-min);
	}


	Index invRotPix(const Index& ind)
	{
		int yb = input_img_.width()-1;
		int x = (ind[0]-ind[1]+1+yb)/2;
		int y = (ind[1]+ind[0]+1-yb)/2;

		return gen_index(x,y);
	}


public:
	/**
	 * @brief WangTilesAlgo
	 * @param img input image
	 * @param tw tile width
	 * @param ov size of overlay
	 * @param nbcol nb of color per edge (tile algo)
	 */
	WangTilesAlgo(const IMG& img, int tw, int ov, int nbcol=2):
	input_img_(img),tw_(tw),ov_(ov)//,nb_colors_(nbcol)
	{
		error_func_ = ssd_error_pixels<IMG>;
	}


	template <typename ERROR_PIX>
	void set_error_func(const ERROR_PIX& ef)
	{
		error_func_ = ef;
	}


	IMG Rot45(const IMG& im)
	{
		using PIX = typename IMG::PixelType;
		using T = typename IMG::DataType;
		int s = im.width()+im.height();
		rot_img_.initItk(s,s,true);

		rot_img_.for_all_pixels([&] (PIX& P)
		{
			P=PIX(T(0));
		});


		int yb = im.width()-1;
		im.for_all_pixels([&] (const PIX& P, int x, int y)
		{
			int xx = x+y;
			int yy = yb+y-x;
			rot_img_.pixelAbsolute(xx,yy) = P;
		});

		for (int j=1; j<im.height();++j)
			for (int i=1; i<im.width();++i)
			{
				int x = i+j-1;
				int y = yb+j-i;

				RGBd n(rot_img_.pixelAbsolute(x,y-1));
				RGBd s(rot_img_.pixelAbsolute(x,y+1));
				RGBd w(rot_img_.pixelAbsolute(x-1,y));
				RGBd e(rot_img_.pixelAbsolute(x+1,y));

				RGB<T> p(0.25*n + 0.25*s +0.25*w + 0.25*e);

				rot_img_.pixelAbsolute(x,y) = p;
			}

		return rot_img_;
	}

	IMG invRot45(const IMG& im)
	{
		IMG res;

		using PIX = typename IMG::PixelType;
		using T = typename IMG::DataType;
		int s = im.width()+im.height()+1;
		res.initItk(s,s,true);

		res.for_all_pixels([&] (PIX& P)
		{
			P=PIX(T(255));
		});


		for(int j=1;j<s;j+=2)
		{
			for(int i=0;i<j+1;++i)
			{
				int y =j/2;
				int x = s/2-j/2+i;
				res.pixelAbsolute(x,y) = im.pixelAbsolute(i,j-i);
			}
		}
		return res;
	}



	void random_index_patches45(Size sz, int nb)
	{
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



//	double computeErrorOverlap(const Region& rA, const Region& rB)
//	{
//		using ITER = typename IMG::ConstIterator;
//		double total=0.0;

//		ITER it2 = rot_img_.beginConstIterator(rB);
//		for (ITER it = rot_img_.beginConstIterator(rA); !it.IsAtEnd(); ++it, ++it2)
//		{
//			total += error_func_(it.Get(), it2.Get());
//		}

//		return total;
//	}


//	std::array<double,NBC> error_overlay_(const std::array<Region,NBC>& regions, )


//	std::array<Index,NBC> best_matching_h(const std::array<Index,NBC>& patch_pos, const std::vector<Index>& random_pos)
//	{
//		Size ho = gen_size(tw_,ov_);
//		std::array<Region,NBC> overlay_A;

//		for (std::size_t i=0; i<NBC; ++i)
//			overlay_A[i] = {gen_index(patch_pos[i][0],patch_pos[i][1]+tw_-ov_),ho};

//		// compute the errors
//		std::vector<double> errors(random_pos.size(),0.0);
//		for (int j=0; j<random_pos.size(); ++i)
//			for (int i=0; i<NBC; ++i)
//				errors[j] += computeErrorOverlap(overlay_A[i], {random_pos[j], ho});

//		// find the 5*NBC min
//		const uint32_t nb_best = NBC*5;
//		std::vector<uint32_t> best(nb_best);
//		for(uint32_t j=0; j<nb_best; ++j)
//		{
//			double m_i = j;
//			for(uint32_t i=j+1; i<errors.size(); ++i)
//				if ( errors[i] < errors[m_i] )
//					m_i = i;
//			best[j] = m_i;
//			std::swap(errors[m_i],errors[j]);
//		}


//		// place a random choice first
//		for(int i=0; i<NBC; ++i)
//			std::swap(best[i],best[i+std::rand()%(best.size()-i)]);

//		// and create positions
//		std::array<Index,NBC> res;
//		for (uint32_t i=0; i<NBC; ++i)
//			res[i] = best[i];

//		return res;
//	}


//	template<bool HORIZ>
//	std::array<Index,NBC> best_matching(const std::array<Index,NBC>& patch_pos)
//	{
//		Size hvo;
//		if (HORIZ)
//			hvo = gen_size(tw_,ov_);
//		else
//			hvo = gen_size(ov_,tw_);

//		auto mk_reg = [&] (const Index& ind) -> Region
//		{
//			if (HORIZ)
//				return {gen_index(ind[0],ind[1]+tw_-ov_),hvo};
//			return {gen_index(ind[0]+tw_-ov_,ind[1]),hvo};
//		};

//		std::array<Region,NBC> overlay_A;

//		for (std::size_t i=0; i<NBC; ++i)
//			overlay_A[i] = mk_reg(patch_pos[i]);

//		// compute the errors
//		std::vector<double> errors(random_pos.size(),0.0);
//		for (int j=0; j<random_pos.size(); ++i)
//			for (int i=0; i<NBC; ++i)
//				errors[j] += computeErrorOverlap(overlay_A[i], {random_pos[j], hvo});

//		// find the 5*NBC min
//		const uint32_t nb_best = NBC*5;
//		std::vector<uint32_t> best(nb_best);
//		for(uint32_t j=0; j<nb_best; ++j)
//		{
//			double m_i = j;
//			for(uint32_t i=j+1; i<errors.size(); ++i)
//				if ( errors[i] < errors[m_i] )
//					m_i = i;
//			best[j] = m_i;
//			std::swap(errors[m_i],errors[j]);
//		}


//		// place a random choice first
//		for(int i=0; i<NBC; ++i)
//			std::swap(best[i],best[i+std::rand()%(best.size()-i)]);

//		// and create positions
//		std::array<Index,NBC> res;
//		for (uint32_t i=0; i<NBC; ++i)
//			res[i] = best[i];

//		return res;
//	}

//	std::array<Index,NBC> best_matching_hv(const std::array<Index,NBC>& patch_posH, const std::array<Index,NBC>& patch_posV)
//	{
//		Size ho = gen_size(tw_,ov_);
//		Size vo = gen_size(ov_,tw_);

//		std::array<Region,NBC> overlay_H;
//		for (uint32_t i=0; i<NBC; ++i)
//			overlay_H[i] = {gen_index(patch_posH[i][0],patch_posH[i][1]+tw_-ov_),ho};

//		std::array<Region,NBC> overlay_V;
//		for (uint32_t i=0; i<NBC; ++i)
//			overlay_V[i] = {gen_index(patch_posV[i][0]+tw_-ov_,patch_posH[i][1]),vo};


//		// compute the errors
//		std::vector<double> errors(random_pos.size(),0.0);
//		for (int j=0; j<random_pos.size(); ++i)
//			for (int i=0; i<NBC; ++i)
//				errors[j] += computeErrorOverlap(overlay_H[i], {random_pos[j], ho}) + computeErrorOverlap(overlay_V[i], {random_pos[j], vo});


//		// find the 5*NBC min
//		const uint32_t nb_best = NBC*5;
//		std::vector<uint32_t> best(nb_best);
//		for(uint32_t j=0; j<nb_best; ++j)
//		{
//			double m_i = j;
//			for(uint32_t i=j+1; i<errors.size(); ++i)
//				if ( errors[i] < errors[m_i] )
//					m_i = i;
//			best[j] = m_i;
//			std::swap(errors[m_i],errors[j]);
//		}


//		// place a random choice first
//		for(int i=0; i<NBC; ++i)
//			std::swap(best[i],best[i+std::rand()%(best.size()-i)]);

//		// and create positions
//		std::array<Index,NBC> res;
//		for (uint32_t i=0; i<NBC; ++i)
//			res[i] = best[i];

//		return res;
//	}


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


		std::vector<double> errors(random_pos.size(),0.0);
		for (int j=0; j<random_pos.size(); ++j)
		{
			const Index& iB = random_pos[j];
			Region rBn = {iB,ho};
			Region rBs = {gen_index(iB[0],iB[1]+tw_-ov_),ho};
			Region rBw = {iB,vo};
			Region rBe = {gen_index(iB[0]+tw_-ov_,iB[1]),vo};

			for (int i=0; i<NBC; ++i)
			{
				errors[j] += computeErrorOverlap(rAn[i],rBs);
				errors[j] += computeErrorOverlap(rAs[i],rBn);
				errors[j] += computeErrorOverlap(rAw[i],rBe);
				errors[j] += computeErrorOverlap(rAe[i],rBw);
			}
		}

		// find the 5*NBC min an place them first
		const uint32_t nb_best = NBC*5;
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
		std::array<Index,NBC> res;
		for (uint32_t i=0; i<NBC; ++i)
		{
			res[i] = random_pos[best[i]];
			err += errors[i];
		}

		return err;
	}


//	void init_A(std::array<Index,NBC>& A)
//	{
//		Size ov = gen_size(tw_,tw_);

//		std::vector<double> errors(random_pos.size(),0.0);
//		for (int j=0; j<random_pos.size(); ++j)
//		{
//			const Index& iB = random_pos[j];
//			errors[j] = computeErrorOverlap({A,ov}, {iB,ov });
//		}

//		// find the NBC min an place them first
//		std::vector<uint32_t> best(NBC);
//		for(uint32_t j=0; j<NBC; ++j)
//		{
//			uint32_t m_i = j;
//			for(uint32_t i=j+1; i<errors.size(); ++i)
//				if ( errors[i] > errors[m_i] )
//					m_i = i;
//			best[j] = m_i;
//			std::swap(errors[m_i],errors[j]);
//		}


//		// place a random choice first
//		for(int i=0; i<NBC; ++i)
//			std::swap(best[i],best[i+std::rand()%(best.size()-i)]);

//		double err = 0.0;
//		// and create positions
//		std::array<Index,NBC> res;
//		for (uint32_t i=0; i<NBC; ++i)
//		{
//			res[i] = random_pos[best[i]];
//			err += errors[i];
//		}

//		return err;
//	}


	void create_tiles()
	{
		std::array<Index,NBC> A;

		// choose the A, randomy
		for(uint32_t j=0; j<NBC; ++j)
		{
			A[j] = random_pos.back();
			random_pos.pop_back();
		}

		// find the best B
		std::array<Index,NBC> B;

		double err = best_matching(A,B);

		// try several time and keep those which has minimal error ?

//		copy_pixels(const Index& dst, const Self& img_from, const Region& reg_from)

		std::vector<IMG> tiles;
		tiles.reserve(NBC*NBC*NBC*NBC);
		for(Index nw : A)
			for(Index se : A)
				for(Index ne : B)
					for(Index sw : B)
						tiles.push_back(create_tile(nw,ne,se,sw));


	}

	IMG create_tile(const Index& nw, const Index& ne, const Index& se, const Index& sw)
	{
	MinCutBuffer<IMG,0> mcb(rot_img_,rot_img_,tw_,ov_);
	}


//	double computeErrorFourPatchesOverlaps(const Index* pa)
//	{
//		std::array<Region,16> overlays;
//		Size ho = gen_size(tw_,ov_);
//		Size vo = gen_size(ov_,tw_);

//		for (int i=0; i<4; i++)
//		{
//			overlays[4*i]   = {pa[i],ho};									// N
//			overlays[4*i+1] = {gen_index(pa[i][0],pa[i][1]+tw_-ov_),ho};	// S
//			overlays[4*i+2] = {pa[i],vo};									// E
//			overlays[4*i+3] = {gen_index(pa[i][0]+tw_-ov_,pa[i][1]),ho};	// W
//		}
//		double nb = tw_*ov_;
//		double err4 = 0.0;
//		for (int i=0; i<2; i++)
//		{
//			for (int j=0; j<2; j++)
//			{
//				err4 +=c / nb;  // E-W
//				err4 += computeErrorOverlap( overlays[4*i+3], overlays[4*j+2]) / nb;  // W-E
//			}
//			for (int j=2; j<4; j++)
//			{
//				err4 += computeErrorOverlap( overlays[4*i], overlays[4*j+1]) / nb;  // N-S
//				err4 += computeErrorOverlap( overlays[4*i+1], overlays[4*j]) / nb;  // S-N
//			}
//		}
//		for (int i=2; i<4; i++)
//		{
//			for (int j=0; j<2; j++)
//			{
//				err4 += computeErrorOverlap( overlays[4*i],   overlays[4*j+1]) / nb;  // N-S
//				err4 += computeErrorOverlap( overlays[4*i+1], overlays[4*j]) / nb;  // S-N
//			}
//			for (int j=2; j<4; j++)
//			{
//				err4 += computeErrorOverlap( overlays[4*i+2], overlays[4*j+3]) / nb;  // E-W
//				err4 += computeErrorOverlap( overlays[4*i+3], overlays[4*j+2]) / nb;  // W-E
//			}
//		}

//		return err4;
//	}


//	double computeErrorPatchesOverlaps(/*const Index* pa*/ const std::array<Index,4>& pa)
//	{
//		int nbp = nb_colors_ * 2;

//		std::array<Region,16> overlays;
//		Size ho = gen_size(tw_,ov_);
//		Size vo = gen_size(ov_,tw_);

//		for (int i=0; i<nbp; i++)
//		{
//			overlays[4*i]   = {pa[i],ho};									// N
//			overlays[4*i+1] = {gen_index(pa[i][0],pa[i][1]+tw_-ov_),ho};	// S
//			overlays[4*i+2] = {pa[i],vo};									// E
//			overlays[4*i+3] = {gen_index(pa[i][0]+tw_-ov_,pa[i][1]),vo};	// W
//		}
//		double nb = tw_*ov_;
//		double err4 = 0.0;
//		for (int i=0; i<nb_colors_; i++)
//		{
//			for (int j=0; j<nb_colors_; j++)
//			{
//				err4 += computeErrorOverlap( overlays[4*i+2], overlays[4*j+3]) / nb;  // E-W
//				err4 += computeErrorOverlap( overlays[4*i+3], overlays[4*j+2]) / nb;  // W-E
//			}
//			for (int j=nb_colors_; j<nbp; j++)
//			{
//				err4 += computeErrorOverlap( overlays[4*i], overlays[4*j+1]) / nb;  // N-S
//				err4 += computeErrorOverlap( overlays[4*i+1], overlays[4*j]) / nb;  // S-N
//			}
//		}
//		for (int i=nb_colors_; i<nbp; i++)
//		{
//			for (int j=0; j<nb_colors_; j++)
//			{
//				err4 += computeErrorOverlap( overlays[4*i],   overlays[4*j+1]) / nb;  // N-S
//				err4 += computeErrorOverlap( overlays[4*i+1], overlays[4*j]) / nb;  // S-N
//			}
//			for (int j=nb_colors_; j<nbp; j++)
//			{
//				err4 += computeErrorOverlap( overlays[4*i+2], overlays[4*j+3]) / nb;  // E-W
//				err4 += computeErrorOverlap( overlays[4*i+3], overlays[4*j+2]) / nb;  // W-E
//			}
//		}

//		return err4;
//	}




//	IMG createTile(std::array<Index,4> pa)
//	{
//		ImagePixelIndex imi;
////		IMG img;
//		imi.initItk(2*tw_,2*tw_,false);
//		imi.for_all_pixels([] (uint32_t& i) { i = 0;});

//		using ITER = typename IMG::Iterator;
//		using CITER = typename IMG::ConstIterator;


//		Region src = gen_region(pa[0][0],pa[0][1],tw_-ov_,tw_-ov_);
//		Region dst = gen_region(0,0,tw_-ov_,tw_-ov_);

//		for_indices(src,dst,[&] (int is, int js, int id, int jd)
//		{
//			imi.iset(id,jd,invRotPix(gen_index(is,js)));
//		});


//		src = gen_region(pa[1][0]+ov_,pa[1][1],tw_-ov_,tw_-ov_);
//		dst = gen_region(tw_,0,tw_-ov_,tw_-ov_);

//		for_indices(src,dst,[&] (int is, int js, int id, int jd)
//		{
//			imi.iset(id,jd,invRotPix(gen_index(is,js)));
//		});

//		src = gen_region(pa[2][0],pa[2][1]+ov_,tw_-ov_,tw_-ov_);
//		dst = gen_region(0,tw_,tw_-ov_,tw_-ov_);

//		for_indices(src,dst,[&] (int is, int js, int id, int jd)
//		{
//			imi.iset(id,jd,invRotPix(gen_index(is,js)));
//		});

//		src = gen_region(pa[3][0]+ov_,pa[3][1]+ov_,tw_-ov_,tw_-ov_);
//		dst = gen_region(tw_,tw_,tw_-ov_,tw_-ov_);

//		for_indices(src,dst,[&] (int is, int js, int id, int jd)
//		{
//			imi.iset(id,jd,invRotPix(gen_index(is,js)));
//		});

//		IMG img;
//		img.initItk(2*tw_,2*tw_,true);
//		apply_all_pixels(imi,img, [&] ( uint32_t i) -> PIX
//		{
//			Index ind = ImagePixelIndex::idx_to_pos(i);
//			return input_img_.pixelAbsolute(ind);
//		});

//		return img;
//	}


};


}
using namespace ASTex;

std::vector<Index>* vi_ptr;
ImageViewer* imgv_ptr;
ImageRGBu8* imr_ptr;


void app_key_pressed(int /*code*/, char /*key*/,int /*id*/)
{
////	int j=0;
//	for (auto i: *vi_ptr)
//	{
//		for_indices(0,50,0,50,[&](int x, int y)
//		{
//			imr_ptr->pixelAbsolute(i[0]+x,i[1]+y) = RGBu8(0,0,255);
//		});
//	}

//	for_indices(0,50,0,50,[&](int x, int y)
//	{
//		imr_ptr->pixelAbsolute(x,y) = RGBu8(0,0,255);
//	});


//	imgv_ptr->set_rgb(imr_ptr->getDataPtr(),imr_ptr->width(),imr_ptr->height(),2);
//	imgv_ptr->update();

}




bool test_overlap()
{
	ImageRGBu8 im;
	im.load("/tmp/ASTex_data/quilting_input8.png");
	im.copy_pixels(gen_index(100,10), im, gen_region(10,10,20,200));
	im.copy_pixels(gen_index(130,10), im, gen_region(10,10,20,200));
	im.copy_pixels(gen_index(160,10), im, gen_region(10,10,20,200));
	im.copy_pixels(gen_index(190,10), im, gen_region(10,10,20,200));

	im.for_region_pixels(gen_region(160,10,20,200),[] (ImageRGBu8::PixelType &P)
	{
		P[0] = 0;
	});
	im.for_region_pixels(gen_region(190,10,20,200),[] (ImageRGBu8::PixelType &P)
	{
		P[0] = 0;
		P[1] = 0;
		P[2] = 0;
	});

	ImageRGBu8::PixelType P = im.pixelAbsolute(195,15);
	std::cout << "PIX:" << P << std::endl;

	im.save("/tmp/mel3.png");


	double err = computeErrorOverlap(im, gen_region(10,10,20,200), im, gen_region(10,10,20,200),
									 [] (const ImageRGBu8::PixelType& P, const ImageRGBu8::PixelType& Q)
	{
		double err0 = ImageRGBu8::normalized_value(P[0]) - ImageRGBu8::normalized_value(Q[0]);
		double err1 = ImageRGBu8::normalized_value(P[1]) - ImageRGBu8::normalized_value(Q[1]);
		double err2 = ImageRGBu8::normalized_value(P[2]) - ImageRGBu8::normalized_value(Q[2]);
		return err0*err0 + err1*err1 + err2*err2;
	 });

	std::cout << err<< std::endl;
	err = computeErrorOverlap(im, gen_region(10,10,20,200), im, gen_region(130,10,20,200), ssd_error_pixels<ImageRGBu8>);
	std::cout << err<< std::endl;
	err = computeErrorOverlap(im, gen_region(10,10,20,200), im, gen_region(160,10,20,200), ssd_error_pixels<ImageRGBu8>);
	std::cout << err<< std::endl;
	err = computeErrorOverlap(im, gen_region(10,10,20,200), im, gen_region(190,10,20,200), ssd_error_pixels<ImageRGBu8>);
	std::cout << err<< std::endl;


	return true;
}










int main(int argc, char** argv)
{
	test_overlap();
	return 0;


	QApplication app(argc, argv);

	ImageRGBu8 im;
	im.load("/tmp/ASTex_data/quilting_input8.png");


	auto start_chrono = std::chrono::system_clock::now();

	WangTilesAlgo<ImageRGBu8,2> wta(im,60,10);
	ImageRGBu8 imr = wta.Rot45(im);

	int nb_samples = 8*1024;

	Size lsz = gen_size(im.width(),im.height());
	wta.random_index_patches45(lsz, nb_samples);

	std::vector<int> errs(nb_samples);

//#pragma omp parallel for
//	for(int j=0;j<nb_samples-4;++j)
//	{
//		errs[j] =  wta.computeErrorPatchesOverlaps(vi.data()+4*j);
//	}




	int min_pos=0;
	for(int j=1;j<nb_samples/4;++j)
	{
		if (errs[j] < errs[min_pos])
			min_pos = j;
	}


//	std::array<Index,4> ind;
//	for (int i=0;i<4;++i)
//		ind[i] = vi[4*min_pos+i];
//	ImageRGBu8 imt = wta.createTile(ind);

//	ImageRGBu8 imt2 = wta.invRot45(imt);

	ImageRGBu8 imt;
	ImageRGBu8 imt2;


	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "wang : " << elapsed_seconds.count() << " s." << std::endl;




	ImageViewer imgv("Rot45", &app, 0);
	imgv.set_rgb(imr.getDataPtr(), imr.width(),imr.height(),2);
	imgv.show();

	ImageViewer imgvt("Tile", &app, 1);
	imgvt.set_rgb(imt.getDataPtr(), imt.width(),imt.height(),2);
	imgvt.show();


//	ImageViewer imgvt2("Tile", &app, 1);
//	imgvt2.set_rgb(imt2.getDataPtr(), imt2.width(),imt2.height(),1);
//	imgvt2.show();


	imgv_ptr = &imgv;
	imr_ptr = & imr;
//	vi_ptr = &vi;

	//imr.save("melon_r45.png");

	return app.exec();
}
