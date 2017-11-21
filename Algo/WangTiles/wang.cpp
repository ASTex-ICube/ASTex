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
#include <ASTex/easy_io.h>

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

	imgA.parallel_for_region_pixels(rA, [&] (const typename IMG::PixelType& P,int x, int y, uint16_t t)
	{
		const auto& Q = imgB.pixelAbsolute(x+dx,y+dy);
		totals[t] += error_func(P,Q);
	});

	double total=0.0;
	for(double t: totals)
		total+= t;

//	std::cout << "go" << std::endl;

//	double total=0.0;
//	imgA.for_region_pixels(rA, [&] (const typename IMG::PixelType& P,int x, int y)
//	{
//		const auto& Q = imgB.pixelAbsolute(x+dx,y+dy);
//		total += error_func(P,Q);
//	});

//	std::cout << "done" << std::endl;


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

	/// rotated orginal image
	IMG rot_img_;

	/// buffer for tile creation (min-cut)
	IMG buff_img_;

	std::function<double(const IMG&, const Index&, const IMG&, const Index&)> error_func_;

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
		error_func_ = ssd_error_pixel<IMG>;
	}


	template <typename ERROR_PIX>
	void set_error_func(const ERROR_PIX& ef)
	{
		error_func_ = ef;
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

	IMG invRot45(const IMG& im, /*int xt, int yt,*/ int w)
	{
		IMG res;

		using PIX = typename IMG::PixelType;
		using T = typename IMG::DataType;
		res.initItk(w,w,true);

		int xb = im.width()/2 + 1 - w;
		int yb = im.width()/2 /*yt - 1 + w*/;
		res.for_all_pixels([&] (PIX& P, int x, int y)
		{
			int xx = xb+ x + y;
			int yy = yb + y - x;
			P = im.pixelAbsolute(xx, yy);
//			im.pixelAbsolute(xx, yy)[1] = 255;
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

		std::vector<double> errors(random_pos.size(),0.0);
		for (int j=0; j<random_pos.size(); ++j)
		{
			const Index& X = random_pos[j];
			Region rXnw = gen_region(X[0],X[1],ov_,ov_);
			Region rXse = gen_region(X[0]+tw_-ov_,X[1]+tw_-ov_,ov_,ov_);

			errors[j] += computeErrorOverlap(rot_img_, rAnw, rot_img_, rXse, [&] (const PIX& P,const PIX& Q) {return (IMG::template eigenPixel<double>(Q)-IMG::template eigenPixel<double>(P)).squaredNorm();});
			errors[j] += computeErrorOverlap(rot_img_, rAse,rot_img_, rXnw, [&] (const PIX& P,const PIX& Q) {return (IMG::template eigenPixel<double>(Q)-IMG::template eigenPixel<double>(P)).squaredNorm();});
		}

		// find the 5*NBC min an place them first
		std::vector<uint32_t> best(5);
		for(uint32_t j=0; j<5; ++j)
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
		std::cout << "best_matching" << std::endl;

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
				errors[j] += computeErrorOverlap(rot_img_, rAn[i], rot_img_, rBs, [&] (const PIX& P,const PIX& Q) {return (IMG::template eigenPixel<double>(Q)-IMG::template eigenPixel<double>(P)).squaredNorm();});
				errors[j] += computeErrorOverlap(rot_img_, rAs[i],rot_img_, rBn, [&] (const PIX& P,const PIX& Q) {return (IMG::template eigenPixel<double>(Q)-IMG::template eigenPixel<double>(P)).squaredNorm();});
				errors[j] += computeErrorOverlap(rot_img_, rAw[i],rot_img_, rBe, [&] (const PIX& P,const PIX& Q) {return (IMG::template eigenPixel<double>(Q)-IMG::template eigenPixel<double>(P)).squaredNorm();});
				errors[j] += computeErrorOverlap(rot_img_, rAe[i],rot_img_, rBw, [&] (const PIX& P,const PIX& Q) {return (IMG::template eigenPixel<double>(Q)-IMG::template eigenPixel<double>(P)).squaredNorm();}
				);
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

//		return img;
		return invRot45(img,(img.width()-ov_)/2);
	}


	double choose_tiles_pos(std::array<Index,NBC>& A, std::array<Index,NBC>& B)
	{
		A[0] = random_pos.back();
		random_pos.pop_back();
		for(uint32_t j=1; j<NBC; ++j)
			best_matching_corner(A[0],A[j]);

		random_index_patches45(2000);

		// find the best B
		double err = best_matching(A,B);
		return err;
	}


	std::vector<IMG> create_tiles()
	{
		random_index_patches45(2000);

		std::cout << "create_tiles" << std::endl;
		std::array<Index,NBC> A,B;
		double err = choose_tiles_pos(A,B);

		buff_img_.initItk(2*tw_,2*tw_,true);

		std::vector<IMG> tiles;
		tiles.reserve(NBC*NBC*NBC*NBC);

		for(Index nw : A)
			for(Index se : A)
				for(Index ne : B)
					for(Index sw : B)
						tiles.push_back(create_tile(nw,ne,se,sw));


		return tiles;

	}
};

template <uint32_t N>
uint16_t nord(uint32_t v)
{
	return v/(N*N*N);
}

template <uint32_t N>
uint16_t sud(uint32_t v)
{
	return (v/(N*N))%N;
}

template <uint32_t N>
uint16_t east(uint32_t v)
{
	return (v/N)%N;
}

template <uint32_t N>
uint16_t west(uint32_t v)
{
	return v%N;
}


template <uint32_t N>
uint16_t rand_west(uint32_t e)
{
	uint32_t r;
	do
	{
		r = (rand()%(N*N*N)) * N + east<N>(e);
	} while (r == e);
	return r;
}

template <uint32_t N>
uint16_t rand_nord(uint32_t s)
{
	uint32_t r;
	do
	{
		r = sud<N>(s)*N*N*N + rand()%(N*N*N);
	} while (r == s);
	return r;
}

template <uint32_t N>
uint16_t rand_nw(uint32_t s, uint32_t e)
{
	uint32_t r;
	do
	{
		r = sud<N>(s)*N*N*N + (rand()%(N*N))*N + east<N>(e);
	} while ((r == s)||(r==e));
	return r;
}



template <typename IMG, int NBC>
IMG compose(const std::vector<IMG> tiles, int nbw, int nbh)
{
	std::vector<uint16_t> compo_idx(nbw*nbh);
	auto f = [&compo_idx, nbw] (int i, int j) -> uint16_t& { return compo_idx[i+ j*nbw];};

	int32_t nbc = rand()%tiles.size();

	f(0,0) = rand()%(NBC*NBC*NBC*NBC);
	for(int i=1;i<nbw;++i)
		f(i,0) = rand_west<NBC>(f(i-1,0));

	for(int j=1;j<nbh;++j)
	{
		f(0,j) = rand_nord<NBC>(f(0,j-1));
		for(int i=1;i<nbw;++i)
			f(i,j) = rand_nw<NBC>(f(i,j-1), f(i-1,j));
	}

	std::cout << "create compo" << std::endl;

	int32_t w = tiles[0].width();
	int32_t h = tiles[0].height();

	IMG compo(w*nbw, h*nbh);
	for(int j=0;j<nbh;++j)
		for(int i=0;i<nbw;++i)
		{
			int32_t k = f(i,j);
			if  (i==0)
				std::cout << std::endl;
			std::cout << k << "  ";
			compo.copy_pixels(gen_index(i*w,j*h), tiles[k], gen_region(0,0,w,h));
		}
	return compo;
//	IMG compo(w*4+4, h*4+4);
//	for(int j=0;j<4;++j)
//		for(int i=0;i<4;++i)
//		{
//			int k = i +j*4;
//			compo.copy_pixels(gen_index(i*(w+1), j*(h+1)), tiles[k], gen_region(0,0,w,h));
//		}
//	return compo;

}



}

using namespace ASTex;

//std::vector<Index>* vi_ptr;
//ImageViewer* imgv_ptr;
ImageRGBu8* imr_ptr;


void app_key_pressed(int /*code*/, char /*key*/,int /*id*/)
{
}


int main(int argc, char** argv)
{
//	std::string fn = "C:/Users/thery/Desktop/blue_rust.png";
//	std::string fn = "/Users/thery/Desktop/blue_rust.png";

//	std::string fn = "/tmp/blue_rust.png";
	std::string fn = "/tmp/ASTex_data/quilting_input8.png";

	QApplication app(argc, argv);

	ImageRGBu8 im;
	im.load(fn);

	auto start_chrono = std::chrono::system_clock::now();

	WangTilesAlgo<ImageRGBu8,2> wta(im,120,20);
	wta.rotate45();
	std::vector<ImageRGBu8> tiles = wta.create_tiles();


	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "wang : " << elapsed_seconds.count() << " s." << std::endl;


	ImageRGBu8 gen = compose<ImageRGBu8,2>(tiles, 10, 8);

	ImageViewer imgv("tiled", &app, 0);
	imgv.set_rgb(gen.getDataPtr(), gen.width(),gen.height(),1);
	imgv.show();


	return app.exec();
}
