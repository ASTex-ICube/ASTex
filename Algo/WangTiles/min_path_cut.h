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


namespace ASTex
{



//template <typename IMG>
//inline double ssd_error_pixel(const IMG& imA, const Index& iA, const IMG& imB, const Index& iB)
//{
//	typename IMG::DoublePixelEigen pA = imA.pixelNormEigenAbsolute(iA);
//	typename IMG::DoublePixelEigen pB = imB.pixelNormEigenAbsolute(iB);
//	return (pB-pA).squaredNorm();
//}

template <typename IMG>
inline double ssd_error_pixel(const typename IMG::PixelType& A, const typename IMG::PixelType& B)
{
	auto pA = IMG::template normalized<double>(A);
	auto pB = IMG::template normalized<double>(B);
	return (pB-pA).squaredNorm();
}



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

//	imgA.parallel_for_region_pixels(rA, [&imgB, &totals, dx, dy, &error_func] (const typename IMG::PixelType& P,int x, int y, uint16_t t)
//	{
//		const auto& Q = imgB.pixelAbsolute(x+dx,y+dy);
//		totals[t] += error_func(P,Q);
//	});

//	double total=0.0;
//	for(double t: totals)
//		total+= t;

	double total=0.0;
	imgA.for_region_pixels(rA, [&] (const typename IMG::PixelType& P,int x, int y)
	{
		const auto& Q = imgB.pixelAbsolute(x+dx,y+dy);
		total += error_func(P,Q);
	});

	return total/(rA.GetSize()[0]*rA.GetSize()[1]);
}

template <typename IMG, typename EF>
auto computeErrorOverlap(const IMG& imgA, const Region& rA, const IMG& imgB, const Region& rB, const EF& error_func)
-> typename std::enable_if<(function_traits<EF>::arity==4), double>::type
{
	assert_msg(rA.GetSize() == rB.GetSize(),"computeErrorOverlap: regions must have same size");

	int dx = rB.GetIndex()[0] - rA.GetIndex()[0];
	int dy = rB.GetIndex()[1] - rA.GetIndex()[1];

//	std::vector<double> totals(nb_launched_threads(),0.0);

//	imgA.parallel_for_region_pixels(rA, [&] (int x, int y, uint16_t t)
//	{
//		totals[t] += error_func(imgA,gen_index(x,y),imgB,gen_index(x+dx,y+dy));
//	});
//	double total=0.0;
//	for(double t: totals)
//		total+= t;

	double total=0.0;
	imgA.for_region_pixels(rA, [&] (int x, int y)
	{
		total += error_func(imgA,gen_index(x,y),imgB,gen_index(x+dx,y+dy));
	});

	return total/(rA.GetSize()[0]*rA.GetSize()[1]);
}



/**
 *
 *  DIR: 0=Horizontal  -
 *       1=Vertical    |
 */
template<typename IMG, int DIR>
class MinCutBuffer
{
	using PIX = typename IMG::PixelType;
	using T = typename IMG::DataType;

	const IMG& imA_;
	const IMG& imB_;

	int length_;
	int overlay_;

	std::vector<double> data_err_;
	std::vector<double> data_err_cum_;

	std::vector<int> minPos_;

	std::function<double(const PIX&, const PIX&)> error_func_;

public:
	/**
	 * @brief constructor
	 * @param imA input image 1
	 * @param imB input image 2
	 * @param tw tile width
	 * @param to tile overlay
	 */
	MinCutBuffer(const IMG& imA, const IMG& imB, int tw, int to):
		imA_(imA), imB_(imB),
		length_(tw),overlay_(to),
		data_err_(tw*to),
		data_err_cum_(tw*to),
		minPos_(tw)
	{
		error_func_ = ssd_error_pixel<IMG>;
	}

	~MinCutBuffer()
	{
	}

	template <typename ERROR_PIX>
	void set_error_func(const ERROR_PIX& ef)
	{
		error_func_ = ef;
	}



protected:

	inline double& error_local(int i, int j)
	{
		return data_err_[i+j*overlay_];
	}

	inline double&  error_cumul(int i, int j)
	{
		return data_err_cum_[i+j*overlay_];
	}

	/**
	 * @brief get minimum of a column or row
	 * @param c indice of column or row
	 * @return
	 */
	int minOf(int c)
	{
		double m = std::numeric_limits<double>::max();
		int mi=0;
		for (int i=0;i<overlay_;++i)
		{
			double e = error_cumul(i,c);
			if (e < m)
			{
				m = e;
				mi = i;
			}
		}
		return mi;
	}

	/**
	 * @brief get minimum Of 3 value (i,j) / (i-1,j) / (i+1,j)
	 * @param i col (row)
	 * @param j row (col)
	 * @return
	 */
	int minOf3(int i, int j)
	{
		double m = error_cumul(i,j);
		int mi = i;
		if (i>0)
		{
			double v = error_cumul(i-1,j);
			if (v<m)
			{
				m = v;
				mi = i-1;
			}
		}
		i++;
		if (i<overlay_)
		{
			double v = error_cumul(i,j);
			if (v<m)
			{
				mi = i;
			}
		}
		return mi;
	}


	/**
	 * @brief compute the path of cut
	 * @param posA position in A
	 * @param posB position in B
	 */
	void pathCut(const Index& posA, const Index& posB)
	{
//		using PIX = typename IMG::PixelType;
		// compute initial error
		auto dir_index = [] (int i,int j) {if (DIR==1) return gen_index(i,j); else return gen_index(j,i);};

		// compute initial error

		for (int j=0; j<length_; ++j)
			for (int i=0; i<overlay_; ++i)
			{
				Index inc = dir_index(i,j);
				error_local(i,j) = error_func_(imA_.pixelAbsolute(posA+inc), imB_.pixelAbsolute(posB+inc));
			}

		// compute cumulative error

		// first row
		for(int i=0;i<overlay_;++i)
			error_cumul(i,0) = error_local(i,0);

		// others rows
		for(int j=1;j<length_;++j)
		{
			error_cumul(0,j) = error_local(0,j) + std::min(error_cumul(0,j-1),error_cumul(1,j-1));
			int i=1;
			while (i<overlay_-1)
			{
				error_cumul(i,j) = error_local(i,j) + std::min( std::min(error_cumul(i-1,j-1),error_cumul(i,j-1)), error_cumul(i+1,j-1)) ;
				++i;
			}
			error_cumul(i,j) = error_local(i,j) + std::min(error_cumul(i-1,j-1),error_cumul(i,j-1));
		}

		// up to store local position of min error cut
		int i=length_-1;

		minPos_[i] = minOf(i);
		while (i>0)
		{
			minPos_[i-1] = minOf3(minPos_[i],i-1);
			--i;
		}
	}

public:
	void fusion(const Index& posA, const Index& posB, IMG& dst, const Index& pos)
	{
		pathCut(posA, posB);

		auto dir_index = [] (int i,int j) {if (DIR==1) return gen_index(i,j); return gen_index(j,i);};

		for(int j=0; j<length_;++j)
		{
			int i=0;
			for(; i< minPos_[j];++i)
				dst.pixelAbsolute(pos+dir_index(i,j)) = imA_.pixelAbsolute(posA+dir_index(i,j));

			typename IMG::DoublePixelEigen p = imA_.pixelEigenAbsolute(posA+dir_index(i,j));
			typename IMG::DoublePixelEigen q = imB_.pixelEigenAbsolute(posB+dir_index(i,j));
			dst.pixelEigenAbsolute(pos+dir_index(i,j)) = (p+q)/2;

			i++;
			for(; i<overlay_;++i)
				dst.pixelAbsolute(pos+dir_index(i,j)) = imB_.pixelAbsolute(posB+dir_index(i,j));
		}
	}

};




} // namespace ASTex
