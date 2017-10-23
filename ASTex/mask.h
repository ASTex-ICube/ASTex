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



#ifndef __ASTEX_MASKS__
#define __ASTEX_MASKS__

#include <ASTex/image_gray.h>
#include <iostream>

namespace ASTex
{

/**
 * @brief The Mask class
 */
class ASTEX_API Mask
{
public :
    virtual bool operator()(int i, int j) const=0;

	virtual void set_value(int i, int j, bool b) = 0;

	/**
	 * @brief export binary image
	 * @param im output image (must be initialized/allocated)
	 */
    void export_binary_image (ImageGrayu8& im);
};


/**
 * @brief Mask class that always return true
 */
class ASTEX_API MaskTrue : public Mask
{
public :
	inline bool operator()(int,int) const { return true; }
	inline void set_value(int, int, bool) {	std::cerr << "not possible to set a value in MaskTrue" << std::endl;}
};


template <typename T>
class MaskImage : public Mask
{
public :
	using MaskPixelType = T;
	using MaskImageType = ImageGray<T>;

protected:
	inline MaskImage(typename MaskImageType::ItkImg::Pointer im)
	{
		img_mask_.update_itk(im);
	}

public :

	MaskImage()
	{}

	inline MaskImage(int width, int height)
	{
		img_mask_.initItk(width,height);
	}

	inline const MaskImageType& get_image() const { return img_mask_; }

//	inline MaskImageType& get_image() { return img_mask_; }

	inline int width() const {return img_mask_.width();}

	inline int height() const {return img_mask_.height();}


protected :
	MaskImageType img_mask_;
};




class ASTEX_API MaskBool : public MaskImage<uint8_t>
{
	inline MaskBool(const MaskBool& /*mb*/):MaskImage<uint8_t>() {}
	inline MaskBool& operator=(const MaskBool& /*mb*/) {return *this;}

protected:
	inline MaskBool(MaskImageType::ItkImg::Pointer im):MaskImage(im) {}

public :
	using MaskPixelType = uint8_t;

	using MaskImageType = ImageGray<uint8_t>;

	MaskBool(int width, int height);

	MaskBool(int width, int height, bool val);

	void clone(const MaskBool& mb);

	inline void swap(MaskBool& mb)
	{
		img_mask_.swap(mb.img_mask_);
	}

	inline bool operator()(int i,int j) const
	{
		return img_mask_.pixelAbsolute(i,j) != 0;
	}

	inline void set_value(int i, int j, bool b)
	{
		img_mask_.pixelAbsolute(i,j) = 255*MaskPixelType(b);
	}

	void clear(bool val);

	template <typename MASK>
	MaskBool& operator&=(const MASK& m);

	template <typename MASK>
	MaskBool& operator|=(const MASK& m);

	MaskBool operator~() const;

	void neg();


	/**
	 * @brief fill randomly random
	 * @param proportion proportion of true
	 */
	void random(double proportion);

	/**
	 * @brief address of pixel (i,j)
	 * @param i
	 * @param j
	 * @return
	 */
	inline MaskPixelType* address(int i, int j) { return img_mask_.getPixelsPtr() + j*img_mask_.width()+i; }

	inline const MaskPixelType* address(int i, int j) const { return img_mask_.getPixelsPtr() + j*img_mask_.width()+i; }

	/**
	 * @brief compute the erosion of this
	 * @param radius of structuring element (ball)
	 * @return the eroded image
	 */
	MaskBool erosion(int radius) const;

	/**
	 * @brief count number of true pixels
	 * @return
	 */
	int count_true() const;

	/**
	 * @brief check if all pixels are false
	 * @return
	 */
	bool full_false() const;

	/**
	 * @brief apply this op= A
	 * @param m_a Mask a
	 * @param r_a Region of a
	 * @param ind position of region on this
	 * @param op operation uint8 x uint8 -> uint8
	 */
	template<typename OP>
	void apply_op( const Index& ind, const MaskBool& m_a, const Region& r_a, const OP& op);

	/**
	 * @brief operator OR on a region/region
	 * @param ind index of this region (size deduce from r_a)
	 * @param m_a input mask
	 * @param r_a region of a
	 */
	inline void operatorOR(const Index& ind, const MaskBool& m_a, const Region& r_a)
	{
		apply_op(ind, m_a, r_a, [] (MaskPixelType a, MaskPixelType b) -> MaskPixelType
		{
			return a | b;
		});
	}

	/**
	 * @brief operator OR on a region/region
	 * @param ind index of this region (size deduce from r_a)
	 * @param m_a input mask
	 * @param r_a region of a
	 */
	inline void operatorAND(const Index& ind, const MaskBool& m_a, const Region& r_a)
	{
		apply_op(ind, m_a, r_a, [] (MaskPixelType a, MaskPixelType b) -> MaskPixelType
		{
			return a & b;
		});
	}


	/**
	 * @brief operator OR on a region/region
	 * @param ind index of this region (size deduce from r_a)
	 * @param m_a input mask
	 * @param r_a region of a
	 */
	inline void operatorXOR(const Index& ind, const MaskBool& m_a, const Region& r_a)
	{
		apply_op(ind, m_a, r_a, [] (MaskPixelType a, MaskPixelType b) -> MaskPixelType
		{
			return a ^ b;
		});
	}

};





template <typename MASK>
MaskBool& MaskBool::operator&=(const MASK& m)
{
	img_mask_.for_all_pixels([&] (MaskPixelType& p, int i, int j)
	{
		p = (p!=0) && m(i,j);
	});
	return *this;
}

template <typename MASK>
MaskBool& MaskBool::operator|=(const MASK& m)
{
	img_mask_.for_all_pixels([&] (MaskPixelType& p, int i, int j)
	{
		p = (p!=0) || m(i,j);
	});
	return *this;
}

template<typename OP>
void MaskBool::apply_op( const Index& ind, const MaskBool& m_a, const Region& r_a, const OP& op)
{
	Size s = r_a.GetSize();

	const MaskPixelType* in = m_a.address(int(r_a.GetIndex()[0]), int(r_a.GetIndex()[1]));
	MaskPixelType* out = this->address(ind[0], ind[1]);
	int offseti = m_a.get_image().width()-s[0];
	int offseto = this->get_image().width()-s[0];

	for (uint32_t j=0; j<s[1]; ++j)
	{
		for (uint32_t i=0; i<s[0]; ++i)
		{
			*out = op(*out,*in++);
			++out;
		}
		out += offseto;
		in += offseti;
	}
}

/*
template<typename OP>
void MaskBool::affect_oper(const MaskBool& ma, const MaskBool& mb, const OP& op)
{
	assert(ma.img_mask_.is_initialized_as(mb.img_mask_));
	assert(ma.img_mask_.is_initialized_as(img_mask_));

	int nb = img_mask_.width()*img_mask_.height();
	uint8_t* ptr_src1 = ma.img_mask_.getPixelsPtr();
	uint8_t* ptr_src2 = mb.img_mask_.getPixelsPtr();
	uint8_t* ptr_dst = img_mask_.getPixelsPtr();
	for(int i=0;i<nb;++i)
		*ptr_dst++ = op(*ptr_src1++,*ptr_src2++);
}
*/


class ASTEX_API MaskDist : public MaskBool
{
public :
	MaskDist(int width, int height, double _dist, std::vector<std::vector<int>> _coord_seeds );

	inline bool operator()(int i,int j) const
	{
		return img_mask_.pixelAbsolute(i,j) != 0;
	}

	void update_dist(double _dist);
	void update_seeds(std::vector<std::vector<int> > new_seeds);
	double get_ratio_true();

private :
	void update();
	double dist;
	std::vector<std::vector<int>> coord_seeds;
};




template<typename T>
class MaskLargestValue : public MaskImage<T>
{
public:
	MaskLargestValue(const typename MaskImage<T>::MaskImageType& im, double proportion);
	
	inline bool operator()(int i, int j) const 
	{
		return this->img_mask_.pixelAbsolute(i,j) >= m_th ;
	}

    inline int get_number_of_true () const {return m_n_true;}

	inline void set_value(int, int, bool) {	std::cerr << "not possible to set a value in MaskTrue" << std::endl;}

protected:
	int m_n_true; /**< number of pixels inside (test = true) */
    T m_th; /**< value threshold : test(i,j) = ( value(i,j) >= m_th ) */
};

template<typename T>
MaskLargestValue<T>::MaskLargestValue(const typename MaskImage<T>::MaskImageType& im, double proportion):
	MaskImage<T>()
{
	assert (proportion >=0.0 && proportion <=1.0);
	
	this->img_mask_.itk() = im.itk();

	int W = this->img_mask_.width();
	int H = this->img_mask_.height();
	m_n_true = int(W*H*proportion);
	m_n_true = std::min(std::max(m_n_true,1),W*H-1); // Basile : I don't know why itv (below) must : v.begin() < itv < v.end()

	std::vector<T> v;
	v.reserve(W*H);
	
	this->img_mask_.for_all_pixels([&] (ImageGrayd::PixelType& p)
	{
		v.push_back(p);
	});
		
	typename std::vector<T>::iterator itv = v.begin()+(W*H-m_n_true);
	std::nth_element(v.begin(), itv, v.end());
	m_th = *itv;
}

template<typename T>
class MaskSmallestValue : public MaskImage<T>
{
public:
	MaskSmallestValue(const typename MaskImage<T>::MaskImageType& im, double proportion);
	
	inline bool operator()(int i, int j) const 
	{
		return  this->img_mask_.pixelAbsolute(i,j) <= m_th ;
	}

    inline int get_number_of_true () const {return m_n_true;}

	inline void set_value(int, int, bool) {	std::cerr << "not possible to set a value in MaskTrue" << std::endl;}


protected:
	int m_n_true; /**< number of pixels inside (test = true) */
    T m_th; /**< value threshold : test(i,j) = ( value(i,j) >= m_th ) */
};


template<typename T>
MaskSmallestValue<T>::MaskSmallestValue(const typename MaskImage<T>::MaskImageType& im, double proportion):
	MaskImage<T>()
{
	assert (proportion >=0.0 && proportion <=1.0);
	
	this->img_mask_.itk() = im.itk();

	int W = this->img_mask_.width();
	int H = this->img_mask_.height();
	m_n_true = int(W*H*proportion);
	m_n_true = std::min(std::max(m_n_true,1),W*H-1); // Basile : I don't know why itv (below) must : v.begin() < itv < v.end()

	std::vector<T> v;
	v.reserve(W*H);
	
	this->img_mask_.for_all_pixels([&] (ImageGrayd::PixelType& p)
	{
		v.push_back(p);
	});
		
	typename std::vector<T>::iterator itv = v.begin()+m_n_true; // != from MaskLargestValue
	std::nth_element(v.begin(), itv, v.end());
	m_th = *itv;
}



template<typename T>
class MaskAboveThreshold : public MaskImage<T>
{
public:
	using MaskPixelType = typename MaskImage<T>::MaskPixelType;
	using MaskImageType = typename MaskImage<T>::MaskImageType;

	MaskAboveThreshold(const typename MaskImage<T>::MaskImageType& im, T threshold);
	
	inline bool operator()(int i, int j) const 
	{
		return this->img_mask_.pixelAbsolute(i,j) > m_th ;
	}

    inline int get_number_of_true () const {return m_n_true;}

	inline void set_value(int, int, bool) {	std::cerr << "not possible to set a value in MaskTrue" << std::endl;}

protected:
	int m_n_true; /**< number of pixels inside (test = true) */
    double m_th; /**< value threshold : test(i,j) = ( value(i,j) >= m_th ) */
};

template<typename T>
MaskAboveThreshold<T>::MaskAboveThreshold(const MaskImageType& im, T threshold):
	MaskImage<double>(), m_th(threshold)

{
	this->img_mask_.itk() = im.itk();
	m_n_true = 0;
	this->img_mask_.for_all_pixels([&] (MaskPixelType& p)
	{
		if (p > m_th)
			m_n_true++;
	});
}

template<typename T>
class MaskBelowThreshold : public MaskImage<T>
{
public:
	using MaskPixelType = typename MaskImage<T>::MaskPixelType;
	using MaskImageType = typename MaskImage<T>::MaskImageType;

	MaskBelowThreshold(const typename MaskImage<T>::MaskImageType& im, T threshold);
	
	inline bool operator()(int i, int j) const 
	{
		return this->img_mask_.pixelAbsolute(i,j) < m_th;
	}

    inline int get_number_of_true () const {return m_n_true;}

	inline void set_value(int, int, bool) {	std::cerr << "not possible to set a value in MaskTrue" << std::endl;}

protected:
	int m_n_true; /**< number of pixels inside (test = true) */
    double m_th; /**< value threshold : test(i,j) = ( value(i,j) >= m_th ) */
};

template<typename T>
MaskBelowThreshold<T>::MaskBelowThreshold(const MaskImageType& im, T threshold):
	MaskImage<double>(), m_th(threshold)
{
	this->img_mask_.itk() = im.itk();
	m_n_true = 0;
	this->img_mask_.for_all_pixels([&] (MaskPixelType& p)
	{
		if (p < m_th)
			m_n_true++;
	});
}

} // namespace ASTex

#endif
