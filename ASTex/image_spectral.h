#ifndef __ASTEX_IMAGE_SPECTRAL__
#define __ASTEX_IMAGE_SPECTRAL__

#include "itkImage.h"
#include <ASTex/image_gray.h>

namespace ASTex
{

template <typename INHERIT, bool CST>
class CommonSpectral: public ImageCommon<INHERIT,CST>
{
public:

	typedef typename ImageCommon<INHERIT,CST>::ItkImg ItkImg;
	typedef typename ImageCommon<INHERIT,CST>::IteratorIndexed IteratorIndexed;
	typedef typename ImageCommon<INHERIT,CST>::IteratorIndexed Iterator;
	typedef typename ImageCommon<INHERIT,CST>::PixelType PixelType;
	typedef typename ImageCommon<INHERIT,CST>::ASTexPixelType ASTexPixelType;
	typedef typename ImageCommon<INHERIT,CST>::DataType DataType;

	using Self = ImageCommon<INHERIT, CST>;


	CommonSpectral():
		ImageCommon<INHERIT,CST>()
	{
	}


	CommonSpectral(typename itk::Image< PixelType >::Pointer itk_im):
		ImageCommon<INHERIT,CST>(itk_im)
	{
		this->center_= {{this->width()/2,this->height()/2}};
	}

//	CommonSpectral(int w, int h)
//	{
//		this->initItk(w,h);
//		this->center_= {{w/2, h/2}};
//	}

	explicit CommonSpectral(int s)
	{
		this->initItk(s,s);
		this->center_= {{s/2, s/2}};
	}

	explicit CommonSpectral(int s, bool init_to_zero)
	{
		this->initItk(s,s,init_to_zero);
		this->center_= {{s/2, s/2}};
	}



	void update_itk(const typename ItkImg::Pointer& ptr)
	{
		this->itk_img_ = ptr;
		this->center_= {{this->width()/2,this->height()/2}};
	}

	void initItk(int width, int height, bool init_to_zero = false)
	{
		if (width != height)
			std::cerr << "Spectral images must be square"<< std::endl;
		ImageCommon<INHERIT,CST>::initItk(width,height,init_to_zero);
		this->center_= {{this->width()/2,this->height()/2}};
	}

	void initItk(int width, bool init_to_zero = false)
	{
		this->itk_img_ = ItkImg::New();
		itk::ImageRegion<2> region;
		region.SetIndex(0,0);
		region.SetIndex(1,0);
		region.SetSize(0,width);
		region.SetSize(1,width);
		this->itk_img_->SetRegions(region);
		this->itk_img_->Allocate(init_to_zero);
		this->center_= {{this->width()/2,this->height()/2}}; //width,width ??
	}


	PixelType& pixelModulo( int i, int j)
	{
		typename ItkImg::IndexType pixelIndex = {{ (i+this->center_[0])%this->width() ,(j+this->center_[1])%this->height() }};
		return this->itk_img_->GetPixel(pixelIndex);
	}

	const PixelType& pixelModulo(int i, int j) const
	{
		typename ItkImg::IndexType pixelIndex = {{ (i+this->center_[0])%this->width() ,(j+this->center_[1])%this->height() }};
		return  this->itk_img_->GetPixel(pixelIndex);
	}


	void zero()
	{
		this->for_all_pixels([&](PixelType& p)
		{
			 p = PixelType(0);
		});
	}


	// enable_if INHERIT::PixelType ?

	// operator +
	Self operator + (const Self& img) const
	{
		assert((this->width()==img.width())&&(this->height()==img.height()));

		Self res(this->width(),this->height());

		res.for_all_pixels([&](typename INHERIT::PixelType& p, int x, int y)
		{
			 p = this->pixelAbsolute(x,y) + img.pixelAbsolute(x,y);
		});
		return res;
	}


	// operator +=
	Self& operator += (const Self& img)
	{
		assert((this->width()==img.width())&&(this->height()==img.height()));

		this->for_all_pixels([&](typename INHERIT::PixelType& p, int x, int y)
		{
			 p += img.pixelAbsolute(x,y);
		});
		return *this;
	}


	// operator -
	Self operator - (const Self& img) const
	{
		assert((this->width()==img.width())&&(this->height()==img.height()));

		Self res(this->width(),this->height());

		res.for_all_pixels([&](typename INHERIT::PixelType& p, int x, int y)
		{
			 p = this->pixelAbsolute(x,y) - img.pixelAbsolute(x,y);
		});
		return res;
	}

	// operator -=
	Self& operator -= (const Self& img)
	{
		assert((this->width()==img.width())&&(this->height()==img.height()));

		this->for_all_pixels([&](typename INHERIT::PixelType& p, int x, int y)
		{
			 p -= img.pixelAbsolute(x,y);
		});
		return *this;
	}



	// operator * (img)
	Self operator * (const Self& img) const
	{
		assert((this->width()==img.width())&&(this->height()==img.height()));

		Self res(this->width(),this->height());

		res.for_all_pixels([&](typename INHERIT::PixelType& p, int x, int y)
		{
			 p = this->pixelAbsolute(x,y) * img.pixelAbsolute(x,y);
		});
		return res;
	}

	// operator *= (img
	Self& operator *= (const Self& img)
	{
		assert((this->width()==img.width())&&(this->height()==img.height()));

		this->for_all_pixels([&](typename INHERIT::PixelType& p, int x, int y)
		{
			 p *= img.pixelAbsolute(x,y);
		});
		return *this;
	}




	// operator *
	template<typename SCALAR>
	Self operator * (const SCALAR scal) const
	{
		Self res(this->width(),this->height());

		res.for_all_pixels([&](typename INHERIT::PixelType& p, int x, int y)
		{
			 p = this->pixelAbsolute(x,y) * scal;
		});
		return res;
	}

	// operator *=
	template<typename SCALAR>
	Self operator *= (const SCALAR scal) const
	{

		this->for_all_pixels([&](typename INHERIT::PixelType& p)
		{
			 p *= scal;
		});
		return *this;
	}


	// operator /
	template<typename SCALAR>
	Self operator / (const SCALAR scal) const
	{
		Self res(this->width(),this->height());

		res.for_all_pixels([&](typename INHERIT::PixelType& p, int x, int y)
		{
			 p = this->pixelAbsolute(x,y) / scal;
		});
		return res;
	}

	// operator /=
	template<typename SCALAR>
	Self operator /= (const SCALAR scal)
	{
		this->for_all_pixels([&](typename INHERIT::PixelType& p)
		{
			 p /= scal;
		});
		return *this;
	}

};


// extern operator * (scal * img)
template <typename INHERIT, bool CST, typename SCALAR>
ImageCommon<INHERIT, CST> operator * (const SCALAR scal, const ImageCommon<INHERIT, CST>& img)
{
	ImageCommon<INHERIT, CST> res(img.width(),img.height());
	res.for_all_pixels([&](typename INHERIT::PixelType& p, int x, int y)
	{
		 p = img.pixelAbsolute(x,y) * scal;
	});
	return res;
}




template <class T> using ImageSpectralGen = CommonSpectral< ImageGrayBase< T >, false >;
using ImageSpectralf  = ImageSpectralGen< float >;
using ImageSpectrald  = ImageSpectralGen< double >;
using ImageSpectralcf = ImageSpectralGen< std::complex<float> >;
using ImageSpectralcd = ImageSpectralGen< std::complex<double> >;
template <typename REAL> using ImageSpectralc = ImageSpectralGen< std::complex<REAL> >;

}
#endif



