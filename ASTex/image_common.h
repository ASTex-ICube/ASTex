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



#ifndef __ASTEX_IMAGE_COMMON__
#define __ASTEX_IMAGE_COMMON__

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImportImageFilter.h"

#include <ASTex/internal.h>
#include <ASTex/thread_pool.h>
#include <ASTex/dll.h>

#include <thread>
#include <vector>
#include <algorithm>
#include <cstdint>


namespace ASTex
{



extern std::vector<std::string> ASTEX_API file_paths_search;

using Region = itk::ImageRegion<2>;
using Index  = itk::Index<2>;
using Size   = itk::Size<2>;
using Offset = itk::Offset<2>;

/**
 * @brief generate an Index
 * @param x
 * @param y
 * @return
 */
inline Index gen_index(int x, int y)
{
	return Index({{::itk::IndexValueType(x),::itk::IndexValueType(y)}});
}

/**
 * @brief generate a Size
 * @param w
 * @param h
 * @return
 */
inline Size gen_size(int w, int h)
{
	return Size({{::itk::SizeValueType(w),::itk::SizeValueType(h)}});
}


/**
 * @brief generate a region ( Region reg=gen_region(x,y,w,h); )
 * @param x
 * @param y
 * @param w
 * @param h
 * @return th
 */
inline Region gen_region(int x, int y, int w, int h)
{
	return Region({{::itk::IndexValueType(x),::itk::IndexValueType(y)}}, {{::itk::SizeValueType(w),::itk::SizeValueType(h)}});
}

inline Index operator +(const Index& iA, const Index& iB)
{
	return gen_index(iA[0]+iB[0], iA[1]+iB[1]);
}

/**
 * @brief Base class: only for easy type checking
 */
class ImageBase
{};

//template <typename T, bool CST>
//struct CstDcl
//{};

//template <typename T>
//struct CstDcl<T,true>
//{
//	using type = const T&;
//};

//template <typename T>
//struct CstDcl<T,false>
//{
//	using type = T&;
//};


/**
 * template class for overload Gray/RGB/RGBA classes with all common functionnalities.
 * ImageCommon inherit from its first template parameter
 * The second param select the definition of normal class (false) or a constant class (true)
 */
template <typename INHERIT, bool CST >
class ImageCommon: public INHERIT
{
// 2 macros for more readable code
#define NOT_CONST template <bool PIPOBOOL=true>
#define RETURNED_TYPE(TYPE_T) typename std::enable_if<PIPOBOOL && !CST,TYPE_T>::type

public:

	typedef typename INHERIT::ItkImg               ItkImg;
	typedef typename INHERIT::IteratorIndexed      IteratorIndexed;
	typedef typename INHERIT::Iterator             Iterator;
	typedef typename INHERIT::ConstIteratorIndexed ConstIteratorIndexed;
	typedef typename INHERIT::ConstIterator        ConstIterator;

	typedef typename INHERIT::PixelType            PixelType;
	typedef typename INHERIT::DoublePixelType      DoublePixelType;
	typedef typename INHERIT::ASTexPixelType       ASTexPixelType;
	typedef typename INHERIT::DataType             DataType;
	static const uint32_t NB_CHANNELS = INHERIT::NB_CHANNELS;

	using Self = ImageCommon<INHERIT, CST>;


	inline static double normalized_value(DataType v)
	{
		if (std::is_floating_point<DataType>::value) //compile time resolved
			return double(v);

		if (std::is_unsigned<DataType>::value) //compile time resolved
					return double(v) / double(std::numeric_limits<DataType>::max());

		return double(v - std::numeric_limits<DataType>::lowest()) / (double(std::numeric_limits<DataType>::max()) - double(std::numeric_limits<DataType>::lowest()));
	}

	inline static DataType unnormalized_value(double v)
	{
		if (std::is_floating_point<DataType>::value) //compile time resolved
			return DataType(v);

		if (std::is_unsigned<DataType>::value) //compile time resolved
			return DataType(v * std::numeric_limits<DataType>::max());

		return DataType(v*(double(std::numeric_limits<DataType>::max()) - double(std::numeric_limits<DataType>::lowest())) + std::numeric_limits<DataType>::lowest());
	}

	template <bool B=true>
	auto normalized_pixel(const PixelType& p) -> typename std::enable_if<B && std::is_arithmetic<PixelType>::value, DoublePixelType>::type
	{
		return normalized_value(p);
	}

	template <bool B = true>
	auto normalized_pixel(const PixelType& p) -> typename std::enable_if<B && !std::is_arithmetic<PixelType>::value, DoublePixelType>::type
	{
		DoublePixelType q;
		for (uint32_t i = 0; i<INHERIT::NB_CHANNELS; ++i)
			q[i] = Self::normalized_value(p[i]);
		return q;
	}


	template <bool B = true>
	auto unnormalized_pixel(const DoublePixelType& p) -> typename std::enable_if<B && std::is_arithmetic<PixelType>::value, PixelType>::type
	{
		return unnormalized_value(p);
	}

	template <bool B = true>
	auto unnormalized_pixel(const DoublePixelType& p) -> typename std::enable_if<B && !std::is_arithmetic<PixelType>::value, PixelType>::type
	{
		PixelType q;
		for (uint32_t i = 0; i<INHERIT::NB_CHANNELS; ++i)
			q[i] = Self::unnormalized_value(p[i]);
		return q;
	}


	template <bool B=true>
	inline static auto channel(const PixelType& p, uint32_t i)-> typename std::enable_if<!B || std::is_arithmetic<PixelType>::value , DataType>::type
	{
		assert(i==0);
		std::ignore = i;
		return p;
	}

	template <bool B=true>
	inline static auto channel(const PixelType& p, uint32_t i) -> typename std::enable_if<!B || !std::is_arithmetic<PixelType>::value , DataType>::type
	{
		assert(i<NB_CHANNELS);
		return p[i];
	}


protected:
	Offset center_;
	PixelType clampvalue_;

public:

	ImageCommon():
		center_({{0,0}}),
		clampvalue_(DataType(0))
	{}


	/**
	 * constructor for non-constant version of class
	 * second parameter is for SFINAE
	 */
	template <bool PIPOBOOL=true>
	ImageCommon(typename ItkImg::Pointer itk_ptr,typename std::enable_if<PIPOBOOL && !CST>::type* = nullptr):
		INHERIT(itk_ptr), center_({{0,0}}), clampvalue_(DataType(0))
	{}


	/**
	 * constructor for constant version of class
	 * second parameter is for SFINAE
	 */
	template <bool PIPOBOOL=true>
	ImageCommon(typename ItkImg::ConstPointer itk_ptr, typename std::enable_if<PIPOBOOL && CST>::type* = nullptr):
		INHERIT(typename ItkImg::Pointer(const_cast<ItkImg*>(itk_ptr.GetPointer()))), center_({{0,0}})
	{}

	template <bool PIPOBOOL = true>
	ImageCommon(typename ItkImg::Pointer itk_ptr, typename std::enable_if<PIPOBOOL && CST>::type* = nullptr) :
		INHERIT(itk_ptr), center_({ { 0,0 } }), clampvalue_(DataType(0))
	{}


	ImageCommon(int w, int h, bool init_to_zero = false):
		center_({{0,0}}), clampvalue_(DataType(0))
	{
		initItk(w,h,init_to_zero);
	}

	NOT_CONST inline auto getPixelsPtr() -> RETURNED_TYPE(ASTexPixelType*)
	{
		return INHERIT::getPixelsPtr();
	}

	NOT_CONST inline auto getDataPtr() -> RETURNED_TYPE(DataType*)
	{
		return INHERIT::getDataPtr();
	}


	const ASTexPixelType* getPixelsPtr() const
	{
		return INHERIT::getPixelsPtr();
	}

	inline const DataType* getDataPtr() const
	{
		return INHERIT::getDataPtr();
	}



	NOT_CONST auto initItk() -> RETURNED_TYPE(void)
	{
		this->itk_img_ = ItkImg::New();
	}

	bool is_initialized() const
	{
		return this->itk_img_ != nullptr;
	}

	template <typename IMG>
	auto is_initialized_as(const IMG& img) const -> typename std::enable_if<std::is_base_of<ImageBase,IMG>::value, bool>::type
	{
		if (this->itk_img_ == nullptr)
			return false;
		if (img.itk() == nullptr)
			return false;
		return (this->width() == img.width()) && (this->height() == img.height());
	}


	NOT_CONST auto initItk(int width, int height, bool init_to_zero = false) -> RETURNED_TYPE(void)
	{
		this->itk_img_ = ItkImg::New();
		Region region;
		region.SetIndex(0,0);
		region.SetIndex(1,0);
		region.SetSize(0,width);
		region.SetSize(1,height);
		this->itk_img_->SetRegions(region);
		this->itk_img_->Allocate(init_to_zero);
	}

	NOT_CONST auto initItkValue(int width, int height, const PixelType& value ) -> RETURNED_TYPE(void)
	{
		this->itk_img_ = ItkImg::New();
		Region region;
		region.SetIndex(0,0);
		region.SetIndex(1,0);
		region.SetSize(0,width);
		region.SetSize(1,height);
		this->itk_img_->SetRegions(region);
		this->itk_img_->Allocate(false);
		for (auto it = this->beginIterator(); !it.IsAtEnd(); ++it)
			*it = value;
	}


	NOT_CONST auto swap(Self& im) -> RETURNED_TYPE(void)
	{
		assert(is_initialized_as(im));
		std::swap(this->itk_img_,im.itk_img_);
		std::swap(this->center_,im.center_);

	}

	NOT_CONST auto itk() -> RETURNED_TYPE(typename ItkImg::Pointer&)
	{
		return this->itk_img_;
	}

	typename ItkImg::Pointer itk() const
	{
		return this->itk_img_;
	}

	NOT_CONST auto share(const Self& img) -> RETURNED_TYPE(void)
	{
		this->itk_img_ = img.itk();
		this->setCenter(img.getCenter());
	}


	NOT_CONST auto update_itk(const typename ItkImg::Pointer& ptr) -> RETURNED_TYPE(void)
	{
		this->itk_img_ = ptr;
	}


	NOT_CONST auto setCenter(int x, int y) -> RETURNED_TYPE(void)
	{
		center_ = {{x,y}};
	}

	NOT_CONST auto setCenter(const Offset& c) -> RETURNED_TYPE(void)
	{
		center_ = c;
	}


	const Offset& getCenter() const
	{
		return center_;
	}


	NOT_CONST auto pixelAbsolute(const Index& pos) -> RETURNED_TYPE(PixelType&)
	{
		return this->itk_img_->GetPixel(pos);
	}

	const PixelType& pixelAbsolute(const Index& pos) const
	{
		return this->itk_img_->GetPixel(pos);
	}


	NOT_CONST auto pixelAbsolute( int i, int j) -> RETURNED_TYPE(PixelType&)
	{
		return this->itk_img_->GetPixel({{ i ,j }});
	}

	const PixelType& pixelAbsolute(int i, int j) const
	{
		return this->itk_img_->GetPixel({{i,j}});
	}


	NOT_CONST auto pixelRegion( int i, int j) -> RETURNED_TYPE(PixelType&)
	{
		const Index& imgIdx = this->itk_img_->GetLargestPossibleRegion().GetIndex();
		Index pixelIndex = {{ i+imgIdx[0] ,j+imgIdx[1] }};

		return this->itk_img_->GetPixel(pixelIndex);
	}

	const PixelType& pixelRegion(int i, int j) const
	{
		const Index& imgIdx = this->itk_img_->GetLargestPossibleRegion().GetIndex();
		Index pixelIndex = {{ i+imgIdx[0] ,j+imgIdx[1] }};

		return this->itk_img_->GetPixel(pixelIndex);
	}


	NOT_CONST auto pixelRelative( int i, int j) -> RETURNED_TYPE(PixelType&)
	{
		Index pixelIndex = {{ i+center_[0] ,j+center_[1] }};
		return this->itk_img_->GetPixel(pixelIndex);
	}

	const PixelType& pixelRelative(int i, int j) const
	{
		Index pixelIndex = {{ i+center_[0] ,j+center_[1] }};
		return  this->itk_img_->GetPixel(pixelIndex);
	}


	const PixelType& pixelChecked(int i, int j) const
	{
		if( i < 0 || i >= width()  ||  j < 0 || j >= height() )
			return clampvalue_;

		Index pixelIndex = {{ i , j }};
		return  this->itk_img_->GetPixel(pixelIndex);
	}


	const PixelType& pixelRelativeChecked(int i, int j) const
	{
		int ii = i+center_[0];
		int jj = j+center_[1];

		if( ii < 0 || ii >= width()  ||  jj < 0 || jj >= height() )
			return clampvalue_;

		Index pixelIndex = {{ ii , jj }};
		return  this->itk_img_->GetPixel(pixelIndex);
	}

	const PixelType& get_clamp_value() const
	{
		return clampvalue_;
	}

	void set_clamp_value(const PixelType& cv)
	{
		clampvalue_ = cv;
	}



	inline static std::string real_file_name(const std::string& filename)
	{
		{
			std::ifstream fs(filename.c_str());
			if (fs.good())
			{
				fs.close();
				return filename;
			}
		}
		for(const std::string& path: file_paths_search)
		{
			std::string fn = path+filename;
			std::ifstream fs(fn.c_str());
			if (fs.good())
			{
				fs.close();
				return fn;
			}
		}
		std::cerr << "could not open file " << filename << std::endl;
		return std::string("");
	}


	/**
	 * @brief loadImage
	 * @param filename
	 */
	NOT_CONST auto load(const std::string& filename) -> RETURNED_TYPE(bool)
	{
		std::string rfn = real_file_name(filename);

		if (rfn.empty())
			return false;

		typedef itk::ImageFileReader< ItkImg > ReaderLocalType;
		typename ReaderLocalType::Pointer reader = ReaderLocalType::New();
		reader->SetFileName(rfn);

		try
		{
			reader->Update();
		}
		catch (itk::ExceptionObject& /*error*/)
		{
//			error.Print(std::cerr);
			std::cerr << "Error while trying to read file " << filename << std::endl;
			return false;
		}

		this->itk_img_ = reader->GetOutput();
		return true;
	}


	/**
	 * @brief save an image in file
	 * @param filename
	 */
	bool save(const std::string& filename) const
	{
		typedef itk::ImageFileWriter< ItkImg >	WriterLocalType;
		typename WriterLocalType::Pointer writer = WriterLocalType::New();
		writer->SetFileName(filename);
		writer->SetInput(this->itk_img_);

		try
		{
			writer->Update();
		}
		catch (itk::ExceptionObject & /*error*/)
		{
//			error.Print(std::cerr);
			std::cerr << "Error while trying to write file " << filename << std::endl;
			return false;
		}

		return true;
	}


	inline const Size&  size() const
	{
		return this->itk_img_->GetLargestPossibleRegion().GetSize();
	}

	/**
	 * @brief easy access to width of image (warning O(1) but not so fast)
	 * @return width
	 */
	inline int width() const
	{
		return this->itk_img_->GetLargestPossibleRegion().GetSize()[0];
	}

	/**
	 * @brief easy access to height of image (warning O(1) but not so fast)
	 * @return height
	 */
	inline int height() const
	{
		return this->itk_img_->GetLargestPossibleRegion().GetSize()[1];
	}


	///////////////////////////////////////////////
	//////////////    ITERATORS    ////////////////
	///////////////////////////////////////////////


	/**
	 * @brief create a region IteratorIndexed (set to begin)
	 * @param reg the region on which iterate
	 * @return the iterator
	 */
	NOT_CONST inline auto beginIteratorIndexed(const Region& reg) -> RETURNED_TYPE(IteratorIndexed)
	{
		IteratorIndexed iter( this->itk_img_, reg );
		iter.GoToBegin();
		return iter;
	}

	/**
	 * @brief create a region Iterator (set to begin)
	 * @param reg the region on which iterate
	 * @return the iterator
	 */
	NOT_CONST inline auto beginIterator(const Region& reg) -> RETURNED_TYPE(Iterator)
	{
		Iterator iter( this->itk_img_, reg );
		iter.GoToBegin();
		return iter;
	}

	/**
	 * @brief create a region IteratorIndexed (set to begin)
	 * @param x x pos of begin of region
	 * @param y y pos of begin of region
	 * @param w width of region
	 * @param h height of region
	 * @return the iterator
	 */
	NOT_CONST inline auto beginIteratorIndexed(int x, int y, int w, int h) -> RETURNED_TYPE(IteratorIndexed)
	{
		return beginIteratorIndexed(gen_region(x,y,w,h));
	}

	/**
	 * @brief create a region Iterator (set to begin)
	 * @param x x pos of begin of region
	 * @param y y pos of begin of region
	 * @param w width of region
	 * @param h height of region
	 * @return the iterator
	 */
	NOT_CONST inline auto beginIterator(int x, int y, int w, int h) -> RETURNED_TYPE(Iterator)
	{
	   return beginIterator(gen_region(x,y,w,h));
	}


	/**
	 * @brief create a image IteratorIndexed (set to begin)
	 * @return the iterator
	 */
	NOT_CONST inline auto beginIteratorIndexed() -> RETURNED_TYPE(IteratorIndexed)
	{
		IteratorIndexed iter( this->itk_img_, this->itk_img_->GetBufferedRegion() );
		iter.GoToBegin();
		return iter;
	}


	/**
	 * @brief create a image Iterator (set to begin)
	 * @return the iterator
	 */
	NOT_CONST inline auto beginIterator() -> RETURNED_TYPE(Iterator)
	{
		Iterator iter( this->itk_img_, this->itk_img_->GetBufferedRegion() );
		iter.GoToBegin();
		return iter;
	}


	/**
	 * @brief create a region Iterator (set to end)
	 * @param reg the region on which iterate
	 * @return the iterator
	 */
	NOT_CONST inline auto endIterator(const Region& reg) -> RETURNED_TYPE(Iterator)
	{
		Iterator iter( this->itk_img_, reg );
		iter.GoToEnd();
		return iter;
	}

	/**
	 * @brief create a region Iterator (set to end)
	 * @param x x pos of end of region
	 * @param y y pos of end of region
	 * @param w width of region
	 * @param h height of region
	 * @return the iterator
	 */
	NOT_CONST inline auto endIterator(int x, int y, int w, int h) -> RETURNED_TYPE(Iterator)
	{
	   return endIterator(gen_region(x,y,w,h));
	}



	/**
	 * @brief create a image Iterator (set to end)
	 * @return the iterator
	 */
	NOT_CONST inline auto endIterator() -> RETURNED_TYPE(Iterator)
	{
		Iterator iter( this->itk_img_, this->itk_img_->GetBufferedRegion() );
		iter.GoToEnd();
		return iter;
	}



	///////////////////////////////////////////////
	////////////    CONST ITERATOR    /////////////
	///////////////////////////////////////////////


	inline ConstIteratorIndexed beginConstIteratorIndexed(const Region& reg) const
	{
		ConstIteratorIndexed iter( this->itk_img_, reg );
		iter.GoToBegin();
		return iter;
	}

	inline ConstIterator beginConstIterator(const Region& reg) const
	{
		ConstIterator iter( this->itk_img_, reg );
		iter.GoToBegin();
		return iter;
	}

	inline ConstIteratorIndexed beginConstIteratorIndexed(int x, int y, int w, int h) const
	{
		return beginConstIteratorIndexed(gen_region(x,y,w,h));
	}

	inline ConstIterator beginConstIterator(int x, int y, int w, int h) const
	{
		return beginConstIterator(gen_region(x,y,w,h));
	}

	inline ConstIteratorIndexed beginConstIteratorIndexed() const
	{
		ConstIteratorIndexed iter( this->itk_img_, this->itk_img_->GetBufferedRegion() );
		iter.GoToBegin();
		return iter;
	}

	inline ConstIterator beginConstIterator() const
	{
		ConstIterator iter( this->itk_img_, this->itk_img_->GetBufferedRegion() );
		iter.GoToBegin();
		return iter;
	}


	inline ConstIterator endConstIterator(const Region& reg) const
	{
		ConstIterator iter( this->itk_img_, reg );
		iter.GoToEnd();
		return iter;
	}


	inline ConstIterator endConstIterator(int x, int y, int w, int h) const
	{
		return endConstIterator(gen_region(x,y,w,h));
	}


	inline ConstIterator endConstIterator() const
	{
		ConstIterator iter( this->itk_img_, this->itk_img_->GetBufferedRegion() );
		iter.GoToEnd();
		return iter;
	}



	NOT_CONST auto copy_pixels(const Self& img) -> RETURNED_TYPE(void)
	{
		assert(compatible_size(*this,img));
		auto it_in = img.beginConstIterator();
		for (auto it_out = this->beginIterator(); !it_out.IsAtEnd(); ++it_out,++it_in)
			it_out.Value() = it_in.Value();
	}



	NOT_CONST auto copy_pixels(const Index& dst, const Self& img_src, const Region& reg_src) -> RETURNED_TYPE(void)
	{
		Region reg_dst = {dst, reg_src.GetSize()};

		auto itf = img_src.beginConstIterator(reg_src);
		for (auto it = this->beginIterator(reg_dst); !it.IsAtEnd(); ++it,++itf)
			it.Value() = itf.Value();
	}


	///////////////////////////////////////////////
	//////////////    TRAVERSAL    ////////////////
	///////////////////////////////////////////////

	/**
	 * traversal of all pixels
	 * @param f fonction/lamba to apply with 1 param: PixelType p&
	 */
	template <typename FUNC>
	inline auto for_all_pixels(const FUNC& f) ->  typename std::enable_if<!CST && function_traits<FUNC>::arity==1, void>::type
	{
		for (auto it = this->beginIterator(); !it.IsAtEnd(); ++it)
			f(it.Value());
	}



	/**
	 * traversal of all pixels
	 * @param f fonction/lamba to apply with 3 params: PixelType p&, int x, int y
	 */
	template <typename FUNC>
	inline auto for_all_pixels(const FUNC& f) -> typename std::enable_if<!CST && function_traits<FUNC>::arity==3, void>::type
	{
		for (int y=0,lh=this->height(); y<lh; ++y)
			for (int x=0,lw=this->width(); x<lw; ++x)
			{
				f(itk()->GetPixel({{x,y}}), x, y);
			}
	}


	/**
	 * traversal of pixels that belong to a mask
	 * @param f fonction/lamba to apply with with 3 params: PixelType p&, int x, int y
	 * @param mask Mask/lamba for filtering
	 */
	template <typename FUNC, typename MASK>
	auto for_all_pixels(const FUNC& f, const MASK& mask) -> typename std::enable_if<!CST && ((function_traits<FUNC>::arity==1)||(function_traits<FUNC>::arity==3)), void>::type
	{
		for (int y=0,lh=this->height(); y<lh; ++y)
			for (int x=0,lw=this->width(); x<lw; ++x)
			{
				if (mask(x,y))
				{
					PixelType& P = itk()->GetPixel({{x,y}});
					four_param_binder(f,P,x,y,uint16_t(0));
				}
			}
	}



	/**
	 * traversal of region
	 * @param x,y,w,h region
	 * @param f fonction/lamba to apply that can have: 3 params: PixelType p&, int x, int y or 1 (p)
	 */
	template <typename FUNC>
	auto for_region_pixels(int x, int y, int w, int h, const FUNC& f) -> typename std::enable_if<!CST && ((function_traits<FUNC>::arity==3) || (function_traits<FUNC>::arity==1)), void>::type
	{
		int ye = y+h;
		int xe = x+w;

		for (int j=y; j<ye; ++j)
			for (int i=x; i<xe; ++i)
			{
				PixelType& P = itk()->GetPixel({{i,j}});
				four_param_binder(f,P,i,j,uint16_t(0));
			}
	}

	/**
	 * traversal of region
	 * @param region region
	 * @param f fonction/lamba to apply that can have: 3 params: PixelType p&, int x, int y or 1 (p)
	 */
	template <typename FUNC>
	auto for_region_pixels(const Region& reg, const FUNC& f) -> typename std::enable_if<!CST && ((function_traits<FUNC>::arity==3) || (function_traits<FUNC>::arity==1)), void>::type
	{
		const Index& I = reg.GetIndex();
		const Size& S = reg.GetSize();
		for_region_pixels(I[0], I[1], S[0], S[1], f);
	}


	///////////////////////////////////////////////
	////////////// CONST TRAVERSAL ////////////////
	///////////////////////////////////////////////

	/**
	 * traversal of all pixels
	 * @param f fonction/lamba to apply with 2 params:  int x, int y
	 */
	template <typename FUNC>
	inline auto for_all_pixels(const FUNC& f) const ->  typename std::enable_if<function_traits<FUNC>::arity==2, void>::type
	{
		for (int y=0,lh=this->height(); y<lh; ++y)
			for (int x=0,lw=this->width(); x<lw; ++x)
			{
				f(x, y);
			}
	}

	template <typename FUNC>
	inline auto for_all_pixels(const FUNC& f) const ->  typename std::enable_if<function_traits<FUNC>::arity==1, void>::type
	{
		for (auto it = this->beginConstIterator(); !it.IsAtEnd(); ++it)
		{
			f(it.Value());
		}
	}


	template <typename FUNC>
	inline auto for_all_pixels(const FUNC& f) const -> typename std::enable_if<function_traits<FUNC>::arity==3, void>::type
	{
		for (int y=0,lh=this->height(); y<lh; ++y)
			for (int x=0,lw=this->width(); x<lw; ++x)
			{
				const PixelType& P = itk()->GetPixel({{x,y}});
				f(P, x, y);
			}
	}


	template <typename FUNC, typename MASK>
	inline auto for_all_pixels(const FUNC& f, const MASK& mask) const -> typename std::enable_if<(function_traits<FUNC>::arity==1)||(function_traits<FUNC>::arity==3), void>::type
	{
		for (int y=0,lh=this->height(); y<lh; ++y)
			for (int x=0,lw=this->width(); x<lw; ++x)
			{
				if (mask(x,y))
				{
					const PixelType& P = itk()->GetPixel({{x,y}});
					four_param_binder(f,P,x,y,uint16_t(0));
				}
			}
	}


	/**
	 * const traversal of region
	 * @param x,y,w,h region
	 * @param f fonction/lamba to apply that can have: 3 params: PixelType p&, int x, int y or 1 (p)
	 */
	template <typename FUNC>
	auto for_region_pixels(int x, int y, int w, int h, const FUNC& f) const -> typename std::enable_if<((function_traits<FUNC>::arity==3) || (function_traits<FUNC>::arity==1)), void>::type
	{
		int ye = y+h;
		int xe = x+w;

		for (int j=y; j<ye; ++j)
			for (int i=x; i<xe; ++i)
			{
				const PixelType& P = itk()->GetPixel({{i,j}});
				four_param_binder(f,P,i,j,uint16_t(0));
			}
	}

	/**
	 * const traversal of region
	 * @param region region
	 * @param f fonction/lamba to apply that can have: 3 params: PixelType p&, int x, int y or 1 (p)
	 */
	template <typename FUNC>
	auto for_region_pixels(const Region& reg, const FUNC& f) const -> typename std::enable_if<((function_traits<FUNC>::arity==3) || (function_traits<FUNC>::arity==1)), void>::type
	{
		const Index& I = reg.GetIndex();
		const Size& S = reg.GetSize();
		for_region_pixels(I[0], I[1], S[0], S[1], f);
	}







	///////////////////////////////////////////////
	///////////// PARALLEL TRAVERSAL //////////////
	///////////////////////////////////////////////


	template <typename FUNC>
	void parallel_for_region_pixels(int xb, int yb, int w, int h, const FUNC& f)// -> typename std::enable_if<!CST, void>::type
	{
		ThreadPool* thread_pool = internal_thread_pool();

		uint16_t nbt = nb_launched_threads();//std::thread::hardware_concurrency() * 8;

		std::vector<ThreadPool::Future> futures;
		futures.reserve(nbt);
		std::atomic_bool finished(false);

		int lh = yb+h;
		int lw = xb+w;

		for (uint16_t i = 0; i<nbt; ++i)
		{
			futures.push_back(thread_pool->enqueue([i, nbt, &finished, lw, lh, xb, yb, &f, this]()
			{
				for (int y = yb+i; y<lh; y += nbt)
					for (int x = xb; x<lw; ++x)
					{
						PixelType& P = itk()->GetPixel({ { x,y } });
						if (!four_param_binder(f, P, x, y, i))
						{
							finished = true;
							return;
						}
						if (return_bool<FUNC>::value && finished) return;
					}
			}));
		}

		for (auto& f : futures)
			f.wait();

	}


	template <typename FUNC>
	void parallel_for_region_pixels(const Region& reg, const FUNC& f)
	{
		const Index& I = reg.GetIndex();
		const Size& S = reg.GetSize();
		parallel_for_region_pixels(I[0], I[1], S[0], S[1], f);
	}





	//template <typename FUNC,  typename MASK,
	//		  typename std::enable_if<!CST && ((function_traits<FUNC>::arity>=1)||(function_traits<FUNC>::arity<=4))>::type* = nullptr>
	//void parallel_for_all_pixels(const FUNC& f, const MASK& mask)
	//{
	//	int nbt = std::thread::hardware_concurrency();
	//	if (nbt > this->height())
	//		nbt = this->height();
	//	std::vector<std::thread*> ths(nbt);

	//	for (int i=0; i<nbt;++i)
	//	{
	//		ths[i] = new std::thread([i,nbt,&f,&mask,this]()
	//		{
	//			for (int y=i,lh=this->height(); y<lh; y+=nbt)
	//				for (int x=0,lw=this->width(); x<lw; ++x)
	//				{
	//					if (mask(x,y))
	//					{
	//						PixelType& P = itk()->GetPixel({{x,y}});
	//						four_param_binder(f,P,x,y,i);
	//					}
	//				}
	//		});
	//	}

	//	for (int i=0;i<nbt;++i)
	//	{
	//		ths[i]->join();
	//		delete ths[i];
	//	}
	//}

	template <typename FUNC, typename MASK,
		typename std::enable_if<!CST && (function_traits<FUNC>::arity >= 1) && (function_traits<FUNC>::arity <= 4)>::type* = nullptr>
		void parallel_for_all_pixels(const FUNC& f, const MASK& mask)
	{
		ThreadPool* thread_pool = internal_thread_pool();
		
		uint16_t nbt = nb_launched_threads();
		
		std::vector<ThreadPool::Future> futures;
		futures.reserve(nbt);
		std::atomic_bool finished(false);
		
		int lw = this->width();
		int lh = this->height();

		for (uint16_t i = 0; i<nbt; ++i)
		{
			futures.push_back(thread_pool->enqueue( [i, nbt, lw, lh, &finished, &f, &mask, this]()
			{
				for (int y = i; y<lh; y += nbt)
					for (int x = 0; x<lw; ++x)
					{
						if (mask(x, y))
						{
							PixelType& P = itk()->GetPixel({ { x,y } });
							if (!four_param_binder(f, P , x, y, i))
							{
								finished = true;
								return;
							}
						}
						if (return_bool<FUNC>::value && finished) return;
					}
			}));
		}

		for (auto& f : futures)
			f.wait();
	}


	template <typename FUNC,
			  typename std::enable_if<!CST && (function_traits<FUNC>::arity>=1) && (function_traits<FUNC>::arity<=4)>::type* = nullptr>
	inline void parallel_for_all_pixels(const FUNC& f)
	{
		parallel_for_all_pixels(f,[](int,int){return true;});
	}



	///////////////////////////////////////////////
	////////// CONST PARALLEL TRAVERSAL ///////////
	///////////////////////////////////////////////

	template <typename FUNC>
	void parallel_for_region_pixels(int xb, int yb, int w, int h, const FUNC& f) const
	{
		ThreadPool* thread_pool = internal_thread_pool();

		uint16_t nbt = nb_launched_threads();

		std::vector<ThreadPool::Future> futures;
		futures.reserve(nbt);
		std::atomic_bool finished(false);

		int lh = yb+h;
		int lw = xb+w;

		for (uint16_t i = 0; i<nbt; ++i)
		{
			futures.push_back(thread_pool->enqueue([i, nbt, &finished, xb, yb, lw, lh, &f, this]()
			{
				for (int y = yb+i; y<lh; y += nbt)
					for (int x = xb; x<lw; ++x)
					{
						const PixelType& P = itk()->GetPixel({ { x,y } });
						if (!four_param_binder(f, P, x, y, i))
						{
							finished = true;
							return;
						}
						if (return_bool<FUNC>::value && finished) return;
					}
			}));
		}

		for (auto& f : futures)
			f.wait();

	}

	template <typename FUNC>
	void parallel_for_region_pixels(const Region& reg, const FUNC& f) const
	{
		const Index& I = reg.GetIndex();
		const Size& S = reg.GetSize();
		parallel_for_region_pixels(I[0], I[1], S[0], S[1], f);
	}



	//template <typename FUNC,  typename MASK,
	//		  typename std::enable_if<(function_traits<FUNC>::arity>=1)&&(function_traits<FUNC>::arity<=4)>::type* = nullptr>
	//void parallel_for_all_pixels(const FUNC& f, const MASK& mask) const
	//{
	//	int nbt = std::thread::hardware_concurrency();
	//	if (nbt > this->height())
	//		nbt = this->height();
	//	std::vector<std::thread*> ths(nbt);

	//	for (int i=0; i<nbt;++i)
	//	{
	//		ths[i] = new std::thread([i,nbt,&f,&mask,this]()
	//		{
	//			for (int y=i,lh=this->height(); y<lh; y+=nbt)
	//				for (int x=0,lw=this->width(); x<lw; ++x)
	//				{
	//					if (mask(x,y))
	//					{
	//						const PixelType& P = itk()->GetPixel({{x,y}});
	//						four_param_binder(f,P,x,y,i);
	//					}
	//				}
	//		});
	//	}

	//	for (int i=0;i<nbt;++i)
	//	{
	//		ths[i]->join();
	//		delete ths[i];
	//	}
	//}

	
	template <typename FUNC, typename MASK,
		typename std::enable_if<(function_traits<FUNC>::arity >= 1) && (function_traits<FUNC>::arity <= 4)>::type* = nullptr>
		void parallel_for_all_pixels(const FUNC& f, const MASK& mask) const
	{
		ThreadPool* thread_pool = internal_thread_pool();

		uint16_t nbt = nb_launched_threads();

		std::vector<ThreadPool::Future> futures;
		futures.reserve(nbt);
		std::atomic_bool finished(false);

		int lw = this->width();
		int lh = this->height();

		for (uint16_t i = 0; i<nbt; ++i)
		{
			futures.push_back(thread_pool->enqueue([i, nbt, &finished, lw, lh, &f, &mask, this]()
			{
				for (int y = i; y<lh; y += nbt)
					for (int x = 0; x<lw; ++x)
					{
						if (mask(x, y))
						{
							const PixelType& P = itk()->GetPixel({ { x,y } });
							if (!four_param_binder(f, P, x, y, i))
							{
								finished = true;
								return;
							}
						}
						if (return_bool<FUNC>::value && finished) return;
					}
			}));
		}

		for (auto& f : futures)
			f.wait();
	}



	template <typename FUNC,
		typename std::enable_if<(function_traits<FUNC>::arity >= 1) && (function_traits<FUNC>::arity <= 4)>::type* = nullptr>
	inline void parallel_for_all_pixels(const FUNC& f) const
	{
		parallel_for_all_pixels(f,[](int,int){return true;});
	}

#undef NOT_CONST
#undef RETURNED_TYPE

};




template <typename T_IN, typename T_OUT>
inline bool compatible_size(const T_IN& in, const T_OUT& out)
{
	return (in.width() == out.width()) && (in.height() == out.height());
}

/**
 * traversal of all pixels of two images
 * @param in input image
 * @param out output image
 * @param f fonction/lamba to apply
 */
template <typename T_IN, typename T_OUT, typename FUNC>
auto for_all_pixels(const T_IN& in, T_OUT& out, const FUNC& f)
	-> typename std::enable_if<std::is_base_of<ImageBase,T_IN>::value && std::is_base_of<ImageBase,T_OUT>::value && function_traits<FUNC>::arity==2, void>::type
{
	assert(compatible_size(in,out));

	auto it_in = in.beginConstIterator();
	for (auto it_out = out.beginIterator(); !it_out.IsAtEnd(); ++it_out,++it_in)
		f(it_in.Value(),it_out.Value());
}




template <typename T_IN, typename T_OUT, typename FUNC>
auto apply_all_pixels(const T_IN& in, T_OUT& out, const FUNC& f)
	-> typename std::enable_if<std::is_base_of<ImageBase,T_IN>::value && std::is_base_of<ImageBase,T_OUT>::value && function_traits<FUNC>::arity==1, void>::type
{
	assert(compatible_size(in,out));

	auto it_in = in.beginConstIterator();
	for (auto it_out = out.beginIterator(); !it_out.IsAtEnd(); ++it_out,++it_in)
		it_out.Value() = f(it_in.Value());
}




template <typename FUNC,
		  typename std::enable_if<function_traits<FUNC>::arity == 2>::type* = nullptr>
inline void for_indices(int begin_x, int end_x, int begin_y, int end_y, const FUNC& f)
{
	for (int y=begin_y; y<end_y; ++y)
		for (int x=begin_x; x<end_x; ++x)
			f(x,y);
}



//template <typename FUNC,
//	typename std::enable_if<(function_traits<FUNC>::arity >= 2) && (function_traits<FUNC>::arity <= 3)>::type* = nullptr>
//	void parallel_for_indices(int begin_x, int end_x, int begin_y, int end_y, const FUNC& f)
//{
//	int nbt = std::thread::hardware_concurrency();
//	int height = end_y - begin_y;
//	if (nbt > height)
//		nbt = height;
//	std::vector<std::thread*> ths(nbt);

//	for (int i = 0; i<nbt; ++i)
//	{
//		ths[i] = new std::thread([i, nbt, begin_x, end_x, begin_y, end_y, &f]()
//		{
//			for (int y = begin_y; y < end_y; y += nbt)
//				for (int x = begin_x; x < end_x; ++x)
//				{
//					four_param_binder(f, x, y, i, 0);
//				}
//		});
//	}

//	for (int i = 0; i<nbt; ++i)
//	{
//		ths[i]->join();
//		delete ths[i];
//	}
//}


template <typename FUNC,
			  typename std::enable_if<function_traits<FUNC>::arity == 2>::type* = nullptr>
inline void for_region(const Region& region, const FUNC& f)
{
	int yb = region.GetIndex()[1];
	int xb = region.GetIndex()[0];

	int ye = xb+region.GetSize()[1];
	int xe = yb+region.GetSize()[0];

	for (int y=yb; y<ye; ++y)
		for (int x=xb; x<xe; ++x)
			f(x,y);
}

template <typename FUNC,
		  typename std::enable_if<function_traits<FUNC>::arity == 4>::type* = nullptr>
void for_regions(const Region& r1,const Region& r2, const FUNC &f)
{
	int i0 = r1.GetIndex()[0];
	int j0 = r1.GetIndex()[1];
	int ie = i0+r1.GetSize()[0];
	int je = j0+r1.GetSize()[1];
	int ii0 = r2.GetIndex()[0];
	int jj0 = r2.GetIndex()[1];

	for(int j=j0, jj=jj0; j<je; ++j,++jj)
		for(int i=i0, ii=ii0; i<ie; ++i,++ii)
		{
			f(i,j,ii,jj);
		}
}




// A TESTER

//template <typename T_IN, typename T_OUT, typename FUNC>
//auto parallel_for_all_pixels(const T_IN& in, T_OUT& out, const FUNC& f)
//  -> typename std::enable_if<std::is_base_of<ImageBase,T_IN>::value && std::is_base_of<ImageBase,T_OUT>::value && function_traits<FUNC>::arity==2, void>::type
//{
//	assert(compatible_size(in,out));

//	int nbt = std::thread::hardware_concurrency();
//	std::vector<std::thread*> ths(nbt);

//	std::vector<Region> regions;
//	regions.reserve(nbt);
//	for(int i=0;i<nbt;++i)
//		regions.emplace_back(gen_region(0, i, in.width(), 1));

//	for (int i=0; i<nbt;++i)
//	{
//		ths[i] = new std::thread([i,nbt,&in,&out,&regions,&f]()
//		{
//			while (regions[i].GetIndex()[1] < in.height())
//			{
//				auto it_in = in.beginIterator();
//				for (auto it_out = out.beginIterator(); !it_out.IsAtEnd(); ++it_out,++it_in)
//					f(it_in.Value(),it_out.Value());
//				regions[i].GetModifiableIndex()[1] += nbt;
//			}
//		});
//	}

//	for (int i=0;i<nbt;++i)
//	{
//		ths[i]->join();
//		delete ths[i];
//	}
//}



/**
 * @brief create an image from a buffer (no copy)
 * @param w image width
 * @param h image height
 * @param buffer pixel buffer
 * @param manage_dest delegate buffer destruction to image
 * @return the image
 */
template <typename IMG >
IMG create_from_buffer(int w, int h, char * buffer, bool manage_dest)
{
	using PixelType = typename IMG::PixelType;
	using ImportFilterType = itk::ImportImageFilter<PixelType>;

	typename ImportFilterType::Pointer importer = ImportFilterType::New();

	Region region = gen_region(0,0,w,h);
	importer->SetRegion(region);
//	importer->SetSpacing({0.0,0.0});
//	importer->SetOrigin({0.0,0.0});

	PixelType* pixelData = reinterpret_cast<PixelType*>(buffer);

	importer->SetImportPointer(pixelData, w*h, manage_dest);
	importer->Update();

	return IMG(importer->GetOutput());
}



}

#endif



