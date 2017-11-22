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



#ifndef __ASTEX_IMG_MERGING__
#define __ASTEX_IMG_MERGING__

#include <ASTex/image_gray.h>

namespace ASTex
{



/**
 * \class image_collector io.h
 * \brief collects images and builds a patchwork
 */

namespace Assembler1D
{
const int HorizontalFlush = -1;
const int VerticalFlush = -2;

template<typename IMG>
class ASTEX_API Assembler1D
{	
public:
	std::vector<const IMG*> images_;
	std::vector<Region> regions_;
	IMG& result_;

	Assembler1D(IMG& im):
		result_(im)
	{}

	Assembler1D(Assembler1D&& x):
		result_(x.result_)
	{

		images_ = std::move(x.images_);
		regions_= std::move(x.regions_);
	}

	void clear()
	{
		images_.clear();
		regions_.clear();
	}
};

template<typename IMG>
inline Assembler1D<IMG> into(IMG& im)
{
	return Assembler1D<IMG>(im);
}

template<typename IMG>
Assembler1D<IMG>& operator <<(Assembler1D<IMG>& asmb, const IMG& im)
{
	asmb.images_.push_back(&im);
	asmb.regions_.push_back(gen_region(0,0,im.width(), im.height()));
	return (asmb);
}

template<typename IMG>
Assembler1D<IMG>& operator <<(Assembler1D<IMG>& asmb, const Region& reg)
{
	asmb.regions_.back() = reg;
	return (asmb);
}

template<typename IMG>
Assembler1D<IMG>&& operator <<(Assembler1D<IMG>&& asmb, const IMG& im)
{
	asmb.images_.push_back(&im);
	asmb.regions_.push_back(gen_region(0,0,im.width(), im.height()));
	return std::move(asmb);
}

template<typename IMG>
Assembler1D<IMG>&& operator <<(Assembler1D<IMG>&& asmb, const Region& reg)
{
	asmb.regions_.back() = reg;
	return std::move(asmb);
}


template<typename T>
void internal_intermediate(T asmb, int sp)
{
	if (sp==HorizontalFlush)
	{
		// size of new image
		int w=0;
		int h=0;
		for (auto& r: asmb.regions_)
		{
			w += r.GetSize()[0];
			h = std::max(h,int(r.GetSize()[1]));
		}
		asmb.result_.initItk(w,h);

		std::size_t nb = asmb.images_.size();
		int x=0;
		for(std::size_t i=0; i<nb; ++i)
		{
			if (asmb.images_[i]!=nullptr)
				asmb.result_.copy_pixels(gen_index(x,0),*(asmb.images_[i]), asmb.regions_[i]);
			x += asmb.regions_[i].GetSize()[0];
		}
		asmb.clear();
	}
	else if (sp==VerticalFlush)
	{
		// size of new image
		int w=0;
		int h=0;
		for (auto& r: asmb.regions_)
		{
			h += r.GetSize()[1];
			w = std::max(w,int(r.GetSize()[0]));
		}
		asmb.result_.initItk(w,h);

		std::size_t nb = asmb.images_.size();
		int y=0;
		for(std::size_t i=0; i<nb; ++i)
		{
			if (asmb.images_[i]!=nullptr)
				asmb.result_.copy_pixels(gen_index(0,y),*(asmb.images_[i]), asmb.regions_[i]);
			y += asmb.regions_[i].GetSize()[1];
		}
		asmb.clear();
	}
	else
	{
		asmb.images_.push_back(nullptr);
		asmb.regions_.push_back(gen_region(0,0,sp,sp));
	}
}


template<typename IMG>
Assembler1D<IMG>& operator << (Assembler1D<IMG>& asmb, int sp)
{
	internal_intermediate<Assembler1D<IMG>&>(std::forward<Assembler1D<IMG>&>(asmb),sp);
	return asmb;
}

template<typename IMG>
Assembler1D<IMG>&& operator << (Assembler1D<IMG>&& asmb, int sp)
{
	internal_intermediate<Assembler1D<IMG>&&>(std::forward<Assembler1D<IMG>&&>(asmb),sp);
	return std::move(asmb);
}

}



namespace Assembler2D
{
const int FinalFlush = -1;

inline int EndLine(int sp) { return -100-sp;}


template<typename IMG>
class ASTEX_API Assembler2D
{
public:
	std::vector<const IMG*> images_;
	std::vector<Region> regions_;
	std::vector<int> lines_;
	std::size_t next_;
	IMG& result_;

	Assembler2D(IMG& im):
		next_(0),
		result_(im)
	{}

	Assembler2D(Assembler2D&& x):
		next_(x.next_),
		result_(x.result_)
	{
		images_ = std::move(x.images_);
		regions_= std::move(x.regions_);
		lines_= std::move(x.lines_);
	}

	void clear()
	{
		images_.clear();
		regions_.clear();
		lines_.clear();
	}
};

template<typename IMG>
inline Assembler2D<IMG> into(IMG& im)
{
	return Assembler2D<IMG>(im);
}

template<typename IMG>
Assembler2D<IMG>& operator <<(Assembler2D<IMG>& asmb, const IMG& im)
{
	asmb.images_.push_back(&im);
	asmb.regions_.push_back(gen_region(0,0,im.width(), im.height()));
	return (asmb);
}

template<typename IMG>
Assembler2D<IMG>& operator <<(Assembler2D<IMG>& asmb, const Region& reg)
{
	asmb.regions_.back() = reg;
	return (asmb);
}

template<typename IMG>
Assembler2D<IMG>&& operator <<(Assembler2D<IMG>&& asmb, const IMG& im)
{
	asmb.images_.push_back(&im);
	asmb.regions_.push_back(gen_region(0,0,im.width(), im.height()));
	return std::move(asmb);
}

template<typename IMG>
Assembler2D<IMG>&& operator <<(Assembler2D<IMG>&& asmb, const Region& reg)
{
	asmb.regions_.back() = reg;
	return std::move(asmb);
}


template<typename T>
void internal_intermediate(T asmb, int sp)
{
	if ((sp<=EndLine(0)) || (sp==FinalFlush))
	{
		// size of new image
		int w=0;
		int h=0;
		std::size_t nb = asmb.images_.size();
		std::cout << "de "<< asmb.next_ << "  a "<< nb << std::endl;

		for (std::size_t i=asmb.next_; i<nb;++i)
		{
			auto& r = asmb.regions_[i];
			w += r.GetSize()[0];
			h = std::max(h,int(r.GetSize()[1]));
		}
		asmb.next_ = nb;
		asmb.lines_.push_back(nb);
		asmb.lines_.push_back(w);
		if (sp==FinalFlush)
			asmb.lines_.push_back(h);
		else
			asmb.lines_.push_back(h-sp+EndLine(0));
	}
	if (sp==FinalFlush)
	{
		// size of new image
		int w=0;
		int h=0;

		for(auto x:asmb.lines_)
			std::cout<< x <<" ";
		std::cout << std::endl;

		for (std::size_t i=0;i< asmb.lines_.size();i+=3)
		{
			w = std::max(w,asmb.lines_[i+1]);
			h += asmb.lines_[i+2];
		}
		std::cout <<"SIZE="<< w <<","<<h<< std::endl;
		asmb.result_.initItk(w,h);

		std::size_t nb = asmb.images_.size();
		int x=0;
		int y=0;
		auto lit = asmb.lines_.begin();
		for(std::size_t i=0; i<nb; ++i)
		{
			if (int(i)==*lit)
			{
				++lit; ++lit;
				y += *lit++;
				x=0;
			}
			std::cout << i <<" -> "<<x <<" , "<< y << std::endl;
			if (asmb.images_[i]!=nullptr)
				asmb.result_.copy_pixels(gen_index(x,y),*(asmb.images_[i]), asmb.regions_[i]);
			x += asmb.regions_[i].GetSize()[0];
		}
		asmb.clear();
	}
	else
	{
		asmb.images_.push_back(nullptr);
		asmb.regions_.push_back(gen_region(0,0,sp,0));
	}
}


template<typename IMG>
Assembler2D<IMG>& operator << (Assembler2D<IMG>& asmb, int sp)
{
	internal_intermediate<Assembler2D<IMG>&>(std::forward<Assembler2D<IMG>&>(asmb),sp);
	return asmb;
}

template<typename IMG>
Assembler2D<IMG>&& operator << (Assembler2D<IMG>&& asmb, int sp)
{
	internal_intermediate<Assembler2D<IMG>&&>(std::forward<Assembler2D<IMG>&&>(asmb),sp);
	return std::move(asmb);
}

}

}
#endif
