#ifndef _HV_ASTEX_IO_H_
#define _HV_ASTEX_IO_H_

#include "hv/hvPictRGB.h"
#include "hv/hvPict.h"
#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>

namespace ASTex
{
template<typename T>
void load_to_hv(const std::string& filename, hview::hvPictRGB<T>& pict)
{
	ImageRGB<T> img;
	img.load(filename);
	pict.reset(img.width(),img.height(), hview::hvColRGB<T>(0));
	auto src = img.getDataPtr();
	auto dst = pict.data();
	memcpy(reinterpret_cast<uint8_t*>(dst),reinterpret_cast<uint8_t*>(src),img.width()*img.height()*3*sizeof(T));
}

template<typename T>
void save_from_hv(const std::string& filename, hview::hvPictRGB<T>& pict)
{
	ImageRGB<T> img = create_from_buffer<ImageRGB<T>>(pict.sizeX(), pict.sizeY(), reinterpret_cast<char*>(pict.data()), false);
	img.save(filename);
}

template<typename T>
void load_to_hv(const std::string& filename, hview::hvPict<T>& pict)
{
	ImageGray<T> img;
	img.load(filename);
	pict.reset(img.width(),img.height());
	auto src = img.getDataPtr();
	auto dst = pict.data();
	memcpy(reinterpret_cast<uint8_t*>(dst),reinterpret_cast<uint8_t*>(src),img.width()*img.height()*3*sizeof(T));
}

template<typename T>
void save_from_hv(const std::string& filename, hview::hvPict<T>& pict)
{
	ImageGray<T> img = create_from_buffer<ImageGray<T>>(pict.sizeX(), pict.sizeX(), reinterpret_cast<char*>(pict.data()), false);
	img.save(filename);
}

}

#endif
