#ifndef __ASTEX_IMAGE_EASY_IO__
#define __ASTEX_IMAGE_EASY_IO__

#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include <ASTex/image_rgba.h>

namespace ASTex
{
namespace IO
{

void ASTEX_API save01_in_u8(const ImageGrayd& img, const std::string& filename);

void ASTEX_API save01_in_u8(const ImageGrayf& img, const std::string& filename);

void ASTEX_API save01_in_u8(const ImageRGBd& img, const std::string& filename);

void ASTEX_API save01_in_u8(const ImageRGBf& img, const std::string& filename);

void ASTEX_API save01_in_u8(const ImageRGBAd& img, const std::string& filename);

void ASTEX_API save01_in_u8(const ImageRGBAf& img, const std::string& filename);



bool ASTEX_API loadu8_in_01(ImageGrayd& img, const std::string& filename);

bool ASTEX_API loadu8_in_01(ImageGrayf& img, const std::string& filename);

bool ASTEX_API loadu8_in_01(ImageRGBd& img, const std::string& filename);

bool ASTEX_API loadu8_in_01(ImageRGBf& img, const std::string& filename);

bool ASTEX_API loadu8_in_01(ImageRGBAd& img, const std::string& filename);

bool ASTEX_API loadu8_in_01(ImageRGBAf& img, const std::string& filename);

}
}

#endif
