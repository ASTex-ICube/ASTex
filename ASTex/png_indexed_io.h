#ifndef __ASTEX_PNG_INDEXED_IO__
#define __ASTEX_PNG_INDEXED_IO__

#include "dll.h"

namespace ASTex
{

/**
 * @brief read indexed png image with palette
 * @param filename
 * @param width
 * @param height
 * @param nb_colors
 * @param palette
 * @param data
 * @return error code 0 = ok
 */
int ASTEX_API read_png_palette(const char *filename, int& width, int& height, int& nb_colors, unsigned char*& palette, unsigned char*& data);

/**
 * @brief write_png_palette
 * @param filename
 * @param width
 * @param height
 * @param nb_colors
 * @param palette
 * @param data
 * @return error code 0 = ok
 */
int ASTEX_API write_png_palette(const char *filename, int width, int height, int nb_colors, unsigned char* palette, unsigned char* data);

}

#endif
