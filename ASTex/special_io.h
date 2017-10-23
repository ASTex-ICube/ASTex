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



#ifndef __ASTEX_SPECIAL_IO__
#define __ASTEX_SPECIAL_IO__

#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include <ASTex/image_spectral.h>

/** \file
 * \brief inputs, outputs, and monitoring tools for scalar images
*/

/**
 * \brief namespace ASTex
 */
namespace ASTex
{

/**
 * \brief mean value of a scalar image
 * \param [in] input input scalar image
 * \return mean pixel value
 */
double ASTEX_API getMean (const ImageGrayd& input);

/**
 * \brief standard deviation of a scalar image
 * \param [in] input input scalar image
 * \return standard deviation of pixel value
 */
double ASTEX_API getStDev (const ImageGrayd& input);


namespace IO
{

/**
 * \brief prints first order image statistics (mean, standard deviation, min, max) on std::cout
 * \param [in] input input scalar image
 * \param [in] intro a string that printing will start by.
 */
void ASTEX_API  monitorStats (const ImageGrayd& input, std::string intro ="image stats");

/**
 * \brief load scalar image from file
 * \param [out] image scalar image
 * \param [in] filename
 * \param [in] min min value to be loaded (mapped to 0)
 * \param [in] max max value to be loaded (mapped to 255)
 * \pre the file contains grayscale on 1 byte (values 0 to 255)
 * \post image values v in [min;max] : v = g / 255 * (max-min) + min
 */
void ASTEX_API  load (ImageGrayd& image, const std::string& filename, double min, double max);

/**
 * \brief save scalar image to file
 * \param [in] img scalar image
 * \param [in] filename
 * \param [in] min min value to be saved (mapped to 0)
 * \param [in] max max value to be saved (mapped to 255)
 * \pre image values v in [min;max], otherwise values are clamped
 * \post the file contains grayscale on 1 byte (values 0 to 255) : g = (v - min) / (max-min) * 255
 */
void ASTEX_API  save (const ImageGrayd& img, const std::string& filename, double min, double max);

/**
 * \brief load scalar image from file
 * \param [out] image scalar image
 * \param [in] filename
 * \param [in] shift_mean_to_zero (default = false) activates shift to the mean value
 * \pre the file contains grayscale on 1 byte (values 0 to 255)
 * \post if (shift_mean_to_zero = false) image values in [0;1]
 * \post if (shift_mean_to_zero = true) image values in [-m;1-m] with mean value = 0
 */
double ASTEX_API load (ImageGrayd& image, const std::string& filename, bool shift_mean_to_zero = false);

/**
 * \brief save scalar image to file
 * \param [in] img scalar image
 * \param [in] filename
 * \param [in] mean_shift (default = 0.0) shifts the mean value
 * \pre values (image + mean_shift) are in [0;1], otherwise values are clamped
 * \post the file contains grayscale on 1 byte (values 0 to 255)
 */
void ASTEX_API  save (const ImageGrayd& img, const std::string& filename, double mean_shift = 0.0);

/**
 * \brief load spectrum (encoded as scalar image) from file (encoded as gray level png)
 * \param [out] image scalar image of size T*T
 * \param [in] filename
 * \pre the file contains grayscale on 1 byte (values g from 0 to 255)
 * \post image values v in [0;T] : v = (g/255)^2 * T
 * \post image is centered on (T/2,T/2), i.e. pixelRelative(0,0) = pixelAbsolute(T/2,T/2)
 */
void ASTEX_API  load_spectrum (ImageSpectrald& image, const std::string& filename);

/**
 * \brief save spectrum (encoded as scalar image) to file (encoded as gray level png)
 * \param [in] img scalar image of size T*T
 * \param [in] filename
 * \pre img has values v in [0;T]
 * \post the file contains grayscale on 1 byte (values 0 to 255) : g = sqrt (v/T) * 255
 */
void ASTEX_API  save_spectrum (const ImageSpectrald& img, const std::string& filename);

/**
 * \brief load phase (encoded as scalar image) from file (encoded as gray level png)
 * \param [out] image scalar image of size T*T
 * \param [in] filename
 * \pre the file contains grayscale on 1 byte (values g from 0 to 255)
 * \post image values v in [-Pi;Pi] : v = (g/255 - 0.5) * 2 * Pi
 * \post image is centered on (T/2,T/2), i.e. pixelRelative(0,0) = pixelAbsolute(T/2,T/2)
 */
void ASTEX_API  load_phase (ImageSpectrald& image, const std::string& filename);

/**
 * \brief save phase (encoded as scalar image) to file (encoded as gray level png)
 * \param [in] img scalar image
 * \param [in] filename
 * \pre img has values v in [-Pi;Pi]
 * \post the file contains grayscale on 1 byte (values 0 to 255) : g = (v + Pi) * 255 / (2*Pi)
 */
void ASTEX_API  save_phase (const ImageSpectrald& img, const std::string& filename);


/**
 * \brief load coordinates (encoded as scalar image) from file (encoded as gray level png)
 * \param [out] image scalar image
 * \param [in] filename
 * \pre the file contains grayscale on 1 byte (values g from 0 to 255)
 * \post image values v in [-1;2] : v = g/255 * 3 -1
 */
void ASTEX_API  load_coordinates (ImageGrayd& image, const std::string& filename);

/**
 * \brief save coordinates (encoded as scalar image) to file (encoded as gray level png)
 * \param [in] img scalar image
 * \param [in] filename
 * \pre img has values v in [-1;2]
 * \post the file contains grayscale on 1 byte (values 0 to 255) : g = (v + 1) * 255 / 3
 */
void ASTEX_API  save_coordinates (const ImageGrayd& img, const std::string& filename);

/**
 * \brief load PCA coordinates (encoded as scalar image) from file (encoded as gray level png)
 * \param [out] image scalar image
 * \param [in] filename
 * \pre the file contains grayscale on 1 byte (values g from 0 to 255)
 * \post image values v in [-10/7;10/7]
 */
void ASTEX_API  load_PCA_coordinates (ImageGrayd& image, const std::string& filename);

/**
 * \brief save PCA coordinates (encoded as scalar image) to file (encoded as gray level png)
 * \param [in] img scalar image
 * \param [in] filename
 * \pre img has values v in [-10/7;10/7]
 * \post the file contains grayscale on 1 byte (values 0 to 255)
 */
void ASTEX_API  save_PCA_coordinates (const ImageGrayd& img, const std::string& filename);




/**
 * @brief load an RGB image into gray of double
 * @param output
 * @param filename
 */
void ASTEX_API  load_RGB_2_gray(ImageGrayd& output, const std::string& filename);

/**
 * @brief load an RGB image into luminance image encoded in gray of double
 * @param output
 * @param filename
 */
void ASTEX_API  load_RGB_2_luminance(ImageGrayd& output, const std::string& filename);

/**
 * @brief load an RGB image into lightness image encoded in gray of double
 * @param output
 * @param filename
 */
void ASTEX_API  load_RGB_2_lightness(ImageGrayd& output, const std::string& filename);

/**
 * @brief save RGBd image
 * @param input input image
 * @param name filename
 * @param n ???
 * @param mean_shift_0 mean_shift for red
 * @param mean_shift_1 mean_shift for green
 * @param mean_shift_2 mean_shift for blue
 */
void ASTEX_API  save(const ImageRGBd& input, const std::string&  name, double n,  double mean_shift_0,  double mean_shift_1,  double mean_shift_2);


} // namespace IO

} // namespace ASTex
#endif
