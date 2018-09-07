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



#ifndef __ASTEX_FOURIER__
#define __ASTEX_FOURIER__

#include <ASTex/image_gray.h>
#include <ASTex/image_spectral.h>
#include <ASTex/image_rgb.h>

/** \file
 * \brief Fourier transforms and functions in frequency domain
*/



/**
 * \brief namespace ASTex
 */
namespace ASTex
{

namespace Fourier
{

/**
 * \brief routine for circular shifting in [-PI; PI]
 */
void ASTEX_API shiftDouble (double& d, double s);

/**
 * \brief inplace random shifting of the phases; random is uniform in [-PI; PI]
 * \param [in,out] phase image of phases
 * \param [in] mask the shifting applies only where the mask is true (default : true everywhere)
 * \pre phases in [-PI; PI]
 * \post phases in [-PI; PI]
 * \pre phase has central antisymmetry
 * \post phase has central antisymmetry
 */
template<typename MASK>
void ASTEX_API randomPhaseShift(ImageSpectrald& phase, const MASK& mask);

inline void ASTEX_API randomPhaseShift(ImageSpectrald& phase)
{
	randomPhaseShift(phase,[](int,int){return true;});
}

/**
 * \brief sets random values to phase; random is uniform in [-PI; PI]
 * \param [in,out] phase image of phases
 * \param [in] mask the shifting applies only where the mask is true (default : true everywhere)
 * \param [in] call_srand the pseudo-random generator is initialized (call of srand) iff true (default : true everywhere)
 * \post phases in [-PI; PI]
 * \post phase has central antisymmetry
 */
template<typename MASK>
void randomPhase(ImageSpectrald& phase, const MASK& mask, bool call_srand = true);

inline void ASTEX_API randomPhase(ImageSpectrald& phase)
{
	randomPhase(phase,[](int,int){return true;},true);
}


/**
 * \brief phase shift so as to create a space shift [dx,dy]
 * \param [in,out] phase image of phases
 * \param [in] dx space shift on X axis
 * \param [in] dy space shift on Y axis
 * \pre la fréquence nulle est en [0,0], i.e. FFTShiftImageFilter n'a pas été appelé
 * \pre phase->GetLargestPossibleRegion().GetIndex() vaut [0,0]
 * \warning depending on the shift value (e.g. odd or even), this is not strictly equivalent to the space shift using itk::CyclicShiftImageFilter, I don't know why
 */

void ASTEX_API phaseShiftFromSpaceShift(ImageSpectrald& phase, double dx, double dy);

/**
 * \brief DFT (Discrete Fourier Transform)
 * \param [in] input grayscale input image
 * \param [out] modulus DFT modulus
 * \param [out] phase DFT phase
 * \param [in] preserve energy : preserve the global energy (sqrt of #pixels factor) if true (default : true)
 */
void ASTEX_API fftForwardModulusAndPhase(const ImageGrayd& input, ImageSpectrald& modulus, ImageSpectrald& phase, bool preserve_energy = true );

/**
 * \brief inverse DFT (Discrete Fourier Transform)
 * \param [in] modulus DFT modulus
 * \param [in] phase DFT phase
 * \param [out] output grayscale output image
 * \param [in] preserve energy : preserve the global energy (sqrt of #pixels factor) if true (default : true)
 */
void ASTEX_API fftInverseModulusAndPhase(const ImageSpectrald& modulus, const ImageSpectrald& phase, ImageGrayd& output, bool preserve_energy = true);

/**
 * \brief compute square-root of the PSD (power spectrum density) by Welch's method
 * \param [in] input grayscale input image
 * \param [out] modulus sqrt of the PSD (size = larger power of two < image size)
 * \param [in] step shift step (>= 1 pixel) of the window. Default = 10 (saves 99% time w.r.t. step = 1).
 */
void ASTEX_API welch(const ImageGrayd& input, ImageSpectrald& modulus, uint32_t step = 10 );

void ASTEX_API welch(const ImageGrayd& input, ImageSpectrald& modulus, uint32_t step, int32_t sp_size);

/**
 * \brief compute crossCorrelation of the full image
 * \param [in] input RGB input image
 * \param [out] acorr cross-correlation
 * \param [in] mask pairs of pixels are considered iff mask.test(pix) == true for both pixels of the pair.
 * \pre acorr must be allocated
 * \pre input and acorr are images of the same size
 * \return number of pixels used (mask size)
 */
template<typename MASK>
int crossCorrelation_full_size(const ImageRGBd& input, ImageGrayd& acorr, const MASK& mask, int channel1, int channel2);

/**
 * \brief compute autoCorrelation of the full image
 * \param [in] input grayscale input image
 * \param [out] acorr auto-correlation
 * \param [in] mask pairs of pixels are considered iff mask.test(pix) == true for both pixels of the pair.
 * \pre acorr must be allocated
 * \pre input and acorr are square images of the same size
 * \return number of pixels used (mask size)
 */
template<typename MASK>
int autoCorrelation_full_size(const ImageGrayd& input, ImageGrayd& acorr, const MASK& mask);

inline int ASTEX_API autoCorrelation_full_size(const ImageGrayd& input, ImageGrayd& acorr)
{
	return autoCorrelation_full_size(input,acorr,[](int,int){return true;});
}

/**
 * \brief compute autoCorrelation of a block begining at (bloc_x,bloc_y)
 * \param [in] input grayscale input image
 * \param [out] acorr auto-correlation
 * \param [in] block_x begin abscissa of the block
 * \param [in] block_y begin ordinate of the block
 * \param [in] mask pairs of pixels are considered iff mask.test(pix) == true for both pixels of the pair.
 * \pre acorr must be allocated
 * \pre acorr is a square TxT image
 * \pre block_x+T <= input.width() and block_y+T <= input.height()
 * \return number of pixels used (mask size)
 */
template<typename MASK>
int autoCorrelation_block(const ImageGrayd& input, ImageGrayd& acorr, int block_x, int block_y, const MASK& mask);

inline int ASTEX_API autoCorrelation_block(const ImageGrayd& input, ImageGrayd& acorr, int block_x, int block_y)
{
	return autoCorrelation_block(input,acorr,block_x,block_y,[](int,int){return true;});
}


template<typename MASK>
int autoCorrelation_block_3chan(const ImageGrayd& input_1,const ImageGrayd& input_2,const ImageGrayd& input_3, ImageGrayd& acorr_1, ImageGrayd& acorr_2, ImageGrayd& acorr_3, int block_x, int block_y, const MASK& mask);

inline int autoCorrelation_block_3chan(const ImageGrayd& input_1,const ImageGrayd& input_2,const ImageGrayd& input_3, ImageGrayd& acorr_1, ImageGrayd& acorr_2, ImageGrayd& acorr_3, int block_x, int block_y)
{
	return autoCorrelation_block_3chan(input_1,input_2,input_3,acorr_1,acorr_2,acorr_3,block_x,block_y,[](int,int){return true;});
}

void ASTEX_API autoCorrelation_normalize(ImageGrayd& acorr, double v);

void ASTEX_API autoCorrelation_normalize_3chan(ImageGrayd& acorr_1,ImageGrayd& acorr_2,ImageGrayd& acorr_3, double v);



/**
 * \brief compute square-root of the PSD (power spectrum density) by fft of auto-correlation
 * \param [in] input grayscale input image
 * \param [out] modulus sqrt of the PSD
 * \param [in] mask select the pixels that are considered for computation
 * \pre input and modulus are square images of the same size
 */
template<typename MASK>
void spectrum_by_autocorrelation_full_size(const ImageGrayd& input, ImageSpectrald& modulus, const MASK& mask);

inline void ASTEX_API spectrum_by_autocorrelation_full_size(const ImageGrayd& input, ImageSpectrald& modulus)
{
	return spectrum_by_autocorrelation_full_size(input,modulus,[](int,int){return true;});
}

/**
 * \brief compute square-root of the PSD (power spectrum density) by fft of auto-correlation
 * \param [in] input grayscale input image
 * \param [out] modulus sqrt of the PSD
 * \param [in] mask select the pixels that are considered for computation
 * \param [in] proportion_threshold each block is used iff the proportion of pixels in the block AND in the mask is above the threshold
 * \pre modulus is a square image, smaller than input
 * \return the number of blocks used in the computation
 */
template<typename MASK>
int spectrum_by_autocorrelation_small_size(const ImageGrayd& input, ImageSpectrald& modulus, const MASK& mask, double proportion_threshold, int step);

inline int ASTEX_API spectrum_by_autocorrelation_small_size(const ImageGrayd& input, ImageSpectrald& modulus, double proportion_threshold, int step)
{
	return spectrum_by_autocorrelation_small_size(input,modulus,[](int,int){return true;},proportion_threshold,step);
}

/**
 * \brief power of a scalar image
 * \param [in] input input scalar image
 * \param [in] mask select the pixels that are considered for computation
 * \return power
 */
template<typename MASK>
double getPower(const ImageGrayd& input, const MASK& mask);

inline double ASTEX_API getPower(const ImageGrayd& input)
{
	return getPower(input,[](int,int){return true;});
}


/**
 * \brief scales an image so as to reach a target power
 * \param [in,out] image scalar image
 * \param [in] power the target power
 */
void ASTEX_API setPower(ImageGrayd& image, double power);


/**
 * \brief spectrum downsampling, factor 2x2. Work in progress (Basile).
 */
void ASTEX_API down_sampling(const ImageSpectrald& modulus, ImageSpectrald& output);

/**
 * \brief work in progress (Basile).
 */
void ASTEX_API phaseOfRegion(const ImageGrayd& input, const ImageGrayd::ItkImg::RegionType& region, ImageSpectrald& phase );

/**
 * \brief work in progress (Basile) : average phase from shifted (step) blocks (size = power of 2 lower than input size) with or without countering space shift (activate).
 */
void ASTEX_API phaseAverage_UsingNormalizedComplexAverage(const ImageGrayd& input, ImageSpectrald& phase, bool activate = true, uint32_t step = 10 );

/**
 * \brief work in progress (Basile) : average phase from shifted (step) blocks (size = power of 2 lower than input size) with or without countering space shift (activate).
 */
void ASTEX_API phaseAverage_Naive(const ImageGrayd& input, ImageSpectrald& phase, bool activate = true, uint32_t step = 10 );

/**
 * \brief work in progress (Basile), some problems due to bug in phase shift : average phase from shifted (step) blocks (size = power of 2 lower than input size) with or without countering space shift (activate).
 */
void ASTEX_API phaseAverage_UsingPhaseShift(const ImageGrayd& input, ImageSpectrald& phase, bool activate = true, uint32_t step = 10 );


/**
 * @brief RPnoise_mosaic random phase noise
 * @param [in] modulus a spectrum modulus (encoded as grayscale image)
 * @param [out] result is a mosaic of local RP-noise with linear blending (bands of width = blending_size)
 * @param [in] blending_size is the width of bands for blending
 * @param [in] call_srand the pseudo-random generator is initialized (call of srand) iff true (default : true everywhere)
 */
void ASTEX_API RPnoise_mosaic(const ImageSpectrald& modulus, ImageGrayd& result, int blending_size, bool call_srand = true);


/**
 * @brief RPnoise_mosaic random phase noise
 * @param [in] modulus a spectrum modulus (encoded as grayscale image)
 * @param [out] result is a mosaic of local RP-noise with linear blending (bands of width = blending_size)
 */
void ASTEX_API RPnoise_mosaic_periodique(const ImageSpectrald& modulus, ImageGrayd& result);

/**
 * @brief RPnoise_mosaic random phase noise
 * @param [in] modulus a spectrum modulus (encoded as grayscale image)
 * @param [out] result is a mosaic of local RP-noise with linear blending (bands of width = blending_size)
 */
void ASTEX_API RPnoise_mosaic_bandes(const ImageSpectrald& modulus, ImageGrayd& result);

/**
 * @brief distance_spectrum_to_spectrum_linear_weights computes the distance between 2 spectrums (frequencies are wheighted linearly)
 * @param [in] sp1 spectrum 1
 * @param [in] sp2 spectrum 2
 * \pre spectrums are shifted and centered : pixelRelative(fx,fy) corresponds to frequency (fx,fy) in [-T/2 ; T/2-1] ^ 2
 * @return the distance from sp1 to sp2
 */
double ASTEX_API distance_spectrum_to_spectrum_linear_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2);

/**
 * @brief distance_spectrum_to_spectrum_uniform_weights computes the distance between 2 spectrums (frequencies are wheighted uniformly)
 * @param [in] sp1 spectrum 1
 * @param [in] sp2 spectrum 2
 * \pre spectrums are shifted and centered : pixelRelative(fx,fy) corresponds to frequency (fx,fy) in [-T/2 ; T/2-1] ^ 2
 * @return the distance from sp1 to sp2
 */
double ASTEX_API distance_spectrum_to_spectrum_uniform_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2);

/**
 * @brief distance_spectrum_to_spectrum_linear_weights computes the distance squared between 2 spectrums (frequencies are wheighted linearly)
 * @param [in] sp1 spectrum 1
 * @param [in] sp2 spectrum 2
 * \pre spectrums are shifted and centered : pixelRelative(fx,fy) corresponds to frequency (fx,fy) in [-T/2 ; T/2-1] ^ 2
 * @return the distance squared from sp1 to sp2
 */
double ASTEX_API distance_squared_spectrum_to_spectrum_linear_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2);

/**
 * @brief distance_spectrum_to_spectrum_uniform_weights computes the distance squared between 2 spectrums (frequencies are wheighted uniformly)
 * @param [in] sp1 spectrum 1
 * @param [in] sp2 spectrum 2
 * \pre spectrums are shifted and centered : pixelRelative(fx,fy) corresponds to frequency (fx,fy) in [-T/2 ; T/2-1] ^ 2
 * @return the distance squared from sp1 to sp2
 */
double ASTEX_API distance_squared_spectrum_to_spectrum_uniform_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2);

/**
 * @brief dot_product_spectrum_linear_weight computes the dot product between 2 spectrums (frequencies are wheighted linearly)
 * @param [in] sp1 spectrum 1
 * @param [in] sp2 spectrum 2
 * \pre spectrums are shifted and centered : pixelRelative(fx,fy) corresponds to frequency (fx,fy) in [-T/2 ; T/2-1] ^ 2
 * @return the dot-product between sp1 and sp2
 */
double ASTEX_API dot_product_spectrum_linear_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2);

/**
 * @brief dot_product_spectrum_uniform_weight computes the dot product between 2 spectrums (frequencies are wheighted uniformly)
 * @param [in] sp1 spectrum 1
 * @param [in] sp2 spectrum 2
 * \pre spectrums are shifted and centered : pixelRelative(fx,fy) corresponds to frequency (fx,fy) in [-T/2 ; T/2-1] ^ 2
 * @return the dot-product between sp1 and sp2
 */
double ASTEX_API dot_product_spectrum_uniform_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2);

/**
 * \class spectrum_projector_2D fourier.h
 * \brief projects a spectrum onto a subspace of dimension 2
 */
class ASTEX_API spectrum_projector_2D
{
public :
	/**
	 * \brief add a spectrum in the basis
	 * \param [in] sp spectrum
	 */
	void add (const ImageSpectrald& sp);
	/**
	 * \brief pre-compute the left-hand side of the linear system
	 * \pre all spectra in the basis have been added
	 * \post frequencies are wheighted linearly
	 */
	void precompute_system ();
	/**
	 * \brief computes the coordinates of the projection of sp in the sub-space spanned by the basis
	 * \param [in] sp spectrum
	 * \param [out] coord the coordinates
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post frequencies are wheighted linearly
	 */
	void project (const ImageSpectrald& sp, std::vector<double>& coord);

	/**
	 * \brief pre-compute the left-hand side of the linear system
	 * \pre all spectra in the basis have been added
	 * \post frequencies are wheighted uniformly
	 */
	void precompute_system_uniform_weights ();
	/**
	 * \brief computes the coordinates of the projection of sp in the sub-space spanned by the basis
	 * \param [in] sp spectrum
	 * \param [out] coord the coordinates
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post frequencies are wheighted uniformly
	 */
	void project_uniform_weights (const ImageSpectrald& sp, std::vector<double>& coord);

	/**
	 * \brief projection in "the half-cone of positive coordinates" from the signed coordinates in the the basis (sp0,sp1)
	 * \param [in,out] coord the coordinates
	 * \pre coord contains the coordinates in sp0,sp1 (may be negative)
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post the coordinates are all positive
	 */
	void reproject_coordinates_to_positive (std::vector<double>& coord);

	/**
	 * \brief projection to [0,1] from the positive coordinates in the the basis (sp0,sp1)
	 * \param [in,out] coord the coordinates
	 * \pre coord contains the coordinates in sp0,sp1
	 * \pre the coordinates are all positive (may be above 1)
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post the coordinates are all in [0,1]
	 */
	void reproject_coordinates_to_0_1 (std::vector<double>& coord);

	/**
	 * \brief computes the coordinates of the projection of sp in "the half-cone of positive coordinates" in the sub-space spanned by the basis
	 * \param [in] sp spectrum
	 * \param [out] coord the coordinates
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post the coordinates are all positive
	 */
	void project_with_positive_coordinates (const ImageSpectrald& sp, std::vector<double>& coord);

	/**
	 * \brief dot product expressed in the sub-space spanned by the basis
	 * \param [in] c1 coordinates of the spectrum 1
	 * \param [in] c2 coordinates of the spectrum 2
	 * \param [out] scalar product
	 * \pre the system has been pre-computed, see precompute_system()
	 */
	double dot_product (std::vector<double>& c1 , std::vector<double>& c2);

	/**
	 * \brief distance between spectra expressed in the sub-space spanned by the basis
	 * \param [in] c1 coordinates of the spectrum 1
	 * \param [in] c2 coordinates of the spectrum 2
	 * \param [out] distance between the spectra
	 * \pre the system has been pre-computed, see precompute_system()
	 */
	double distance_spectrum_to_spectrum (std::vector<double>& c1, std::vector<double>& c2);

private :
	std::vector<ImageSpectrald> m_S;
	double m_A [2][2];
	double m_A_inverse [2][2];

	/**
	 * \brief re-projection on the line spanned by spx
	 * \param [in,out] coord the coordinates
	 * \param [in] x
	 * \pre coord contains the coordinates in sp0,sp1
	 * \post the other coordinates (y=1-x) is set to zero
	 */
	void reproject (std::vector<double>& coord, const uint32_t x);
};

/**
 * \class spectrum_projector_3D fourier.h
 * \brief projects a spectrum onto a subspace of dimension 3
 */
class ASTEX_API spectrum_projector_3D
{
public :
	/**
	 * \brief add a spectrum in the basis
	 * \param [in] sp spectrum
	 */
	void add (const ImageSpectrald& sp);
	/**
	 * \brief pre-compute the left-hand side of the linear system
	 * \pre all spectra in the basis have been added
	 * \post frequencies are wheighted linearly
	 */
	void precompute_system ();
	/**
	 * \brief computes the coordinates of the projection of sp in the sub-space spanned by the basis
	 * \param [in] sp spectrum
	 * \param [out] coord the coordinates
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post frequencies are wheighted linearly
	 */
	void project (const ImageSpectrald& sp, std::vector<double>& coord);

	/**
	 * \brief pre-compute the left-hand side of the linear system
	 * \pre all spectra in the basis have been added
	 * \post frequencies are wheighted uniformly
	 */
	void precompute_system_uniform_weights ();
	/**
	 * \brief computes the coordinates of the projection of sp in the sub-space spanned by the basis
	 * \param [in] sp spectrum
	 * \param [out] coord the coordinates
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post frequencies are wheighted uniformly
	 */
	void project_uniform_weights (const ImageSpectrald& sp, std::vector<double>& coord);


	/**
	 * \brief projection in "the half-cone of positive coordinates" from the signed coordinates in the the basis (sp0,sp1,sp2)
	 * \param [in,out] coord the coordinates
	 * \pre coord contains the coordinates in sp0,sp1,sp2 (may be negative)
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post the coordinates are all positive
	 */
	void reproject_coordinates_to_positive (std::vector<double>& coord);

	/**
	 * \brief projection to [0,1] from the positive coordinates in the the basis (sp0,sp1,sp2)
	 * \param [in,out] coord the coordinates
	 * \pre coord contains the coordinates in sp0,sp1,sp2
	 * \pre the coordinates are all positive (may be above 1)
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post the coordinates are all in [0,1]
	 */
	void reproject_coordinates_to_0_1 (std::vector<double>& coord);

	/**
	 * \brief computes the coordinates of the projection of sp in "the half-cone of positive coordinates" in the sub-space spanned by the basis
	 * \param [in] sp spectrum
	 * \param [out] coord the coordinates
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post the coordinates are all positive
	 */
	void project_with_positive_coordinates (const ImageSpectrald& sp, std::vector<double>& coord);

	/**
	  * \brief dot product expressed in the sub-space spanned by the basis
	  * \param [in] c1 coordinates of the spectrum 1
	  * \param [in] c2 coordinates of the spectrum 2
	  * \param [out] scalar product
	  * \pre the system has been pre-computed, see precompute_system()
	  */
	double dot_product (std::vector<double>& c1 , std::vector<double>& c2);

	/**
	  * \brief distance between spectra expressed in the sub-space spanned by the basis
	  * \param [in] c1 coordinates of the spectrum 1
	  * \param [in] c2 coordinates of the spectrum 2
	  * \param [out] distance between the spectra
	  * \pre the system has been pre-computed, see precompute_system()
	  */
	double distance_spectrum_to_spectrum (std::vector<double>& c1 , std::vector<double>& c2);

	//    void test(double x, double y, double z);

private :
	std::vector<ImageSpectrald> m_S;
	double m_A [3][3];
	double m_A_inverse [3][3];

	/**
	 * \brief re-projection on the plane spanned by spx and spy
	 * \param [in,out] coord the coordinates
	 * \param [in] x
	 * \param [in] y
	 * \pre coord contains the coordinates in sp0,sp1,sp2
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post the third coordinate (z = 3-x-y) is set to zero
	 */
	void reproject (std::vector<double>& coord, const uint32_t x, const uint32_t y);

	/**
	 * \brief re-projection on the line spanned by spx
	 * \param [in,out] coord the coordinates
	 * \param [in] x
	 * \pre coord contains the coordinates in sp0,sp1,sp2
	 * \pre the system has been pre-computed, see precompute_system()
	 * \post the two other coordinates are set to zero
	 */
	void reproject (std::vector<double>& coord, const uint32_t x);
};

}

}

#include <ASTex/fourier.hpp>


#endif
