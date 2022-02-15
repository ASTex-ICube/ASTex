
#ifndef __STAMP__H__
#define __STAMP__H__

#include "vector"

#include <vector>
#include <iostream>
#include <fstream>
#include <memory.h>

#include <cmath>
#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>
#include <ASTex/fourier.h>
#include <ASTex/local_spectrum.h>
#include <ASTex/utils.h>
#include <ASTex/distances_maps.h>
#include <Eigen/Core>
#include <ASTex/mask.h>

namespace ASTex
{

namespace Stamping
{

//Declaration of StampBase
//TODO: for every class, constrain I to be an ImageCommon
template<typename I>
/**
 * @brief The StampBase class is an interface for every Stamps.
 * The main idea is to provide an interface for images and functions
 * and allow the user to access floating point x and y coordinates through pixel(x, y).
 */
class StampBase
{
public:
	/**
	 * @brief StampBase constructor for StampBase.
	 */
	StampBase() {}

	/**
	 * @brief StampBase copy constructor for StampBase.
	 * @param other object to be copied.
	 */
	StampBase(const StampBase &other) : StampBase() {}

	typedef typename I::PixelType PixelType;

	/**
	 * @brief pixel access function to the underlying functional or image.
	 * @param x width coordinate.
	 * @param y height coordinate.
	 * @return the value of the pixel p(x, y).
	 */
	virtual PixelType   pixel(double x, double y) const = 0;

	/**
	 * @brief size
	 * @return the size of the underlying functional or image.
	 */
	virtual Size        size() const = 0;

	/**
	 * @brief width
	 * @return the width of the underlying functional or image.
	 */
	virtual int         width() const = 0;

	/**
	 * @brief height
	 * @return the height of the underlying functional or image.
	 */
	virtual int         height() const = 0;

private:
};

//Declaration of StampDiscrete
//TODO: add mask?
template<typename I>
/**
 * @brief The StampDiscrete class is a class which contains a discrete image.
 */
class StampDiscrete : public StampBase<I>
{
public:

	/**
	 * @brief StampDiscrete copy constructor for StampDiscrete.
	 * @param other the object to be copied.
	 */
	StampDiscrete(const StampDiscrete &other);

	/**
	 * @brief StampDiscrete special constructor for StampDiscrete.
	 * @param image an image to be used by default. It is copied into the class.
	 */
	StampDiscrete(const I &image);

	//nested types

	typedef typename StampBase<I>::PixelType PixelType;
	/**
	 * @brief The Dimensions class is a small nested class which contains the
	 * ratio between pixel size and resolution.
	 * =1.0 means they have the same dimensions (say 1 pixel = 1 meter),
	 * >1.0 means the pixels are bigger (for 2, 1 pixel = 2 meters),
	 * <1.0 means the pixels are smaller (for 0.5, 2 pixels = 1 meter).
	 */
	class Dimensions
	{
	public:
		Dimensions(double x=1.0, double y=1.0): dimX(x), dimY(y) {}
		double dimX;
		double dimY;
	};
	/**
	  * @enum interpolation_rule_t is a type specifying the rule for interpolation,
	  * @value BILINEAR linear interpolation (in two directions).
	  * @value NEAREST chooses nearest pixel.
	  * @value TRUNCATE casts the coordinates to int.
	  */
	typedef enum{BILINEAR, NEAREST, TRUNCATE, BILINEAR_PERIODIC} interpolation_rule_t;

	//get

	/**
	 * @return the size of the underlying image.
	 */
	Size        size() const                        {return m_image.size();}
	/**
	 * @return the width of the underlying image.
	 */
	int         width() const                       {return m_image.width();}
	/**
	 * @return the height of the underlying image.
	 */
	int         height() const                      {return m_image.height();}

	/**
	 * @return the dimensions stored by the class.
	 * they contain the ratio between pixel size and resolution.
	 * =1.0 means they have the same dimensions (say 1 pixel = 1 meter),
	 * >1.0 means the pixels are bigger (for 2, 1 pixel = 2 meters),
	 * <1.0 means the pixels are smaller (for 0.5, 2 pixels = 1 meter).
	 * 1.0 for both dimensions by default (see dimX() and dimY())
	 */
	Dimensions  dimensions() const                  {return m_dimensions;}
	/**
	 * @return the x dimension stored by the class.
	 */
	double dimX() const                             {return m_dimensions.dimX;}
	/**
	 * @return the y dimension stored by the class.
	 */
	double dimY() const                             {return m_dimensions.dimY;}
	/**
	 * @return the interpolation rule stored by the class.
	 * See @enum interpolation_rule_t for details.
	 */
	interpolation_rule_t interpolationRule() const  {return m_rule;}

	/**
	 * @brief image
	 * @return the underlying image.
	 */
	const I& image() const                          {return m_image;}

	/**
	 * @brief pixel access function of the underlying image.
	 * Utilizes the interpolation rule to interpolate if x or y are not integers.
	 * @param x width coordinate.
	 * @param y height coordinate.
	 * @return the value of the pixel p(x, y).
	 */
	virtual PixelType   pixel(double x, double y) const;

	//set

	/**
	 * @param rule the interpolation rule stored by the class.
	 * See @enum interpolation_rule_t for details.
	 */
	void setInterpolationRule(interpolation_rule_t rule);
	/**
	 * @param dimensions the dimensions stored by the class.
	 * they contain the ratio between pixel size and resolution.
	 * =1.0 means they have the same dimensions (say 1 pixel = 1 meter),
	 * >1.0 means the pixels are bigger (for 2, 1 pixel = 2 meters),
	 * <1.0 means the pixels are smaller (for 0.5, 2 pixels = 1 meter).
	 * 1.0 for both dimensions by default (see dimX() and dimY())
	 * @pre dimensions.dimX > 0 and dimensions.dimY > 0.
	 */
	void setDimensions(Dimensions dimensions);

	/**
	 * @param dimX the x dimension stored by the class.
	 * @pre dimX>0
	 */
	void setDimX(double dimX);
	/**
	 * @param dimY the y dimension stored by the class.
	 * @pre dimY>0
	 */
	void setDimY(double dimY);

	/**
	 * @param image the underlying image to be set.
	 * @pre image is initialized.
	 */
	void setImage(const I& image); //TODO: might be sensitive, check for itk documentation

protected:

	virtual bool is_in_range(int x, int y) const;

private:

	/**
	 * @brief bilinear_interpolation_ computes a bilinear interpolation
	 * utilizing the values of the four nearest pixel.
	 * @param q11 top left pixel value.
	 * @param q12 bottom left pixel value.
	 * @param q21 top right pixel value.
	 * @param q22 bottom right pixel value.
	 * @param x1 leftmost coordinate.
	 * @param y1 bottom coordinate.
	 * @param x width coordinate.
	 * @param y height coordinate.
	 * @return the interpolated pixel value.
	 */
	inline typename I::PixelType bilinear_interpolation_(typename I::PixelType q11,
														 typename I::PixelType q12,
														 typename I::PixelType q21,
														 typename I::PixelType q22,
														 int x1, int y1,
														 double x, double y) const;

	const I&                m_image;
	interpolation_rule_t    m_rule;
	Dimensions              m_dimensions;

	/**
	 * @brief ms_zero a "hack" variable to have a pixel set to zero by default.
	 * It could be removed if we ensure that
	 * every pixel value is set to zero with its default constructor.
	 */
	static PixelType ms_zero;
};

//Declaration of StampDiscreteWithMask

template<typename I>
class StampDiscreteWithMask : public StampDiscrete<I>
{
public:
	StampDiscreteWithMask(const StampDiscreteWithMask<I> &other);
	StampDiscreteWithMask(const I& image);

	void setMask(const MaskBool &mask);

	const MaskBool& mask() const {return m_mask;}
	MaskBool &mask() {return m_mask;}

protected:

	bool is_in_range(int x, int y) const;

private:

	MaskBool m_mask;
};

//Declaration of StampContinuous?

//
//Implementation of StampDiscrete
//

template<typename I>
StampDiscrete<I>::StampDiscrete(const StampDiscrete &other) :
	m_image(other.image()),
	m_rule(other.interpolationRule()),
	m_dimensions(other.dimensions())
{}

template<typename I>
typename StampDiscrete<I>::PixelType StampDiscrete<I>::ms_zero;

template<typename I>
StampDiscrete<I>::StampDiscrete(const I &image) :
	StampBase<I>::StampBase(),
	m_image(image),
	m_rule(BILINEAR),
	m_dimensions(Dimensions(1.0, 1.0))
{}

template<typename I>
typename StampDiscrete<I>::PixelType StampDiscrete<I>::pixel(double x, double y) const
{
	int dx, dy, tx, ty;
	PixelType pixel = ms_zero;
	double dimX, dimY;
	typename I::PixelType q11, q12, q21, q22;

	dimX = x / m_dimensions.dimX;
	dimY = y / m_dimensions.dimY; //set x and y to scale;
	tx = (int) dimX;
	ty = (int) dimY;

	if(m_rule == BILINEAR)
	{
		q11 = is_in_range(tx, ty) ? m_image.pixelAbsolute(tx, ty) : ms_zero;
		q12 = is_in_range(tx, ty+1) ? m_image.pixelAbsolute(tx, ty+1) : ms_zero;
		q21 = is_in_range(tx+1, ty) ? m_image.pixelAbsolute(tx+1, ty) : ms_zero;
		q22 = is_in_range(tx+1, ty+1) ? m_image.pixelAbsolute(tx+1, ty+1) : ms_zero;

		pixel = bilinear_interpolation_(q11, q12, q21, q22, tx, ty, dimX, dimY);
	}
	else if(m_rule == BILINEAR_PERIODIC)
	{
		pixel = bilinear_interpolation(m_image, dimX, dimY, true);
	}
	else if(m_rule == NEAREST) //TODO: that's not actually a real nearest is it
	{
		dx = dimX - tx;
		dy = dimY - ty;
		if(dx >= 0.5)
			++tx;
		if(dy >= 0.5)
			++ty;

		pixel = is_in_range(tx, ty) ? m_image.pixelAbsolute(tx, ty) : ms_zero;
	}
	else //m_rule == SD_TRUNC
		pixel = is_in_range(tx, ty) ? m_image.pixelAbsolute(tx, ty) : ms_zero;

	return pixel;
}

template<typename I>
bool StampDiscrete<I>::is_in_range(int x, int y) const
{
	return x >= 0 && x < m_image.width() && y >= 0 && y < m_image.height();
}


template<typename I>
typename I::PixelType StampDiscrete<I>::bilinear_interpolation_(typename I::PixelType q11,
																typename I::PixelType q12,
																typename I::PixelType q21,
																typename I::PixelType q22,
																int x1, int y1,
																double x, double y) const
{
	double x2x, y2y, yy1, xx1;
	yy1 = y - y1;
	xx1 = x - x1;
	x2x = 1.0 - xx1;
	y2y = 1.0 - yy1;
	return
		q11 * (x2x * y2y) +
		q21 * (xx1 * y2y) +
		q12 * (x2x * yy1) +
		q22 * (xx1 * yy1);
}

//set

template<typename I>
void StampDiscrete<I>::setInterpolationRule(typename StampDiscrete<I>::interpolation_rule_t rule)
{
	m_rule = rule;
}

template<typename I>
void StampDiscrete<I>::setDimensions(typename StampDiscrete<I>::Dimensions dimensions)
{
	assert(dimensions.dimX>0 && dimensions.dimY>0
		   && "StampDiscrete::setDimensions: dimensions must be > 0");
	m_dimensions = dimensions;
}

template<typename I>
void StampDiscrete<I>::setDimX(double dimX)
{
	assert(dimX>0
		   && "setDimX: dimX must be > 0");
	m_dimensions.dimX = dimX;
}

template<typename I>
void StampDiscrete<I>::setDimY(double dimY)
{
	assert(dimY>0
		   && "setDimY: dimY must be > 0");
	m_dimensions.dimY = dimY;
}

template<typename I>
void StampDiscrete<I>::setImage(const I& image)
{
	assert(image.isInitialized());
	m_image = image;
}

//
//Implementation of StampDiscreteWithMask
//

template<typename I>
StampDiscreteWithMask<I>::StampDiscreteWithMask(const StampDiscreteWithMask<I> &other) :
	StampDiscrete<I>(other),
	m_mask(other.mask())
{}

template<typename I>
StampDiscreteWithMask<I>::StampDiscreteWithMask(const I& image) :
	StampDiscrete<I>(image),
	m_mask(image.width(), image.height(), true)
{}

template<typename I>
void StampDiscreteWithMask<I>::setMask(const MaskBool &mask)
{
	assert(mask.width() == this->width() && mask.height() == this->height() &&
		   "StampDiscreteWithMask::setMask(mask): mask must have the same dimensions as the underlying image");
	m_mask=mask;
}

template<typename I>
bool StampDiscreteWithMask<I>::is_in_range(int x, int y) const
{
	return StampDiscrete<I>::is_in_range(x, y) && m_mask(x, y);
}

} //namespace Stamping

} //namespace ASTex


#endif //__STAMP__H__
