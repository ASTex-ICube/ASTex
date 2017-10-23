
#ifndef __ASTEX_COLORSPACE_FILTERS__
#define __ASTEX_COLORSPACE_FILTERS__

#include <cmath>


#include "itkUnaryFunctorImageFilter.h"

namespace ASTex
{

namespace ColorSpace
{

// #define  avec undef a lafin du fichier ?
static const double WP_x   = 0.3127;
static const double WP_y   = 0.3290;
static const double WP_Xn  = 0.950456;
static const double WP_Yn  = 1.0;
static const double WP_Zn  = 1.088754;
static const double WP_upn = 0.197832;
static const double WP_vpn = 0.468340;

/**
 * @brief function that transform a RGB real pixel color into XYZ real pixel
 * @param p input color pixel
 * @return corresponding XYZ color
 */
template< typename TInputPixel, typename TOutputPixel >
inline TOutputPixel colorRGBtoXYZ(const TInputPixel& p)
{
	TOutputPixel q;
	q[0] = 0.412453*p[0] + 0.357580*p[1] + 0.180423*p[2];
	q[1] = 0.212671*p[0] + 0.715160*p[1] + 0.072169*p[2];
	q[2] = 0.019334*p[0] + 0.119193*p[1] + 0.950227*p[2];
	return q;
}

/// \brief functor version of function for filter usage
template< class TInput, class TOutput>
class fonctorRGBtoXYZ
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		return colorRGBtoXYZ<TInput,TOutput>(p);
	}
};




/**
 * @brief function that transform a XYZ real pixel color into RGB real pixel
 * @param p input color pixel
 * @return corresponding RGB color
 */
template< typename TInputPixel, typename TOutputPixel >
inline TOutputPixel colorXYZtoRGB(const TInputPixel& p)
{
	TOutputPixel q;
	q[0] =  3.240479*p[0] - 1.537150*p[1] - 0.498535*p[2];
	q[1] = -0.969256*p[0] + 1.875992*p[1] + 0.041556*p[2];
	q[2] =  0.055648*p[0] - 0.204043*p[1] + 1.057311*p[2];
	return q;
}

/// \brief functor version of function for filter usage
template< class TInput, class TOutput>
class fonctorXYZtoRGB
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		return colorXYZtoRGB<TInput,TOutput>(p);
	}
};




/**
 * @brief function that transform a XYZ real pixel color into LUV real pixel
 * @param p input color pixel
 * @return corresponding LUV color
 */
template< typename TInputPixel, typename TOutputPixel >
inline TOutputPixel colorXYZtoLUV(const TInputPixel& p)
{
	TOutputPixel q; // L*,u*,v*

	double YY = p[1]/WP_Yn;

	// compute L*
	if ( YY > 0.008856)
		q[0] = 116.0 * std::pow(YY,1.0/3.0) - 16.0;
	else
		q[0] = 903.3 * YY;

	double sum = (p[0] + 15.0*p[1] + 3.0*p[2]);
	double up = 4.0*p[0] / sum;
	double vp = 9.0*p[1] / sum;

	// compute u*
	q[1] = 13.0*q[0]*(up-WP_upn);

	// compute v*
	q[2] = 13.0*q[0]*(vp-WP_vpn);

	return q;
}



/// \brief functor version of function for filter usage
template< class TInput, class TOutput>
class fonctorXYZtoLUV
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		return colorXYZtoLUV<TInput,TOutput>(p);
	}
};


#define min_f(a, b, c)  (fminf(a, fminf(b, c)))
#define max_f(a, b, c)  (fmaxf(a, fmaxf(b, c)))
/**
 * @brief function that transform a RGB real pixel color into HSV real pixel
 * @param p input color pixel
 * @return corresponding HSV color
 */
template< typename TInputPixel, typename TOutputPixel >
inline TOutputPixel colorRGBtoHSV(const TInputPixel& p)
{
	TOutputPixel q;

	float r = p[0] / 255.0f;
    float g = p[1] / 255.0f;
    float b = p[2] / 255.0f;

    float h, s, v;

    float max = max_f(r, g, b);
    float min = min_f(r, g, b);

    v = max;

    if (max == 0.0f) {
        s = 0;
        h = 0;
    }
    else if (max - min == 0.0f) {
        s = 0;
        h = 0;
    }
    else {
        s = (max - min) / max;

        if (max == r) {
            h = 60 * ((g - b) / (max - min)) + 0;
        }
        else if (max == g) {
            h = 60 * ((b - r) / (max - min)) + 120;
        }
        else {
            h = 60 * ((r - g) / (max - min)) + 240;
        }
    }

    if (h < 0) h += 360.0f;

    q[0] = (uint8_t)(h / 2);   // 
    q[1] = (uint8_t)(s * 255); // 
    q[2] = (uint8_t)(v * 255); // 

    return q;


}

#undef min_f
#undef max_f


/// \brief functor version of function for filter usage
template< class TInput, class TOutput>
class fonctorRGBtoHSV
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		return colorRGBtoHSV<TInput,TOutput>(p);
	}
};




/**
 * @brief function that transform a LUV real pixel color into XYZ real pixel
 * @param p input color pixel
 * @return corresponding XYZ color
 */
template< typename TInputPixel, typename TOutputPixel >
inline TOutputPixel colorHSVtoRGB(const TInputPixel& p)
{
	TOutputPixel ret; // X,Y,Z


    float h = p[0] *   2.0f; // 0-360
    float s = p[1] / 255.0f; // 0.0-1.0
    float v = p[2] / 255.0f; // 0.0-1.0

    float r, g, b; // 0.0-1.0

    int   hi = (int)(h / 60.0f) % 6;
    float f  = (h / 60.0f) - hi;
    float p2  = v * (1.0f - s);
    float q  = v * (1.0f - s * f);
    float t  = v * (1.0f - s * (1.0f - f));

    switch(hi) {
        case 0: r = v, g = t, b = p2; break;
        case 1: r = q, g = v, b = p2; break;
        case 2: r = p2, g = v, b = t; break;
        case 3: r = p2, g = q, b = v; break;
        case 4: r = t, g = p2, b = v; break;
        case 5: r = v, g = p2, b = q; break;
    }

    ret[0] = (uint8_t)(r * 255); // dst_r : 0-255
    ret[1] = (uint8_t)(g * 255); // dst_r : 0-255
    ret[2] = (uint8_t)(b * 255); // dst_r : 0-255

    return ret;

}

/// \brief functor version of function for filter usage
template< class TInput, class TOutput>
class fonctorHSVtoRGB
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		return colorHSVtoRGB<TInput,TOutput>(p);
	}
};




/**
 * @brief function that transform a LUV real pixel color into XYZ real pixel
 * @param p input color pixel
 * @return corresponding XYZ color
 */
template< typename TInputPixel, typename TOutputPixel >
inline TOutputPixel colorLUVtoXYZ(const TInputPixel& p)
{
	TOutputPixel q; // X,Y,Z

	// compute Y
	if (p[0]> 8.0)
	{
		double d = (p[0] + 16.0)/116.0;
		q[1] = WP_Yn * d*d*d;
	}
	else
	{
		q[1] = WP_Yn * p[0] / 903.3;
	}

	double up = double(p[1]) / (13.0*p[0]) + WP_upn;
	double vp = double(p[2]) / (13.0*p[0]) + WP_vpn;



	//compute X
	q[0] = (9.0*up)/(4.0*vp)*q[1];

	// compute Z
	q[2] = (12.0 - 3.0*up - 20.0*vp)/(4.0*vp)*q[1];

	return q;
}

/// \brief functor version of function for filter usage
template< class TInput, class TOutput>
class fonctorLUVtoXYZ
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		return colorLUVtoXYZ<TInput,TOutput>(p);
	}
};



/// \brief small convenient computing function
inline double func_LAB( double t)
{
	if (t > 0.008856)
		return std::pow(t,1.0/3.0);
	return 7.787*t + 16.0/116.0;
}


/**
 * @brief function that transform a XYZ real pixel color into LAB real pixel
 * @param p input color pixel
 * @return corresponding LAB color
 */
template< typename TInputPixel, typename TOutputPixel >
inline TOutputPixel colorXYZtoLAB(const TInputPixel& p)
{
	TOutputPixel q; // L*,a*,b*

	double YY = p[1]/WP_Yn;

	// compute L*
	if ( YY > 0.008856)
		q[0] = 116.0 * std::pow(YY,1.0/3.0) - 16.0;
	else
		q[0] = 903.3 * YY;

	// compute a*
	q[1] = 500.0*(func_LAB(p[0]/WP_Xn) - func_LAB(p[1]/WP_Yn));

	// compute b*
	q[2] = 200.0*(func_LAB(p[1]/WP_Yn) - func_LAB(p[2]/WP_Zn));

	return q;
}

/// \brief functor version of function for filter usage
template< class TInput, class TOutput>
class fonctorXYZtoLAB
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		return colorXYZtoLAB<TInput,TOutput>(p);
	}
};




/// \brief small convenient computing function
inline double func_inv_LAB( double s)
{
	if (s > 0.206893)
		return s*s*s;
	return s/7.787 - 16.0/903.3;
}

/**
 * @brief function that transform a LAB real pixel color into XYZ real pixel
 * @param p input color pixel
 * @return corresponding XYZ color
 */
template< typename TInputPixel, typename TOutputPixel >
inline TOutputPixel colorLABtoXYZ(const TInputPixel& p)
{
	TOutputPixel q; // X,Y,Z

	// compute Y
	if (p[0]> 8)
	{
		double d = (p[0] + 16.0)/116.0;
		q[1] = WP_Yn * d*d*d;
	}
	else
	{
		q[1] = WP_Yn * p[0] / 903.3;
	}

	//compute X
	q[0] = WP_Xn * func_inv_LAB((p[0]+16.0)/116.0 + p[1]/500.0);

	// compute Z
	q[2] = WP_Zn * func_inv_LAB((p[0]+16.0)/116.0 - p[2]/200.0);

	return q;
}

/// \brief functor version of function for filter usage
template< class TInput, class TOutput>
class fonctorLABtoXYZ
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		return colorLABtoXYZ<TInput,TOutput>(p);
	}
};



/// \brief type definition for filter RGB -> XYZ
template< class TInput, class TOutput>
using FilterRGBtoXYZ = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorRGBtoXYZ< typename TInput::PixelType, typename TOutput::PixelType> >;


/// \brief type definition for filter XYZ -> RGB
template< class TInput, class TOutput>
using FilterXYZtoRGB = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorXYZtoRGB< typename TInput::PixelType, typename TOutput::PixelType> >;


/// \brief type definition for filter XYZ -> LUV
template< class TInput, class TOutput>
using FilterXYZtoLUV = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorXYZtoLUV< typename TInput::PixelType, typename TOutput::PixelType> >;

/// \brief type definition for filter RGB -> HSV
template< class TInput, class TOutput>
using FilterRGBtoHSV = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorRGBtoHSV< typename TInput::PixelType, typename TOutput::PixelType> >;

/// \brief type definition for filter HSV -> RGB
template< class TInput, class TOutput>
using FilterHSVtoRGB = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorHSVtoRGB< typename TInput::PixelType, typename TOutput::PixelType> >;

/// \brief type definition for filter LUV -> XYZ
template< class TInput, class TOutput>
using FilterLUVtoXYZ = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorLUVtoXYZ< typename TInput::PixelType, typename TOutput::PixelType> >;


/// \brief type definition for filter XYZ -> LAB
template< class TInput, class TOutput>
using FilterXYZtoLAB = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorXYZtoLAB< typename TInput::PixelType, typename TOutput::PixelType> >;


/// \brief type definition for filter LAB -> XYZ
template< class TInput, class TOutput>
using FilterLABtoXYZ = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorLABtoXYZ< typename TInput::PixelType, typename TOutput::PixelType> >;



// 255 -> 01


template< class TInput, class TOutput>
class fonctorGray255To01
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		return 1.0/255.0*p;
	}
};



template< class TInput, class TOutput>
using FilterGray255To01 = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorGray255To01< typename TInput::PixelType, typename TOutput::PixelType> >;


template< class TInput, class TOutput>
class fonctorGray01To255
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		return 255.0*p;
	}
};


/// \brief type definition for filter LAB -> XYZ
template< class TInput, class TOutput>
using FilterGray01To255 = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorGray01To255< typename TInput::PixelType, typename TOutput::PixelType> >;





template< class TInput, class TOutput>
class fonctorRGB255To01
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		TOutput q;
		for (uint32_t i=0; i<TInput::Length; ++i)
			q[i] = 1.0/255.0*p[i];
		return q;
	}
};


/// \brief type definition for filter LAB -> XYZ
template< class TInput, class TOutput>
using FilterRGB255To01 = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorRGB255To01< typename TInput::PixelType, typename TOutput::PixelType> >;


template< class TInput, class TOutput>
class fonctorRGB01To255
{
public:
	inline TOutput operator()( const TInput & p ) const
	{
		TOutput q;
		for (uint32_t i=0; i<TInput::Length; ++i)
			q[i] = 255.0*p[i];
		return q;
	}
};



template< class TInput, class TOutput>
using FilterRGB01To255 = itk::UnaryFunctorImageFilter< TInput, TOutput,
						fonctorRGB01To255< typename TInput::PixelType, typename TOutput::PixelType> >;





} //end namespace

} //end namespace





#endif
