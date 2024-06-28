#ifndef _MP_INPUTS_
#define _MP_INPUTS_

#include <Eigen/Dense>
#include <ASTex/image_gray.h>


class MicroPatternsInput
{
public:
	/**
	* @param uv coordinate to fetch
	* @param level level of mipmap in which to fetch
	* @return Eigen::Vector4d(noiseA-average(NoiseA),noiseB-average(noiseB),var(noiseA),var(noiseB)
	*/
	virtual Eigen::Vector4d fetch(const Eigen::Vector2d& uv, double level) const = 0;

	/**
	* @return average value (noiseA,noiseB)
	*/
	virtual Eigen::Vector2d fetch_average() const = 0;

	/**
	* @return width of based squared textures of mipmap
	*/
	virtual int width() const = 0;

	/**
	* @return return number of level in the mipmap
	*/
	virtual int depth() const = 0;
};

#endif