#ifndef __ASTEX_QUILTING_ERROR_BUFFER__
#define __ASTEX_QUILTING_ERROR_BUFFER__

#include <iostream>
#include <map>
#include <vector>
#include <ASTex/image_rgb.h>


using namespace ASTex;

/**
 * @brief The ErrorCutPathBuffer class
 */
class ErrorCutPathBuffer
{
	double* data_;
	int width_;
	int height_;

public:
	ErrorCutPathBuffer(int w, int h);

	~ErrorCutPathBuffer();

	inline double  operator()(int i, int j) const { return data_[i+j*width_]; }
	inline double& operator()(int i, int j)       { return data_[i+j*width_]; }
	inline double  operator()(const itk::Offset<2>& off) const { return data_[off[0]+off[1]*width_]; }
	inline double& operator()(const itk::Offset<2>& off)       { return data_[off[0]+off[1]*width_]; }

	int minOfRow(int r);

	int minOfColumn(int c);

	int minOfRow3(int i, int j);

	int minOfColumn3(int i, int j);
};



class ContentLabel
{
	std::vector<ImageRGBu8> label_images_;

public:

	/**
	 * @brief load nb image (bnameX  with X in [1,nb])
	 * @param bname base name
	 * @param nb
	 */
	void load(const std::string& bname);

	/**
	 * @brief do pixels have same label
	 * @param p position pixel 1
	 * @param q position pixel 1
	 * @param l level of content
	 * @return true if same level
	 */
	inline bool has_same_label(const Index& p, const Index& q, int l)
	{
		if (l<0)
		{
			int ll = label_images_.size() + l;
			return label_images_[ll].pixelAbsolute(p) == label_images_[ll].pixelAbsolute(q);
		}
		return label_images_[l].pixelAbsolute(p) == label_images_[l].pixelAbsolute(q);
	}

	/**
	 * @brief compute number of level in which labels are equals between 2 pixels
	 * @param p first pixel position
	 * @param q second pixel position
	 * @return num of level where label are equal
	 */
	int levels_equality(const Index& p, const Index& q);

	/**
	 * @brief get number of levels
	 * @return
	 */
	inline int get_number_of_levels() {	return int(label_images_.size()); }
};




//template<uint32_t NB>
//class GaussianValues
//{
//	double gaussian_[NB+1];
//	double max_;
//public:
//	GaussianValues():max_(0.0)
//	{}

//	GaussianValues(double max_val, double sigma)
//	{
//		init(max_val,sigma);
//	}

//	inline void init(double max_val, double sigma)
//	{
//		max_ = max_val;

//		double inc = max_val/(NB-1);
//		double x=0.0;
//		for (uint32_t i=0; i<NB;++i)
//		{
//			gaussian_[i] = (std::exp(-(x*x)/(2.0*sigma*sigma)) / (sigma*std::sqrt(2.0*M_PI)));
//			x+=inc;
//		}
//		// to avoid rounding pb
//		gaussian_[NB] = gaussian_[NB-1];
//	}

//	inline double operator[](double d)
//	{
//		return gaussian_[size_t((d*(NB-1))/max_)];
//	}
//};


#endif

