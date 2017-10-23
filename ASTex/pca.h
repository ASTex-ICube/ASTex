#ifndef _ASTEX_PCA_H_
#define _ASTEX_PCA_H_

#include <ASTex/mask.h>
#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

namespace ASTex
{

class ASTEX_API PCA
{
protected:
	ImageRGBd m_im;
	// eigenvectors
	double m_v1 [3];
	double m_v2 [3];
	double m_v3 [3];
	// mean values
	double m_meancolor [3];
	// covariance matrix
	double m_covar [3][3];


	double dot(double v[3], double w[3]);

	double normalize(double v[3]);

	void cross(double v[3], double w[3], double res[3]);


public :
	PCA (const ImageRGBd& input);

	void computeMeanColor(const Mask& mask);

	void computeCovariance(const Mask& mask);

	void mat_mult(double vin[3], double vout[3]);

	void computeEigenVectors();

	void computePCA (const Mask& m);

	void exportToTXT(std::string savefile);

	void project(ImageGrayd& res1, ImageGrayd& res2, ImageGrayd& res3) const;

	void back_project(const ImageGrayd& coord1, const ImageGrayd& coord2, const ImageGrayd& coord3, ImageRGBd& res) const;
};

} //ASTex

#endif
