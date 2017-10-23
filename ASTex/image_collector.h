#ifndef __ASTEX_IMG__COLLECTOR__
#define __ASTEX_IMG__COLLECTOR__

#include <ASTex/image_gray.h>

namespace ASTex
{

/**
 * \class image_collector io.h
 * \brief collects images and builds a patchwork
 */
class ASTEX_API ImageCollector
{
public :
	/**
	 * \brief add an image in the collection
     * \param [in] im scalar image
	 * \param [in] mean_shift (default = 0.0) value for the mean value shifting
	 */
	void add (const ImageGrayd& im, double mean_shift = 0.0);
	/**
	 * \brief computes a patchwork image from the collection
	 * \return patchwork image
	 */
	ImageGrayd collect ();

private :
	std::vector<ImageGrayd> m_images;
    std::vector<double> m_shifts;
};


}
#endif
