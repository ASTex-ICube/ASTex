#include <ASTex/image_collector.h>
#include <iostream>

namespace ASTex
{

void ImageCollector::add (const ImageGrayd& im, double mean_shift)
{
    m_images.push_back(im);
    m_shifts.push_back(mean_shift);
}

ImageGrayd ImageCollector::collect ()
{
	int W = 0, H = 0;

    for (uint32_t i=0; i<m_images.size(); ++i)
    {
		W += m_images[i].width();
		H = std::max(H,m_images[i].height());
    }

	ImageGrayd image(W,H,true);

    uint32_t dx=0;
    for (uint32_t i=0; i<m_images.size(); ++i)
    {
		uint32_t Wi = m_images[i].width();
		uint32_t Hi = m_images[i].height();

		uint32_t dxi = m_images[i].itk()->GetOrigin()[0];
		uint32_t dyi = m_images[i].itk()->GetOrigin()[1];

        double shift = m_shifts[i];

        for (uint32_t x=0 ; x < Wi ; ++x)
        {
            for (uint32_t y=0 ; y < Hi ; ++y)
            {
				image.pixelAbsolute(x+dx,y) = m_images[i].pixelAbsolute(x+dxi,y+dyi) + shift;
            }
        }
        dx += Wi;
    }
    return image;
}

} // namespace
