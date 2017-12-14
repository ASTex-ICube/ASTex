#include "texton_cpu.h"

namespace ASTex
{

Stamper::Stamper(const std::vector<vec2> &pointArray, const ImageRGBd &tampon) :
    m_pointArray(pointArray), m_stamp(tampon)
{}

BombingStamper::BombingStamper(const std::vector<vec2> &pointArray, const ImageRGBd &tampon) :
    Stamper(pointArray, tampon)
{}


ImageRGBd BombingStamper::generate(int imageWidth, int imageHeight)
{
    ImageRGBd im_out;

//    std::vector<vec2> points;
//    points = m_pointArray->Generate(nb_points);

    im_out.initItk(imageWidth, imageHeight, true);

    for(std::vector<vec2>::const_iterator it=m_pointArray.begin(); it!=m_pointArray.end(); ++it)
    {
        int i = im_out.width() * (*it)[0]; //i & j: single point coordinates in im_out
        int j = im_out.height() * (*it)[1];

        int textonWidth = m_stamp.width();
        int textonHeight = m_stamp.height();
        Region reg = gen_region(i-textonWidth/2, j-textonHeight/2, textonWidth, textonHeight); //the region we stamp
        im_out.for_region_pixels(reg, [&] (ImageRGBd::PixelType& pix, int x, int y)
        {
            if(x>0 && x<im_out.width() && y>0 && y<im_out.height())
            {
                int tx=x-i+textonWidth/2; //texton coordinate
                int ty=y-j+textonHeight/2; //texton coordinate
                pix = m_stamp.pixelAbsolute(tx, ty);
            }
        });
    }

    return im_out;
}

TextonStamper::TextonStamper(const std::vector<vec2> &pointArray, const ImageRGBd &tampon) :
    Stamper(pointArray, tampon), m_ratioX(1.0), m_ratioY(1.0), m_periodicity(false), m_bilinearInterpolation(false)
{}

ImageRGBd TextonStamper::generate(int imageWidth, int imageHeight)
{
    ImageRGBd im_out;

    ImageRGBd::PixelType mean, sum;
    for(int i=0; i<3; ++i)
    {
        sum[i]=0;
        mean[i]=0;
    }
    m_stamp.for_all_pixels([&] (const ImageRGBd::PixelType &pix)
    {
        for(int i=0; i<3; ++i)
            mean[i] += pix[i];
    });
    //sum=mean;
    for(int i=0; i<3; ++i)
        mean[i]/=m_stamp.width()*m_stamp.height();

    m_stamp.for_all_pixels([&] (ImageRGBd::PixelType &pix)
    {
        for(int i=0; i<3; ++i)
        {
            pix[i] -= mean[i];
            pix[i] /= std::sqrt(m_stamp.width() * m_stamp.height());
            sum[i] += pix[i];
        }
    });

//    std::vector<vec2> points;
//    points = m_pointArray->Generate(nb_points);

    im_out.initItk(imageWidth, imageHeight, true);

    double textonWidth = m_stamp.width();
    double textonHeight = m_stamp.height();

    for(std::vector<vec2>::const_iterator it=m_pointArray.begin(); it!=m_pointArray.end(); ++it)
    {
//        int i = im_out.width() * (*it)[0]; //i & j: single point coordinates in im_out (double)
//        int j = im_out.height() * (*it)[1];

        double i, j;

        if(m_periodicity)
        {
            i = imageWidth * (*it)[0];
            j = imageHeight * (*it)[1];
        }

        if(!m_periodicity)
        {
            //we need to introduce margins here, to shoot textons outside of the domain
            i = (imageWidth + textonWidth * 2 ) * (*it)[0] - textonWidth;
            j = (imageHeight + textonHeight * 2 ) * (*it)[1] - textonHeight;
        }

        int otx=int(i-textonWidth/2.0), oty=int(j-textonHeight/2.0); //texton origin in texture space (top left)
        double dx= m_bilinearInterpolation ? (i-textonWidth/2.0)-otx : 0; //gap between texton and image pixels
        double dy= m_bilinearInterpolation ? (j-textonHeight/2.0)-oty : 0; //0 if no bilinear interpolation : equivalent to clamping dx and dy

        Region reg = gen_region(otx, oty, textonWidth, textonHeight); //the region we stamp
        if(m_periodicity)
            im_out.for_region_pixels(reg, [&] (ImageRGBd::PixelType& pix, int x, int y) //with periodicity
            {

                int tx=x-otx; //texton coordinate in texton space
                int ty=y-oty; //texton coordinate

                for(int c=0; c<3; ++c) //c: channel
                {
                    double interpolatedValue = 0.0;
                    if(tx-1 > 0)
                        if(ty-1 > 0)
                            interpolatedValue += (dx*dy)*m_stamp.pixelAbsolute(tx-1, ty-1)[c];
                        if(ty < textonHeight-1)
                            interpolatedValue += (dx*(1-dy))*m_stamp.pixelAbsolute(tx-1, ty)[c];
                    if(tx < textonWidth-1)
                        if(ty-1 > 0)
                            interpolatedValue += ((1-dx)*dy)*m_stamp.pixelAbsolute(tx-1, ty-1)[c];
                        if(ty < textonHeight-1)
                            interpolatedValue += ((1-dx)*(1-dy))*m_stamp.pixelAbsolute(tx, ty)[c];

                    if(x<0)
                        x+=im_out.width();
                    if(y<0)
                        y+=im_out.height();
                    im_out.pixelAbsolute(x%imageWidth, y%imageHeight)[c] += interpolatedValue;
                }
            });
        else
            im_out.for_region_pixels(reg, [&] (ImageRGBd::PixelType& pix, int x, int y) //without periodicity
            {
                if(x>0 && x<imageWidth && y>0 && y<imageHeight)
                {
                    int tx=x-otx; //texton coordinate
                    int ty=y-oty; //texton coordinate

                    for(int c=0; c<3; ++c) //c: channel
                    {
                        double interpolatedValue = 0.0;
                        if(tx-1 > 0)
                            if(ty-1 > 0)
                                interpolatedValue += (dx*dy)*m_stamp.pixelAbsolute(tx-1, ty-1)[c];
                            if(ty < textonHeight-1)
                                interpolatedValue += (dx*(1-dy))*m_stamp.pixelAbsolute(tx-1, ty)[c];
                        if(tx < textonWidth-1)
                            if(ty-1 > 0)
                                interpolatedValue += ((1-dx)*dy)*m_stamp.pixelAbsolute(tx-1, ty-1)[c];
                            if(ty < textonHeight-1)
                                interpolatedValue += ((1-dx)*(1-dy))*m_stamp.pixelAbsolute(tx, ty)[c];

                        pix[c] += interpolatedValue;

                        //std::cout << interpolatedValue << std::endl;
                    }
                }
            });
    }

    int size2Normalize = m_periodicity ? imageWidth*imageHeight : (imageWidth+textonWidth * 2)*(imageHeight+textonHeight * 2);

    float mean_nb_of_impacts = float(m_pointArray.size()) * (m_stamp.width() * m_stamp.height()) / size2Normalize;
    float lambda = mean_nb_of_impacts/(m_stamp.width() * m_stamp.height());

    im_out.for_all_pixels([&] (ImageRGBd::PixelType &pix) {
        for(int i=0; i<3; ++i)
        {
            pix[i] = 1.0/sqrt(lambda) * (pix[i] - lambda * sum[i]) + mean[i];
            pix[i] = pix[i] > 1.0 ? 1.0 : (pix[i] < 0.0 ? 0.0 : pix[i]);
        }
    });

    return im_out;
}

//vérifier paramètres
//le tampon s'échantillonne
//mip-map => taille vs resolution

} //namespace ASTex
