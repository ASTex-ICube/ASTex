#include <ASTex/easy_io.h>
#include <ASTex/fourier.h>
#include <ASTex/utils.h>
#include <cmath>

#include <itkForwardFFTImageFilter.h>
#include <itkInverseFFTImageFilter.h>

#include <iostream>

using namespace ASTex;
/**
 * @brief decomposition periodic + smooth
 * @param input input grayd image
 * @return a pair of image (periodic,smooth)
 */
std::pair<ImageGrayd,ImageGrayd> decompo(const ImageGrayd& input)
{
    int W = input.width();
    int H = input.height();

    // 1 - compute boundary intensity (v)
    ImageGrayd bound(W,H,true);
    for (int x=0;x<W;x++)
    {
        double v = input.pixelAbsolute(x,0) - input.pixelAbsolute(x,H-1);
        bound.pixelAbsolute(x,0) += v;
        bound.pixelAbsolute(x,H-1) -= v;
    }
    for (int y=0;y<H;y++)
    {
        double v = input.pixelAbsolute(0,y) - input.pixelAbsolute(W-1,y);
        bound.pixelAbsolute(0,y) += v;
        bound.pixelAbsolute(W-1,y) -= v;
    }

    // 2 - DFT of v
    using  FFTType = itk::ForwardFFTImageFilter< ImageGrayd::ItkImg, ImageSpectralcd::ItkImg >;
    FFTType::Pointer fftFilter = FFTType::New();
    fftFilter->SetInput(bound.itk());
    fftFilter->Update();
    ImageSpectralcd ft_bound(fftFilter->GetOutput());

    // 3 - divide by cosin....
    //
    // compute coef by row
    double cx = 2.0*M_PI/double(W);
    std::vector<double> coef_x(W);
    for (int x=0; x<W; ++x)
        coef_x[x] = 2*std::cos(cx*x);
    // by column
    double cy = 2.0*M_PI/double(H);
    std::vector<double> coef_y(H);
    for (int y=0; y<H; ++y)
        coef_y[y] = 2*std::cos(cy*y);
    // apply computations on all pixels
    ft_bound.for_all_pixels([&] (ImageSpectralcd::PixelType& P, int x, int y)
    {
        P /= 4 - coef_x[x] - coef_y[y]; // WARNING OPPOSITE SIGN COMPARE TO ARTICLE
    });
    // except for (0,0) -> 0,0
    ft_bound.pixelAbsolute(0,0) = ImageSpectralcd::PixelType(0,0);


    // 4 - FFT inverse => smooth
    //
    using  IFFTType = itk::InverseFFTImageFilter< ImageSpectralcd::ItkImg, ImageGrayd::ItkImg >;
    IFFTType::Pointer ifftFilter = IFFTType::New();
    ifftFilter->SetInput(ft_bound.itk());
    ifftFilter->Update();
    ImageGrayd smooth(ifftFilter->GetOutput());


    // 4' - compute periodic = input - smooth
    //
    ImageGrayd perio(W,H,true);
	for(auto its = std::make_tuple(perio.beginIterator(),input.beginConstIterator(),smooth.beginConstIterator());
        !std::get<0>(its).IsAtEnd(); ++std::get<0>(its), ++std::get<1>(its), ++std::get<2>(its))
    {
        std::get<0>(its).Value() = std::get<1>(its).Value() - std::get<2>(its).Value();
    }

    return std::make_pair(perio,smooth);
}

/**
 * @brief correct images ?
 * @param in input
 * @param per periodic (clamped to [0,1])
 * @param smo smooth (shifted +0.5 & clamped to [0,1])
 */
void correct(const ImageGrayd& in, ImageGrayd& per, ImageGrayd& smo)
{
    double min0 = compute_min(in);
    double max0 = compute_max(in);
    std::cout << "input " << min0 << " / " << max0 << std::endl;
    double min1 = compute_min(per);
    double max1 = compute_max(per);
    std::cout << "perio " << min1 << " / " << max1 << std::endl;
    double min2 = compute_min(smo);
    double max2 = compute_max(smo);
    std::cout << "smooth " << min2 << " / " << max2<< std::endl;

    per.for_all_pixels([&] (double& p)
    {
        if (p>1.0)
            p=1.0;
        if (p<0.0)
            p=0.0;
    });

    smo.for_all_pixels([&] (double& p)
    {
        p += 0.5;
        if (p>1.0)
            p=1.0;
        if (p<0.0)
            p=0.0;
    });
}

ImageGrayd x4(const ImageGrayd& input)
{
    int W = input.width();
    int H = input.height();

	ImageGrayd quad(W*2,H*2);

    for_indices(0,W,0,H,[&] (int x, int y)
    {
        auto v = input.pixelAbsolute(x,y);
        quad.pixelAbsolute(x,y) = v;
        quad.pixelAbsolute(W+x,y) = v;
        quad.pixelAbsolute(x,H+y) = v;
        quad.pixelAbsolute(W+x,H+y) = v;
    });

    return quad;
}


int main(int argc, char** argv)
{
    if( argc < 2 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " <Source_file>" << std::endl;
        return EXIT_FAILURE;
    }

    // Get the arguments
    std::string source_dir = argv[1];

    // Get the name of the file without extention en creation the folder for the res
//    std::string name_file = ASTex::IO::remove_ext(ASTex::IO::remove_path(source_dir));
	std::string out_dir = ASTex::IO::remove_ext(source_dir) +"_perio_smooth/";

    ASTex::create_directory(out_dir);

    // LOAD INPUT
    ASTex::ImageGrayd input;

    ASTex::IO::loadu8_in_01(input,source_dir);
   
    //Split  into Smooth + Periodic using "Periodic plus smooth image decomposition , Lionel Moisan" and keep the periodic part as rectified input
    auto outs = decompo(input);

    ASTex::IO::save01_in_u8(outs.first,out_dir+"input_periodic.png");
    ASTex::IO::save01_in_u8(x4(outs.first),out_dir+"input_smooth_test_X4.png");
    ASTex::IO::save01_in_u8(x4(input),out_dir+"input_X4.png");

}
