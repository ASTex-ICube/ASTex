#include <ASTex/saliency/contrast_based_saliency.h>
#include <ASTex/easy_io.h>
#include <stdlib.h>



ASTex::ImageGrayf saliencyFiltersTest( ASTex::ImageRGBu8 &img )
{
    const float sizeRatio     = 0.02f;
    const int compactness     = 15;
    const int iterations      = 10;
    const double uniqueness   = 0.25;
    const double distribution = 20.0;
    const double k            = 3.0;
    const double alpha        = 1.0 / 30.0;
    const double beta         = 1.0 / 30.0;

    const int superPixelSize = std::ceil( img.width()*sizeRatio * img.height()*sizeRatio );

    typedef ASTex::ContrastBasedSaliencyFilter< ASTex::ImageGrayf::ItkImg > SaliencyFilterType;
    SaliencyFilterType::Pointer saliencyFilter = SaliencyFilterType::New();
    saliencyFilter->SetInput( img.itk() );
    saliencyFilter->setSuperPixelSize( superPixelSize );
    saliencyFilter->setIterationNumber( iterations );
    saliencyFilter->setCompactness( compactness );
    saliencyFilter->setUniquenessVariance( uniqueness );
    saliencyFilter->setDistributionVariance( distribution );
    saliencyFilter->setExponentialScaling( k );
    saliencyFilter->setAlphaBeta( alpha, beta );

    saliencyFilter->Update();

	ASTex::ImageGrayf im_out(saliencyFilter->GetOutput());

	return im_out;
}


int main( int argc, char **argv )
{
	ASTex::ImageRGBu8 img;
	bool ok = false;

	if( argc < 2 )
		ok = img.load(TEMPO_PATH+"single-leaf-photo-03.jpg");
	else
		ok = img.load( argv[1] );

	if (ok)
	{
		ASTex::ImageGrayf out = saliencyFiltersTest( img );
		std::string out_name = ASTex::IO::remove_ext(argv[1])+"saliency_out.png";
		ASTex::IO::save01_in_u8(out,out_name);
		std::cout << out_name << " generated"<< std::endl;
	}

    return EXIT_SUCCESS;
}
