#include <ASTex/easy_io.h>
#include <ASTex/slic.h>

using namespace ASTex;


void slicTest(const ASTex::ImageRGBu8 &img, ASTex::ImageRGBu8 &imgSLIC, ASTex::ImageRGBu8& imgSLIC_false )
{
	const float sizeRatio = 0.02f;
	const int compactness = 15;
	const int iterations  = 10;

	const int superPixelSize = std::ceil( img.width()*sizeRatio * img.height()*sizeRatio );

	typedef ASTex::SLICSuperPixelFilter< ASTex::NeighborhoodPeriodic > SLICFilterType;
	SLICFilterType::Pointer slicFilter = SLICFilterType::New();
	slicFilter->SetInput( img.itk() );
	slicFilter->setSuperPixelSize( superPixelSize );
	slicFilter->setIterationNumber( iterations );
	slicFilter->setCompactness( compactness );

	slicFilter->Update();


	// compute 2 color image
	// one with mean color
	// one with random color
	ASTex::ImageGrayu32 labelMap( slicFilter->GetOutput() );

	uint32_t labelMax = 0;
	labelMap.for_all_pixels( [&]( ASTex::ImageGrayu32::PixelType p )
	{
		if( p > labelMax )
			labelMax = p;
	});

	std::vector< itk::RGBPixel<float> > labelColors( labelMax+1 );

	for( size_t i=0; i<=labelMax; ++i )
		labelColors[i].Set( 0.0f, 0.0f, 0.0f );

	for( unsigned int i=0; i<slicFilter->getSuperPixelCount(); ++i )
	{
		auto &sp = slicFilter->getSuperPixel( i );
		for( auto &p : sp.pixels )
			labelColors[i] += img.pixelAbsolute( p[0], p[1] );

		for( int k=0; k<3; ++k )
			labelColors[i][k] /= sp.pixels.size();
	}

	imgSLIC.initItk( img.width(), img.height() );
	imgSLIC.for_all_pixels( [&](ASTex::ImageRGBu8::PixelType &p, int x, int y)
	{
		auto &c = labelColors[ labelMap.pixelAbsolute(x,y) ];
		p.Set( (uint8_t) c.GetRed(), (uint8_t) c.GetGreen(), (uint8_t) c.GetBlue() );
	});


	for( size_t i=0; i<=labelMax; ++i )
	{
		labelColors[i].SetRed  ( 64+(rand()%192) );
		labelColors[i].SetGreen( 64+(rand()%192) );
		labelColors[i].SetBlue ( 64+(rand()%192) );
	}

	imgSLIC_false.initItk( img.width(), img.height() );
	imgSLIC_false.for_all_pixels( [&](ASTex::ImageRGBu8::PixelType &p, int x, int y)
	{
		auto &c = labelColors[ labelMap.pixelAbsolute(x,y) ];
		p.Set( (uint8_t) c.GetRed(), (uint8_t) c.GetGreen(), (uint8_t) c.GetBlue() );
	});
//#endif
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
		ASTex::ImageRGBu8 imgSLIC;
		ASTex::ImageRGBu8 imgSLIC_false;
		slicTest(img, imgSLIC, imgSLIC_false);

		std::string name1 = ASTex::IO::remove_ext(argv[1])+"_SLIC_MEAN.png";
		imgSLIC.save(name1);
		std::cout << name1 << " generated"<< std::endl;

		std::string name2 = ASTex::IO::remove_ext(argv[1])+"_SLIC_FALSE.png";
		imgSLIC_false.save(name2);
		std::cout << name2 << " generated"<< std::endl;
	}

    return EXIT_SUCCESS;
}
