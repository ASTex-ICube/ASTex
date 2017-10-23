#include "image_gray.h"
#include "masks.h"
#include "fourier.h"
#include "io.h"

#include <itkRegionOfInterestImageFilter.h>
#include <itkCyclicShiftImageFilter.h>

/** \file
 * \brief test functions for dev of noise synthesis
*/

using namespace ASTex;

/**
 * \brief tests extraction of spectrum by auto-correlation
*/
int testSpectrumExtraction3(std::string inputfile, std::string weights0_file, std::string weights1_file, std::string outputfile){
    // LOAD INPUT
    ImageGrayd input;
    double mean = load(input,inputfile, true);
    monitorStats (input, "input");

    ImageGrayd weights0;
    load(weights0,weights0_file, false);

    ImageGrayd weights1;
    load(weights1,weights1_file, false);

    // COLLECTION OF IMAGES
    image_collector collec;
    collec.add(input,mean);
    collec.add(weights0);
    collec.add(weights1);

    // AUTO-CORRELATION

    double mask_threshold = 0.4;
    mask_largest_values mask0 (weights0,mask_threshold);
    mask_largest_values mask1 (weights1,mask_threshold);

    const int Ntest = 4;
    double prop_th [Ntest] = {0.0, 0.5, 0.8, 0.9};
    int sp_size [Ntest] = {128, 64, 32, 16};
//    const int Ntest = 5;
//    double prop_th [Ntest] = {0.0, 0.5, 0.8, 0.9, 0.95};
//    int sp_size [Ntest] = {32, 32, 32, 32, 32};

    for (int i = 0; i < Ntest; ++i) {
        int T = sp_size[i];
        double proportion_threshold = prop_th[i] ;

        ImageGrayd modulus;
        modulus.initItk(T,T);
        spectrum_by_autocorrelation_small_size(input,modulus,mask0,proportion_threshold);

        ImageGrayd synthesis0;
//        synthesis0.initItk(input.width(),input.height());
        synthesis0.initItk(128,128);
        RPnoise_mosaic(modulus,synthesis0,4);

        std::cout << "spectrum size =" << T << "\t proportion threshold = " << proportion_threshold << std::endl;

        collec.add(modulus, 0.0);
        collec.add(synthesis0, mean);
        monitorStats (modulus, "modulus 0");
        monitorStats (synthesis0, "synthesis 0");

        spectrum_by_autocorrelation_small_size(input,modulus,mask1,proportion_threshold);

        ImageGrayd synthesis1;
//        synthesis1.initItk(input.width(),input.height());
        synthesis1.initItk(128,128);
        RPnoise_mosaic(modulus,synthesis1,4);

        collec.add(modulus, 0.0);
        collec.add(synthesis1, mean);
        monitorStats (modulus, "modulus 1");
        monitorStats (synthesis1, "synthesis 1");
    }

    save(collec.collect(),outputfile, 0.0);

    return 0;
}


/**
 * \brief tests extraction of spectrum by auto-correlation
*/
int testSpectrumExtraction2(std::string inputfile, std::string outputfile){
    // LOAD INPUT
    ImageGrayd input;
    double mean = load(input,inputfile, true);
    monitorStats (input, "input");

    // COLLECTION OF IMAGES
    image_collector collec;

    // FFT modulus
    ImageGrayd fft_modulus, fft_phase;
    fftForwardModulusAndPhase(input, fft_modulus, fft_phase);
    collec.add(fft_modulus);
    collec.add(input, mean); // add input to collection

    // AUTO-CORRELATION

    double prop = 0.8;
    mask_random m (input.width(),input.height(),prop);
//    for (double prop = 1.0; prop >= 0.2; prop -= 0.2 )
    for (int T = input.width(); T >= 16; T/=2 )
    {
//        mask_random m (input.width(),input.height(),prop);

        ImageGrayd modulus;
        modulus.initItk(T,T);
        spectrum_by_autocorrelation_small_size(input,modulus,m);
//        setPower(modulus,getPower(input,m));

        // SYNTHESIS
        ImageGrayd synthesis;
        ImageGrayd phase;
        phase.initItk(modulus.width(),modulus.height());
        randomPhase(phase);
        fftInverseModulusAndPhase(modulus,phase,synthesis);

        collec.add(modulus, 0.0);
        collec.add(synthesis, mean);
        std::cout << T << std::endl;
        monitorStats (modulus, "modulus");
        monitorStats (synthesis, "synthesis");
    }

    save(collec.collect(),outputfile, 0.0);

    return 0;
}


/**
 * \brief tests extraction of spectrum by auto-correlation
*/
int testSpectrumExtraction1(std::string inputfile, std::string outputfile){
    // LOAD INPUT
    ImageGrayd input;
    double mean = load(input,inputfile, true);
    monitorStats (input, "input");

    // COLLECTION OF IMAGES
    image_collector collec;
    collec.add(input, mean); // add input to collection

    // WELCH
//    ImageGrayd welch_modulus;
//    welch(input, welch_modulus, 1);
//    collec.add(welch_modulus);
//    monitorStats (welch_modulus, "welch_modulus");

    // FFT modulus
    ImageGrayd fft_modulus, fft_phase;
    fftForwardModulusAndPhase(input, fft_modulus, fft_phase);
    collec.add(fft_modulus);
//    monitorStats (welch_modulus, "welch_modulus");

    // AUTO-CORRELATION

//    mask_true m;
    for (double prop = 1.0; prop >= 0.2; prop -= 0.2 )
    {
        mask_random m (input.width(),input.height(),prop);

        ImageGrayd modulus;
        spectrum_by_autocorrelation_full_size(input,modulus,m);

        // SYNTHESIS
        ImageGrayd synthesis;
        ImageGrayd phase;
        phase.initItk(modulus.width(),modulus.height());
        randomPhase(phase);
        fftInverseModulusAndPhase(modulus,phase,synthesis);


        // DOWN SAMPLING
        ImageGrayd modulus_DS;
        modulus_DS.initItk(modulus.width()/2,modulus.height()/2);
        modulus_DS.setCenter(modulus.width()/4,modulus.height()/4);
        down_sampling(modulus,modulus_DS);

        ImageGrayd synthesis_DS;
        ImageGrayd phase_DS;
        phase_DS.initItk(modulus_DS.width(),modulus_DS.height());
        randomPhase(phase_DS);
        fftInverseModulusAndPhase(modulus_DS,phase_DS,synthesis_DS);

        // DOWN SAMPLING
        ImageGrayd modulus_DDS;
        modulus_DDS.initItk(modulus.width()/4,modulus.height()/4);
        modulus_DDS.setCenter(modulus.width()/8,modulus.height()/8);
        down_sampling(modulus_DS,modulus_DDS);

        ImageGrayd synthesis_DDS;
        ImageGrayd phase_DDS;
        phase_DDS.initItk(modulus_DDS.width(),modulus_DDS.height());
        randomPhase(phase_DDS);
        fftInverseModulusAndPhase(modulus_DDS,phase_DDS,synthesis_DDS);


        collec.add(modulus, 0.0);
        collec.add(synthesis, mean);
        collec.add(modulus_DS, 0.0);
        collec.add(synthesis_DS, mean);
        collec.add(modulus_DDS, 0.0);
        collec.add(synthesis_DDS, mean);

        std::cout << prop << std::endl;
        monitorStats (modulus, "modulus");
//        monitorStats (synthesis, "synthesis");
//        monitorStats (modulus_DS, "modulus down sampled");
//        monitorStats (synthesis_DS, "synthesis down sampled");
    }

    save(collec.collect(),outputfile, 0.0);

    return 0;
}

/**
 * \brief tests masks : results for various proportions of fixed phase
*/
int testMask(std::string inputfile, std::string outputfile){

	// LOAD INPUT
	ImageGrayd input;
	double mean = load(input,inputfile, true);
	monitorStats(input, "Input (space)"); // look stats

	// COLLECTION OF IMAGES
	image_collector collec;
	collec.add(input, mean); // add input to collection

	// COMPUTE FFT
	ImageGrayd modulus;
	welch(input, modulus, 10);
	collec.add(modulus); // add modulus to collection
//	monitorStats(modulus, "modulus"); // look stats

	double prop [8] = {0,0.01,0.02, 0.05, 0.1, 0.2, 0.5, 1};
	for (int i=0; i<8; ++i)
	{
		// RANDOM PHASE
		ImageGrayd phase;
		phase.initItk(modulus.width(),modulus.height(),true);
//		randomPhase(phase,mask_smallest_values(modulus,1-prop[i]));
		randomPhase(phase,mask_largest_values(modulus,prop[i]));
//		monitorStats(phase, "phase"); // look stats
		collec.add(phase,0.0); // add result to collection

		// INVERSE FFT
		ImageGrayd result;
		fftInverseModulusAndPhase(modulus, phase, result);
		monitorStats(result, "Output (space)"); // look stats
		collec.add(result,mean); // add result to collection
	}
	for (int i=0; i<8; ++i)
	{
		// RANDOM PHASE
		ImageGrayd phase;
		phase.initItk(modulus.width(),modulus.height(),true);
		randomPhase(phase,mask_smallest_values(modulus,1-prop[i]));
//		randomPhase(phase,mask_largest_values(modulus,prop[i]));
//		monitorStats(phase, "phase"); // look stats
		collec.add(phase,0.0); // add result to collection

		// INVERSE FFT
		ImageGrayd result;
		fftInverseModulusAndPhase(modulus, phase, result);
//		monitorStats(result, "Output (space)"); // look stats
		collec.add(result,mean); // add result to collection
	}
	// SAVE
	save(collec.collect(),outputfile, 0.0);

/*
	ImageGrayd::ItkImg::RegionType region;
	ImageGrayd::ItkImg::IndexType regionIndex;
	ImageGrayd::ItkImg::SizeType regionSize;
	regionSize[0] = 4;
	regionSize[1] = 4;
	region.SetSize(regionSize);
	regionIndex[0] = 0;
	regionIndex[1] = 0;
	region.SetIndex(regionIndex);

	ImageGrayd image;
	image.initItk(4,4,true);

	image.pixelAbsolute(2,0) = 0.0;
	image.pixelAbsolute(2,1) = 0.1;
	image.pixelAbsolute(2,2) = 0.2;
	image.pixelAbsolute(2,3) = 0.3;
	image.pixelAbsolute(1,0) = 1.0;
	image.pixelAbsolute(1,1) = 1.1;
	image.pixelAbsolute(1,2) = 1.2;
	image.pixelAbsolute(1,3) = 1.3;
	image.pixelAbsolute(0,0) = 2.0;
	image.pixelAbsolute(0,1) = 2.1;
	image.pixelAbsolute(0,2) = 2.2;
	image.pixelAbsolute(0,3) = 2.3;
	image.pixelAbsolute(3,0) = 3.0;
	image.pixelAbsolute(3,1) = 3.1;
	image.pixelAbsolute(3,2) = 3.2;
	image.pixelAbsolute(3,3) = 3.3;

	for (int i=0; i<=20; ++i)
		mask_largest_values m (image,0.05*i);
*/
	return EXIT_SUCCESS;
}

/**
 * \brief tests animation : use command line "convert -delay 20 result_*.png anim.gif"
*/
int testAnimation(std::string inputfile, std::string output_prefix)
{
    // LOAD
    ImageGrayd input;
    double mean = load(input,inputfile, true);

    // COMPUTE FFT
    ImageGrayd modulus, phase;

    welch(input, modulus, 10);

    ImageGrayd::ItkImg::RegionType phaseRegion (modulus.itk()->GetLargestPossibleRegion());
    phaseOfRegion(input, phaseRegion,phase);

    // INVERSE FFT
    ImageGrayd result;
//    fftInverseModulusAndPhase(modulus, phase, result);

	for (int i=100; i<228; ++i)
	{

		phaseShiftFromSpaceShift(phase, 1, -1);

		fftInverseModulusAndPhase(modulus, phase, result);


		std::stringstream outputfile;
		outputfile << output_prefix << i << ".png";
		save(result,outputfile.str(), 0.0);
	}

	/*
	// TEST
	typedef itk::RegionOfInterestImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > ROIFilterType;
	ROIFilterType::Pointer roiFilter = ROIFilterType::New();
	roiFilter->SetInput(input.itk());
	roiFilter->SetRegionOfInterest(phaseRegion);
	roiFilter->Update();
	typedef itk::CyclicShiftImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > CyclicShiftFilterType;
	CyclicShiftFilterType::Pointer cyclicShiftFilter = CyclicShiftFilterType::New();
	cyclicShiftFilter->SetInput(roiFilter->GetOutput());
	ImageGrayd::ItkImg::OffsetType offset;
	offset[0]=0;
	offset[1]=0;

    for (int i=100; i<228; ++i)
    {
		offset[0]+=1;
		offset[1]+=-1;
		cyclicShiftFilter->SetShift(offset);


		std::stringstream outputfile;
        outputfile << output_prefix << i << ".png";
		save(ImageGrayd(cyclicShiftFilter->GetOutput()),outputfile.str(), 0.0);

    }
	*/

    return EXIT_SUCCESS;
}

int testAnimation2(std::string inputfile, std::string output_prefix)
{
    // LOAD
    ImageGrayd input;
    double mean = load(input,inputfile, true);

    // COMPUTE FFT
    ImageGrayd modulus, phase;

    welch(input, modulus, 10);

    ImageGrayd::ItkImg::RegionType phaseRegion (modulus.itk()->GetLargestPossibleRegion());
    phaseOfRegion(input, phaseRegion,phase);

    // INVERSE FFT
    ImageGrayd result;
    fftInverseModulusAndPhase(modulus, phase, result);

    for (int i=100; i<228; ++i)
    {
        ImageGrayd::ItkImg::OffsetType offset;
        offset[0]=1*(i-99);
        offset[1]=-1*(i-99);

        typedef itk::CyclicShiftImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > CyclicShiftFilterType;
        CyclicShiftFilterType::Pointer cyclicShiftFilter = CyclicShiftFilterType::New();
        cyclicShiftFilter->SetInput(result.itk());
        cyclicShiftFilter->SetShift(offset);


        std::stringstream outputfile;
        outputfile << output_prefix << i << ".png";

        save(ImageGrayd(cyclicShiftFilter->GetOutput()),outputfile.str(), 0.0);
    }

    return EXIT_SUCCESS;
}

/**
 * \brief tests phase :
*/
int testPhase(std::string inputfile, std::string outputfile)
{
	// LOAD
	ImageGrayd input;
	double mean = load(input,inputfile, true);

    image_collector collec;
    collec.add(input, mean);

	// COMPUTE IMAGE STATS
//	monitorStats(input, "Input (space)");

    // COMPUTE FFT
	ImageGrayd modulus, phase;


	welch(input, modulus, 10);
    collec.add(modulus);

	ImageGrayd::ItkImg::RegionType phaseRegion (modulus.itk()->GetLargestPossibleRegion());
	phaseOfRegion(input, phaseRegion,phase);
//	monitorStats(phase, "phase of input top left block");
//	std::cout << phaseRegion << std::endl;

	// INVERSE FFT
	ImageGrayd result;
	fftInverseModulusAndPhase(modulus, phase, result);
    collec.add(result,mean);

//    monitorStats(phase, "phase 1");
//    monitorStats(modulus, "modulus 1");
//    monitorStats(result, "result 1");

    // RANDOM PHASE
/*    randomPhaseShift(phase); */
/*    randomPhase(phase); */

//    const int nTest = 7;
//    const bool activateShift [nTest] = {true,true,true,true,true,true,true};
////	const int stepSizes [nTest] = {4,8,16,32,64,128,256};
//    const int stepSizes [nTest] = {256,128,64,32,16,8,4};

    const int nTest = 14;
    const bool activateShift [nTest] = {false,true,false,true,false,true,false,true,false,true,false,true,false,true,};
    const int stepSizes [nTest] = {256,256,128,128,64,64,32,32,16,16,8,8,4,4};


    for(int t= 0; t<nTest; ++t)
	{
        phaseAverage_UsingNormalizedComplexAverage(input,phase,activateShift[t],stepSizes[t]);
		std::stringstream title;
		title << "phaseAverage (activate shift = " << activateShift[t] << ", step = " << stepSizes[t] << ")";
        monitorStats(phase, title.str());
		fftInverseModulusAndPhase(modulus, phase, result);
        collec.add(result,mean);
        collec.add(phase,0.0);
//        monitorStats(phase, "phase 2");
//        monitorStats(modulus, "modulus 2");
//        monitorStats(result, "result 2");
    }

	// SAVE
	save(collec.collect(),outputfile, 0.0);

	return EXIT_SUCCESS;
}

/**
 * \brief tests welch's method : different parameters.
*/
int testWelch(std::string inputfile, std::string outputfile)
{
	// LOAD
	ImageGrayd input;
	double mean = load(input,inputfile, true);

	image_collector collec;
	collec.add(input, mean);

	// COMPUTE IMAGE STATS
	monitorStats(input, "Input (space)");

	// COMPUTE FFT
	ImageGrayd modulus, phase;


	welch(input, modulus, 5);
	collec.add(modulus);

	ImageGrayd::ItkImg::RegionType phaseRegion (modulus.itk()->GetLargestPossibleRegion());
	phaseOfRegion(input, phaseRegion,phase);
//    phaseAverage(input,phase);

	// INVERSE FFT
	ImageGrayd result;
	fftInverseModulusAndPhase(modulus, phase, result);
	collec.add(result,mean);

	// RANDOM PHASE
//    randomPhaseShift(phase);
//    randomPhase(phase);

	welch(input, modulus, 10);
	collec.add(modulus);
	// INVERSE FFT
	fftInverseModulusAndPhase(modulus, phase, result);
	collec.add(result,mean);

	welch(input, modulus, 20);
	collec.add(modulus);
	// INVERSE FFT
	fftInverseModulusAndPhase(modulus, phase, result);
	collec.add(result,mean);

	welch(input, modulus, 50);
	collec.add(modulus);
	// INVERSE FFT
	fftInverseModulusAndPhase(modulus, phase, result);
	collec.add(result,mean);

	monitorStats(result, "Output (space) random");

	// SAVE
	save(collec.collect(),outputfile, 0.0);

	return EXIT_SUCCESS;
}

/**
 * \brief image x linear function
*/
int test3(std::string inputfile, std::string outputfile)
{
    // LOAD
	ImageGrayd input;
	double mean = load(input,inputfile, true);

	int W = input.width();
	int H = input.height();

    for (int i=0; i<W; ++i)
    {
        for (int j=0; j<H; ++j)
        {
			input.pixelAbsolute(i,j) *= (1 - i / double(W));
        }
    }

    // SAVE
	save(input,outputfile, mean);

    return EXIT_SUCCESS;
}

/**
 * \brief basic test with random phase
*/
int test2(std::string inputfile, std::string outputfile)
{
    // LOAD
	ImageGrayd input;
	load(input,inputfile);

    // COMPUTE FFT
	ImageGrayd modulus, phase;

    welch(input, modulus);

//    phaseAverage(input,phase, true, 5);
//    monitorStats(phase, "true, 5");

//    phaseAverage(input,phase, true, 10);
//    monitorStats(phase, "true, 10");

    phaseAverage_UsingNormalizedComplexAverage(input,phase, true, 20);
    monitorStats(phase, "true, 10");

    // INVERSE FFT
	ImageGrayd result;
    fftInverseModulusAndPhase(modulus, phase, result);

    // SAVE
	save(result, outputfile);

    return EXIT_SUCCESS;
}

/**
 * \brief hand made 4x4 image
*/
int test1()
{
//    const unsigned int W = 4;
//    const unsigned int H = 4;

	ImageGrayd::ItkImg::RegionType region;
	ImageGrayd::ItkImg::IndexType regionIndex;
	ImageGrayd::ItkImg::SizeType regionSize;
    regionSize[0] = 4;
    regionSize[1] = 4;
    region.SetSize(regionSize);
    regionIndex[0] = 0;
    regionIndex[1] = 0;
    region.SetIndex(regionIndex);

	ImageGrayd image;
	image.initItk(4,4,true);
//    ImageGrayd::Pointer image = ImageGrayd::New();
//    image->SetRegions(region);
//    image->Allocate(true); // param to true initializes the pixels to 0

	image.pixelAbsolute(0,0) = 0.0;
	image.pixelAbsolute(0,1) = 0.1;
	image.pixelAbsolute(0,2) = 0.2;
	image.pixelAbsolute(0,3) = 0.3;
	image.pixelAbsolute(1,0) = 1.0;
	image.pixelAbsolute(1,1) = 1.1;
	image.pixelAbsolute(1,2) = 1.2;
	image.pixelAbsolute(1,3) = 1.3;
	image.pixelAbsolute(2,0) = 2.0;
	image.pixelAbsolute(2,1) = 2.1;
	image.pixelAbsolute(2,2) = 2.2;
	image.pixelAbsolute(2,3) = 2.3;
	image.pixelAbsolute(3,0) = 3.0;
	image.pixelAbsolute(3,1) = 3.1;
	image.pixelAbsolute(3,2) = 3.2;
	image.pixelAbsolute(3,3) = 3.3;

//    std::cout << image->GetLargestPossibleRegion() << std::endl;
//    std::cout << image.pixelAbsolute(0,0) << std::endl;

//    regionSize[0] = 2;
//    regionSize[1] = 2;
//    region.SetSize(regionSize);

//    ImageGrayd::Pointer accu = ImageGrayd::New();
//    accu->SetRegions(region);
//    accu->Allocate(true); // param to true initializes the pixels to 0

//    typedef itk::RegionOfInterestImageFilter< ImageGrayd, ImageGrayd > ROIFilterType;
//    ROIFilterType::Pointer roiFilter = ROIFilterType::New();
//    roiFilter->SetInput(image);

//    typedef itk::AddImageFilter<ImageGrayd> AddFilterType;
//    AddFilterType::Pointer addFilter = AddFilterType::New();

//    addFilter->SetInput1(accu);
//    addFilter->SetInput2(roiFilter->GetOutput());


//    regionIndex[0] = 1;
//    regionIndex[1] = 2;
//    region.SetIndex(regionIndex);
//    roiFilter->SetRegionOfInterest(region);

//    roiFilter->Update();
//    roiFilter->GetOutput()->SetOrigin(accu->GetOrigin());

//    addFilter->Update();
//    accu=addFilter->GetOutput();
//    std::cout << accu.pixelAbsolute(0,0) << std::endl;

//    addFilter->SetInput1(accu);
//    addFilter->SetInput2(roiFilter->GetOutput());

//    regionIndex[0] = 2;
//    regionIndex[1] = 0;
//    region.SetIndex(regionIndex);
//    roiFilter->SetRegionOfInterest(region);

//    roiFilter->Update();
//    roiFilter->GetOutput()->SetOrigin(accu->GetOrigin());

//    addFilter->Update();
//    accu=addFilter->GetOutput();
//    std::cout << accu.pixelAbsolute(0,0) << std::endl;

    //    std::cout << roiFilter->GetOutput()->GetLargestPossibleRegion() << std::endl;
    //    std::cout << roiFilter->GetOutput()->GetOrigin() << std::endl;

//    typedef itk::ExtractImageFilter< ImageGrayd, ImageGrayd > ExtractFilterType;
//    ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
//    extractFilter->SetInput(image);
//    extractFilter->SetExtractionRegion(region);

//    extractFilter->Update();
//    std::cout << extractFilter->GetOutput()->GetLargestPossibleRegion() << std::endl;
//    std::cout << extractFilter->GetOutput()->GetOrigin() << std::endl;

    //    roiFilter->Update();
//    std::cout << roiFilter->GetOutput()->GetLargestPossibleRegion() << std::endl;
//    std::cout << roiFilter->GetOutput()->pixel(0,0) << std::endl;

    return EXIT_SUCCESS;
}

