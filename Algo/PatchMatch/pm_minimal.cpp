
//-------------------------------------------------------------------------
//Minimal (unoptimized) example of PatchMatch. Requires that ImageMagick be installed.

//To improve generality you can:
// - Use whichever distance function you want in dist(), e.g. compare SIFT descriptors computed densely.
// - Search over a larger search space, such as rotating+scaling patches (see MATLAB mex for examples of both)
//
//To improve speed you can:
// - Turn on optimizations (/Ox /Oi /Oy /fp:fast or -O6 -s -ffast-math -fomit-frame-pointer -fstrength-reduce -msse2 -funroll-loops)
// - Use the MATLAB mex which is already tuned for speed
// - Use multiple cores, tiling the input. See our publication "The Generalized PatchMatch Correspondence Algorithm"
// - Tune the distance computation: manually unroll loops for each patch size, use SSE instructions (see readme)
// - Precompute random search samples (to avoid using rand, and mod)
// - Move to the GPU
//--------------------------------------------------------------------------


#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include <ASTex/image_rgb.h>
#include <ASTex/colorspace_filters.h>
#include <itkCastImageFilter.h>

using namespace ASTex;

using itkRGBu8 = itk::RGBPixel<uint8_t> ;

typedef ImageRGBd::ItkImg IMG_DBL;
typedef ImageRGBu8::ItkImg IMG_U8;

// For conversion to LAB
void getLABfFromRGBu8(const ImageRGBu8& image_rgbu8, ImageRGBd& image_lab)
{
	// first transform [0,255] double -> [0,1] double
	ColorSpace::FilterRGB255To01<IMG_U8, IMG_DBL>::Pointer filter0 = ColorSpace::FilterRGB255To01<IMG_U8, IMG_DBL>::New();
	filter0->SetInput(image_rgbu8.itk());
	
	// RGB double -> XYZ double
	ColorSpace::FilterRGBtoXYZ<IMG_DBL, IMG_DBL>::Pointer filter1 =
	ColorSpace::FilterRGBtoXYZ<IMG_DBL, IMG_DBL>::New();
	filter1->SetInput(filter0->GetOutput());
	
	// XYZ DBL -> LAB DBL
	ColorSpace::FilterXYZtoLAB<IMG_DBL, IMG_DBL>::Pointer filter2 =
	ColorSpace::FilterXYZtoLAB<IMG_DBL, IMG_DBL>::New();
	filter2->SetInput(filter1->GetOutput());

	ImageRGBd output(filter2->GetOutput());
	image_lab = output;
}

// RGB to L*a*b*
// input : r,g,b in range [0,255] (sRGB colorspace)
// ouput : double L,A,B in range [0,100],[-150,+150],[-150,+150] (CIELAB colorspace)
void RGBu8toLAB(const ImageRGBu8& image, int32_t x, int32_t y, double& L, double& A, double& B)
{
	itkRGBu8 c = image.pixelAbsolute(x, y);
	double r = c.GetRed() / 255.0;
	double g = c.GetGreen() / 255.0;
	double b = c.GetBlue() / 255.0;

	// Convert sRGB (r,g,b) to linear-rgb (lr,lg,lb)
	float lr = (r <= 0.04045) ? r / 12.92 : pow((r + 0.055) / 1.055, 2.4);
	float lg = (g <= 0.04045) ? g / 12.92 : pow((g + 0.055) / 1.055, 2.4);
	float lb = (b <= 0.04045) ? b / 12.92 : pow((b + 0.055) / 1.055, 2.4);

	// Convert to XYZ (assuming sRGB was D65)
	double xx = lr*0.4124564 + lg*0.3575761 + lb*0.1804375;
	double yy = lr*0.2126729 + lg*0.7151522 + lb*0.0721750;
	double zz = lr*0.0193339 + lg*0.1191920 + lb*0.9503041;

	// Rescale X/Y/Z relative to white point D65
	double xr = xx / 0.95047;
	double yr = yy / 1.0;
	double zr = zz / 1.08883;

	// Tristimulus function
	double eps = 216.0 / 24389.0, k = 24389.0 / 27.0;
	double fx = (xr <= eps) ? (k * xr + 16.0) / 116.0 : pow(xr, 1 / 3.0);
	double fy = (yr <= eps) ? (k * yr + 16.0) / 116.0 : pow(yr, 1 / 3.0);
	double fz = (zr <= eps) ? (k * zr + 16.0) / 116.0 : pow(zr, 1 / 3.0);

	// Tranform to LAB
	L = (116.0*fy) - 16.0;
	A = 500.0*(fx - fy);
	B = 200.0*(fy - fz);
}

// -------------------------------------------------------------------------
// PatchMatch, using L2 distance between upright patches that translate only
// -------------------------------------------------------------------------

//const int32_t patch_w = 7;
//const int32_t patch_w = 5;
const int32_t patch_w = 5;
//const int32_t patch_w = 30;
const double sqrt3 = 1.7320508075688772935274463415059;
const double max_dist = patch_w * patch_w * sqrt3;

//int32_t pm_iters = 2;
//int32_t pm_iters = 5;
//int32_t pm_iters = 10;
int32_t pm_iters = 20;

// Measure distance between 2 patches with upper left corners (ax, ay) and (bx, by), terminating early if we exceed a cutoff distance.
// You could implement your own descriptor here.
double dist(const ImageRGBu8& a, const ImageRGBu8& b, int32_t ax, int32_t ay, int32_t bx, int32_t by, double cutoff = max_dist) {

	double ans = 0.0;

	for (int32_t dy = 0; dy < patch_w; dy++) {
		for (int32_t dx = 0; dx < patch_w; dx++) {
			itkRGBu8 ac = a.pixelAbsolute(ax + dx, ay + dy);
			double ar = ac.GetRed() / 255.0;
			double ag = ac.GetGreen() / 255.0;
			double ab = ac.GetBlue() / 255.0;
			itkRGBu8 bc = b.pixelAbsolute(bx + dx, by + dy);
			double br = bc.GetRed() / 255.0;
			double bg = bc.GetGreen() / 255.0;
			double bb = bc.GetBlue() / 255.0;
			double dr = br - ar;
			double dg = bg - ag;
			double db = bb - ab;
			ans += sqrt(dr*dr + dg*dg + db*db);
		}

		if (ans >= cutoff) { return cutoff; }
	}

	return ans;
}

double distPatchSize(const ImageRGBu8& a, const ImageRGBu8& b, int32_t ax, int32_t ay, int32_t bx, int32_t by, int32_t sx, int32_t sy, double cutoff = max_dist) {

	double ans = 0.0;

	for (int32_t dy = 0; dy < sy; dy++) {
		for (int32_t dx = 0; dx < sx; dx++) {
			itkRGBu8 ac = a.pixelAbsolute(ax + dx, ay + dy);
			double ar = ac.GetRed() / 255.0;
			double ag = ac.GetGreen() / 255.0;
			double ab = ac.GetBlue() / 255.0;
			itkRGBu8 bc = b.pixelAbsolute(bx + dx, by + dy);
			double br = bc.GetRed() / 255.0;
			double bg = bc.GetGreen() / 255.0;
			double bb = bc.GetBlue() / 255.0;
			double dr = br - ar;
			double dg = bg - ag;
			double db = bb - ab;
			ans += sqrt(dr*dr + dg*dg + db*db);
		}

		if (ans >= cutoff) { return cutoff; }
	}

	return ans;
}

RGBu8 avcolor(const ImageRGBu8& a, int32_t ax, int32_t ay) {

	double avr = 0.0;
	double avg = 0.0;
	double avb = 0.0;

	for (int32_t dy = 0; dy < patch_w; dy++) {
		for (int32_t dx = 0; dx < patch_w; dx++) {
			itkRGBu8 ac = a.pixelAbsolute(ax + dx, ay + dy);
			double ar = ac.GetRed();
			double ag = ac.GetGreen();
			double ab = ac.GetBlue();
			avr += ar;
			avg += ag;
			avb += ab;
		}
	}

	avr /= patch_w * patch_w;
	avg /= patch_w * patch_w;
	avb /= patch_w * patch_w;

	RGBu8 av((uint8_t)avr, (uint8_t)avg, (uint8_t)avb);

	return av;
}

void improve_guess(const ImageRGBu8& a, const ImageRGBu8& b, int32_t ax, int32_t ay, uint32_t& xbest, uint32_t& ybest, double& dbest, uint32_t bx, uint32_t by) {
	double d = dist(a, b, ax, ay, bx, by, dbest);
	if (d < dbest) {
		dbest = d;
		xbest = bx;
		ybest = by;
	}
}

// Match image a to image b, returning the nearest neighbor field mapping a => b coords, stored in an RGB 24-bit image.
void patchmatch(const ImageRGBu8& a, const ImageRGBu8& b, ImageRGBu8& ann, ImageRGBu8& annd,
	                                                      ImageRGBu8& subst, const std::string& synthName) {

	uint32_t* p = new uint32_t[2*a.width()*a.height()];
	double* d = new double[a.width()*a.height()];
	
	// Initialize with random nearest neighbor field (NNF).

	// Effective width and height (possible upper left corners of patches).

	uint32_t aew = a.width() - patch_w + 1, aeh = a.height() - patch_w + 1;
	uint32_t bew = b.width() - patch_w + 1, beh = b.height() - patch_w + 1;

	for (uint32_t ay = 0; ay < aeh; ay++) {
		for (uint32_t ax = 0; ax < aew; ax++) {
			uint32_t bx = std::rand() % bew;
			uint32_t by = std::rand() % beh;
			p[2 * (a.width()*ay + ax)    ] = bx;
			p[2 * (a.width()*ay + ax) + 1] = by;
			d[a.width()*ay + ax] = dist(a, b, ax, ay, bx, by);
			uint8_t ds = (uint8_t)((d[a.width()*ay + ax] / max_dist) * 255.0);
			annd.pixelAbsolute(ax, ay) = RGBu8(ds, 255 - ds, 0);
		}
	}

	for (uint32_t ay = 0; ay < aeh; ay++) {
		for (uint32_t ax = 0; ax < aew; ax++) {
			uint32_t bx = p[2 * (a.width()*ay + ax)];
			uint32_t by = p[2 * (a.width()*ay + ax) + 1];
			uint8_t cr = (uint8_t)((bx * 255.0) / a.width());
			uint8_t cg = (uint8_t)((by * 255.0) / a.height());
			uint8_t cb = 0.0;
			ann.pixelAbsolute(ax, ay) = RGBu8(cr, cg, cb);
		}
	}
	
	std::ostringstream initNameANN, initNameANND;
	initNameANN << synthName.c_str() << "_ann_init.png";
	initNameANND << synthName.c_str() << "_annd_init.png";
	std::cout << initNameANN.str() << std::endl;
	ann.save(initNameANN.str());
	std::cout << initNameANND.str() << std::endl;
	annd.save(initNameANND.str());
	
	for (int32_t iter = 0; iter < pm_iters; iter++) {

		std::cout << "Iteration " << iter << std::endl;

		// In each iteration, improve the NNF, by looping in scanline or reverse-scanline order.

		int32_t ystart = 0, yend = aeh, ychange = 1;
		int32_t xstart = 0, xend = aew, xchange = 1;

		if (iter % 2 == 1) {
			xstart = xend - 1; xend = -1; xchange = -1;
			ystart = yend - 1; yend = -1; ychange = -1;
		}

		for (int32_t ay = ystart; ay != yend; ay += ychange) {
			for (int32_t ax = xstart; ax != xend; ax += xchange) {

				// Current (best) guess.

				uint32_t xbest = p[2 * (a.width()*ay + ax)    ];
				uint32_t ybest = p[2 * (a.width()*ay + ax) + 1];
				double dbest = d[a.width()*ay + ax];

				// Propagation: Improve current guess by trying instead correspondences from left and above (below and right on odd iterations).

				if ((uint32_t)(ax - xchange) < (uint32_t)aew) {

					uint32_t xp = p[2 * (a.width()*ay + (ax - xchange))    ] + xchange;
					uint32_t yp = p[2 * (a.width()*ay + (ax - xchange)) + 1];
					if (xp < bew) {
						improve_guess(a, b, ax, ay, xbest, ybest, dbest, xp, yp);
					}
				}

				if ((uint32_t)(ay - ychange) < (uint32_t)aeh) {

					uint32_t xp = p[2 * (a.width()*(ay - ychange) + ax)];
					uint32_t yp = p[2 * (a.width()*(ay - ychange) + ax) + 1] + ychange;
					if (yp < beh) {
						improve_guess(a, b, ax, ay, xbest, ybest, dbest, xp, yp);
					}
				}

				// Random search: Improve current guess by searching in boxes of exponentially decreasing size around the current best guess.

				int32_t rs_start = std::max(b.width(), b.height());
				for (int32_t mag = rs_start; mag >= 1; mag /= 2) {

					// Sampling window

					int32_t xmin = std::max((xbest - mag), 0u);
					int32_t xmax = std::min(xbest + mag + 1, bew);
					int32_t ymin = std::max((ybest - mag), 0u);
					int32_t ymax = std::min(ybest + mag + 1, beh);
					int32_t xp = xmin + std::rand() % (xmax - xmin);
					int32_t yp = ymin + std::rand() % (ymax - ymin);
					improve_guess(a, b, ax, ay, xbest, ybest, dbest, xp, yp);
				}

				p[2 * (a.width()*ay + ax)    ] = xbest;
				p[2 * (a.width()*ay + ax) + 1] = ybest;
				d[a.width()*ay + ax] = dbest;
				uint8_t ds = (uint8_t)((dbest / max_dist) * 255.0);
				annd.pixelAbsolute(ax, ay) = RGBu8(ds, 255 - ds, 0);
			}
		}

		for (uint32_t ay = 0; ay < aeh; ay++) {
			for (uint32_t ax = 0; ax < aew; ax++) {
				uint32_t bx = p[2 * (a.width()*ay + ax)    ];
				uint32_t by = p[2 * (a.width()*ay + ax) + 1];
				uint8_t cr = (uint8_t)((bx * 255.0) / a.width());
				uint8_t cg = (uint8_t)((by * 255.0) / a.height());
				uint8_t cb = 0.0;
				ann.pixelAbsolute(ax, ay) = RGBu8(cr, cg, cb);
			}
		}

		std::ostringstream resultNameANN, resultNameANND;
		resultNameANN << synthName.c_str() << "_ann_" << std::setfill('0') << std::setw(3) << iter << ".png";
		resultNameANND << synthName.c_str() << "_annd_" << std::setfill('0') << std::setw(3) << iter << ".png";
		std::cout << resultNameANN.str() << std::endl;
		ann.save(resultNameANN.str());
		std::cout << resultNameANND.str() << std::endl;
		annd.save(resultNameANND.str());
	}

	for (uint32_t ay = 0; ay < aeh; ay++) {
		for (uint32_t ax = 0; ax < aew; ax++) {
			uint32_t bx = p[2 * (a.width()*ay + ax)];
			uint32_t by = p[2 * (a.width()*ay + ax) + 1];
			//subst.pixelAbsolute(ax, ay) = b.pixelAbsolute(bx, by);
			subst.pixelAbsolute(ax, ay) = avcolor(b, bx, by);
		}
	}

	std::ostringstream substName;
	substName << synthName.c_str() << "_subst.png";
	subst.save(substName.str());

	delete[] p;
	delete[] d;
}

int32_t main(int32_t argc, char *argv[]) {

	time_t t = time(NULL);
	std::srand(t);
	std::cout << "Random number generator seed: " << t << std::endl;
	
	ImageRGBu8 inputImage;
	ImageRGBu8 synthImage;

	if (argc < 3)
	{
		std::cerr << "Usage: test_patchmatch input_image synth_image\n";
		return EXIT_FAILURE;
	}

	std::string inputPath(argv[1]);
	std::string synthPath(argv[2]);

	typedef itk::ImageFileReader< ImageRGBu8::ItkImg > ReaderLocalType;

	// Load input image

	//inputImage.load(inputPath);

	ReaderLocalType::Pointer ireader = ReaderLocalType::New();
	ireader->SetFileName(inputPath);
	try
	{
		ireader->Update();
	}
	catch (itk::ExceptionObject & error)
	{
		std::cerr << "Failed to load file " << inputPath << std::endl;
		std::cerr << error << std::endl;
		return EXIT_FAILURE;
	}
	inputImage.itk() = ireader->GetOutput();

	int32_t input_w = inputImage.width();
	int32_t input_h = inputImage.height();

	std::cout << "Input image resolution: " << input_w << "*" << input_h << std::endl;

	size_t sep = inputPath.find_last_of("\\/");
	if (sep != std::string::npos)
		inputPath = inputPath.substr(sep + 1, inputPath.size() - sep - 1);

	std::string inputName;
	std::string inputExt;

	size_t dot = inputPath.find_last_of(".");
	if (dot != std::string::npos)
	{
		inputName = inputPath.substr(0, dot);
		inputExt = inputPath.substr(dot, inputPath.size() - dot);
	}
	else
	{
		inputName = inputPath;
		inputExt = "";
	}

	// Load synthesized image

	//synthImage.load(synthPath);

	ReaderLocalType::Pointer sreader = ReaderLocalType::New();
	sreader->SetFileName(synthPath);
	try
	{
		sreader->Update();
	}
	catch (itk::ExceptionObject & error)
	{
		std::cerr << "Failed to load file " << synthPath << std::endl;
		std::cerr << error << std::endl;
		return EXIT_FAILURE;
	}
	synthImage.itk() = sreader->GetOutput();

	int32_t synth_w = synthImage.width();
	int32_t synth_h = synthImage.height();

	std::cout << "Synthesized image resolution: " << synth_w << "*" << synth_h << std::endl;

	sep = synthPath.find_last_of("\\/");
	if (sep != std::string::npos)
		synthPath = synthPath.substr(sep + 1, synthPath.size() - sep - 1);

	std::string synthName;
	std::string synthExt;

	dot = synthPath.find_last_of(".");
	if (dot != std::string::npos)
	{
		synthName = synthPath.substr(0, dot);
		synthExt = synthPath.substr(dot, synthPath.size() - dot);
	}
	else
	{
		synthName = synthPath;
		synthExt = "";
	}

	std::cout << "Input image name (B): " << inputName << std::endl;
	std::cout << "Synthesized image name (A): " << synthName << std::endl;

	//std::string outputNameANN(synthName);
	//outputNameANN.append("_ann_");
	//std::string outputNameANND(synthName);
	//outputNameANND.append("_annd_");

	ImageRGBu8 annImage;
	annImage.initItk(synthImage.width() - patch_w + 1, synthImage.height() - patch_w + 1, true);
	//annImage.initItk(synthImage.width(), synthImage.height(), true);

	ImageRGBu8 anndImage;
	anndImage.initItk(synthImage.width() - patch_w + 1, synthImage.height() - patch_w + 1, true);
	//anndImage.initItk(synthImage.width(), synthImage.height(), true);

	ImageRGBu8 substImage;
	substImage.initItk(synthImage.width() - patch_w + 1, synthImage.height() - patch_w + 1, true);

	std::cout << "Running PatchMatch\n";

	patchmatch(synthImage, inputImage, annImage, anndImage, substImage, synthName);

	//annImage.save(outputNameANN);

	//anndImage.save(outputNameANND);

	/*
	uint32_t bx = std::rand() % 1600;
	uint32_t by = std::rand() % 1600;

	std::cout << bx << " " << by << std::endl;

	uint8_t cr = ((bx << 20) >> 20) >> 4;
	uint8_t cg = (((bx << 28) >> 28) << 4) | (((by << 20) >> 20) >> 8);
	uint8_t cb = (by << 24) >> 24;

	std::cout << (uint32_t)cr << " " << (uint32_t)cg << " " << (uint32_t)cb << std::endl;

	bx = (cr << 4) | (cg >> 4);
	by = (((cg << 28) >> 28) << 8) | cb;

	std::cout << bx << " " << by << std::endl;
	*/

	/*
	argc--;
	argv++;
	if (argc != 4) { fprintf(stderr, "pm_minimal a b ann annd\n"
	"Given input images a, b outputs nearest neighbor field 'ann' mapping a => b coords, and the squared L2 distance 'annd'\n"
	"These are stored as RGB 24-bit images, with a 24-bit int32_t at every pixel. For the NNF we store (by<<12)|bx."); exit(1); }
	printf("Loading input images\n");
	BITMAP *a = load_bitmap(argv[0]);
	BITMAP *b = load_bitmap(argv[1]);
	BITMAP *ann = NULL, *annd = NULL;
	printf("Running PatchMatch\n");
	patchmatch(a, b, ann, annd);
	printf("Saving output images\n");
	save_bitmap(ann, argv[2]);
	save_bitmap(annd, argv[3]);
	*/

	return 0;
}
