// hvProcWPhase.h: interface for the procedural texture class with Phase.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PROCEDURALWPHASE_H__B6AC0A32_75EF_428E_BC10_6219F619FA29__INCLUDED_)
#define AFX_PROCEDURALWPHASE_H__B6AC0A32_75EF_428E_BC10_6219F619FA29__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

////////////////////////////////////////////
/*

*/ 
////////////////////////////////////////////
////////////////////////////////////////////

#include <math.h>
#include "hvNoise.h"
#include "hvBitmap.h"
#include "hvPicture.h"
#include "../hv_astex_io.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif 

namespace hview {

// precision: size = 2^pow_2
//const int pow_2=3;			// size of spatial subdivision (a square) for subFFT, from 5-> 32x32 (precise) to 7-> 128x128 (low precision)
//const int size = (2<<pow_2)/2;	    // resolution of spatial subdivision (a square) in pixels
const int NTURB_STRAT=3;
const int NTURB_COS = 4;
const bool mirrored=true;
const double pcafactor = 0.7;  // for PCA
const int maxstrates = 4;	// amount of noise bands for phase random part
const int NRNDCOS = 10;
//const double lumweight=1.0; // for color quantization between 0.1 and 1.0
const int MAXCOLORS = 16;		// max number of color features
const int WELCH=50;
const hvColRGB<unsigned char> stratcoltab[8] = 
{ hvColRGB<double>(255,255,140), hvColRGB<double>(255,200,60), 
hvColRGB<double>(255,100,100), hvColRGB<double>(120,255,120), hvColRGB<double>(60,200,200), 
hvColRGB<double>(100,100,255), hvColRGB<double>(0,200,200), hvColRGB<double>(30,30,150) }; 
const double fftdrawlog = 50000.0;
const double BILAT=60000.0;

class hvProcWPhase {

char tname[512];
bool exportc1, exportc2, exportc3;
int NCOLORS;
int NREGCOS;
int pow_2;		// size of spatial subdivision (a square) for subFFT, from 5-> 32x32 (precise) to 7-> 128x128 (low precision)
int size;	    // resolution of spatial subdivision (a square) in pixelsint ssx, ssy;
int ssx, ssy;
int maxsize;
int maxpow;

hvPictRGB<unsigned char> pexample;
hvQPictRGB<unsigned char,64> piquantized;

hvPictRGB<unsigned char> pimlr, pvar, piquant, exfeature, colfeature;
hvPictRGB<unsigned char> fexample, fexfeature, fpvar;
//hvPictRGB<unsigned char> pstratif, pstratifc;
hvPict<double> meanenrg, meanenrgv, meanenrgv2;

int histofeature[256], histovar[256], histovar2[256];
int shistofeature[512], shistovar[512], shistovar2[512];
hvArray1<unsigned char> transfeature, transvar, transvar2;

hvPictRGB<unsigned char> gausstrat[3][32];
hvPict<double> gaussapprox[3];

std::vector<hvFrame3<double> > lfr;
hvVec3<double>  meancol[MAXCOLORS], principaldir[MAXCOLORS], principal2dir[MAXCOLORS];
double meancolvalue[MAXCOLORS];
int featcount[MAXCOLORS];
int featappcount[1024];
double indexfeat[MAXCOLORS+5];
double savg, smin, smax, savgmin, savgmax; 
double savgv, sminv, smaxv, savgminv, savgmaxv; 
double savgv2, sminv2, smaxv2, savgminv2, savgmaxv2; 
unsigned char ominv, omaxv, ominv2, omaxv2;
//double appsavg, appsmin, appsmax, appsavgmin, appsavgmax; 
//double appsavgv, appsminv, appsmaxv, appsavgminv, appsavgmaxv; 
//double appsavgv2, appsminv2, appsmaxv2, appsavgminv2, appsavgmaxv2; 

class FF {
	int fx, fy;
	double amp, ph;
public:
	FF() { fx=0; fy=0; amp=0.0; ph=0.0; }
	FF(int f, int g, double a, double ff=0.0) { fx=f; fy=g; amp=a; ph=ff; }
	int getfx() const { return fx; }
	int getfy() const { return fy; }
	double geta() const { return amp; }
	double getphi() const { return ph; }
	bool operator<(const FF &f) { return amp<f.amp; }
	bool operator==(const FF &f) { return amp==f.amp; }
};

hvArray1<hvPair<double,double> > *exFT, *exFTvar, *exFTvar2;
hvSortedList<FF> glfreq, glfreqv, glfreqv2;
double genrgtot, genrgtotv, genrgtotv2;
double gamptot, gamptotv, gamptotv2;

hvSortedList<FF> lfreq[4096], lfreqv[4096], lfreqv2[4096];
double enrgtot[4096], enrgtotv[4096], enrgtotv2[4096];
// regular part (array of set of frequencies with phase
double aaf0[4096],aavf0[4096],aav2f0[4096];
double *aa[4096], *phi[4096];	
double *aav[4096], *phiv[4096];
double *aav2[4096], *phiv2[4096];
int *ffx[4096], *ffy[4096];
int *ffxv[4096], *ffyv[4096];
int *ffxv2[4096], *ffyv2[4096];
// random part
double gaaf0, gaavf0, gaav2f0;
double *raa, *rphi;	
double *raav, *rphiv;
double *raav2, *rphiv2;
int *rffx, *rffy;
int *rffxv, *rffyv;
int *rffxv2, *rffyv2;

int maxcosinesi, maxcosinesiv, maxcosinesiv2; 
int maxrndcosines, maxrndcosinesv, maxrndcosinesv2;
int nstrates, nstratesv, nstratesv2;
int ceilind, ceilindv, ceilindv2;
double gauss[32], gaussv[32], gaussv2[32];
double thvar[32], thvarv[32], thvarv2[32];

// Turbulence data
double gausst[NTURB_STRAT];
double taa[NTURB_STRAT][NTURB_COS];
int tffx[NTURB_STRAT][NTURB_COS], tffy[NTURB_STRAT][NTURB_COS];

public:

hvProcWPhase(char *name, int ncol, bool e1=false, bool e2=false, bool e3=false):lfr(100)
{
	char buff[512];
	strcpy(tname, name);
	NCOLORS = ncol;
	exportc1=e1; exportc2=e2; exportc3=e3;
	lfr.clear();

//	printf("Loading file : %s\n", name);
	sprintf(buff,"%s.ppm",name);
	FILE *fdes = fopen(buff, "rb");
	if (fdes==0) perror("cannot load image");
	pexample.loadPPM(fdes, 1);
	fclose(fdes);
//	printf("Resolution: %d, %d\n", pexample.sizeX(), pexample.sizeY());

//	printf("Doing color space...\n");
	this->doColorSpace();
	//this->computeColExtr();
//	printf("Spectral analysis...\n");
	this->doFFT();

	for (int i=0; i<4096; i++)
	{
	aa[i]=0; phi[i]=0;
	aav[i]=0; phiv[i]=0;
	aav2[i]=0; phiv2[i]=0;
	ffx[i]=0; ffy[i]=0;
	ffxv[i]=0; ffyv[i]=0;
	ffxv2[i]=0; ffyv2[i]=0;
	}
	raa=0; rphi=0;	raav=0; rphiv=0; raav2=0; rphiv2=0;
	rffx=0; rffy=0; rffxv=0; rffyv=0; rffxv2=0; rffyv2=0;

	for (int idist=0; idist<NTURB_STRAT; idist++)
	{
		gausst[idist]=1.0;
		for (int nncos=0; nncos<NTURB_COS; nncos++)
		{
			taa[idist][nncos]=0.0;
			tffx[idist][nncos]=0;
			tffy[idist][nncos]=0;
		}
	}

	transfeature.reset(512, 0);
	transvar.reset(512, 0);
	transvar2.reset(512, 0);
}


hvProcWPhase(const std::string& fname, int ncol, bool e1=false, bool e2=false, bool e3=false):lfr(100)
{
//	char buff[512];
//	strcpy(tname, name);
	NCOLORS = ncol;
	exportc1=e1; exportc2=e2; exportc3=e3;
	lfr.clear();

	ASTex::load_to_hv(fname, pexample);
//	printf("Loading file : %s\n", name);
//	sprintf(buff,"%s.ppm",name);
//	FILE *fdes = fopen(buff, "rb");
//	if (fdes==0) perror("cannot load image");
//	pexample.loadPPM(fdes, 1);
//	fclose(fdes);
	//printf("Resolution: %d, %d\n", pexample.sizeX(), pexample.sizeY());

	//printf("Doing color space...\n");
	this->doColorSpace();
	//this->computeColExtr();
	//printf("Spectral analysis...\n");
	this->doFFT();

	for (int i=0; i<4096; i++)
	{
	aa[i]=0; phi[i]=0;
	aav[i]=0; phiv[i]=0;
	aav2[i]=0; phiv2[i]=0;
	ffx[i]=0; ffy[i]=0;
	ffxv[i]=0; ffyv[i]=0;
	ffxv2[i]=0; ffyv2[i]=0;
	}
	raa=0; rphi=0;	raav=0; rphiv=0; raav2=0; rphiv2=0;
	rffx=0; rffy=0; rffxv=0; rffyv=0; rffxv2=0; rffyv2=0;

	for (int idist=0; idist<NTURB_STRAT; idist++)
	{
		gausst[idist]=1.0;
		for (int nncos=0; nncos<NTURB_COS; nncos++)
		{
			taa[idist][nncos]=0.0;
			tffx[idist][nncos]=0;
			tffy[idist][nncos]=0;
		}
	}

	transfeature.reset(512, 0);
	transvar.reset(512, 0);
	transvar2.reset(512, 0);
}


hvProcWPhase(char *name, int nclusters, int mask, 
	bool e1 = false, bool e2 = false, bool e3 = false) :lfr(100)
{
	char buff[512];
	strcpy(tname, name);
	NCOLORS = 1;
	exportc1 = e1; exportc2 = e2; exportc3 = e3;

	//printf("Loading file : %s\n", name);
	sprintf(buff, "%s.ppm", name);
	FILE *fdes = fopen(buff, "rb");
	if (fdes == 0) perror("cannot load image");
	pexample.loadPPM(fdes, 1);
	fclose(fdes);
	//printf("Resolution: %d, %d\n", pexample.sizeX(), pexample.sizeY());

	hvPictRGB<unsigned char> bmask;
	sprintf(buff, "%s_nclusters_%d_binary_mask_%d.ppm", name, nclusters, mask);
	fdes = fopen(buff, "rb");
	if (fdes == 0) perror("cannot load image");
	bmask.loadPPM(fdes, 1);
	fclose(fdes);
	hvPict<unsigned char> bmaskp(bmask, 1, 0, 0, bmask.sizeX() - 1, bmask.sizeY() - 1);
	hvBitmap bm(bmaskp, hvBitmap::GREATER, 1);

	//printf("Doing color space...\n");
	this->doColorSpace(&bm);
	//this->computeColExtr();

	printf("Spectral analysis...loading FT files\n");
	hvPictRGB<unsigned char> fftv;
	sprintf(buff, "%s_spectrum_nclusters_%d_mask_%d_coord0.ppm", name, nclusters, mask);
	fdes = fopen(buff, "rb");
	if (fdes == 0) perror("cannot load image");
	fftv.loadPPM(fdes, 1);
	fclose(fdes);
	hvPict<unsigned char> fftvp(fftv, 1, 0, 0, fftv.sizeX() - 1, fftv.sizeY() - 1);
	hvPictRGB<unsigned char> fftv2;
	sprintf(buff, "%s_spectrum_nclusters_%d_mask_%d_coord1.ppm", name, nclusters, mask);
	fdes = fopen(buff, "rb");
	if (fdes == 0) perror("cannot load image");
	fftv2.loadPPM(fdes, 1);
	fclose(fdes);
	hvPict<unsigned char> fftv2p(fftv2, 1, 0, 0, fftv2.sizeX() - 1, fftv2.sizeY() - 1);

	this->doFFT(&fftvp,&fftv2p);

	for (int i = 0; i<4096; i++)
	{
		aa[i] = 0; phi[i] = 0;
		aav[i] = 0; phiv[i] = 0;
		aav2[i] = 0; phiv2[i] = 0;
		ffx[i] = 0; ffy[i] = 0;
		ffxv[i] = 0; ffyv[i] = 0;
		ffxv2[i] = 0; ffyv2[i] = 0;
	}
	raa = 0; rphi = 0;	raav = 0; rphiv = 0; raav2 = 0; rphiv2 = 0;
	rffx = 0; rffy = 0; rffxv = 0; rffyv = 0; rffxv2 = 0; rffyv2 = 0;

	for (int idist = 0; idist<NTURB_STRAT; idist++)
	{
		gausst[idist] = 1.0;
		for (int nncos = 0; nncos<NTURB_COS; nncos++)
		{
			taa[idist][nncos] = 0.0;
			tffx[idist][nncos] = 0;
			tffy[idist][nncos] = 0;
		}
	}

	transfeature.reset(512, 0);
	transvar.reset(512, 0);
	transvar2.reset(512, 0);
}

int nColors() const { return piquantized.ncolors(); }
int getResolution() const { return maxsize; }
int getPowResolution() const { return maxpow; }
int getStrates() const { return nstrates; }
int getStratesv() const { return nstratesv; }
int getStratesv2() const { return nstratesv2; }

void setThresh(int ncos, double percentage)
{
	NREGCOS = (int)(percentage*(double)ncos);
	this->doStructureSeparation(percentage);
	this->doRegFreqSel(ncos, percentage);
	this->doRndFreqSel(ncos, percentage);
	if (this->nColors()>1) this->doNonLinearIndex(); // color indices
}

static double noise(double px, double py, int sstart, int send, double gauss[], int maxsize, int maxrndcosines, double raa[], int rffx[], int rffy[])
{
	int i, di,dj;
	double vv=0.0; 
	double coeff=0.0;
	for (int idist=sstart; idist<send; idist++)
	{
			double gg = gauss[idist]*1.0; if (gg<1.0/(double)maxsize) gg=1.0/(double)maxsize;
			//for (int di=0; di<=0; di++) for (int dj=0; dj<=0; dj++)
			for (int di=-1; di<=1; di++) for (int dj=-1; dj<=1; dj++)
				{
					double x=(double)((int)floor(px)%int((double)maxsize*gg))-(double)(di*maxsize)*gg + px-floor(px);
					double y=(double)((int)floor(py)%int((double)maxsize*gg))-(double)(dj*maxsize)*gg + py-floor(py);
					hvNoise::seeding((int)(((double)floor(px)+(double)(di*maxsize)*gg)/((double)maxsize*gg)), (int)(((double)floor(py)+(double)(dj*maxsize)*gg)/((double)maxsize*gg)), 0);
					double sum=0.0;
					//printf("%g %g (%d,%d) -> %g,%g\n", px,py,di,dj,x,y);
					for (i=idist*maxrndcosines*NRNDCOS; i<(idist+1)*maxrndcosines*NRNDCOS; i+=NRNDCOS)
					{
						int k=(int)((double)NRNDCOS*(hvNoise::next()+1.0)/2.0);
						if (k<0) k=0; else if (k>=NRNDCOS) k=NRNDCOS-1;
						double fx=(double)(rffx[i+k]-maxsize/2)/(double)maxsize;
						double fy=(double)(rffy[i+k]-maxsize/2)/(double)maxsize;
						sum += raa[i+k]*cos(2.0*M_PI*(x*fx+y*fy)+M_PI*hvNoise::next());
					}
					x=(double)(((int)floor(px))%int(gg*(double)maxsize))-(double)(di*maxsize)*gg + px - floor(px);
					y=(double)(((int)floor(py))%int(gg*(double)maxsize))-(double)(dj*maxsize)*gg + py - floor(py);
					double cx = (x-(double)(maxsize/2)*gg)/((double)(maxsize)*gg);
					double cy = (y-(double)(maxsize/2)*gg)/((double)(maxsize)*gg);
					double factor = hvNoise::Kaiser((cx*cx+cy*cy)/1.5/1.5);
					vv += factor*sum; //*pow(1.0-(double)idist/(double)(send),0.05);
					//printf("%g %g (%d,%d) -> %g,%g %g,%g : %g %g\n", px,py,di,dj,x,y, cx,cy, sum, vv);
				}
	}
	return vv;
}

static double noiseprecomp(double px, double py, int sstart, int send, double gauss[], int maxsize, hvPictRGB<unsigned char> stratpi[])
{
	int i, di, dj;
	double vv = 0.0;
	double coeff = 0.0;
	for (int idist = sstart; idist<send; idist++)
	{
		double gg = gauss[idist] * 2.0; if (gg<1.0 / (double)maxsize) gg = 1.0 / (double)maxsize;
		//for (int di=0; di<=0; di++) for (int dj=0; dj<=0; dj++)
		for (int di = -1; di <= 1; di++) for (int dj = -1; dj <= 1; dj++)
		{
			double x = (double)((int)floor(px) % int((double)maxsize*gg)) - (double)(di*maxsize)*gg + px - floor(px);
			double y = (double)((int)floor(py) % int((double)maxsize*gg)) - (double)(dj*maxsize)*gg + py - floor(py);
			hvNoise::seeding((int)(((double)floor(px) + (double)(di*maxsize)*gg) / ((double)maxsize*gg)), (int)(((double)floor(py) + (double)(dj*maxsize)*gg) / ((double)maxsize*gg)), 0);
			double sum = 0.0;
			x = (double)(((int)floor(px)) % int(gg*(double)maxsize)) - (double)(di*maxsize)*gg + px - floor(px);
			y = (double)(((int)floor(py)) % int(gg*(double)maxsize)) - (double)(dj*maxsize)*gg + py - floor(py);
			double cx = (x - (double)(maxsize / 2)*gg) / ((double)(maxsize)*gg);
			double cy = (y - (double)(maxsize / 2)*gg) / ((double)(maxsize)*gg);
			double factor = hvNoise::Kaiser((cx*cx + cy*cy) / 1.5 / 1.5);

			int ix = int(cx*(double)(maxsize)*gg+(double)stratpi[idist].sizeX()*(hvNoise::next() + 1.0) / 2.0);
			while (ix < 0) ix += stratpi[idist].sizeX(); while (ix >= stratpi[idist].sizeX()) ix -= stratpi[idist].sizeX();
			int iy = int(cy*(double)(maxsize)*gg + (double)stratpi[idist].sizeY()*(hvNoise::next() + 1.0) / 2.0);
			while (iy < 0) iy += stratpi[idist].sizeY(); while (iy >= stratpi[idist].sizeY()) iy -= stratpi[idist].sizeY();
			sum = ((double)(stratpi[idist].get(ix, iy).RED()) - 128.0) / 128.0;
			vv += factor*sum; //*pow(1.0-(double)idist/(double)(send),0.05);
							  //printf("%g %g (%d,%d) -> %g,%g %g,%g : %g %g\n", px,py,di,dj,x,y, cx,cy, sum, vv);
		}
	}
	return vv;
}
static double periodic(double px, double py, int ssx, int ssy, int size, int maxcosinesi, double aaf0[], double *aa[], double *phi[], int *ffx[], int *ffy[])
{
	int i, di,dj;
	double sum=0.0, coeff=0.0;

	for (di=0; di<=1; di++) for (dj=0; dj<=1; dj++)
//	for (di=0; di<=0; di++) for (dj=0; dj<=0; dj++)
	{
		int signx = ((int)floor(px))%size<=size/2 ? -1:1;
		int signy = ((int)floor(py))%size<=size/2 ? -1:1;
		double x=(double)(((int)floor(px))%size)+px-floor(px)-(double)signx*(double)(di*size);
		double y=(double)(((int)floor(py))%size)+py-floor(py)-(double)signy*(double)(dj*size);
		int ix = ((int)floor(px))/size+signx*di;
		ix = (ix%ssx+ssx)%ssx;
		int iy = ((int)floor(py))/size+signy*dj;
		iy = (iy%ssy+ssy)%ssy;
		int fid = ix+(iy)*ssx;
		double vv=aaf0[fid]; //*(1.0+noiseamp*hvNoise::next()/(double)piquantized.ncolors());
		for (i=0; i<maxcosinesi; i++)
		{
			double fx=(double)(ffx[fid][i]-size)/(double)size/2.0;
			double fy=(double)(ffy[fid][i]-size)/(double)size/2.0;
			vv += aa[fid][i]*cos(2.0*M_PI*(x*fx+y*fy)+phi[fid][i]);
		}
		double cx = (x-(double)(size/2))/(double)size;
		double cy = (y-(double)(size/2))/(double)size;
		double factor = hvNoise::Kaiser((cx*cx+cy*cy));
		sum += factor*vv;
		coeff += factor;
	}
	return sum/coeff;
}

void evalReg(double px, double py, double &vv, double &vvar, double &vvar2, bool tiling=false)
{
	const int TSIZE=maxsize/2;
	const double MARGIN=16.0; 
	if (tiling)
	{
		int x = ((int)floor(px))/(TSIZE);
		int y = ((int)floor(py))/(TSIZE);
		hvNoise::seeding(x,y,10);
		double rx = px-double(x*TSIZE);
		double ry = py-double(y*TSIZE);
		int nx = (int)((hvNoise::next()+1.0)/2.0*(double)(maxsize-TSIZE*0-16));
		int ny = (int)((hvNoise::next()+1.0)/2.0*(double)(maxsize-TSIZE*0-16));
		double vx = (double)nx+rx, vy = (double)ny+ry;
		vv = hvProcWPhase::periodic(vx,vy,ssx,ssy, size, maxcosinesi, aaf0, aa, phi, ffx, ffy);
		vvar = hvProcWPhase::periodic(vx,vy,ssx,ssy, size, maxcosinesiv, aavf0, aav, phiv, ffxv, ffyv);
		vvar2 = hvProcWPhase::periodic(vx,vy,ssx,ssy, size, maxcosinesiv2, aav2f0, aav2, phiv2, ffxv2, ffyv2);
		if (rx<(int)MARGIN)
		{
		hvNoise::seeding(x-1,y,10);
		nx = (int)((hvNoise::next()+1.0)/2.0*(double)(maxsize-TSIZE*0-16));
		ny = (int)((hvNoise::next()+1.0)/2.0*(double)(maxsize-TSIZE*0-16));
		vx = (double)nx+rx+(double)(TSIZE); vy=(double)ny+ry;
		vv = vv*(rx/MARGIN)+(1.0-rx/MARGIN)*hvProcWPhase::periodic(vx,vy,ssx,ssy, size, maxcosinesi, aaf0, aa, phi, ffx, ffy);
		vvar = vvar*(rx/MARGIN)+(1.0-rx/MARGIN)*hvProcWPhase::periodic(vx,vy,ssx,ssy, size, maxcosinesiv, aavf0, aav, phiv, ffxv, ffyv);
		vvar2 = vvar2*(rx/MARGIN)+(1.0-rx/MARGIN)*hvProcWPhase::periodic(vx,vy,ssx,ssy, size, maxcosinesiv2, aav2f0, aav2, phiv2, ffxv2, ffyv2);
		}
		if (ry<(int)MARGIN)
		{
		hvNoise::seeding(x,y-1,10);
		nx = (int)((hvNoise::next()+1.0)/2.0*(double)(maxsize-TSIZE*0-16));
		ny = (int)((hvNoise::next()+1.0)/2.0*(double)(maxsize-TSIZE*0-16));
		vx = (double)nx+rx; vy=(double)ny+ry+(double)(TSIZE);
		vv = vv*(ry/MARGIN)+(1.0-ry/MARGIN)*hvProcWPhase::periodic(vx,vy,ssx,ssy, size, maxcosinesi, aaf0, aa, phi, ffx, ffy);
		vvar = vvar*(ry/MARGIN)+(1.0-ry/MARGIN)*hvProcWPhase::periodic(vx,vy,ssx,ssy, size, maxcosinesiv, aavf0, aav, phiv, ffxv, ffyv);
		vvar2 = vvar2*(ry/MARGIN)+(1.0-ry/MARGIN)*hvProcWPhase::periodic(vx,vy,ssx,ssy, size, maxcosinesiv2, aav2f0, aav2, phiv2, ffxv2, ffyv2);
		}

		if (rx<(int)MARGIN && ry<(int)MARGIN)
		{
			hvNoise::seeding(x - 1, y - 1, 10);
			nx = (int)((hvNoise::next() + 1.0) / 2.0*(double)(maxsize - TSIZE * 0 - 16));
			ny = (int)((hvNoise::next() + 1.0) / 2.0*(double)(maxsize - TSIZE * 0 - 16));
			vx = (double)nx + rx + (double)(TSIZE); vy = (double)ny + ry + (double)(TSIZE);
			double rmax = rx; if (ry > rx) rmax = ry;
			vv = vv*(rmax / MARGIN) + (1.0 - rmax / MARGIN)*hvProcWPhase::periodic(vx, vy, ssx, ssy, size, maxcosinesi, aaf0, aa, phi, ffx, ffy);
			vvar = vvar*(rmax / MARGIN) + (1.0 - rmax / MARGIN)*hvProcWPhase::periodic(vx, vy, ssx, ssy, size, maxcosinesiv, aavf0, aav, phiv, ffxv, ffyv);
			vvar2 = vvar2*(rmax / MARGIN) + (1.0 - rmax / MARGIN)*hvProcWPhase::periodic(vx, vy, ssx, ssy, size, maxcosinesiv2, aav2f0, aav2, phiv2, ffxv2, ffyv2);
		}
	}
	else
	{
		vv = hvProcWPhase::periodic(px,py,ssx,ssy, size, maxcosinesi, aaf0, aa, phi, ffx, ffy);
		vvar = hvProcWPhase::periodic(px,py,ssx,ssy, size, maxcosinesiv, aavf0, aav, phiv, ffxv, ffyv);
		vvar2 = hvProcWPhase::periodic(px,py,ssx,ssy, size, maxcosinesiv2, aav2f0, aav2, phiv2, ffxv2, ffyv2);
	}
}

void evalRnd(double px, double py, int sstart, int send, double &vv, bool pre = false)
{
	if (!pre) vv = hvProcWPhase::noise(px,py,sstart, send, gauss, maxsize, maxrndcosines, raa, rffx, rffy);
	else vv = hvProcWPhase::noiseprecomp(px, py, sstart, send, gauss, maxsize, gausstrat[0]);
}
void evalRndv(double px, double py, int sstart, int send, double &vvar, bool pre = false)
{
	if (!pre) vvar = hvProcWPhase::noise(px,py,sstart, send, gaussv, maxsize, maxrndcosinesv, raav, rffxv, rffyv);
	else vvar = hvProcWPhase::noiseprecomp(px, py, sstart, send, gaussv, maxsize, gausstrat[1]);
}

void evalRndv2(double px, double py, int sstart, int send, double &vvar2, bool pre = false)
{
	if (!pre) vvar2 = hvProcWPhase::noise(px,py,sstart, send, gaussv2, maxsize, maxrndcosinesv2, raav2, rffxv2, rffyv2);
	else vvar2 = hvProcWPhase::noiseprecomp(px, py, sstart, send, gaussv2, maxsize, gausstrat[2]);
}


hvColRGB<unsigned char> eval(double px, double py, bool linear=false)
{
	double vv = hvProcWPhase::periodic(px,py,ssx,ssy, size, maxcosinesi, aaf0, aa, phi, ffx, ffy);
	double vvar = hvProcWPhase::periodic(px,py,ssx,ssy, size, maxcosinesiv, aavf0, aav, phiv, ffxv, ffyv);
	double vvar2 = hvProcWPhase::periodic(px,py,ssx,ssy, size, maxcosinesiv2, aav2f0, aav2, phiv2, ffxv2, ffyv2);
	double rvv = hvProcWPhase::noise(px,py, 0, nstrates, gauss, maxsize, maxrndcosines, raa, rffx, rffy);
	double rvvar = hvProcWPhase::noise(px,py,0, nstratesv, gaussv, maxsize, maxrndcosinesv, raav, rffxv, rffyv);
	double rvvar2 = hvProcWPhase::noise(px,py,0, nstratesv2, gaussv2, maxsize, maxrndcosinesv2, raav2, rffxv2, rffyv2);
	//return 	 hvColRGB<unsigned char>(  (unsigned char)((rvvar*0.5+0.5)*255.0), 
	//	(unsigned char)((rvvar*0.5+0.5)*255.0), 
	//	(unsigned char)((rvvar*0.5+0.5)*255.0) );
	return this->toRGB(vv+rvv,vvar+rvvar,vvar2+rvvar2, linear);
}
hvColRGB<unsigned char> evallinear(double px, double py)
{
	return this->eval(px,py,true);
}

double transfervar(double vvar)
{
	vvar = (vvar - sminv) / (smaxv - sminv)*512.0;
	if (vvar < 0.0) vvar = 0.0; else if (vvar >= 511.0) vvar = 511.0;
	return ((((double)transvar[(int)vvar])) - 128.0) / 255.0 / pcafactor;
}
double transfervar2(double vvar)
{
	if (smaxv2 == sminv2) return 0.0;
	vvar = (vvar - sminv2) / (smaxv2 - sminv2)*512.0;
	if (vvar < 0.0) vvar = 0.0; else if (vvar >= 511.0) vvar = 511.0;
	return ((((double)transvar2[(int)vvar])) - 128.0) / 255.0 / pcafactor;
}
hvColRGB<unsigned char> toRGB(double vv, double vvar, double vvar2, bool linear)
{
	int ind; double rest;
	this->convertColIndex(vv,linear, ind,rest, this->nColors(), indexfeat);
	return this->toRGB(ind,rest,vvar,vvar2, meancol,principaldir,principal2dir);
}

static void convertColIndex(double vv, bool linear, int &ind, double &rest, int ncolors, double indexfeat[])
{
	ind=0; rest=0.0;
	if (ncolors>=2)
	{
		if (linear)
		{
			if (vv<-1.0) vv=-1.0;
			if (vv>1.0) vv=1.0;
			vv=(double)(ncolors)*(vv*128.0+128.0)/255.0+0.5;
			ind = (int)(vv);
			//rest = vv-(double)ind;
			// contrast
			//if (rest<0.5) rest=0.5*(rest*rest*rest*rest/0.5/0.5); else rest = 0.5+0.5*sqrt(sqrt((rest-0.5)/0.5));
			if (ind>=ncolors-1) { ind=ncolors-1; rest=0.0; }
		}
		else
		{
			while (vv>indexfeat[ind]) ind++;
			ind-=1;
			//rest=(vv-indexfeat[ind])/(indexfeat[ind+1]-indexfeat[ind]);
			if (ind>=ncolors-1) { ind=ncolors-1; rest=0.0; }
			if (ind<0) { ind=0; rest=0.0; }
			//printf("vv=%g : %d (%g)\n", vv, ind, rest);
			// contrast
			//if (rest<0.5) rest=0.5*(rest*rest*rest*rest/0.5/0.5); else rest = 0.5+0.5*sqrt(sqrt((rest-0.5)/0.5));
			//rest=rest*rest*rest*rest;
		}
	}
}


static hvColRGB<unsigned char> toRGB(int ind, double rest, double vvar, double vvar2, 
	hvVec3<double> meancol[], hvVec3<double> principaldir[], hvVec3<double> principal2dir[])
{
	if (vvar<-1.0) vvar=-1.0;
	if (vvar>1.0) vvar=1.0;
	
	if (vvar2<-1.0) vvar2=-1.0;
	if (vvar2>1.0) vvar2=1.0;

	double mx = meancol[ind].X()*(1.0-rest)+meancol[ind+1].X()*rest;
	double my = meancol[ind].Y()*(1.0-rest)+meancol[ind+1].Y()*rest;
	double mz = meancol[ind].Z()*(1.0-rest)+meancol[ind+1].Z()*rest;
	double ddx = principaldir[ind].X()*(1.0-rest)+principaldir[ind+1].X()*rest;
	double ddy = principaldir[ind].Y()*(1.0-rest)+principaldir[ind+1].Y()*rest;
	double ddz = principaldir[ind].Z()*(1.0-rest)+principaldir[ind+1].Z()*rest;
	double ddx2 = principal2dir[ind].X()*(1.0-rest)+principal2dir[ind+1].X()*rest;
	double ddy2 = principal2dir[ind].Y()*(1.0-rest)+principal2dir[ind+1].Y()*rest;
	double ddz2 = principal2dir[ind].Z()*(1.0-rest)+principal2dir[ind+1].Z()*rest;

	double mmx=(mx+ddx*vvar+ddx2*vvar2)*255.0; if (mmx<0.0) mmx=0.0; else if (mmx>255.0) mmx=255.0;
	double mmy=(my+ddy*vvar+ddy2*vvar2)*255.0; if (mmy<0.0) mmy=0.0; else if (mmy>255.0) mmy=255.0;
	double mmz=(mz+ddz*vvar+ddz2*vvar2)*255.0; if (mmz<0.0) mmz=0.0; else if (mmz>255.0) mmz=255.0;

	return hvColRGB<unsigned char>(  (unsigned char)(mmx), (unsigned char)(mmy), (unsigned char)(mmz)  );
	//return hvColRGB<unsigned char>(  (unsigned char)((float)ind/(float)(piquantized.ncolors())*255.0), 
	//	(unsigned char)((float)ind/(float)(piquantized.ncolors())*255.0), 
	//	(unsigned char)((float)ind/(float)(piquantized.ncolors())*255.0) );
}

void doNonLinearIndex()
{
int i,j;
for (i=0; i<1024; i++) featappcount[i]=0;
for (i=0; i<128; i++) for (j=0; j<128; j++)
	{
		double vv=0.0, vvar=0.0, vvar2=0.0, rvv, rvvar, rvvar2;
		double ddx=8.0*(double)rand()/(double)RAND_MAX;
		double ddy=8.0*(double)rand()/(double)RAND_MAX;
		this->evalReg((double)i*8.0+ddx, (double)j*8.0+ddy, vv,vvar,vvar2);
		this->evalRnd((double)i*8.0+ddx, (double)j*8.0+ddy, 0, nstrates, rvv);
		vv += rvv*0.75;
		vv=(vv+1.0)/2.0*1023.0+0.5;
		int ind = (int)(vv);
		if (ind<0) ind=0; else if (ind>=1024) ind=1023;
		featappcount[ind]++;
	}
//for (i=0; i<255; i++) printf("t %3d: %g\n", i,  (double)featappcount[i]/64.0/64.0);
int indcurr=0;
indexfeat[0]=-2.0;
for (i=0; i<piquantized.ncolors(); i++)
{
	double sum=0.0;
	int indstart=indcurr;
	while (sum<(double)featcount[i]/(double)(exfeature.sizeX()*exfeature.sizeY()) && indcurr<1024) 
	{ 
		if (indcurr >= 1024) hvFatal("indcurr of of 1024 range");
			sum+=(double)featappcount[indcurr]/(double)(128*128);
			indcurr++;
	}
	indcurr--; sum-=(double)featappcount[indcurr]/(double)(128*128);
	indexfeat[i+1]=(double)(indcurr-1)/1024.0*2.0-1.0;
//	printf("feat %d from %g - %g cc=%g vs %g\n", i, (double)indstart/1024.0*2.0-1.0, (double)(indcurr-1)/1024.0*2.0-1.0, sum,  (double)featcount[i]/(double)(exfeature.sizeX()*exfeature.sizeY()));
}
indexfeat[piquantized.ncolors()+1]=2.0;
}

void makeSHisto()
{
	int i, j, k;
	hvPict<double> piv1(exfeature.sizeX(), exfeature.sizeY(), 0.0);
	hvPict<double> piv2(exfeature.sizeX(), exfeature.sizeY(), 0.0);
	hvPict<double> piv3(exfeature.sizeX(), exfeature.sizeY(), 0.0);
	for (i = 0; i < exfeature.sizeX() * 8; i += 8) for (j = 0; j < exfeature.sizeY() * 8; j += 8)
	{
		double vv, vvar, vvar2, xx, rvv = 0.0, rvvar = 0.0, rvvar2 = 0.0;
		double ppx = (double)i + 8.0*(double)rand() / (double)RAND_MAX, ppy = (double)j + 8.0*(double)rand() / (double)RAND_MAX;
		this->evalReg(ppx, ppy, vv, vvar, vvar2, false);
		for (k = 0; k < getStrates(); k++)
		{
			this->evalRnd(ppx, ppy, k, k + 1, xx, true);
			rvv += xx;
		}
		for (k = 0; k < getStratesv(); k++)
		{
			this->evalRndv(ppx, ppy, k, k + 1, xx, true);
			rvvar += xx;
		}
		for (k = 0; k < getStratesv2(); k++)
		{
			this->evalRndv2(ppx, ppy, k, k + 1, xx, true);
			rvvar2 += xx;
		}
		piv1.update(i / 8, j / 8, vv+rvv);
		piv2.update(i / 8, j / 8, vvar+rvvar);
		piv3.update(i / 8, j / 8, vvar2+rvvar2);
	}
	piv1.minmax(smin, smax);
	piv2.minmax(sminv, smaxv);
	piv3.minmax(sminv2, smaxv2);
	//smin = -1.0; smax = 1.0;
	//sminv = -1.0; smaxv = 1.0;
	//sminv2 = -1.0; smaxv2 = 1.0;
	for (i = 0; i < piv1.sizeX(); i ++) for (j = 0; j < piv1.sizeY(); j ++)
	{
		int bin = (int)((piv1.get(i,j)-smin)/(smax-smin)*512.0);
		if (bin < 0) bin = 0; else if (bin >= 512) bin = 511;
		shistofeature[bin] += 1;
		bin = (int)((piv2.get(i, j)-sminv)/(smaxv-sminv)*512.0);
		if (bin < 0) bin = 0; else if (bin >= 512) bin = 511;
		shistovar[bin] += 1;
		if (sminv2!=smaxv2) bin = (int)((piv3.get(i, j)-sminv2)/(smaxv2-sminv2)*512.0);
		else bin = 0;
		if (bin < 0) bin = 0; else if (bin >= 512) bin = 511;
		shistovar2[bin] += 1;
	}
	int hcount = 0;
	for (i = 0; i < 256; i++) hcount += histofeature[i];
	computeTransfer(transfeature, histofeature, shistofeature, 0, 255, hcount, 0, 511, piv1.sizeX()*piv1.sizeY());
	hcount = 0;
	for (i = 0; i < 256; i++) hcount += histovar[i];
	computeTransfer(transvar, histovar, shistovar, 0, 255, hcount, 0, 511, piv1.sizeX()*piv1.sizeY());
	hcount = 0;
	for (i = 0; i < 256; i++) hcount += histovar2[i];
	computeTransfer(transvar2, histovar2, shistovar2, 0, 255, hcount, 0, 511, piv1.sizeX()*piv1.sizeY());
//	for (k = 0; k < 512; k+=2) printf("transf %d = %d (%d - %d)\n", k, (int)transvar[k], histovar[k/2], shistovar[k]);
}

private:
///////////////////////////////////////
void doColorSpace(hvBitmap *mask=0)
{
int i,j,k;
pimlr.clone(pexample, 0,0, pexample.sizeX()-1, pexample.sizeY()-1);

double choixa=1.0, choixb=1.0;
bool first=true;
double minerr;
if (mask == 0)
{
	for (double aa = 0.0; aa <= 1.0; aa += 0.2) for (double bb = 0.0; bb <= 1.0; bb += 0.2)
	{
		hvQPictRGB<unsigned char, 64> qpiq;
		qpiq.quantizeHue(pimlr, NCOLORS, hvColRGB<double>(1.0 / (1.0 + 2.0*bb + 2.0*aa), bb, aa), 0);
		std::vector<hvFrame3<double> > fr(100);
		hvPictRGB<unsigned char> plr;
		fr.clear();
		qpiq.plsr(255, pimlr, fr);
		plr.clone(pimlr, 0, 0, pimlr.sizeX() - 1, pimlr.sizeY() - 1);
		qpiq.apply(255, plr, fr, 0.5, pcafactor);
		double err = 0.0;
		for (i = 0; i < pimlr.sizeX(); i++) for (j = 0; j < pimlr.sizeY(); j++)
		{
			double dd = ((double)plr.get(i, j).RED() - 128.0) / 128.0;
			err += dd*dd;
		}
		//hvPictRGB<unsigned char> ppr(plr,hvPictComponent::HV_RED);
		//ppr.extract(true,false, false, 0, 128,128);
		//qpiq.applyInverse(255,ppr,fr,0.5,pcafactor);
		//hvPictRGB<unsigned char> dd; dd.difference(ppr,pimlr);
		//double err = sqrt(dd.avg().normSquaredDouble());
		if (first) { first = false;  minerr = err; choixa = aa; choixb = bb; }
		else if (err < minerr) { minerr = err; choixa = aa; choixb = bb; }
//		printf("."); fflush(stdout);
	}
//	printf("\nbest color weights=(%g,%g)\n", choixa, choixb);
}
else { choixa = 0.5; choixb = 0.5; }
piquantized.quantizeHue(pimlr, NCOLORS, hvColRGB<double>(1.0/(1.0+2.0*choixa+2.0*choixb),choixa,choixb), 0);
//piquantized.quantizeHue(pimlr, NCOLORS, hvColRGB<double>(lumweight, 0.5, 0.5));
//piquantized.quantize(pimlr, NCOLORS);

int colindex[MAXCOLORS];
double cval[MAXCOLORS];
for (i=0; i<piquantized.ncolors(); i++) { colindex[i]=i; cval[i]=(double)piquantized.std::vector<hvColRGB<unsigned char> >::at(i).luminance()/255.0; }
for (i=0; i<piquantized.ncolors(); i++) 
	{
		int jmax=i, idmax=colindex[i];
		double vmax=cval[i];
		for (j=i+1; j<piquantized.ncolors(); j++)
			{
				if (cval[j]<vmax) { vmax=cval[j]; jmax=j; idmax=colindex[j]; }
			}
		cval[jmax]=cval[i]; colindex[jmax]=colindex[i];
		cval[i]=vmax; colindex[i]=idmax;
	}
//printf("sorted color indices: ");
//for (i=0; i<piquantized.ncolors(); i++) printf("%g (%d) ", cval[i], colindex[i]);
//printf("\n");
for (i=0; i<piquantized.sizeX(); i++) for (j=0; j<piquantized.sizeY(); j++)
	{
	int ind = (int)piquantized.hvArray2<unsigned char>::get(i,j);
	for (k=0; k<piquantized.ncolors(); k++) if (colindex[k]==ind) break;
	piquantized.hvArray2<unsigned char>::update(i,j,k);
	}
piquantized.updateTable(pexample);
piquantized.convert(piquant,1.0);
if (mask==0) piquantized.plsr(255, pexample, lfr);
else {	
	lfr.clear();
	hvFrame3<double> fr = pexample.pca(255, mask);
	lfr.push_back(fr);
}
piquantized.apply(255,pimlr,lfr,0.5,pcafactor);

hvLinearTransform3<double> t; 
for (i=0; i<piquantized.ncolors(); i++)
	{
	hvFrame3<double> fr = lfr.at(i);
	hvLinearTransform3<double> t; t = hvLinearTransform3<double>(fr);
	meancol[i]=hvVec3<double>(t.hvMat4<double>::L().X(), t.hvMat4<double>::L().Y(), t.hvMat4<double>::L().Z());
	meancolvalue[i]=(meancol[i].X()+meancol[i].Y()+meancol[i].Z())/3.0;
	principaldir[i]=hvVec3<double>(t.hvMat4<double>::I().X(), t.hvMat4<double>::I().Y(), t.hvMat4<double>::I().Z());
	principal2dir[i]=hvVec3<double>(t.hvMat4<double>::J().X(), t.hvMat4<double>::J().Y(), t.hvMat4<double>::J().Z());
	}
pvar.clone(pimlr, 0,0,pimlr.sizeX()-1, pimlr.sizeY()-1);
exfeature.clone(pexample, 0, 0, pexample.sizeX()-1, pexample.sizeY()-1);
colfeature.reset(pexample.sizeX(), pexample.sizeY(), hvColRGB<unsigned char>());
for (i=0; i<piquantized.ncolors(); i++) featcount[i]=0;
for (i = 0; i < 256; i++) { histofeature[i] = 0; histovar[i] = 0; histovar2[i] = 0; }
ominv = 255; omaxv = 0; ominv2 = 255; omaxv2 = 0;
for (i=0; i<exfeature.sizeX(); i++) for (j=0; j<exfeature.sizeY(); j++)
	{
	int ind = piquantized.hvArray2< unsigned char >::get(i,j);
	featcount[ind]++;
	int vv = (int)((float)ind/(float)(piquantized.ncolors())*255.0);
	if (vv > 255) vv = 255;
	histofeature[vv] += 1;
	exfeature.update(i,j,hvColRGB<unsigned char>((unsigned char)vv,(unsigned char)vv,(unsigned char)vv));
	colfeature.update(i,j,hvColRGB<unsigned char>((unsigned char)(meancol[ind].X()*255.0),(unsigned char)(meancol[ind].Y()*255.0),(unsigned char)(meancol[ind].Z()*255.0)));
	hvColRGB<unsigned char> col = pvar.get(i, j);
	if (mask == 0)
	{
		histovar[col.RED()] += 1;
		if (col.RED() < ominv) ominv = col.RED();
		if (col.RED() > omaxv) omaxv = col.RED();
		histovar2[col.GREEN()] += 1;
		if (col.GREEN() < ominv2) ominv2 = col.GREEN();
		if (col.GREEN() > omaxv2) omaxv2 = col.GREEN();
	}
	else if (mask->get(i,j))
	{
		histovar[col.RED()] += 1;
		if (col.RED() < ominv) ominv = col.RED();
		if (col.RED() > omaxv) omaxv = col.RED();
		histovar2[col.GREEN()] += 1;
		if (col.GREEN() < ominv2) ominv2 = col.GREEN();
		if (col.GREEN() > omaxv2) omaxv2 = col.GREEN();
	}
}
//double coeff = (double)(exfeature.sizeX()*exfeature.sizeY());
//for (i = 0; i < 256; i++) { histofeature[i] /= coeff; histovar[i] /= coeff; histovar2[i] /= coeff; }
//for (i=0; i<piquantized.ncolors(); i++)
//	{
//	printf("PCA color %d (%g):  mean %g,%g,%g \n dir %g, %g, %g\n dir2 %g, %g, %g\n", i, (double)featcount[i]/(double)(exfeature.sizeX()*exfeature.sizeY()), meancol[i].X(), meancol[i].Y(), meancol[i].Z(), principaldir[i].X(), principaldir[i].Y(), principaldir[i].Z(), principal2dir[i].X(), principal2dir[i].Y(), principal2dir[i].Z());
//	}
if (exportc1)
{
	char buff[512];
	sprintf(buff, "%s_C1.ppm",tname);
	FILE *fd = fopen(buff, "wb");
	colfeature.savePPM(fd,1);
	fclose(fd);
}
if (exportc2)
{
	hvPictRGB<unsigned char> pp(pvar, HV_RED);
	char buff[512];
	sprintf(buff, "%s_C2.ppm",tname);
	FILE *fd = fopen(buff, "wb");
	pp.savePPM(fd,1);
	fclose(fd);
}
if (exportc3)
{
	hvPictRGB<unsigned char> pp(pvar, HV_GREEN);
	char buff[512];
	sprintf(buff, "%s_C3.ppm",tname);
	FILE *fd = fopen(buff, "wb");
	pp.savePPM(fd,1);
	fclose(fd);
}
}

///////////////////////////////////////
/*
void computeColExtr()
{
int ii,jj;
for (ii=0; ii<pexample.sizeX(); ii++) for (jj=0; jj<pexample.sizeY(); jj++)
{
	if (ii==0 && jj==0) { savg=(double)exfeature.get(ii,jj).RED(); smin=(double)exfeature.get(ii,jj).RED(); smax=(double)exfeature.get(ii,jj).RED(); }
	else { savg+=(double)exfeature.get(ii,jj).RED(); if ((double)exfeature.get(ii,jj).RED()<smin) smin=(double)exfeature.get(ii,jj).RED(); if ((double)exfeature.get(ii,jj).RED()>smax) smax=(double)exfeature.get(ii,jj).RED(); }
	if (ii==0 && jj==0) { savgv=(double)pvar.get(ii,jj).RED(); sminv=(double)pvar.get(ii,jj).RED(); smaxv=(double)pvar.get(ii,jj).RED(); }
	else { savgv+=(double)pvar.get(ii,jj).RED(); if ((double)pvar.get(ii,jj).RED()<sminv) sminv=(double)pvar.get(ii,jj).RED(); if ((double)pvar.get(ii,jj).RED()>smaxv) smaxv=(double)pvar.get(ii,jj).RED(); }
	if (ii==0 && jj==0) { savgv2=(double)pvar.get(ii,jj).GREEN(); sminv2=(double)pvar.get(ii,jj).GREEN(); smaxv2=(double)pvar.get(ii,jj).GREEN(); }
	else { savgv2+=(double)pvar.get(ii,jj).GREEN(); if ((double)pvar.get(ii,jj).GREEN()<sminv2) sminv2=(double)pvar.get(ii,jj).GREEN(); if ((double)pvar.get(ii,jj).GREEN()>smaxv2) smaxv2=(double)pvar.get(ii,jj).GREEN(); }
}
savg=(savg/(double)(pexample.sizeX()*pexample.sizeY())-128.0)/128.0; smin=(smin-128.0)/128.0; smax=(smax-128.0)/128.0;
savgv=(savgv/(double)(pexample.sizeX()*pexample.sizeY())-128.0)/128.0; sminv=(sminv-128.0)/128.0; smaxv=(smaxv-128.0)/128.0;
savgv2=(savgv2/(double)(pexample.sizeX()*pexample.sizeY())-128.0)/128.0; sminv2=(sminv2-128.0)/128.0; smaxv2=(smaxv2-128.0)/128.0;
printf("avg=%g, min=%g, max=%g, (%g,%g,%g), (%g,%g,%g)\n", savg,smin,smax, savgv,sminv,smaxv,savgv2,sminv2,smaxv2);
savgmin=0.0; savgmax=0.0;
int avgminc=0, avgmaxc=0;
savgminv=0.0; savgmaxv=0.0;
int avgmincv=0, avgmaxcv=0;
savgminv2=0.0; savgmaxv2=0.0;
int avgmincv2=0, avgmaxcv2=0;
for (ii=0; ii<pexample.sizeX(); ii++) for (jj=0; jj<pexample.sizeY(); jj++)
{
	if (((double)exfeature.get(ii,jj).RED()-128.0)/128.0<savg) { savgmin+=((double)exfeature.get(ii,jj).RED()-128.0)/128.0; avgminc++; }
	else { savgmax+=((double)exfeature.get(ii,jj).RED()-128.0)/128.0; avgmaxc++; }
	if (((double)pvar.get(ii,jj).RED()-128.0)/128.0<savgv) { savgminv+=((double)pvar.get(ii,jj).RED()-128.0)/128.0; avgmincv++; }
	else { savgmaxv+=((double)pvar.get(ii,jj).RED()-128.0)/128.0; avgmaxcv++; }
	if (((double)pvar.get(ii,jj).GREEN()-128.0)/128.0<savgv2) { savgminv2+=((double)pvar.get(ii,jj).GREEN()-128.0)/128.0; avgmincv2++; }
	else { savgmaxv2+=((double)pvar.get(ii,jj).GREEN()-128.0)/128.0; avgmaxcv2++; }
}
savgmin /= (double)avgminc; savgmax /= (double)avgmaxc;
savgminv /= (double)avgmincv; savgmaxv /= (double)avgmaxcv;
savgminv2 /= (double)avgmincv2; savgmaxv2 /= (double)avgmaxcv2;
}
*/

void computeTransfer(hvArray1<unsigned char> &trans, int histo[], int proc_histo[], int i, int j, int sumh, int a, int b, int sumr)
{
	int k, sumnh = 0, sumnr = 0;
	int ii, aa, count = 0, ind;
	if (i == j) { for (k = a; k <= b; k++) trans[k] = (unsigned char)i; return; }
	if (a == b)
	{
		for (ind = i; histo[ind] == 0; ind++)
			;
		for (k = a; k <= b; k++)
			trans[k] = (unsigned char)ind; return;
	}

	for (k = i; k <= j; k++) if (histo[k] != 0) { count++; ind = k; }
	if (count == 1) { for (k = a; k <= b; k++) trans[k] = (unsigned char)ind; return; }
	count = 0;
	for (k = a; k <= b; k++) if (proc_histo[k] != 0) { count++; ind = k; }
	if (count == 1) { for (k = a; k <= b; k++) trans[k] = (unsigned char)((i + j) / 2); return; }

	for (k = i; k <= j && sumnh<sumh / 2; k++) { sumnh += histo[k]; }
	ii = k - 1; if (sumh - sumnh == 0) { sumnh -= histo[ii]; ii--; }
	for (k = a; k <= b && sumnr<sumr / 2; k++) { sumnr += proc_histo[k]; }
	aa = k - 1; if (sumr - sumnr == 0) { sumnr -= proc_histo[aa]; aa--; }
	computeTransfer(trans, histo, proc_histo, i, ii, sumnh, a, aa, sumnr);
	computeTransfer(trans, histo, proc_histo, ii + 1, j, sumh - sumnh, aa + 1, b, sumr - sumnr);
}

///////////////////////////////////////
void doFFT(hvPict<unsigned char> *ftpivar=0, hvPict<unsigned char> *ftpivar2=0)
{
int i,j,k,ii,jj;
int count=0;
exFT=0; exFTvar=0; exFTvar2=0;
// FFT of input example
if (ftpivar == 0)
{
	maxsize = 2;
	maxpow = 1;
	while (maxsize < pexample.sizeX() && maxsize < pexample.sizeY()) { maxpow++; maxsize *= 2; }
	maxsize /= 2; maxpow--;
}
else
{
	maxsize = ftpivar->sizeX();
	maxpow = 0;
	int s = maxsize; while (s != 1) {
		maxpow += 1; s /= 2;
	}
}
//printf("FFT size = %d\n",maxsize);

meanenrg.reset(maxsize,maxsize,0.0);
//for (i=0; i<pexample.sizeX()-maxsize-1; i+=1) for (j=0; j<pexample.sizeY()-maxsize-1; j+=1)
for (k=0; k<WELCH; k++)
{
	i = (int)((double)rand()/(double)RAND_MAX*(double)(pexample.sizeX()-maxsize-1));
	j = (int)((double)rand()/(double)RAND_MAX*(double)(pexample.sizeY()-maxsize-1));
count++;
hvPict<unsigned char> exfft(maxsize, maxsize, 0);
for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++) exfft.update(ii,jj,exfeature.get(ii+i, jj+j).RED());
if (exFT!=0) delete exFT;
exFT = exfft.fft(true, 128.0, 128.0, maxpow);
for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++) meanenrg.update(ii,jj,exFT->get(ii+maxsize*jj).energy()+meanenrg.get(ii,jj));
}
for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++) meanenrg.update(ii,jj,meanenrg.get(ii,jj)/(double)count);
//meanenrg.bilateral(BILAT,1.0,15,15);
if (exportc1)
{
	char buff[512];
	hvPict<unsigned char> expen(meanenrg, meanenrg.maxValue(), fftdrawlog, 255.0,0,0,meanenrg.sizeX()-1, meanenrg.sizeY()-1);
	expen.update(expen.sizeX()/2, expen.sizeY()/2, 0);
	hvPictRGB<unsigned char> expergb(expen, 1);
	sprintf(buff, "%s_C1_fft.ppm",tname);
	FILE *fdes = fopen(buff, "wb");
	expergb.savePPM(fdes, 1);
	fclose(fdes);
}

meanenrgv.reset(maxsize,maxsize,0.0);
if (ftpivar == 0)
{
	count = 0;
	//for (i=0; i<pexample.sizeX()-maxsize-1; i+=1) for (j=0; j<pexample.sizeY()-maxsize-1; j+=1)
	for (k = 0; k < WELCH; k++)
	{
		i = (int)((double)rand() / (double)RAND_MAX*(double)(pexample.sizeX() - maxsize - 1));
		j = (int)((double)rand() / (double)RAND_MAX*(double)(pexample.sizeY() - maxsize - 1));
		count++;
		hvPict<unsigned char> exfftvar(maxsize, maxsize, 0);
		for (ii = 0; ii < maxsize; ii++) for (jj = 0; jj < maxsize; jj++) exfftvar.update(ii, jj, pvar.get(ii + i, jj + j).RED());
		if (exFTvar != 0) delete exFTvar;
		exFTvar = exfftvar.fft(true, 128.0, 128.0, maxpow);
		for (ii = 0; ii < maxsize; ii++) for (jj = 0; jj < maxsize; jj++) meanenrgv.update(ii, jj, exFTvar->get(ii + maxsize*jj).energy() + meanenrgv.get(ii, jj));
	}
	for (ii = 0; ii < maxsize; ii++) for (jj = 0; jj < maxsize; jj++) meanenrgv.update(ii, jj, meanenrgv.get(ii, jj) / (double)count);
	//meanenrgv.bilateral(BILAT,1.0,15,15);
}
else
{
	exFTvar = new hvArray1<hvPair<double, double> >(maxsize*maxsize, hvPair<double, double>(0.0,0.0));
	for (ii = 0; ii < maxsize; ii++) for (jj = 0; jj < maxsize; jj++)
	{
		double vv = pow((double)ftpivar->get(ii, jj) / 255.0 / sqrt(2.0), 2.0)*3.0;
		exFTvar->update(ii + maxsize*jj, hvPair<double, double>(vv, vv));
		meanenrgv.update(ii, jj, exFTvar->get(ii + maxsize*jj).energy());
	}
}
if (exportc2)
{
	char buff[512];
	hvPict<unsigned char> expen(meanenrgv, meanenrgv.maxValue(), fftdrawlog, 255.0,0,0,meanenrgv.sizeX()-1, meanenrgv.sizeY()-1);
	expen.update(expen.sizeX()/2, expen.sizeY()/2, 0);
	hvPictRGB<unsigned char> expergb(expen, 1);
	sprintf(buff, "%s_C2_fft.ppm",tname);
	FILE *fdes = fopen(buff, "wb");
	expergb.savePPM(fdes, 1);
	fclose(fdes);
}
meanenrgv2.reset(maxsize,maxsize,0.0);
if (ftpivar == 0)
{
	count = 0;
	//for (i=0; i<pexample.sizeX()-maxsize-1; i+=1) for (j=0; j<pexample.sizeY()-maxsize-1; j+=1)
	for (k = 0; k < WELCH; k++)
	{
		i = (int)((double)rand() / (double)RAND_MAX*(double)(pexample.sizeX() - maxsize - 1));
		j = (int)((double)rand() / (double)RAND_MAX*(double)(pexample.sizeY() - maxsize - 1));
		count++;
		hvPict<unsigned char> exfftvar2(maxsize, maxsize, 0);
		for (ii = 0; ii < maxsize; ii++) for (jj = 0; jj < maxsize; jj++) exfftvar2.update(ii, jj, pvar.get(ii + i, jj + j).GREEN());
		if (exFTvar2 != 0) delete exFTvar2;
		exFTvar2 = exfftvar2.fft(true, 128.0, 128.0, maxpow);
		for (ii = 0; ii < maxsize; ii++) for (jj = 0; jj < maxsize; jj++) meanenrgv2.update(ii, jj, exFTvar2->get(ii + maxsize*jj).energy() + meanenrgv2.get(ii, jj));
	}
	for (ii = 0; ii < maxsize; ii++) for (jj = 0; jj < maxsize; jj++) meanenrgv2.update(ii, jj, meanenrgv2.get(ii, jj) / (double)count);
	//meanenrgv2.bilateral(BILAT,1.0,15,15);
}
else
{
	exFTvar2 = new hvArray1<hvPair<double, double> >(maxsize*maxsize, hvPair<double, double>(0.0, 0.0));
	for (ii = 0; ii < maxsize; ii++) for (jj = 0; jj < maxsize; jj++)
	{
		double vv = pow((double)ftpivar2->get(ii, jj) / 255.0 / sqrt(2.0) , 2.0)*3.0;
		exFTvar2->update(ii + maxsize*jj, hvPair<double, double>(vv, vv));
		meanenrgv2.update(ii, jj, exFTvar2->get(ii + maxsize*jj).energy());
	}
}
if (exportc3)
{
	char buff[512];
	hvPict<unsigned char> expen(meanenrgv2, meanenrgv2.maxValue(), fftdrawlog, 255.0,0,0,meanenrgv2.sizeX()-1, meanenrgv2.sizeY()-1);
	expen.update(expen.sizeX()/2, expen.sizeY()/2, 0);
	hvPictRGB<unsigned char> expergb(expen, 1);
	sprintf(buff, "%s_C3_fft.ppm",tname);
	FILE *fdes = fopen(buff, "wb");
	expergb.savePPM(fdes, 1);
	fclose(fdes);
}

glfreq.resize(maxsize*maxsize / 2); glfreq.clear();
glfreqv.resize(maxsize*maxsize/2); glfreqv.clear();
glfreqv2.resize(maxsize*maxsize/2); glfreqv2.clear();
genrgtot=0.0; genrgtotv=0.0; genrgtotv2=0.0;
gamptot=0.0; gamptotv=0.0; gamptotv2=0.0;
for (ii=0; ii<maxsize; ii++) for (jj=0;jj<maxsize; jj++)
	{
	//if (jj<=ii && (ii>maxsize/2 || ii!=jj)) { genrgtot+=exFT->get(ii+jj*maxsize).energy(); glfreq.pushBack(FF(ii,jj,exFT->get(ii+jj*maxsize).energy(),exFT->get(ii+jj*maxsize).phase())); }
	//if (jj<=ii && (ii>maxsize/2 || ii!=jj)) { genrgtotv+=exFTvar->get(ii+jj*maxsize).energy(); glfreqv.pushBack(FF(ii,jj,exFTvar->get(ii+jj*maxsize).energy(),exFTvar->get(ii+jj*maxsize).phase())); }
	//if (jj<=ii && (ii>maxsize/2 || ii!=jj)) { genrgtotv2+=exFTvar2->get(ii+jj*maxsize).energy(); glfreqv2.pushBack(FF(ii,jj,exFTvar2->get(ii+jj*maxsize).energy(),exFTvar2->get(ii+jj*maxsize).phase())); }
	if (jj<=ii && (ii>maxsize/2 || ii!=jj)) { genrgtot+=meanenrg.get(ii,jj); gamptot+=sqrt(meanenrg.get(ii,jj));  glfreq.push_back(FF(ii,jj,meanenrg.get(ii,jj),exFT->get(ii+jj*maxsize).phase())); }
	if (jj<=ii && (ii>maxsize/2 || ii!=jj)) { genrgtotv+=meanenrgv.get(ii,jj); gamptotv+=sqrt(meanenrgv.get(ii,jj)); glfreqv.push_back(FF(ii,jj,meanenrgv.get(ii,jj),exFTvar->get(ii+jj*maxsize).phase())); }
	if (jj<=ii && (ii>maxsize/2 || ii!=jj)) { genrgtotv2+=meanenrgv2.get(ii,jj); gamptotv2+=sqrt(meanenrgv2.get(ii,jj)); glfreqv2.push_back(FF(ii,jj,meanenrgv2.get(ii,jj),exFTvar2->get(ii+jj*maxsize).phase())); }
	}
glfreq.sort(); glfreqv.sort(); glfreqv2.sort();
glfreq.reverse(); glfreqv.reverse(); glfreqv2.reverse();
//printf("Energies: %g, %g, %g\n", genrgtot,genrgtotv,genrgtotv2);
gaaf0=exFT->get(maxsize/2+maxsize*maxsize/2).mod()*cos(exFT->get(maxsize/2+maxsize*maxsize/2).phase());
gaavf0=exFTvar->get(maxsize/2+maxsize*maxsize/2).mod()*cos(exFTvar->get(maxsize/2+maxsize*maxsize/2).phase());
gaav2f0=exFTvar2->get(maxsize/2+maxsize*maxsize/2).mod()*cos(exFTvar2->get(maxsize/2+maxsize*maxsize/2).phase());
if (exportc1 || exportc2 || exportc3)
{
	hvArray1<hvPair<double, double> > *rft = new hvArray1<hvPair<double, double> >(maxsize*maxsize, hvPair<double,double>(0.0,0.0));
	for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++)
	{
				double phase = M_PI*(double)rand()/(double)RAND_MAX;
				int fx=ii, fy=jj; 
				double c = sqrt(meanenrg.get(ii,jj));
				rft->update((fx>=maxsize/2?fx-maxsize/2:fx+maxsize/2)+(fy>=maxsize/2?fy-maxsize/2:fy+maxsize/2)*maxsize, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
				if (fx>0 && fy>0) rft->update(((maxsize-fx)>=maxsize/2?(maxsize-fx)-maxsize/2:(maxsize-fx)+maxsize/2)+((maxsize-fy)>=maxsize/2?(maxsize-fy)-maxsize/2:(maxsize-fy)+maxsize/2)*maxsize, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
	}
	for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,1,ii*maxsize,false);
	for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,maxsize,ii,false);
	hvPictRGB<unsigned char> pp(maxsize,maxsize, hvColRGB<unsigned char>(0));
	hvPict<double> ppf(maxsize,maxsize, 0.0);
	for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++)
	{
		hvPair<double,double> c = rft->get(ii+jj*maxsize);
		double val=gaaf0+c.mod()*cos(c.phase());
		ppf.update(ii,jj,val);
		pp.update(ii,jj,this->toRGB(val,0.0,0.0, true));
	}
	delete rft;
	char buff[512];
	sprintf(buff, "%s_C1_gaussian.ppm",tname);
	FILE *fdes = fopen(buff, "wb");
	pp.savePPM(fdes, 1);
	fclose(fdes);

	rft = new hvArray1<hvPair<double, double> >(maxsize*maxsize, hvPair<double,double>(0.0,0.0));
	for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++)
	{
				double phase = M_PI*(double)rand()/(double)RAND_MAX;
				int fx=ii, fy=jj; 
				double c = sqrt(meanenrgv.get(ii,jj));
				rft->update((fx>=maxsize/2?fx-maxsize/2:fx+maxsize/2)+(fy>=maxsize/2?fy-maxsize/2:fy+maxsize/2)*maxsize, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
				if (fx>0 && fy>0) rft->update(((maxsize-fx)>=maxsize/2?(maxsize-fx)-maxsize/2:(maxsize-fx)+maxsize/2)+((maxsize-fy)>=maxsize/2?(maxsize-fy)-maxsize/2:(maxsize-fy)+maxsize/2)*maxsize, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
	}
	for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,1,ii*maxsize,false);
	for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,maxsize,ii,false);
	hvPictRGB<unsigned char> ppv(maxsize,maxsize, hvColRGB<unsigned char>(0));
	for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++)
	{
		hvPair<double,double> c = rft->get(ii+jj*maxsize);
		double val=gaavf0+c.mod()*cos(c.phase());
		if (val>0.99) val=0.99; else if (val<-1.0) val=-1.0;
		unsigned char vv = (unsigned char)(val*128.0+128.0);
		ppv.update(ii,jj,hvColRGB<unsigned char>(vv,vv,vv));
	}
	delete rft;
	sprintf(buff, "%s_C2_gaussian.ppm",tname);
	fdes = fopen(buff, "wb");
	ppv.savePPM(fdes, 1);
	fclose(fdes);

	rft = new hvArray1<hvPair<double, double> >(maxsize*maxsize, hvPair<double,double>(0.0,0.0));
	for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++)
	{
				double phase = M_PI*(double)rand()/(double)RAND_MAX;
				int fx=ii, fy=jj; 
				double c = sqrt(meanenrgv2.get(ii,jj));
				rft->update((fx>=maxsize/2?fx-maxsize/2:fx+maxsize/2)+(fy>=maxsize/2?fy-maxsize/2:fy+maxsize/2)*maxsize, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
				if (fx>0 && fy>0) rft->update(((maxsize-fx)>=maxsize/2?(maxsize-fx)-maxsize/2:(maxsize-fx)+maxsize/2)+((maxsize-fy)>=maxsize/2?(maxsize-fy)-maxsize/2:(maxsize-fy)+maxsize/2)*maxsize, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
	}
	for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,1,ii*maxsize,false);
	for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,maxsize,ii,false);
	hvPictRGB<unsigned char> ppv2(maxsize,maxsize, hvColRGB<unsigned char>(0));
	for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++)
	{
		hvPair<double,double> c = rft->get(ii+jj*maxsize);
		double val=gaav2f0+c.mod()*cos(c.phase());
		if (val>0.99) val=0.99; else if (val<-1.0) val=-1.0;
		unsigned char vv = (unsigned char)(val*128.0+128.0);
		ppv2.update(ii,jj,hvColRGB<unsigned char>(vv,vv,vv));
	}
	delete rft;
	sprintf(buff, "%s_C3_gaussian.ppm",tname);
	fdes = fopen(buff, "wb");
	ppv2.savePPM(fdes, 1);
	fclose(fdes);

	for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++)
	{
		pp.update(ii,jj,this->toRGB(ppf.get(ii,jj),((double)ppv.get(ii,jj).RED()-128.0)/128.0,((double)ppv2.get(ii,jj).RED()-128.0)/128.0, true));
	}
	sprintf(buff, "%s_gaussian.ppm",tname);
	fdes = fopen(buff, "wb");
	pp.savePPM(fdes, 1);
	fclose(fdes);
}
}

///////////////////////////////////
void doStructureSeparation(double ratio)
{
int ii,jj;
double cumulate=0.0, cumulatev=0.0, cumulatev2=0.0;
int ind=0, indv=0, indv2=0;

//printf("Separation structure / noise:\n");
//printf("Energies: %g,%g,%g\n", genrgtot, genrgtotv, genrgtotv2);
if (this->nColors()==1) 
{ 
	ind=glfreq.size(); 
	while(indv<glfreqv.size() && cumulatev<gamptotv*ratio) { cumulatev+=sqrt(glfreqv.at(indv).geta()); indv++; }
	while(indv2<glfreqv2.size() && cumulatev2<gamptotv2*ratio) { cumulatev2+=sqrt(glfreqv2.at(indv2).geta()); indv2++; }
}
else 
{
	while(ind<glfreq.size() && cumulate<gamptot*ratio) { cumulate+=sqrt(glfreq.at(ind).geta()); ind++; }
	while(indv<glfreqv.size() && cumulatev<gamptotv*ratio) { cumulatev+=sqrt(glfreqv.at(indv).geta()); indv++; }
	while(indv2<glfreqv2.size() && cumulatev2<gamptotv2*ratio) { cumulatev2+=sqrt(glfreqv2.at(indv2).geta()); indv2++; }
}
ceilind=ind; ceilindv=indv; ceilindv2=indv2;
printf("corresponding ceil indices: %d, %d, %d\n", ind,indv,indv2);

hvArray1<hvPair<double, double> > *rft = new hvArray1<hvPair<double, double> >(maxsize*maxsize, hvPair<double,double>(0.0,0.0));
for (ii=0; ii<ind; ii++)
{
			double phase = glfreq.at(ii).getphi();
			int fx=glfreq.at(ii).getfx(), fy=glfreq.at(ii).getfy(); 
			//double c = 2.0*sqrt(glfreq.get(ii).geta());
			double c = 2.0*sqrt(exFT->get(fx+fy*maxsize).energy());
			rft->update((fx>=maxsize/2?fx-maxsize/2:fx+maxsize/2)+(fy>=maxsize/2?fy-maxsize/2:fy+maxsize/2)*maxsize, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
}
for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,1,ii*maxsize,false);
for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,maxsize,ii,false);
fexfeature.reset(maxsize,maxsize, hvColRGB<unsigned char>(0));
hvPictRGB<unsigned char> fcolfeature(maxsize,maxsize, hvColRGB<unsigned char>(0));
for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++)
{
	hvPair<double,double> c = rft->get(ii+jj*maxsize);
	double val=gaaf0+c.mod()*cos(c.phase());
	if (val>0.999) val=0.999; else if (val<-1.0) val=-1.0;
	unsigned char vv = (unsigned char)(val*128.0+128.0);
	fexfeature.update(ii,jj,hvColRGB<unsigned char>(vv,vv,vv));
	fcolfeature.update(ii,jj,this->toRGB(val,0.0,0.0, true));
}
delete rft;
rft = new hvArray1<hvPair<double, double> >(maxsize*maxsize, hvPair<double,double>(0.0,0.0));
hvArray1<hvPair<double, double> > *rft2 = new hvArray1<hvPair<double, double> >(maxsize*maxsize, hvPair<double,double>(0.0,0.0));
for (ii=0; ii<indv; ii++)
{
			double phase = glfreqv.at(ii).getphi();
			int fx=glfreqv.at(ii).getfx(), fy=glfreqv.at(ii).getfy(); 
			//double c = 2.0*sqrt(glfreqv.get(ii).geta());
			double c = 2.0*sqrt(exFTvar->get(fx+fy*maxsize).energy());
			rft->update((fx>=maxsize/2?fx-maxsize/2:fx+maxsize/2)+(fy>=maxsize/2?fy-maxsize/2:fy+maxsize/2)*maxsize, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
}
for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,1,ii*maxsize,false);
for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,maxsize,ii,false);
for (ii=0; ii<indv2; ii++)
{
			double phase = glfreqv2.at(ii).getphi();
			int fx=glfreqv2.at(ii).getfx(), fy=glfreqv2.at(ii).getfy(); 
			//double c = 2.0*sqrt(glfreqv2.get(ii).geta());
			double c = 2.0*sqrt(exFTvar2->get(fx+fy*maxsize).energy());
			rft2->update((fx>=maxsize/2?fx-maxsize/2:fx+maxsize/2)+(fy>=maxsize/2?fy-maxsize/2:fy+maxsize/2)*maxsize, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
}
for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft2, maxpow,1,ii*maxsize,false);
for (ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft2, maxpow,maxsize,ii,false);fpvar.reset(maxsize,maxsize, hvColRGB<unsigned char>(0));
fpvar.reset(maxsize,maxsize, hvColRGB<unsigned char>(0));
for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++)
{
	hvPair<double,double> c = rft->get(ii+jj*maxsize);
	double val=gaavf0+c.mod()*cos(c.phase());
	if (val>0.95) val=0.95; else if (val<-0.95) val=-0.95;
	unsigned char vv = (unsigned char)(val*128.0+128.0);
	c = rft2->get(ii+jj*maxsize);
	val=gaav2f0+c.mod()*cos(c.phase());
	if (val>0.95) val=0.95; else if (val<-0.95) val=-0.95;
	unsigned char vv2 = (unsigned char)(val*128.0+128.0);
	fpvar.update(ii,jj,hvColRGB<unsigned char>(vv,vv2,128));
}
delete rft;
delete rft2;

fexample.reset(maxsize,maxsize, hvColRGB<unsigned char>(0));
for (ii=0; ii<maxsize; ii++) for (jj=0; jj<maxsize; jj++)
{
	fexample.update(ii,jj,this->toRGB(((double)fexfeature.get(ii,jj).RED()-128.0)/128.0, ((double)fpvar.get(ii,jj).RED()-128.0)/128.0, ((double)fpvar.get(ii,jj).GREEN()-128.0)/128.0, true)); 
}

if (exportc1)
{
	char buff[512];
	sprintf(buff,"%s_C1_period.ppm", tname);
	FILE *fd=fopen(buff,"wb");
	fcolfeature.savePPM(fd,1);
	fclose(fd);
//	hvPictRGB<unsigned char> spectrum(maxsize,maxsize, hvColRGB<unsigned char>(0));
	hvPictRGB<unsigned char> spectrum(maxsize,maxsize, hvColRGB<unsigned char>(128,128,255));
	//for (ii=0; ii<glfreq.length(); ii++)
	for (ii=0; ii<ind; ii++)
	{
		int fx=glfreq.at(ii).getfx(), fy=glfreq.at(ii).getfy(); 
		spectrum.update(fx,fy, hvColRGB<unsigned char>(255,255,128));
		if (fx>0 && fy>0) spectrum.update(maxsize-fx,maxsize-fy, hvColRGB<unsigned char>(255,255,128));
//		double vv = log(1.0+fftdrawlog*(double)glfreq.get(ii).geta()/glfreq.get(0).geta())/log(fftdrawlog+1.0);
//		spectrum.update(fx,fy, hvColRGB<unsigned char>((unsigned char)255.0*vv));
//		if (fx>0 && fy>0) spectrum.update(maxsize-fx,maxsize-fy, hvColRGB<unsigned char>((unsigned char)255.0*vv));

		//printf("freq:%d -> e=%g\n", ii, glfreq.get(ii).geta());
	}
	sprintf(buff,"%s_C1_FFT_period.ppm", tname);
	fd=fopen(buff,"wb");
	spectrum.savePPM(fd,1);
	fclose(fd);
}
if (exportc2)
{
	char buff[512];
	sprintf(buff,"%s_C2_period.ppm", tname);
	hvPictRGB<unsigned char> pp(fpvar,HV_RED);
	FILE *fd=fopen(buff,"wb");
	pp.savePPM(fd,1);
	fclose(fd);
	hvPictRGB<unsigned char> spectrum(maxsize,maxsize, hvColRGB<unsigned char>(128,128,255));
	for (ii=0; ii<indv; ii++)
	{
		int fx=glfreqv.at(ii).getfx(), fy=glfreqv.at(ii).getfy(); 
		spectrum.update(fx,fy, hvColRGB<unsigned char>(255,255,128));
		if (fx>0 && fy>0) spectrum.update(maxsize-fx,maxsize-fy, hvColRGB<unsigned char>(255,255,128));
	}
	sprintf(buff,"%s_C2_FFT_period.ppm", tname);
	fd=fopen(buff,"wb");
	spectrum.savePPM(fd,1);
	fclose(fd);
}
if (exportc3)
{
	char buff[512];
	sprintf(buff,"%s_C3_period.ppm", tname);
	hvPictRGB<unsigned char> pp(fpvar,HV_GREEN);
	FILE *fd=fopen(buff,"wb");
	pp.savePPM(fd,1);
	fclose(fd);
	hvPictRGB<unsigned char> spectrum(maxsize,maxsize, hvColRGB<unsigned char>(128,128,255));
	for (ii=0; ii<indv2; ii++)
	{
		int fx=glfreqv2.at(ii).getfx(), fy=glfreqv2.at(ii).getfy(); 
		spectrum.update(fx,fy, hvColRGB<unsigned char>(255,255,128));
		if (fx>0 && fy>0) spectrum.update(maxsize-fx,maxsize-fy, hvColRGB<unsigned char>(255,255,128));
	}
	sprintf(buff,"%s_C3_FFT_period.ppm", tname);
	fd=fopen(buff,"wb");
	spectrum.savePPM(fd,1);
	fclose(fd);
}
if (exportc1 || exportc2 || exportc3)
{
	char buff[512];
	sprintf(buff,"%s_period.ppm", tname);
	FILE *fd=fopen(buff,"wb");
	fexample.savePPM(fd,1);
	fclose(fd);
}
//char buff[128]; gets(buff);
//exit(0);
}

///////////////////////////////////////
void doRegFreqSel(int ncos, double percentage)
{
	int subx, suby, i, j;
	hvPict<double> grid;
	hvPict<double> gridv;
	hvPict<double> gridv2;
	hvBitmap gsel;
	hvBitmap gselv;
	hvBitmap gselv2;
	double error;

	if (this->nColors()>1)
	{
		maxcosinesi=NREGCOS; maxcosinesiv=NREGCOS/2; maxcosinesiv2=NREGCOS/2;
	}
	else
	{
		maxcosinesi=0; maxcosinesiv=NREGCOS; maxcosinesiv2=NREGCOS;
	}
	pow_2=maxpow-1;	
	do { 
		size = (2<<pow_2)/2;	  	
		ssx=fexfeature.sizeX()/size; ssy=fexfeature.sizeY()/size;
		if (ssx<=0 || ssy<=0) hvFatal("pow_2 is too high for input image resolution!");
		if (ssx*ssy>4096) hvFatal("size (pow_2) is too small");
		printf("Windowed FFT size: %d, WFFT Grid dimension: %d,%d\n", size, ssx, ssy);
		printf("cos reg = %d,%d,%d\n", maxcosinesi,maxcosinesiv,maxcosinesiv2);
		grid.reset(ssx*size*2, ssy*size*2, 0.0);
		gridv.reset(ssx*size*2, ssy*size*2, 0.0);
		gridv2.reset(ssx*size*2, ssy*size*2, 0.0);
		gsel.reset(ssx*size*2, ssy*size*2, false);
		gselv.reset(ssx*size*2, ssy*size*2, false);
		gselv2.reset(ssx*size*2, ssy*size*2, false);
		for (suby=0; suby<ssy; suby++) for (subx=0; subx<ssx; subx++)
		{
			this->doFreqListsGrid(subx,suby, &grid, &gridv, &gridv2);
		}
		for (suby=0; suby<ssy; suby++) for (subx=0; subx<ssx; subx++)
		{
			this->regFreqSel(subx,suby, &gsel, &gselv, &gselv2);
		}
		error=0.0;
		for (i=0; i<ssx*size; i++) for (j=0; j<ssy*size; j++)
		{
			double vv, vvar, vvar2;
			this->evalReg((double)i, (double)j, vv, vvar, vvar2, false);
			error += (vv-((double)fexfeature.get(i,j).RED()-128.0)/128.0)*(vv-((double)fexfeature.get(i,j).RED()-128.0)/128.0);
			error += (vvar-((double)fpvar.get(i,j).RED()-128.0)/128.0)*(vvar-((double)fpvar.get(i,j).RED()-128.0)/128.0);
		}
		error /= 2.0*(double)(ssx*size*ssy*size);
		printf("error=%g\n", sqrt(error));
	pow_2--;
	} while(pow_2>2 && sqrt(error)>0.05);
pow_2++;
if (exportc1)
{
	hvPict<unsigned char> expen(grid, grid.maxValue(), fftdrawlog, 255.0,0,0,grid.sizeX()-1, grid.sizeY()-1);
	hvPictRGB<unsigned char> spectrum(ssx*size*2,ssy*size*2,hvColRGB<unsigned char>(0));
	for (int ii=0; ii<ssx*size*2; ii++) for (int jj=0; jj<ssy*size*2; jj++)
	{
		//if (gsel.get(ii,jj)) spectrum.update(ii,jj, hvColRGB<unsigned char>(expen.get(ii,jj), expen.get(ii,jj)/2, expen.get(ii,jj)/2));
		if (gsel.get(ii,jj)) spectrum.update(ii,jj, hvColRGB<unsigned char>(255,50,50));
		else spectrum.update(ii,jj, hvColRGB<unsigned char>(expen.get(ii,jj)));
	}
	char buff[512];
	sprintf(buff,"%s_C1_grid.ppm", tname);
	FILE *fd=fopen(buff,"wb");
	spectrum.savePPM(fd,1);
	fclose(fd);
}
if (exportc2)
{
	hvPict<unsigned char> expen(gridv, gridv.maxValue(), fftdrawlog, 255.0,0,0,gridv.sizeX()-1, gridv.sizeY()-1);
	hvPictRGB<unsigned char> spectrum(ssx*size*2,ssy*size*2,hvColRGB<unsigned char>(0));
	for (int ii=0; ii<ssx*size*2; ii++) for (int jj=0; jj<ssy*size*2; jj++)
	{
		if (gselv.get(ii,jj)) spectrum.update(ii,jj, hvColRGB<unsigned char>(expen.get(ii,jj), expen.get(ii,jj)/2, expen.get(ii,jj)/2));
		else spectrum.update(ii,jj, hvColRGB<unsigned char>(expen.get(ii,jj)));
	}
	char buff[512];
	sprintf(buff,"%s_C2_grid.ppm", tname);
	FILE *fd=fopen(buff,"wb");
	spectrum.savePPM(fd,1);
	fclose(fd);
}
if (exportc3)
{
	hvPict<unsigned char> expen(gridv2, gridv2.maxValue(), fftdrawlog, 255.0,0,0,gridv2.sizeX()-1, gridv2.sizeY()-1);
	hvPictRGB<unsigned char> spectrum(ssx*size*2,ssy*size*2,hvColRGB<unsigned char>(0));
	for (int ii=0; ii<ssx*size*2; ii++) for (int jj=0; jj<ssy*size*2; jj++)
	{
		if (gselv2.get(ii,jj)) spectrum.update(ii,jj, hvColRGB<unsigned char>(expen.get(ii,jj), expen.get(ii,jj)/2, expen.get(ii,jj)/2));
		else spectrum.update(ii,jj, hvColRGB<unsigned char>(expen.get(ii,jj)));
	}
	char buff[512];
	sprintf(buff,"%s_C3_grid.ppm", tname);
	FILE *fd=fopen(buff,"wb");
	spectrum.savePPM(fd,1);
	fclose(fd);
}
}



///////////////////////////////////////
void doFreqListsGrid(int subx, int suby, hvPict<double> *grid=0, hvPict<double> *gridv=0, hvPict<double> *gridv2=0)
{
	int ii,jj;

	int indsub=subx+suby*ssx;
	if (indsub >= 4096) hvFatal("indsub out of4096 limit");
	hvArray1<hvPair<double,double> > *ftfeature, *ftvar2, *ftvar22; 
	int indx = subx*size; if (indx<0) indx=0;
	int indy = suby*size; if (indy<0) indy=0;
	hvPict<unsigned char> pifft(2*size, 2*size, 0);
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) pifft.update(ii,jj,fexfeature.get(mirrored?indx+ii:indx+ii/2, mirrored?indy+jj:indy+jj/2).RED());
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) pifft.update(size+ii,jj,fexfeature.get(mirrored?indx+size-ii-1:indx+size/2+ii/2, mirrored?indy+jj:indy+jj/2).RED());
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) pifft.update(size+ii,size+jj,fexfeature.get(mirrored?indx+size-ii-1:indx+size/2+ii/2, mirrored?indy+size-jj-1:indy+size/2+jj/2).RED());
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) pifft.update(ii,size+jj,fexfeature.get(mirrored?indx+ii:indx+ii/2, mirrored?indy+size-jj-1:indy+size/2+jj/2).RED());
	ftfeature = pifft.fft(true, 128.0, 128.0, pow_2+1);
	if (grid!=0) for (ii=0; ii<2*size; ii++) for (jj=0;jj<2*size; jj++)
	{
		grid->update(2*subx*size+ii, 2*suby*size+jj, ftfeature->get(ii+jj*2*size).energy());
	}
	hvPict<unsigned char> pivfft(2*size, 2*size, 0);
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) pivfft.update(ii,jj,fpvar.get(mirrored?indx+ii:indx+ii/2, mirrored?indy+jj:indy+jj/2).RED());
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) pivfft.update(size+ii,jj,fpvar.get(mirrored?indx+size-ii-1:indx+size/2+ii/2, mirrored?indy+jj:indy+jj/2).RED());
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) pivfft.update(size+ii,size+jj,fpvar.get(mirrored?indx+size-ii-1:indx+size/2+ii/2, mirrored?indy+size-jj-1:indy+size/2+jj/2).RED());
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) pivfft.update(ii,size+jj,fpvar.get(mirrored?indx+ii:indx+ii/2, mirrored?indy+size-jj-1:indy+size/2+jj/2).RED());
	ftvar2 = pivfft.fft(true, 128.0, 128.0, pow_2+1);
	if (gridv!=0) for (ii=0; ii<2*size; ii++) for (jj=0;jj<2*size; jj++)
	{
		gridv->update(2*subx*size+ii, 2*suby*size+jj, ftvar2->get(ii+jj*2*size).energy());
	}
	hvPict<unsigned char> piv2fft(2*size, 2*size, 0);
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) piv2fft.update(ii,jj,fpvar.get(mirrored?indx+ii:indx+ii/2, mirrored?indy+jj:indy+jj/2).GREEN());
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) piv2fft.update(size+ii,jj,fpvar.get(mirrored?indx+size-ii-1:indx+size/2+ii/2, mirrored?indy+jj:indy+jj/2).GREEN());
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) piv2fft.update(size+ii,size+jj,fpvar.get(mirrored?indx+size-ii-1:indx+size/2+ii/2, mirrored?indy+size-jj-1:indy+size/2+jj/2).GREEN());
	for (ii=0; ii<size; ii++) for (jj=0; jj<size; jj++) piv2fft.update(ii,size+jj,fpvar.get(mirrored?indx+ii:indx+ii/2, mirrored?indy+size-jj-1:indy+size/2+jj/2).GREEN());
	ftvar22 = piv2fft.fft(true, 128.0, 128.0, pow_2+1);
	if (gridv2!=0) for (ii=0; ii<2*size; ii++) for (jj=0;jj<2*size; jj++)
	{
		gridv2->update(2*subx*size+ii, 2*suby*size+jj, ftvar22->get(ii+jj*2*size).energy());
	}

	lfreq[indsub].resize(2 * size*size); lfreq[indsub].clear();
	lfreqv[indsub].resize(2 * size*size); lfreqv[indsub].clear();
	lfreqv2[indsub].resize(2 * size*size); lfreqv2[indsub].clear();
	hvBitmap stillfree(2*size, 2*size, true);
	for (ii=0; ii<2*size; ii++) for (jj=0;jj<2*size; jj++) if (jj>ii) stillfree.set(ii,jj, false);
	for (ii=0; ii<=size; ii++) stillfree.set(ii,ii, false);
	hvBitmap stillfreev(2*size, 2*size, true);
	for (ii=0; ii<2*size; ii++) for (jj=0;jj<2*size; jj++) if (jj>ii) stillfreev.set(ii,jj, false);
	for (ii=0; ii<=size; ii++) stillfreev.set(ii,ii, false);
	hvBitmap stillfreev2(2*size, 2*size, true);
	for (ii=0; ii<2*size; ii++) for (jj=0;jj<2*size; jj++) if (jj>ii) stillfreev2.set(ii,jj, false);
	for (ii=0; ii<=size; ii++) stillfreev2.set(ii,ii, false);
	enrgtot[indsub]=0.0; enrgtotv[indsub]=0.0; enrgtotv2[indsub]=0.0;
	for (ii=0; ii<2*size; ii++) for (jj=0;jj<2*size; jj++)
	{
		if (stillfree.get(ii,jj)) { enrgtot[indsub]+=ftfeature->get(ii+jj*2*size).energy(); lfreq[indsub].push_back(FF(ii,jj,ftfeature->get(ii+jj*2*size).energy(),ftfeature->get(ii+jj*2*size).phase())); }
		if (stillfreev.get(ii,jj)) { enrgtotv[indsub]+=ftvar2->get(ii+jj*2*size).energy(); lfreqv[indsub].push_back(FF(ii,jj,ftvar2->get(ii+jj*2*size).energy(),ftvar2->get(ii+jj*2*size).phase())); }
		if (stillfreev2.get(ii,jj)) { enrgtotv2[indsub]+=ftvar22->get(ii+jj*2*size).energy(); lfreqv2[indsub].push_back(FF(ii,jj,ftvar22->get(ii+jj*2*size).energy(),ftvar22->get(ii+jj*2*size).phase())); }
	}
	//printf("Total Energies: %g, %g, %g\n",  enrgtot[indsub], enrgtotv[indsub], enrgtotv2[indsub]);
	lfreq[indsub].sort(); lfreqv[indsub].sort(); lfreqv2[indsub].sort();
	lfreq[indsub].reverse(); lfreqv[indsub].reverse(); lfreqv2[indsub].reverse();
	aaf0[indsub]=ftfeature->get(size+size*2*size).mod()*cos(ftfeature->get(size+size*2*size).phase());
	aavf0[indsub]=ftvar2->get(size+size*2*size).mod()*cos(ftvar2->get(size+size*2*size).phase());
	aav2f0[indsub]=ftvar22->get(size+size*2*size).mod()*cos(ftvar22->get(size+size*2*size).phase());
	delete ftfeature;
	delete ftvar2;
	delete ftvar22;
}


/////////////////////////////////////////////
void regFreqSel(int subx, int suby, hvBitmap *gsel=0, hvBitmap *gselv=0, hvBitmap *gselv2=0)
{
int i,j;

	int fid = subx+suby*ssx;
	if (fid >= 4096) hvFatal("out of range");
	if (aa[fid]!=0) delete aa[fid];
	aa[fid]=new double[maxcosinesi];
	if (phi[fid]!=0) delete phi[fid];
	phi[fid]=new double[maxcosinesi];
	if (aav[fid]!=0) delete aav[fid];
	aav[fid]=new double[maxcosinesiv];
	if (phiv[fid]!=0) delete phiv[fid];
	phiv[fid]=new double[maxcosinesiv];
	if (aav2[fid]!=0) delete aav2[fid];
	aav2[fid]=new double[maxcosinesiv2];
	if (phiv2[fid]!=0) delete phiv2[fid];
	phiv2[fid]=new double[maxcosinesiv2];
	if (ffx[fid]!=0) delete ffx[fid];
	ffx[fid]=new int[maxcosinesi];
	if (ffy[fid]!=0) delete ffy[fid];
	ffy[fid]=new int[maxcosinesi];
	if (ffxv[fid]!=0) delete ffxv[fid];
	ffxv[fid]=new int[maxcosinesiv];
	if (ffyv[fid]!=0) delete ffyv[fid];
	ffyv[fid]=new int[maxcosinesiv];
	if (ffxv2[fid]!=0) delete ffxv2[fid];
	ffxv2[fid]=new int[maxcosinesiv2];
	if (ffyv2[fid]!=0) delete ffyv2[fid];
	ffyv2[fid]=new int[maxcosinesiv2];

	// regular part (phase fixed): simply keep the maxcosinesi first highest
	double engsum=0.0, engsumv=0.0, engsumv2=0.0;
	for (i=0; i<maxcosinesi; i++)
	{
		FF nextfreq = lfreq[fid].at(i);
		aa[fid][i]=2.0*sqrt(nextfreq.geta()); phi[fid][i]=nextfreq.getphi(); 
		ffx[fid][i]=nextfreq.getfx(); ffy[fid][i]=nextfreq.getfy();
		engsum+=nextfreq.geta();
		if (gsel!=0) 
		{
			gsel->set(2*subx*size+nextfreq.getfx(), 2*suby*size+nextfreq.getfy(), true);
			if (nextfreq.getfx()>0 && nextfreq.getfy()>0) gsel->set(2*subx*size+2*size-nextfreq.getfx(), 2*suby*size+2*size-nextfreq.getfy(), true);
		}
	}
	for (i=0; i<maxcosinesiv; i++)
	{
		FF nextfreq = lfreqv[fid].at(i);
		aav[fid][i]=2.0*sqrt(nextfreq.geta()); phiv[fid][i]=nextfreq.getphi(); 
		ffxv[fid][i]=nextfreq.getfx(); ffyv[fid][i]=nextfreq.getfy();
		engsumv+=nextfreq.geta();
		if (gselv!=0) 
		{
			gselv->set(2*subx*size+nextfreq.getfx(), 2*suby*size+nextfreq.getfy(), true);
			if (nextfreq.getfx()>0 && nextfreq.getfy()>0) gselv->set(2*subx*size+2*size-nextfreq.getfx(), 2*suby*size+2*size-nextfreq.getfy(), true);
		}
	}
	for (i=0; i<maxcosinesiv2; i++)
	{
		FF nextfreq = lfreqv2[fid].at(i);
		aav2[fid][i]=2.0*sqrt(nextfreq.geta()); phiv2[fid][i]=nextfreq.getphi(); 
		ffxv2[fid][i]=nextfreq.getfx(); ffyv2[fid][i]=nextfreq.getfy();
		engsumv2+=nextfreq.geta();
		if (gselv2!=0) 
		{
			gselv2->set(2*subx*size+nextfreq.getfx(), 2*suby*size+nextfreq.getfy(), true);
			if (nextfreq.getfx()>0 && nextfreq.getfy()>0) gselv2->set(2*subx*size+2*size-nextfreq.getfx(), 2*suby*size+2*size-nextfreq.getfy(), true);
		}
	}
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void doRndFreqSel(int ncos, double percentage)
{
	if (this->nColors()>1)
	{
		nstrates = (int)((1.0-percentage)*(double)maxstrates);
		if (nstrates>0) maxrndcosines=(ncos-NREGCOS)/nstrates;
		else maxrndcosines=0;
		nstratesv = (int)((1.0-percentage)*(double)maxstrates);
		if (nstratesv>0) maxrndcosinesv=(ncos-NREGCOS/2)/nstratesv;
		else maxrndcosinesv=0;
		nstratesv2 = (int)((1.0-percentage)*(double)maxstrates);
		if (nstratesv2>0) maxrndcosinesv2=(ncos-NREGCOS/2)/nstratesv2;
		else maxrndcosinesv2=0;
	}
	else
	{
		nstrates=0; maxrndcosines=0;
		nstratesv=(int)((1.0-percentage)*(double)maxstrates);
		maxrndcosinesv=(ncos-NREGCOS)/nstratesv;
		nstratesv2=(int)((1.0-percentage)*(double)maxstrates);
		maxrndcosinesv2=(ncos-NREGCOS)/nstratesv2;

		if (ceilindv2 == 0) {
			nstratesv2 = 0; maxrndcosinesv2 = 0;
		}
	}

	//maxrndcosines=(int)((1.0-percentage)*(double)ncos)/(nstrates-1);
	printf("Compute random part...\n");
	printf("cos rnd=%d(%d), %d(%d), %d(%d)\n", maxrndcosines, nstrates, maxrndcosinesv, nstratesv, maxrndcosinesv2, nstratesv2);
	int i,j, ii,jj;
	hvBitmap stratumfreq[32], stratumfreqv[32], stratumfreqv2[32];
	int findex[32], findexv[32], findexv2[32];
	
	this->stratification(maxsize,maxpow,ceilind,glfreq,meanenrg,nstrates,thvar,findex,stratumfreq,exportc1,1,tname);
	this->stratification(maxsize,maxpow,ceilindv,glfreqv,meanenrgv,nstratesv,thvarv,findexv,stratumfreqv,exportc2,2,tname);
	this->stratification(maxsize,maxpow,ceilindv2,glfreqv2,meanenrgv2,nstratesv2,thvarv2,findexv2,stratumfreqv2,exportc3,3,tname);

	if (raa!=0) delete [] raa;
	raa=new double[maxrndcosines*nstrates*NRNDCOS];
	if (raa == 0) hvFatal("not enough memory");
	if (rphi!=0) delete [] rphi;
	rphi=new double[maxrndcosines*nstrates*NRNDCOS];
	if (rphi == 0) hvFatal("not enough memory");
	if (raav!=0) delete[] raav;
	raav=new double[maxrndcosinesv*nstratesv*NRNDCOS];
	if (raav == 0) hvFatal("not enough memory");
	if (rphiv!=0) delete [] rphiv;
	rphiv=new double[maxrndcosinesv*nstratesv*NRNDCOS];
	if (rphiv == 0) hvFatal("not enough memory");
	if (raav2!=0) delete [] raav2;
	raav2=new double[maxrndcosinesv2*nstratesv2*NRNDCOS];
	if (raav2 == 0) hvFatal("not enough memory");
	if (rphiv2!=0) delete [] rphiv2;
	rphiv2=new double[maxrndcosinesv2*nstratesv2*NRNDCOS];
	if (rphiv2 == 0) hvFatal("not enough memory");
	if (rffx!=0) delete[] rffx;
	rffx=new int[maxrndcosines*nstrates*NRNDCOS];
	if (rffx == 0) hvFatal("not enough memory");
	if (rffy!=0) delete[] rffy;
	rffy=new int[maxrndcosines*nstrates*NRNDCOS];
	if (rffy == 0) hvFatal("not enough memory");
	if (rffxv!=0) delete[] rffxv;
	rffxv=new int[maxrndcosinesv*nstratesv*NRNDCOS];
	if (rffxv == 0) hvFatal("not enough memory");
	if (rffyv!=0) delete[] rffyv;
	rffyv=new int[maxrndcosinesv*nstratesv*NRNDCOS];
	if (rffyv == 0) hvFatal("not enough memory");
	if (rffxv2!=0) delete[] rffxv2;
	rffxv2=new int[maxrndcosinesv2*nstratesv2*NRNDCOS];
	if (rffxv2 == 0) hvFatal("not enough memory");
	if (rffyv2!=0) delete[] rffyv2;
	rffyv2=new int[maxrndcosinesv2*nstratesv2*NRNDCOS];
	if (rffyv2 == 0) hvFatal("not enough memory");

//	pstratifc.reset(maxsize,maxsize,hvColRGB<unsigned char>(0));
//	hvList<FF> freqselection, freqselectionv, freqselectionv2;

	if(exportc1)
	{
		hvPictRGB<unsigned char> pp(maxsize,maxsize,hvColRGB<unsigned char>(0));
		hvPictRGB<unsigned char> ppe(maxsize,maxsize,hvColRGB<unsigned char>(0));
		hvPictRGB<unsigned char> ppp(maxsize, maxsize, hvColRGB<unsigned char>(255));
		//printf("tname(%d)=%s\n", strlen(tname), tname);
		hvProcWPhase::doRndCosSel(maxsize, ceilind, nstrates, maxrndcosines, stratumfreq, findex, raa, rphi, rffx, rffy, gauss, meanenrg, &pp, &ppe, &ppp, false);
		char buff[512];
		sprintf(buff,"%s_C1_FFT_strates_erosion.ppm",tname);
		FILE *fd=fopen(buff,"wb");
		pp.savePPM(fd,1);
		fclose(fd);
		sprintf(buff,"%s_C1_FFT_strates_erosion_clusters.ppm",tname);
		fd=fopen(buff,"wb");
		ppe.savePPM(fd,1);
		fclose(fd);
	}
	else hvProcWPhase::doRndCosSel(maxsize, ceilind, nstrates, maxrndcosines, stratumfreq, findex, raa, rphi, rffx, rffy, gauss, meanenrg);
	
	if(exportc2)
	{
		hvPictRGB<unsigned char> pp(maxsize,maxsize,hvColRGB<unsigned char>(0));
		hvPictRGB<unsigned char> ppe(maxsize,maxsize,hvColRGB<unsigned char>(0));
		hvPictRGB<unsigned char> ppp(maxsize,maxsize,hvColRGB<unsigned char>(255));
		hvProcWPhase::doRndCosSel(maxsize, ceilindv, nstratesv, maxrndcosinesv, stratumfreqv, findexv, raav, rphiv, rffxv, rffyv, gaussv, meanenrgv,&pp,&ppe, &ppp, false);
		char buff[512];
		sprintf(buff,"%s_C2_FFT_strates_erosion.ppm",tname);
		FILE *fd=fopen(buff,"wb");
		pp.savePPM(fd,1);
		fclose(fd);
		sprintf(buff,"%s_C2_FFT_strates_erosion_clusters.ppm",tname);
		fd=fopen(buff,"wb");
		ppe.savePPM(fd,1);
		fclose(fd);
		sprintf(buff,"%s_C2_FFT_frequencies.ppm",tname);
		fd=fopen(buff,"wb");
		ppp.savePPM(fd,1);
		fclose(fd);
	}
	else hvProcWPhase::doRndCosSel(maxsize, ceilindv, nstratesv, maxrndcosinesv, stratumfreqv, findexv, raav, rphiv, rffxv, rffyv, gaussv, meanenrgv); //, &pstratifc, true);

	if(exportc3)
	{
		hvPictRGB<unsigned char> pp(maxsize,maxsize,hvColRGB<unsigned char>(0));
		hvPictRGB<unsigned char> ppe(maxsize,maxsize,hvColRGB<unsigned char>(0));
		hvPictRGB<unsigned char> ppp(maxsize, maxsize, hvColRGB<unsigned char>(255));
		hvProcWPhase::doRndCosSel(maxsize, ceilindv2, nstratesv2, maxrndcosinesv2, stratumfreqv2, findexv2, raav2, rphiv2, rffxv2, rffyv2, gaussv2, meanenrgv2, &pp, &ppe, &ppp);
		char buff[512];
		sprintf(buff,"%s_C3_FFT_strates_erosion.ppm",tname);
		FILE *fd=fopen(buff,"wb");
		pp.savePPM(fd,1);
		fclose(fd);
		sprintf(buff,"%s_C3_FFT_strates_erosion_clusters.ppm",tname);
		fd=fopen(buff,"wb");
		ppe.savePPM(fd,1);
		fclose(fd);
	}
	else hvProcWPhase::doRndCosSel(maxsize, ceilindv2, nstratesv2, maxrndcosinesv2, stratumfreqv2, findexv2, raav2, rphiv2, rffxv2, rffyv2, gaussv2, meanenrgv2);

//FILE *fd=fopen("../stratif_erod.ppm","wb");
//pstratifc.savePPM(fd,1);
//fclose(fd);
//char buff[128]; gets(buff);
//exit(0);
/*
for (int idist=0; idist<nstratesv; idist++)
{
	double varcos=0.0;
	for (ii=0; ii<maxsize; ii+=2) for (jj=0;jj<maxsize; jj+=2)
	{
		double xx; this->evalRndv((double)ii, (double)jj, idist, idist+1, xx);
		varcos += 4.0*xx*xx/(float)(maxsize*maxsize);
	}
	printf("theoretical var %d =%g versus num = %g\n", idist, thvarv[idist], varcos);
}
*/
//printf("finish.\n");
}

void stratification(int maxsize, int maxpow, int ceilindv, std::vector<FF> &glfreqv, hvPict<double> &meanenrgv, 
	int nstratesv, double thvarv[], int ccv[], hvBitmap stratumfreqv[], bool exportc=false, int expcind=1, char *tname=0)
{
const double LG=1.0; 
hvPictRGB<unsigned char>pstratif(maxsize,maxsize,hvColRGB<unsigned char>(0));
if (expcind < 1 || expcind>3) hvFatal("out of expcind range, must be 1,2 or 3");
gaussapprox[expcind-1].reset(maxsize, maxsize, 0.0);

	int countv=ceilindv;
	for (int idist=0; idist<nstratesv; idist++)
	{
		printf("stratification %d -> strat %d, tname=%s\n", expcind, idist, tname);
		thvarv[idist]=0.0; 
		gausstrat[expcind - 1][idist].reset(maxsize,maxsize, hvColRGB<unsigned char>(0));
		stratumfreqv[idist].reset(maxsize,maxsize, false);
		hvArray1<hvPair<double, double> > *rft = new hvArray1<hvPair<double, double> >(maxsize*maxsize, hvPair<double,double>(0.0,0.0));
		double genrgtotv=0.0; for (int i=0; i<glfreqv.size(); i++) genrgtotv+=pow(glfreqv.at(i).geta(),LG);
		double enrgregv=0.0; for (int i=0; i<ceilindv; i++) enrgregv+=pow(glfreqv.at(i).geta(),LG);
		double cumulated = 0.0;
		ccv[idist]=countv;
		bool loop=true;
		while (loop && cumulated<(genrgtotv-enrgregv)/(double)nstratesv)
		{
			if (countv<glfreqv.size())
			{
				FF nextfreq = glfreqv.at(countv);
				double engmax=nextfreq.geta();
				thvarv[idist] += meanenrgv.get(nextfreq.getfx(),nextfreq.getfy());
				cumulated+=pow(engmax,LG); countv++; ccv[idist]=countv;
				stratumfreqv[idist].set(nextfreq.getfx(), nextfreq.getfy(), true);
				if (nextfreq.getfx()>0 && nextfreq.getfy()>0) stratumfreqv[idist].set(maxsize-nextfreq.getfx(), maxsize-nextfreq.getfy(), true);
				
				{
					pstratif.update(nextfreq.getfx(), nextfreq.getfy(),stratcoltab[idist]);
					if (nextfreq.getfx()>0 && nextfreq.getfy()>0) pstratif.update(maxsize-nextfreq.getfx(), maxsize-nextfreq.getfy(),stratcoltab[idist]);
					double phase = M_PI*(double)rand()/(double)RAND_MAX;
					int fx=nextfreq.getfx(), fy=nextfreq.getfy(); 
					double c = sqrt(nextfreq.geta());
					rft->update((fx>=maxsize/2?fx-maxsize/2:fx+maxsize/2)+(fy>=maxsize/2?fy-maxsize/2:fy+maxsize/2)*maxsize, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
					if (nextfreq.getfx()>0 && nextfreq.getfy()>0) rft->update(((maxsize-fx)>=maxsize/2?(maxsize-fx)-maxsize/2:(maxsize-fx)+maxsize/2)+((maxsize-fy)>=maxsize/2?(maxsize-fy)-maxsize/2:(maxsize-fy)+maxsize/2)*maxsize, hvPair<double,double>( c*cos(phase), -c*sin(phase) ));
				}
			}
			else loop=false;
		}
		{
			for (int ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,1,ii*maxsize,false);
			for (int ii=0; ii<maxsize; ii++) hvArray1<double>::fft(*rft, maxpow,maxsize,ii,false);
			for (int ii=0; ii<maxsize; ii++) for (int jj=0; jj<maxsize; jj++)
				{
				hvPair<double,double> c = rft->get(ii+jj*maxsize);
				double val=c.mod()*cos(c.phase());
				gaussapprox[expcind - 1].update(ii,jj,gaussapprox[expcind - 1].get(ii,jj)+val);
				if (val>1.0) val=1.0; else if (val<-1.0) val=-1.0;
				gausstrat[expcind - 1][idist].update(ii,jj,hvColRGB<unsigned char>((unsigned char)(255.0*(1.0+val)/2.0)));
				}
			delete rft;
		}
	}
if (exportc)
{
	char buff[512];
	sprintf(buff,"%s_C%d_FFT_strates.ppm",tname,expcind);
	FILE *fd=fopen(buff,"wb");
	pstratif.savePPM(fd,1);
	fclose(fd);
	for (int idist=0; idist<nstratesv; idist++)
	{
		sprintf(buff,"%s_C%d_gauss_rnd_strat%02d.ppm", tname, expcind, idist);
		fd=fopen(buff,"wb");
		gausstrat[expcind - 1][idist].savePPM(fd,1);
		fclose(fd);
	}
	hvPictRGB<unsigned char> gg(maxsize,maxsize,hvColRGB<unsigned char>());
	for (int ii=0; ii<maxsize; ii++) for (int jj=0; jj<maxsize; jj++)
		{
			double vv = gaussapprox[expcind - 1].get(ii,jj);
			if (vv<-1.0) vv=-1.0; else if (vv>1.0) vv=1.0;
			gg.update(ii,jj,hvColRGB<unsigned char>((unsigned char)(255.0*(vv+1.0)/2.0)));
		}
	sprintf(buff,"%s_C%d_gauss_rnd.ppm", tname, expcind);
	fd=fopen(buff,"wb");
	gg.savePPM(fd,1);
	fclose(fd);
}
}

static void doRndCosSel(int maxsize, int ceilind, int nstrates, int maxrndcosines, 
	hvBitmap stratumfreq[], int findex[], double raa[], double rphi[], int rffx[], int rffy[], 
	double gauss[], const hvPict<double> &meanenrg,
	hvPictRGB<unsigned char> *pstratifc=0, hvPictRGB<unsigned char> *pstratifce=0, hvPictRGB<unsigned char> *pstratifcp=0, bool verbose=false)
{
	int i,j,ii,jj;
	hvSortedList<FF> freqselection;
	if (nstrates >= 32) hvFatal("too many strates");
	for (int idist=0; idist<nstrates; idist++)
	{
		int count=findex[idist]-(idist>0?findex[idist-1]:ceilind);
		int radius=(int)(sqrt((double)count/(double)maxrndcosines)+0.5); if (radius<1) radius=1;
		gauss[idist] = (double)maxsize/(double)(radius);
		double gg = gauss[idist] / (double)maxsize;
		double gn = 0.5; while(gn>gg) gn*=0.5;
		gauss[idist]=gn*(radius==1?1.0:0.5);

		// generate list of frequencies for this stratum
		//hvBitmap straterod(maxsize,maxsize,false); 
		//straterod = stratumfreq[idist]; 
		//for (ii=0; ii<radius/8 && straterod.count()>10; ii++) straterod.erosion(3,3);
		//straterod.erosion2(radius/14);
		//straterod.erosion2((int)sqrt((double)radius)/2);
		int counterod = stratumfreq[idist].count();
		std::vector<FF> flist(2*count+1), flisterod(counterod);		
		freqselection.resize(counterod);
		flist.clear(); flisterod.clear(); freqselection.clear();
		printf("Stratum %d, %d freq: %d radius -> %d (%g)\n", idist, counterod, count, radius, gauss[idist]);
		for (ii=0; ii<maxsize; ii++) for (jj=0;jj<maxsize; jj++)
		{
			if (jj<=ii && (ii>maxsize/2 || ii!=jj) && stratumfreq[idist].get(ii,jj)) flisterod.push_back(FF(ii,jj,2.0*sqrt(meanenrg.get(ii,jj))));
			if (jj<=ii && (ii>maxsize/2 || ii!=jj) && stratumfreq[idist].get(ii,jj)) freqselection.push_back(FF(ii,jj,2.0*sqrt(meanenrg.get(ii,jj))));
			if (stratumfreq[idist].get(ii,jj)) flist.push_back(FF(ii,jj,2.0*sqrt(meanenrg.get(ii,jj))));
		}
		//freqselection.sort(); freqselection.reverse();
		counterod = flisterod.size();
		int nbtotfreq=flist.size();
		printf("maxrndcosines=%d\n", maxrndcosines);
		if (maxrndcosines>=freqselection.size())
		{
			int lmax = freqselection.size();
			for (int nncos=0; nncos<maxrndcosines; nncos++)
			{
				FF freq=freqselection.at(nncos%lmax);
				for (i=0; i<NRNDCOS; i++)
				{
				if (nncos*NRNDCOS + i + maxrndcosines*idist*NRNDCOS >= maxrndcosines*nstrates*NRNDCOS) hvFatal("out of freq tab range");
				raa[nncos*NRNDCOS+i+maxrndcosines*idist*NRNDCOS]=freq.geta();
				rffx[nncos*NRNDCOS+i+maxrndcosines*idist*NRNDCOS]=freq.getfx();
				rffy[nncos*NRNDCOS+i+maxrndcosines*idist*NRNDCOS]=freq.getfy();
				}
			}
			for (ii=0; ii<maxsize; ii++) for (jj=0;jj<maxsize; jj++)
			{
				if (pstratifc!=0) if (stratumfreq[idist].get(ii,jj))
					pstratifc->update(ii,jj,hvColRGB<unsigned char>(
						(unsigned char)((double)stratcoltab[idist].RED()),
						(unsigned char)((double)stratcoltab[idist].GREEN()),
						(unsigned char)((double)stratcoltab[idist].BLUE())
						));
				if (pstratifce!=0) if (stratumfreq[idist].get(ii,jj))
					pstratifce->update(ii,jj,hvColRGB<unsigned char>(
						(unsigned char)((double)stratcoltab[idist].RED()*(0.5+0.5*(double)rand()/(double)RAND_MAX)),
						(unsigned char)((double)stratcoltab[idist].GREEN()*(0.5+0.5*(double)rand()/(double)RAND_MAX)),
						(unsigned char)((double)stratcoltab[idist].BLUE()*(0.5+0.5*(double)rand()/(double)RAND_MAX))
						));
			}

		}
		else
		{
			// K means (clustering) : initial step by dart throwing
			int centerx[200], centery[200];
			double FRADIUS=0.5*sqrt((double)flisterod.size()/(double)maxrndcosines)/(double)radius;
			printf("FRADIUS=%g\n", FRADIUS);
			double amax=-1.0;
			for (int nncos=0; nncos<maxrndcosines; nncos++)
			{
				float vaa, fxaa, fyaa;
				std::vector<FF> clusterfreq;
				int clusternbtotfreq=0;
				int indmax;
				int fxmax, fymax;
				double variancecos=0.0;
			
				if (nbtotfreq>0)
				{
					variancecos=0.0;
					amax=-1.0;
					int ccmax=0;
					indmax=-1;
					int sec=0; 
					do { for (i=0; i<60; i++)
					{
						int ind = (int)((double)rand()/(double)RAND_MAX*(double)counterod); if (ind>=counterod) ind=counterod-1;
						double sum=0.0;
						int cc=0;
						FF freqind = flisterod.at(ind);
						for (j=0; j<flisterod.size(); j++) 
						{
							FF freq = flisterod.at(j);
							double dd = sqrt((double)((freqind.getfx()-freq.getfx())*(freqind.getfx()-freq.getfx()))+(double)((freqind.getfy()-freq.getfy())*(freqind.getfy()-freq.getfy())));
							if (dd<=(double)radius) { cc++; sum += freq.geta(); }
						}
						if (cc>ccmax) { amax=sum; indmax=ind; fxmax=freqind.getfx(); fymax=freqind.getfy(); ccmax=cc; }
					}
					sec++; } while (indmax==-1 && sec<10);
					if (sec==10) hvFatal("boucle infinie?");
					clusterfreq.resize(flist.size()); clusterfreq.clear();
					int cc=0;
					FF freqindmax = flisterod.at(indmax);
					for (j=0; j<flist.size(); j++) 
					{
						FF freq = flist.at(j);
						double dd = sqrt((double)((freqindmax.getfx()-freq.getfx())*(freqindmax.getfx()-freq.getfx()))+(double)((freqindmax.getfy()-freq.getfy())*(freqindmax.getfy()-freq.getfy())));
						if (dd<=(double)radius*FRADIUS)
						{ 
							clusterfreq.push_back(FF(freq.getfx(), freq.getfy(), freq.geta())),
							variancecos += freq.geta()*freq.geta()/2.0;
							cc++; 
						}
					}
					vaa=(float)amax/(float)cc; fxaa=(float)(fxmax-maxsize/2)/(float)maxsize; fyaa=(float)(fymax-maxsize/2)/(float)maxsize;
					clusternbtotfreq=cc; amax/=(double)cc;
				}
				else { amax=0.0; fxmax=0; fymax=0; vaa=0.0; fxaa=0.0; fyaa=0.0;  }
				// compute variances / energies of cluster
				if (clusternbtotfreq==0) { variancecos=0.0; }
				else { variancecos=sqrt(variancecos/(2.0*0.1125883*M_PI)); }
				if (nncos >= 200) hvFatal("out of nncos range");
				centerx[nncos]=fxmax; 
				centery[nncos]=fymax;

				if (nbtotfreq>0)
				{
				clusterfreq.resize(nbtotfreq); clusterfreq.clear();
				int cc=0;
				for (j=0; j<flist.size(); j++) 
				{
					FF freq = flist.at(j);
					double dd=sqrt((double)((fxmax-freq.getfx())*(fxmax-freq.getfx()))+(double)((fymax-freq.getfy())*(fymax-freq.getfy())));
					double dds=sqrt((double)((maxsize-fxmax-freq.getfx())*(maxsize-fxmax-freq.getfx()))+(double)((maxsize-fymax-freq.getfy())*(maxsize-fymax-freq.getfy())));
					if (dd>(double)radius*FRADIUS && dds>(double)radius*FRADIUS) clusterfreq.push_back(freq); 
				}
				flist=clusterfreq;
				nbtotfreq=flist.size();

				clusterfreq.resize(counterod); clusterfreq.clear();
				for (j=0; j<flisterod.size(); j++) 
				{
					FF freq = flisterod.at(j);
					double dd=sqrt((double)((fxmax-freq.getfx())*(fxmax-freq.getfx()))+(double)((fymax-freq.getfy())*(fymax-freq.getfy())));
					if (dd>(double)radius*FRADIUS) clusterfreq.push_back(freq); 
				}
				flisterod=clusterfreq;
				counterod=flisterod.size();
				if (counterod==0) nbtotfreq=0;
				}
			}

			// re-centering centroids
			
			double cxx[200], cyy[200];
			int ncc[200];
			if (maxrndcosines > 200) hvFatal("too much maxrndcosines");
			for (int nncos=0; nncos<maxrndcosines; nncos++) { cxx[nncos]=0.0; cyy[nncos]=0.0; ncc[nncos]=0; }
			for (ii=0; ii<maxsize; ii++) for (jj=0;jj<maxsize; jj++)
			{
				double dmin=(double)maxsize*sqrt((double)2.0);
				int indmin=-1;
				bool addit=false;
				for (int nncos=0; nncos<maxrndcosines; nncos++)
				{
					if (nncos >= 200) hvFatal("out of nncos range");
					int cfx=centerx[nncos];
					int cfy=centery[nncos];
					double dd = sqrt((double)((cfx-ii)*(cfx-ii))+(double)((cfy-jj)*(cfy-jj)));
					double dds=sqrt((double)((maxsize-cfx-ii)*(maxsize-cfx-ii))+(double)((maxsize-cfy-jj)*(maxsize-cfy-jj)));
					if (stratumfreq[idist].get(ii,jj))
						if (dd<dmin||dds<dmin) { dmin=(dd<dds?dd:dds); indmin=nncos; addit=dd<dds; }
				}
				if (addit) { if (indmin >= 200) hvFatal("out of indmin range");  cxx[indmin]+=(double)ii; cyy[indmin]+=(double)jj; ncc[indmin]++; }
			}
			for (int nncos=0; nncos<maxrndcosines; nncos++) 
			{ 
				if (nncos >= 200) hvFatal("out of nncos range");
				centerx[nncos]=(int)(cxx[nncos]/(double)ncc[nncos]+0.5);
				centery[nncos]=(int)(cyy[nncos]/(double)ncc[nncos]+0.5);
			}
			
			// finalize clusters, by computing energy
			double varcos[200];
			for (int nncos=0; nncos<maxrndcosines; nncos++) { varcos[nncos]=0.0; }
			for (ii=0; ii<maxsize; ii++) for (jj=0;jj<maxsize; jj++)
			{
				double dmin=(double)maxsize*sqrt((double)2.0);
				int indmin=-1;
				for (int nncos=0; nncos<maxrndcosines; nncos++)
				{
					if (nncos >= 200) hvFatal("out of nncos range");
					int cfx=centerx[nncos];
					int cfy=centery[nncos];
					double dd = sqrt((double)((cfx-ii)*(cfx-ii))+(double)((cfy-jj)*(cfy-jj)));
					double dds=sqrt((double)((maxsize-cfx-ii)*(maxsize-cfx-ii))+(double)((maxsize-cfy-jj)*(maxsize-cfy-jj)));
					if (stratumfreq[idist].get(ii,jj))
						if (dd<dmin||dds<dmin) { dmin=(dd<dds?dd:dds); indmin=nncos; }
				}
				if (pstratifc!=0) if (stratumfreq[idist].get(ii,jj))
					pstratifc->update(ii,jj,hvColRGB<unsigned char>(
						(unsigned char)((double)stratcoltab[idist].RED()),
						(unsigned char)((double)stratcoltab[idist].GREEN()),
						(unsigned char)((double)stratcoltab[idist].BLUE())
						));
				if (pstratifce!=0) if (stratumfreq[idist].get(ii,jj))
					pstratifce->update(ii,jj,hvColRGB<unsigned char>(
						(unsigned char)((double)stratcoltab[idist].RED()*(0.5+0.5*(double)indmin/(double)maxrndcosines)),
						(unsigned char)((double)stratcoltab[idist].GREEN()*(0.5+0.5*(double)indmin/(double)maxrndcosines)),
						(unsigned char)((double)stratcoltab[idist].BLUE()*(0.5+0.5*(double)indmin/(double)maxrndcosines))
						));
				if (jj<=ii && (ii>maxsize/2 || ii!=jj))
				{
					if (indmin >= 200) hvFatal("out of varcos indmin range");
					varcos[indmin] += 4.0*meanenrg.get(ii,jj)/2.0;
				}
			}

			for (int nncos=0; nncos<maxrndcosines; nncos++)
			{
				double enrgf=0.0;
				int count=0;
				std::vector<FF> freqcluster(freqselection.size()); freqcluster.clear();
				for (i=0; i<freqselection.size(); i++)
				{
					FF freq=freqselection.at(i);
					double dmin=(double)maxsize*sqrt((double)2.0);
					int indmin=-1;
					for (j=0; j<maxrndcosines; j++)
					{
						if (j >= 200) hvFatal("out of nncos range");
						int cfx=centerx[j];
						int cfy=centery[j];
						double dd = sqrt((double)((cfx-freq.getfx())*(cfx-freq.getfx()))+(double)((cfy-freq.getfy())*(cfy-freq.getfy())));
						double dds=sqrt((double)((maxsize-cfx-freq.getfx())*(maxsize-cfx-freq.getfx()))+(double)((maxsize-cfy-freq.getfy())*(maxsize-cfy-freq.getfy())));
						if (dd<dmin||dds<dmin) { dmin=(dd<dds?dd:dds); indmin=j; }
					}
					if (indmin==nncos) freqcluster.push_back(freq);
				}
				if (freqcluster.size()==0) 
				{ 
					printf("warning cosinus %d, 0 freq out of %d\n",nncos,(int)freqselection.size());
					for (count=0; count<NRNDCOS; count++)
					{
						if (nncos*NRNDCOS + count + maxrndcosines*idist*NRNDCOS >= maxrndcosines*nstrates*NRNDCOS) hvFatal("out of freq tab range");
						raa[nncos*NRNDCOS+count+maxrndcosines*idist*NRNDCOS]=0.0;
						rffx[nncos*NRNDCOS+count+maxrndcosines*idist*NRNDCOS]=0;
						rffy[nncos*NRNDCOS+count+maxrndcosines*idist*NRNDCOS]=0;
					}
				}
				else for (count=0; count<NRNDCOS; count++)
				{
					int indmin;
					if (freqcluster.size()<NRNDCOS) indmin=count%freqcluster.size();
					else indmin = (int)((double)rand()/(double)RAND_MAX*(double)freqcluster.size());
					if (indmin>=freqcluster.size()) indmin=freqcluster.size()-1;
					FF freq=freqcluster.at(indmin);
					if (nncos*NRNDCOS + count + maxrndcosines*idist*NRNDCOS >= maxrndcosines*nstrates*NRNDCOS) hvFatal("out of freq tab range");
					raa[nncos*NRNDCOS+count+maxrndcosines*idist*NRNDCOS]=freq.geta();
					rffx[nncos*NRNDCOS+count+maxrndcosines*idist*NRNDCOS]=freq.getfx();
					rffy[nncos*NRNDCOS+count+maxrndcosines*idist*NRNDCOS]=freq.getfy();
					enrgf+=freq.geta();
					if (pstratifcp!=0 && count==0) 
					{
							pstratifcp->update(freq.getfx(), freq.getfy(), hvColRGB<unsigned char>(
								(unsigned char)((double)stratcoltab[idist].RED()*(1.0+0.0*(double)nncos/(double)maxrndcosines)),
								(unsigned char)((double)stratcoltab[idist].GREEN()*(1.0+0.0*(double)nncos/(double)maxrndcosines)),
								(unsigned char)((double)stratcoltab[idist].BLUE()*(1.0+0.0*(double)nncos/(double)maxrndcosines))
								));
							if (freq.getfx()>0 && freq.getfy()>0) pstratifcp->update(maxsize-freq.getfx(), maxsize-freq.getfy(), hvColRGB<unsigned char>(
								(unsigned char)((double)stratcoltab[idist].RED()*(1.0+0.0*(double)nncos/(double)maxrndcosines)),
								(unsigned char)((double)stratcoltab[idist].GREEN()*(1.0+0.0*(double)nncos/(double)maxrndcosines)),
								(unsigned char)((double)stratcoltab[idist].BLUE()*(1.0+0.0*(double)nncos/(double)maxrndcosines))
								));
					}

				}
				enrgf /= (double)NRNDCOS;
				for (i=0; i<NRNDCOS; i++) 
				{
					if (nncos*NRNDCOS + i + maxrndcosines*idist*NRNDCOS >= maxrndcosines*nstrates*NRNDCOS) hvFatal("out of freq tab range");
					if (enrgf!=0.0) raa[nncos*NRNDCOS+i+maxrndcosines*idist*NRNDCOS]*=sqrt(0.9*varcos[nncos]/(2.0*0.1125883*M_PI))/enrgf;
					if (verbose) printf("%02d: cos%d (%d,%d) -> A0=%g (%d,%d) (eng=%g)\n", idist,nncos,centerx[nncos]-maxsize/2, centery[nncos]-maxsize/2, 
						raa[nncos*NRNDCOS+i+maxrndcosines*idist*NRNDCOS],
						rffx[nncos*NRNDCOS+i+maxrndcosines*idist*NRNDCOS]-maxsize/2,
						rffy[nncos*NRNDCOS+i+maxrndcosines*idist*NRNDCOS]-maxsize/2, 
						varcos[nncos]);
				}
			}

		} // end if maxrndcosines>=freqselection.length()
	} // end for (int idist=0; idist<nstrates; idist++)
	//printf("tname(%d)=%s\n", strlen(tname), tname);

}

public:
void doTurbulence()
{
	if (this->nColors()==1) hvProcWPhase::doTurbulence(glfreqv,maxsize,gausst,taa,tffx,tffy);
	else hvProcWPhase::doTurbulence(glfreq,maxsize,gausst,taa,tffx,tffy);
}

static void doTurbulence(hvSortedList<FF> &glfreq, int maxsize, double gausst[],
	double taa[NTURB_STRAT][NTURB_COS], int tffx[NTURB_STRAT][NTURB_COS], int tffy[NTURB_STRAT][NTURB_COS])
{
	int ccv[NTURB_STRAT];
	double normalise=0.0;
	int count=0;
	double genrgtot=0.0; 
	for (int i=0; i<glfreq.size()/3; i++) genrgtot+=pow(glfreq.at(i).geta(),1.0);
	for (int idist=0; idist<NTURB_STRAT; idist++)
	{
		double cumulated = 0.0;
		ccv[idist]=0;
		bool loop=true;
		std::vector<FF> sfreq(glfreq.size()); sfreq.clear();
		while (loop && cumulated<(genrgtot)/(double)NTURB_STRAT)
		{
				FF nextfreq = glfreq.at(count);
				double engmax=nextfreq.geta();
				cumulated+=pow(engmax,1.0);
				sfreq.push_back(nextfreq);
				count++; ccv[idist]=count;
				loop = count < glfreq.size();
		}

		int num=ccv[idist]-(idist>0?ccv[idist-1]:0);
		int radius=(int)(sqrt((double)num/(double)NTURB_COS)+0.5); if (radius<1) radius=1;
		gausst[idist] = (double)maxsize/(double)(radius);
		double gg = gausst[idist] / (double)maxsize;
		double gn = 0.5; while(gn>gg) gn*=0.5;
		gausst[idist]=gn*(radius==1?1.0:0.5);
//		printf("Turbulence Stratum %d: %d radius -> %d (%g)\n", idist, (int)sfreq.size(), radius, gausst[idist]);
		
		int centerx[NTURB_COS], centery[NTURB_COS], index[NTURB_COS];
		int centerxmax[NTURB_COS], centerymax[NTURB_COS], indexmax[NTURB_COS];
		double summax=0.0;
		for (int i=0; i<200; i++)
		{
			for (int nncos=0; nncos<NTURB_COS; nncos++)
			{
				int ind = (int)((double)rand()/(double)RAND_MAX*(double)sfreq.size()); if (ind>=sfreq.size()) ind=sfreq.size()-1;
				FF nextfreq = sfreq.at(ind);
				centerx[nncos]=nextfreq.getfx(); centery[nncos]=nextfreq.getfy(); index[nncos]=ind;
			}
			double sum=0.0;
			for (int j=0; j<sfreq.size(); j++)
			{
				FF nextfreq = sfreq.at(j);
				bool add=false;
				for (int nncos=0; nncos<NTURB_COS; nncos++)
				{
					double dd = sqrt((double)((centerx[nncos]-nextfreq.getfx())*(centerx[nncos]-nextfreq.getfx()))+(double)((centery[nncos]-nextfreq.getfy())*(centery[nncos]-nextfreq.getfy())));
					double dds=sqrt((double)((maxsize-centerx[nncos]-nextfreq.getfx())*(maxsize-centerx[nncos]-nextfreq.getfx()))+(double)((maxsize-centery[nncos]-nextfreq.getfy())*(maxsize-centery[nncos]-nextfreq.getfy())));
					if (dd<(double)radius||dds<(double)radius) { add=true; }
				}
				if (add)
				{
					sum+=nextfreq.geta();
				}
			}
			if (sum>summax)
			{
				summax=sum;
				for (int nncos=0; nncos<NTURB_COS; nncos++) { indexmax[nncos]=index[nncos]; centerxmax[nncos]=centerx[nncos]; centerymax[nncos]=centery[nncos]; } 
			}
		}
		//printf("indmax: ");
		for (int nncos=0; nncos<NTURB_COS; nncos++)
		{
			//printf("%d ", indexmax[nncos]);
			taa[idist][nncos]=sqrt(sfreq.at(indexmax[nncos]).geta());
			normalise += taa[idist][nncos];
			tffx[idist][nncos]=sfreq.at(indexmax[nncos]).getfx();
			tffy[idist][nncos]=sfreq.at(indexmax[nncos]).getfy();
		}
		//printf("\n");
	}
	for (int idist=0; idist<NTURB_STRAT; idist++)
	{
		for (int nncos=0; nncos<NTURB_COS; nncos++)
		{
			//taa[idist][nncos]= pow(taa[idist][nncos]/normalise, 1.0/(double)(idist+1));
			taa[idist][nncos] = taa[idist][nncos] / normalise;
			//printf("turb: s %d: cos %d -> A=%g (%d,%d)\n", idist,nncos, taa[idist][nncos], tffx[idist][nncos]-maxsize/2, tffy[idist][nncos]-maxsize/2);
		}
	}

}
void saveTurbulence(char *name)
{
	FILE *fd = fopen(name, "w");
	if (fd == 0) { perror("cannot save turbulence:"); return ; }
	fprintf(fd, "%d\n", maxsize);
	fprintf(fd, "%d\n", NTURB_STRAT);
	int i, j;
	for (i = 0; i < NTURB_STRAT; i++) fprintf(fd, "%g\n", gausst[i]);
	fprintf(fd, "%d\n", NTURB_COS);
	for (i = 0; i < NTURB_STRAT; i++) for (j = 0; j < NTURB_COS; j++)
	{
		fprintf(fd, "%g %d %d\n", taa[i][j], tffx[i][j], tffy[i][j]);
	}
	fclose(fd);
}
double turbulence(double px, double py)
{
	int i, di,dj;
	double vv=0.0; 
	double coeff=0.0;
	for (int idist=0; idist<NTURB_STRAT; idist++)
	{
			double gg = gausst[idist];
			double svv=0.0;
			//for (int di=0; di<=0; di++) for (int dj=0; dj<=0; dj++)
			for (int di=-1; di<=1; di++) for (int dj=-1; dj<=1; dj++)
				{
					double x=(double)((int)floor(px)%int((double)maxsize*gg))-(double)(di*maxsize)*gg + px - floor(px);
					double y=(double)((int)floor(py)%int((double)maxsize*gg))-(double)(dj*maxsize)*gg + py - floor(py);
					hvNoise::seeding((int)(((double)floor(px)+(double)(di*maxsize)*gg)/((double)maxsize*gg)), (int)(((double)floor(py)+(double)(dj*maxsize)*gg)/((double)maxsize*gg)), 0);
					double sum=0.0;
					//printf("%g %g (%d,%d) -> %g,%g\n", px,py,di,dj,x,y);
					for (i=0; i<NTURB_COS; i++)
					{
						double fx=(double)(tffx[idist][i]-maxsize/2)/(double)maxsize;
						double fy=(double)(tffy[idist][i]-maxsize/2)/(double)maxsize;
						sum += taa[idist][i]*cos(2.0*M_PI*(x*fx+y*fy)+M_PI*hvNoise::next());
					}
					//x=(double)(((int)floor(px))%int(gg*(double)maxsize))-(double)(di*maxsize)*gg + px - floor(px);
					//y=(double)(((int)floor(py))%int(gg*(double)maxsize))-(double)(dj*maxsize)*gg + py - floor(py);
					double cx = (x-(double)(maxsize/2)*gg)/((double)(maxsize)*gg);
					double cy = (y-(double)(maxsize/2)*gg)/((double)(maxsize)*gg);
					double factor = hvNoise::Kaiser((cx*cx+cy*cy)/1.5/1.5);
					svv += factor*sum;
				}
			//if (svv<0.0) svv=-svv;
			vv+=svv;
	}
	return vv;
}

};

}

#endif // !defined(AFX_FUNCTION_H__B6AC0A32_75EF_428E_BC10_6219F619FA29__INCLUDED_)
