
// IMAGE
#include <itkShiftScaleImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkSquareImageFilter.h>
#include <itkSqrtImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkStatisticsImageFilter.h>

// FFT
#include <itkForwardFFTImageFilter.h>
#include <itkInverseFFTImageFilter.h>
#include <itkFFTShiftImageFilter.h>

// COMPLEX
#include <itkComplexToRealImageFilter.h>
#include <itkComplexToImaginaryImageFilter.h>
#include <itkComplexToModulusImageFilter.h>
#include <itkComplexToPhaseImageFilter.h>
#include <itkMagnitudeAndPhaseToComplexImageFilter.h>
#include <itkForwardFFTImageFilter.h>
#include <itkInverseFFTImageFilter.h>
#include <itkFFTShiftImageFilter.h>


namespace ASTex
{

namespace Fourier
{


template<typename MASK>
void randomPhaseShift(ImageSpectrald& phase, const MASK& mask)
{
	// random shifting of the phases

	// phase.pixelAbsolute(0,0)=0 and must not be changed (freq 0 is a real constant).
	// we don't know yet what to do with the 3 remaining frequencies that have no symmetric counterpart : (W/2,0) , (0,H/2)  and (W/2,H/2)
	// thay seem to be real too (phase = 0) though I don't know what they are exactly.
	// until new insight we don't shift them

	int W = phase.width();
	int H = phase.height();

	int WD = phase.itk()->GetLargestPossibleRegion().GetIndex()[0];
	int HD = phase.itk()->GetLargestPossibleRegion().GetIndex()[1];

	int WC = W/2 + WD;
	int HC = H/2 + HD;

	srand (time(NULL));
	for (int j=1; j< H/2;++j)
	{
		for (int i=1; i< W/2;++i)
		{
			if (mask(WC+i,HC+j))
			{
				double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
				shiftDouble(phase.pixelAbsolute(WC+i,HC+j), r);
				shiftDouble(phase.pixelAbsolute(WC-i,HC-j), -r);
			}

			if (mask(WC+i,HC-j))
			{
				double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
				shiftDouble(phase.pixelAbsolute(WC+i,HC-j), r);
				shiftDouble(phase.pixelAbsolute(WC-i,HC+j), -r);
			}
		}
	}

	for (int i=1; i< W/2;++i)
	{
		if (mask(WC+i,HC))
		{
			double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
			shiftDouble(phase.pixelAbsolute(WC+i,HC), r);
			shiftDouble(phase.pixelAbsolute(WC-i,HC), -r);
		}

		if (mask(WC+i,HD))
		{
			double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
			shiftDouble(phase.pixelAbsolute(WC+i,HD), r);
			shiftDouble(phase.pixelAbsolute(WC-i,HD), -r);
		}
	}

	for (int j=1; j< H/2;++j)
	{
		if (mask(WC,HC+j))
		{
			double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
			shiftDouble(phase.pixelAbsolute(WC,HC+j), r);
			shiftDouble(phase.pixelAbsolute(WC,HC-j), -r);
		}

		if (mask(WD,HC+j))
		{
			double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
			shiftDouble(phase.pixelAbsolute(WD,HC+j), r);
			shiftDouble(phase.pixelAbsolute(WD,HC-j), -r);
		}
	}
}


template<typename MASK>
void randomPhase(ImageSpectrald& phase, const MASK& mask, bool call_srand)
{
	// random value of the phases

	// phase.pixelAbsolute(0,0)=0 and must not be changed (freq 0 is a real constant).
	// we don't know yet what to do with the 3 remaining frequencies that have no symmetric counterpart : (W/2,0) , (0,H/2)  and (W/2,H/2)
	// thay seem to be real too (phase = 0) though I don't know what they are exactly.
	// until new insight we don't change them

	int W = phase.width();
	int H = phase.height();

	int WD = phase.itk()->GetLargestPossibleRegion().GetIndex()[0];
	int HD = phase.itk()->GetLargestPossibleRegion().GetIndex()[1];

	int WC = W/2 + WD;
	int HC = H/2 + HD;

	if (call_srand) srand (time(NULL));
	for (int j=1; j< H/2;++j)
	{
		for (int i=1; i< W/2;++i)
		{
			if (mask(WC+i,HC+j))
			{
				double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
				phase.pixelAbsolute(WC+i,HC+j) = r;
				phase.pixelAbsolute(WC-i,HC-j) = -r;
			}

			if (mask(WC+i,HC-j))
			{
				double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
				phase.pixelAbsolute(WC+i,HC-j) = r;
				phase.pixelAbsolute(WC-i,HC+j) = -r;
			}
		}
	}

	for (int i=1; i< W/2;++i)
	{
		if (mask(WC+i,HC))
		{
			double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
			phase.pixelAbsolute(WC+i,HC) = r;
			phase.pixelAbsolute(WC-i,HC) = -r;
		}

		if (mask(WC+i,HD))
		{
			double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
			phase.pixelAbsolute(WC+i,HD) = r;
			phase.pixelAbsolute(WC-i,HD) = -r;

		}
	}

	for (int j=1; j< H/2;++j)
	{
		if (mask(WC,HC+j))
		{
			double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
			phase.pixelAbsolute(WC,HC+j) = r;
			phase.pixelAbsolute(WC,HC-j) = -r;
		}

		if (mask(WD,HC+j))
		{
			double r = M_PI*((double)(rand())*2.0/RAND_MAX-1);
			phase.pixelAbsolute(WD,HC+j) = r;
			phase.pixelAbsolute(WD,HC-j) = -r;
		}
	}
}



template<typename MASK>
int autoCorrelation_full_size(const ImageGrayd& input, ImageGrayd& acorr, const MASK& mask)
{
	const int im_w = input.width();
	const int im_h = input.height();
	const int T = acorr.height();
	assert (T==acorr.width() && T==acorr.height() && T==input.width() && T==input.height());

	int64_t mask_size = 0;
	for (int y = 0; y < im_h; ++y)
		for (int x=0; x < im_w; ++x)
			if (mask(x,y))
				++mask_size;

	// compute auto-correlation
	acorr.for_all_pixels([&] (double& P, int dx, int dy)
	{
		double v = 0;

		input.for_all_pixels([&] (int x, int y)
		{
			if(mask(x,y) && mask((x+dx)%im_w,(y+dy)%im_h))
			{
				v += input.pixelAbsolute(x,y) * input.pixelAbsolute((x+dx)%im_w,(y+dy)%im_h);
			}
		});

		P = v / mask_size * T;
	});

	return mask_size;
}


template<typename MASK>
int true_prop(const ImageGrayd& input, ImageGrayd& acorr, int block_x, int block_y, const MASK& mask)
{
	assert (acorr.width() == acorr.height()
			&& block_x+acorr.width()<=input.width() && block_y+acorr.height()<=input.height());

	std::ignore = input; // to avoid warning in release

	// compute auto-correlation
	long int mask_size = 0;
	acorr.for_all_pixels([&](double& /*P*/, int dx, int dy)
	{
		if (mask(block_x+dx,block_y+dy))
		{
			++mask_size;
		}
	});

	return mask_size;
}

template<typename MASK>
int true_prop(int block_x, int block_y, int taille_x, int taille_y, const MASK& mask)
{
	int mask_size = 0;
	for (int x = 0; x < taille_x; ++x) {
		for (int y = 0; y < taille_y; ++y) {
			if (mask(block_x+x,block_y+y))
			{
				++mask_size;
			}
		}
	}
	return mask_size;
}

template<typename MASK>
void true_prop_precompute(int mask_size_x, int mask_size_y, int block_size_x, int block_size_y, const MASK& mask, std::vector<std::vector<int>>& v)
{
	//init v
	v.resize(mask_size_x-block_size_x+1);
	for (int x = 0; x < mask_size_x-block_size_x+1; ++x)
		v[x].resize(mask_size_y-block_size_y+1);

	for (int block_x = 0; block_x + block_size_x <= mask_size_x ; ++block_x)
	{
		v[block_x][0]=true_prop(block_x,0,block_size_x,block_size_y,mask);
		for (int block_y = 1; block_y + block_size_y <= mask_size_y ; ++block_y)
		{
			v[block_x][block_y]= v[block_x][block_y-1]
					-true_prop(block_x,block_y-1,block_size_x,1,mask)
					+true_prop(block_x,block_y+block_size_y-1,block_size_x,1,mask);
		}
	}
}

template<typename MASK>
int autoCorrelation_block_backup(const ImageGrayd& input, ImageGrayd& acorr, int block_x, int block_y, const MASK& mask)
{
	const int T = acorr.height();
	assert (acorr.width() == acorr.height() && block_x+T<=input.width() && block_y+T<=input.height());

	// compute auto-correlation
	int mask_size = true_prop(input,acorr,block_x,block_y,mask);

	if (mask_size==0)
		return 0;

	acorr.for_all_pixels([&](double& P, int dx, int dy)
	{
		double v = 0;
		for (int y = block_y; y < block_y+T; ++y)
		{
			for (int x=block_x; x < block_x+T; ++x)
			{
				const int xx = block_x + (x-block_x+dx)%T;
				const int yy = block_y + (y-block_y+dy)%T;
				if(mask(x,y) && mask(xx,yy))
				{
					v += input.pixelAbsolute(x,y) * input.pixelAbsolute(xx,yy);
				}
			}
		}
		P = v / mask_size * T;
	});

	return mask_size;
}

template<typename MASK>
double autoCorrelation_pointwise(const ImageGrayd& input, int block_x, int block_y, const MASK& mask, int dx, int dy, int Tx, int Ty)
{
	double v = 0;
	for (int y = block_y; y < block_y+Ty; ++y)
	{
		for (int x=block_x; x < block_x+Tx; ++x)
		{
			const int xx = block_x + (x-block_x+dx)%Tx;
			const int yy = block_y + (y-block_y+dy)%Ty;
			if(mask(x,y) && mask(xx,yy))
			{
				v += input.pixelAbsolute(x,y) * input.pixelAbsolute(xx,yy);
			}
		}
	}
	return v;
}

template<typename MASK>
double autoCorrelation_pointwise_column(const ImageGrayd& input, int block_x, int block_y, const MASK& mask, int dx, int dy, int Tx, int Ty, int column_y)
{
	double v = 0;
	int y = column_y;
	for (int x=block_x; x < block_x+Tx; ++x)
	{
		const int xx = block_x + (x-block_x+dx)%Tx;
		const int yy = block_y + (y-block_y+dy)%Ty;
		if(mask(x,y) && mask(xx,yy))
		{
			v += input.pixelAbsolute(x,y) * input.pixelAbsolute(xx,yy);
		}
	}
	return v;
}



template<typename MASK>
int autoCorrelation_block(const ImageGrayd& input, ImageGrayd& acorr, int block_x, int block_y, const MASK& mask)
{
	const int T = acorr.height();
	assert (acorr.width() == acorr.height() && block_x+T<=input.width() && block_y+T<=input.height());

	// compute auto-correlation
	int mask_size = true_prop(input,acorr,block_x,block_y,mask);

	if (mask_size==0)
		return 0;

	acorr.for_all_pixels([&](double& P, int dx, int dy)
	{
		P = autoCorrelation_pointwise(input, block_x, block_y, mask, dx, dy, T, T);
	});

	autoCorrelation_normalize(acorr, double(T) / double(mask_size));

	return mask_size;
}



template<typename MASK>
int autoCorrelation_block_3chan(const ImageGrayd& input_1,const ImageGrayd& input_2,const ImageGrayd& input_3, ImageGrayd& acorr_1, ImageGrayd& acorr_2, ImageGrayd& acorr_3, int block_x, int block_y, const MASK& mask)
{
	const int T = acorr_1.height();
	assert (acorr_1.width() == acorr_1.height() && block_x+T<=input_1.width() && block_y+T<=input_1.height());
	assert (acorr_2.width() == acorr_2.height() && block_x+T<=input_2.width() && block_y+T<=input_2.height());
	assert (acorr_3.width() == acorr_3.height() && block_x+T<=input_3.width() && block_y+T<=input_3.height());

	// compute auto-correlation
	int mask_size = true_prop(input_1,acorr_1,block_x,block_y,mask);

	if (mask_size==0)
		return 0;

	for (int x = 0; x < T; ++x)
	{
		for (int y = 0; y < T; ++y)
		{
			acorr_1.pixelAbsolute(x,y) = autoCorrelation_pointwise(input_1, block_x, block_y, mask, x, y, T, T);
			acorr_2.pixelAbsolute(x,y) = autoCorrelation_pointwise(input_2, block_x, block_y, mask, x, y, T, T);
			acorr_3.pixelAbsolute(x,y) = autoCorrelation_pointwise(input_3, block_x, block_y, mask, x, y, T, T);
		}
	}

	autoCorrelation_normalize_3chan(acorr_1,acorr_2,acorr_3, double(T) / double(mask_size));

	return mask_size;
}



template<typename MASK>
int autoCorrelation_block_update(const ImageGrayd& input, ImageGrayd& acorr, int block_x, int block_y, const MASK& mask)
{ // passage de y-1 à y
	// precond : acorr est remplie avec l'AC du bloc juste à gauche
	const int T = acorr.height();
	assert (acorr.width() == acorr.height() && block_x+T<=input.width() && block_y+T<=input.height());

	// compute auto-correlation
	int mask_size = true_prop(input,acorr,block_x,block_y,mask);

	if (mask_size==0) // TODO warning
		return 0;

	acorr.for_all_pixels([&](double& P, int dx, int dy)
	{
		P -= autoCorrelation_pointwise_column(input, block_x, block_y-1, mask, dx, dy, T, T, block_y-1);
		P += autoCorrelation_pointwise_column(input, block_x, block_y,   mask, dx, dy, T, T, block_y+T-1);

		if (dy != 0)
		{
			P -= autoCorrelation_pointwise_column(input, block_x, block_y-1, mask, dx, dy, T, T, block_y+T-1-dy);
			P += autoCorrelation_pointwise_column(input, block_x, block_y,   mask, dx, dy, T, T, block_y+T-1-dy);
		}
	});

	autoCorrelation_normalize(acorr, double(T) / double(mask_size));

	return mask_size;
}



template<typename MASK>
int autoCorrelation_block_update_3chan(const ImageGrayd& input_1,const ImageGrayd& input_2,const ImageGrayd& input_3, ImageGrayd& acorr_1, ImageGrayd& acorr_2, ImageGrayd& acorr_3, int block_x, int block_y, const MASK& mask)
{ // passage de y-1 à y
	// precond : acorr est remplie avec l'AC du bloc juste à gauche
	const int T = acorr_1.height();
	assert (acorr_1.width() == acorr_1.height() && block_x+T<=input_1.width() && block_y+T<=input_1.height());
	assert (acorr_2.width() == acorr_2.height() && block_x+T<=input_2.width() && block_y+T<=input_2.height());
	assert (acorr_3.width() == acorr_3.height() && block_x+T<=input_3.width() && block_y+T<=input_3.height());

	// compute auto-correlation
	int mask_size = true_prop(input_1,acorr_1,block_x,block_y,mask);

	if (mask_size==0) // TODO warning
		return 0;

	for (int x = 0; x < T; ++x)
	{
		for (int y = 0; y < T; ++y)
		{
			acorr_1.pixelAbsolute(x,y) -= autoCorrelation_pointwise_column(input_1, block_x, block_y-1, mask, x, y, T, T, block_y-1);
			acorr_2.pixelAbsolute(x,y) -= autoCorrelation_pointwise_column(input_2, block_x, block_y-1, mask, x, y, T, T, block_y-1);
			acorr_3.pixelAbsolute(x,y) -= autoCorrelation_pointwise_column(input_3, block_x, block_y-1, mask, x, y, T, T, block_y-1);

			acorr_1.pixelAbsolute(x,y) += autoCorrelation_pointwise_column(input_1, block_x, block_y,   mask, x, y, T, T, block_y+T-1);
			acorr_2.pixelAbsolute(x,y) += autoCorrelation_pointwise_column(input_2, block_x, block_y,   mask, x, y, T, T, block_y+T-1);
			acorr_3.pixelAbsolute(x,y) += autoCorrelation_pointwise_column(input_3, block_x, block_y,   mask, x, y, T, T, block_y+T-1);

			if (y != 0)
			{
				acorr_1.pixelAbsolute(x,y) -= autoCorrelation_pointwise_column(input_1, block_x, block_y-1, mask, x, y, T, T, block_y+T-1-y);
				acorr_2.pixelAbsolute(x,y) -= autoCorrelation_pointwise_column(input_2, block_x, block_y-1, mask, x, y, T, T, block_y+T-1-y);
				acorr_3.pixelAbsolute(x,y) -= autoCorrelation_pointwise_column(input_3, block_x, block_y-1, mask, x, y, T, T, block_y+T-1-y);

				acorr_1.pixelAbsolute(x,y) += autoCorrelation_pointwise_column(input_1, block_x, block_y,   mask, x, y, T, T, block_y+T-1-y);
				acorr_2.pixelAbsolute(x,y) += autoCorrelation_pointwise_column(input_2, block_x, block_y,   mask, x, y, T, T, block_y+T-1-y);
				acorr_3.pixelAbsolute(x,y) += autoCorrelation_pointwise_column(input_3, block_x, block_y,   mask, x, y, T, T, block_y+T-1-y);
			}
		}
	}

	autoCorrelation_normalize_3chan(acorr_1,acorr_2,acorr_3, double(T) / double(mask_size));

	return mask_size;
}





template<typename MASK>
int autoCorrelation_small_size_backup(const ImageGrayd& input, ImageGrayd& acorr, const MASK& mask)
{
	const int im_w = input.width();
	const int im_h = input.height();
	const int T = acorr.height();
	assert (acorr.width() == acorr.height() && acorr.width()<=input.width() && acorr.height()<=input.height());

	// compute auto-correlation
	for (int dy = 0; dy < T; ++dy)
	{
		for (int dx = 0; dx < T; ++dx)
		{
			acorr.pixelAbsolute(dx,dy) = 0.0;
			int nb_blocs = 0;
			for (int bloc_y=0; bloc_y+T<=im_h; bloc_y+=T)
			{
				for (int bloc_x=0; bloc_x+T<=im_w; bloc_x+=T)
				{
					++nb_blocs;
					double v = 0;
					int64_t mask_size = 0;
					for (int y = bloc_y; y < bloc_y+T; ++y)
					{
						for (int x=bloc_x; x < bloc_x+T; ++x)
						{
							if (mask(x,y))
							{
								++mask_size;
							}
							const int xx = bloc_x + (x-bloc_x+dx)%T;
							const int yy = bloc_y + (y-bloc_y+dy)%T;
							if(mask(x,y) && mask(xx,yy))
							{
								v += input.pixelAbsolute(x,y) * input.pixelAbsolute(xx,yy);
							}
						}
					}
					acorr.pixelAbsolute(dx,dy) += v / mask_size;
				}
			}
			acorr.pixelAbsolute(dx,dy) /= nb_blocs;
			acorr.pixelAbsolute(dx,dy) *= T;
		}
	}

	return 0;
}


template<typename MASK>
void spectrum_by_autocorrelation_full_size(const ImageGrayd& input, ImageSpectrald& modulus, const MASK& mask)
{
	ImageGrayd acorr(input.width(),input.height());
	autoCorrelation_full_size(input, acorr, mask);

	ImageSpectrald phase;
	fftForwardModulusAndPhase(acorr, modulus, phase);

	modulus.for_all_pixels([](double& p)
	{
		p = std::sqrt(p);
	});
}



template<typename MASK>
int spectrum_by_autocorrelation_small_size_backup(const ImageGrayd& input, ImageSpectrald& modulus, const MASK& mask, double proportion_threshold, int step)
{
	assert (modulus.width()==modulus.height());
	int T = modulus.width();

	ImageGrayd acorr(T,T);
	ImageSpectrald accu(T,true);
	ImageSpectrald phase;

	int nb_blocks=0;
	int64_t weight = 0;
	// ATTENTION, PENSER A REMPLACER 1 par ZERO une fois terminé et la SQRT et step = T

	std::vector<std::vector<int>> v;
	true_prop_precompute(input.width(),input.height(),T,T,mask,v);

	for (int block_y=0; block_y+T<=input.height(); block_y+=step)
	{
		for (int block_x=0; block_x+T<=input.width(); block_x+=step)
		{
			int mask_size =v[block_x][block_y];
			if (mask_size >= T*T*proportion_threshold && mask_size>0)
			{
				++nb_blocks;
				weight += mask_size;
				autoCorrelation_block_backup(input,acorr,block_x,block_y,mask);
				fftForwardModulusAndPhase(acorr, modulus, phase);
				accu.for_all_pixels([&](double& P, int x, int y)
				{
					P += modulus.pixelAbsolute(x,y) * mask_size;
				});
			}
		}
	}

	modulus.for_all_pixels([&](double& P, int x, int y)
	{
		P = std::sqrt( accu.pixelAbsolute(x,y) / weight) ;
	});
	return nb_blocks;
}


template<typename MASK>
int spectrum_by_autocorrelation_small_size(const ImageGrayd& input, ImageSpectrald& modulus, const MASK& mask, double proportion_threshold, int step)
{
	assert (modulus.width()==modulus.height());
	int T = modulus.width();
	ImageGrayd acorr(T,T);
	ImageSpectrald accu(T,true);
	ImageSpectrald phase;

	int32_t nb_blocks=0;
	int64_t weight = 0;
	// ATTENTION, PENSER A REMPLACER 1 par ZERO une fois terminé et la SQRT et step = T

	std::vector<std::vector<int>> v;
	true_prop_precompute(input.width(),input.height(),T,T,mask,v);

	bool recompute_from_scratch = true;

	for (int32_t block_x=0; block_x+T<=input.width(); block_x+=step)
	{
		for (int32_t block_y=0; block_y+T<=input.height(); block_y+=step)
		{
			int32_t mask_size =v[block_x][block_y];
			if (mask_size >= T*T*proportion_threshold && mask_size>0)
			{
				++nb_blocks;
				weight += mask_size;
				if (recompute_from_scratch)
				{
					autoCorrelation_block(input,acorr,block_x,block_y,mask);
				}
				else
				{
					autoCorrelation_block_update(input,acorr,block_x,block_y,mask);
				}
				fftForwardModulusAndPhase(acorr, modulus, phase);
				accu.for_all_pixels([&](double& P, int x, int y)
				{
					P += modulus.pixelAbsolute(x,y) * mask_size;
				});
				autoCorrelation_normalize(acorr,double(mask_size)/double(T));
				recompute_from_scratch=false;
			}
			else
			{
				recompute_from_scratch = true;
			}
		}
		recompute_from_scratch = true;
	}

	modulus.for_all_pixels([&](double& P, int x, int y)
	{
		P = std::sqrt( accu.pixelAbsolute(x,y) / weight) ;
	});
	return nb_blocks;
}


template<typename MASK>
int spectrum_by_autocorrelation_3chan_Nmasks(const std::vector<std::vector<ImageGrayd>>& coords_per_mask, std::vector<std::vector<ImageSpectrald>>& modulus, const std::vector<MASK>& mask, double proportion_threshold, int step)
{
	float im_w = coords_per_mask[0][0].width();
	float im_h = coords_per_mask[0][0].height();

	int32_t T = modulus[0][0].width();

	std::vector<std::vector<ImageGrayd>> acorr;
	std::vector<std::vector<ImageSpectrald>> accu;

	acorr.resize( 3 , std::vector<ImageGrayd>( mask.size() ) );
	accu.resize( 3 , std::vector<ImageSpectrald>( mask.size() ) );

	for (int32_t i = 0; i < 3; ++i)
	{
		for (uint32_t j = 0; j < mask.size()	; ++j)
		{
			acorr[i][j].initItk(T,T,true);
			accu[i][j].initItk(T,T,true);
		}
	}

	ImageSpectrald phase;

	std::vector<std::vector<std::vector<int>>> v;
	v.resize( mask.size() , std::vector<std::vector<int>>(T) );

	std::vector<int>  weight;
	weight.resize(mask.size());

	std::vector<bool>  recompute_from_scratch;
	recompute_from_scratch.resize(mask.size());

	for (int32_t i = 0; i < mask.size(); ++i)
	{
		true_prop_precompute(im_w,im_h,T,T,mask[i],v[i]);
		weight[i]=0;
		recompute_from_scratch[i]=true;
	}

	for (int32_t block_x=0; block_x+T<=im_w; block_x+=step)
	{
		for (int32_t block_y=0; block_y+T<=im_h; block_y+=step)
		{
			for (int i = 0; i < mask.size(); ++i)
			{
				int32_t mask_size =v[i][block_x][block_y];
				if (mask_size >= T*T*proportion_threshold && mask_size>0)
				{
					weight[i] += mask_size;
					if (recompute_from_scratch[i])
					{
						autoCorrelation_block_3chan(coords_per_mask[0][i],coords_per_mask[1][i],coords_per_mask[2][i],acorr[0][i],acorr[1][i],acorr[2][i],block_x,block_y,mask[i]);
					}
					else
					{
						autoCorrelation_block_update_3chan(coords_per_mask[0][i],coords_per_mask[1][i],coords_per_mask[2][i],acorr[0][i],acorr[1][i],acorr[2][i],block_x,block_y,mask[i]);
					}
					fftForwardModulusAndPhase(acorr[0][i], modulus[0][i], phase);
					fftForwardModulusAndPhase(acorr[1][i], modulus[1][i], phase);
					fftForwardModulusAndPhase(acorr[2][i], modulus[2][i], phase);

					for (int32_t y = 0; y < T; ++y)
					{
						for (int32_t x = 0; x < T; ++x)
						{
							accu[0][i].pixelAbsolute(x,y) += modulus[0][i].pixelAbsolute(x,y) * mask_size;
							accu[1][i].pixelAbsolute(x,y) += modulus[1][i].pixelAbsolute(x,y) * mask_size;
							accu[2][i].pixelAbsolute(x,y) += modulus[2][i].pixelAbsolute(x,y) * mask_size;
						}
					}
					autoCorrelation_normalize_3chan(acorr[0][i],acorr[1][i],acorr[2][i],double(mask_size)/double(T));
					recompute_from_scratch[i]=false;

				}
				else
				{
					recompute_from_scratch[i] = true;
				}
			}
		}
		for (int i = 0; i < mask.size(); ++i)
		{
			recompute_from_scratch[i] = true;
		}
	}

	for (int i = 0; i < mask.size(); ++i)
	{
		for (int y = 0; y < T; ++y)
		{
			for (int x = 0; x < T; ++x)
			{
				modulus[0][i].pixelAbsolute(x,y) = std::sqrt( accu[0][i].pixelAbsolute(x,y) / weight[i]) ;
				modulus[1][i].pixelAbsolute(x,y) = std::sqrt( accu[1][i].pixelAbsolute(x,y) / weight[i]) ;
				modulus[2][i].pixelAbsolute(x,y) = std::sqrt( accu[2][i].pixelAbsolute(x,y) / weight[i]) ;
			}
		}
	}

	return 0;
}







/////////////////////////////////////////////////////////////////

template<typename MASK>
void spectrum_by_autocorrelation_small_size_backup(const ImageGrayd& input, ImageSpectrald& modulus, const MASK& mask)
{
	ImageGrayd acorr(modulus.width(),modulus.height());
	autoCorrelation_small_size_backup(input,acorr, mask);

	ImageSpectrald phase;
	fftForwardModulusAndPhase(acorr, modulus, phase);

	modulus.for_all_pixels([](double& p)
	{
		p = std::sqrt(p);
	});

}


template<typename MASK>
double getPower (const ImageGrayd& input, const MASK& mask)
{
	double power = 0.0;
	int64_t nb_pix = 0;
	input.for_all_pixels([&](const double& p)
	{
		++nb_pix;
		power += p*p;
	},mask);


	power /= double(nb_pix);
	return power;
}



} // end namespace Fourier
} // end namespace ASTex
