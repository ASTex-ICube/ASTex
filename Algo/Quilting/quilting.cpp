#include "quilting.h"
#include <ASTex/colorspace_filters.h>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iomanip>

using namespace ASTex;



int Quilting::random_int(int min, int max)
{
	double d = double(std::rand())/RAND_MAX;
	return  int(d*(max-min) + min);
}



void Quilting::loadInput(const std::string& imgName)
{
	ImageRGBu8 img;
	img.load(imgName);
	ColorSpace::FilterRGB255To01<ImageRGBu8::ItkImg,Quilting::ImageType::ItkImg>::Pointer filter =
			ColorSpace::FilterRGB255To01<ImageRGBu8::ItkImg,Quilting::ImageType::ItkImg>::New();
	filter->SetInput(img.itk());
	filter->Update();
	img_orig_.itk() = filter->GetOutput();
}

void Quilting::saveOutput(const std::string& imgName)
{
	ColorSpace::FilterRGB01To255<Quilting::ImageType::ItkImg,ImageRGBu8::ItkImg>::Pointer filter =
			ColorSpace::FilterRGB01To255<Quilting::ImageType::ItkImg,ImageRGBu8::ItkImg>::New();
	filter->SetInput(img_generated_.itk());
	ImageRGBu8 img(filter->GetOutput());
	img.save(imgName);
}




void Quilting::init(int tw, int ovl, int rc, int nbrt)
{
//	std::srand(std::time(0));

	std::srand(0);

	tile_width_ = tw;
	overlap_width_ = ovl;
	multi_mins_.nbMaxElts_ = rc;

	overlap_h_sz_[0] = tw;
	overlap_h_sz_[1] = ovl;

	overlap_v_sz_[0] = ovl;
	overlap_v_sz_[1] = tw;

	tile_sz_[0] = tw;
	tile_sz_[1] = tw;

	bufferPositions_.resize(nbrt);
	bufferErrors_.resize(nbrt);

}


double Quilting::computeErrorOverlap(const Region& or_reg, const Region& gen_reg)
{
	double total=0.0;

	auto it2 = img_generated_.beginIterator(gen_reg);
	for (auto it = img_orig_.beginIterator(or_reg); !it.IsAtEnd(); ++it, ++it2)
	{
		double dR = it.Get()[0] - it2.Get()[0];
		double dG = it.Get()[1] - it2.Get()[1];
		double dB = it.Get()[2] - it2.Get()[2];
		double dP2 = dR*dR + dG*dG + dB*dB; // sqrt ??
		total += dP2;
	}

	return total/(or_reg.GetSize()[0]*or_reg.GetSize()[1]);
}



double Quilting::computeErrorTile(const Index& or_pos, const Index& gen_pos)
{
	double res = 0.0;
	// Vertical overlapping (except on left)
	if (gen_pos[0]>=overlap_width_)
	{
		Region r1;
		r1.SetIndex(or_pos);
		r1.SetSize(overlap_v_sz_);

		Region r2;
		r2.SetIndex(gen_pos);
		r2.SetSize(overlap_v_sz_);

		res += computeErrorOverlap(r1,r2);
	}

	// horizontal overlapping (except on top)
	if (gen_pos[1]>=overlap_width_)
	{
		Region r3;
		r3.SetIndex(or_pos);
		r3.SetSize(overlap_h_sz_);

		Region r4;
		r4.SetIndex(gen_pos);
		r4.SetSize(overlap_h_sz_);

		res += computeErrorOverlap(r3,r4);
	}


	return res;
}


Index Quilting::bestFittingTile(const Index& gen_pos, double error_max)
{
	if (gen_pos[0]+gen_pos[1] == 0)
	{
		Index best;
		best[0]=random_int(0,img_orig_.width()-tile_width_);
		best[1]=random_int(0,img_orig_.height()-tile_width_);

		return best;
	}

	Region cropped_original;
	cropped_original.SetIndex(0,0);
	cropped_original.SetIndex(1,0);
	cropped_original.SetSize(0,img_orig_.width()-tile_width_);
	cropped_original.SetSize(1,img_orig_.height()-tile_width_);

	multi_mins_.reset();

//	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(cropped_original); !it.IsAtEnd(); ++it)
//	{
//		const Index& localPos =  it.GetIndex();
//		double err = computeErrorTile(localPos, gen_pos);
//		if (err < error_max)
//		{
//			multi_mins_.try_add(localPos,err);
//		}
//	}


	//
	//
	std::size_t nbt = bufferPositions_.size();

	for (std::size_t i=0; i<nbt;++i)
	{
		bufferPositions_[i][0] = random_int(0,img_orig_.width()-tile_width_);
		bufferPositions_[i][1] = random_int(0,img_orig_.height()-tile_width_);
	}



	#pragma omp parallel for
	for (long i=0; i<long(nbt);++i)
		bufferErrors_[i] = computeErrorTile(bufferPositions_[i], gen_pos);

	for (std::size_t i=0; i<nbt;++i)
	{
		if (bufferErrors_[i] < error_max)
		{
			multi_mins_.try_add(bufferPositions_[i],bufferErrors_[i]);
		}
	}

	// take one of the best
//	int r = random() % (multi_mins_.nbMaxElts_) ;
	int r = random_int(0,multi_mins_.nbMaxElts_);


	std::multimap<double,Index>::iterator it = multi_mins_.mpos_.begin();
	for (int i=0;i<r;++i)
		++it;

	return it->second;
}



void Quilting::copy(const Region& from, const Region& to)
{
	auto it2 = img_generated_.beginIterator(to);
	for (auto it = img_orig_.beginIterator(from); !it.IsAtEnd(); ++it, ++it2)
	{
		it2.Set(it.Get());
	}
}

void Quilting::computeRandom(int width, int height, int tw )
{
	// ensure image covered by complete tiles
	width = (width /  tw )* tw ;
	height = (height / tw ) * tw ;


	img_generated_.initItk(width,height);
	init(tw,0,0,0);

	Region from;
	from.SetSize(tile_sz_);
	Region to;
	to.SetSize(tile_sz_);

	for (int i = 0; i<height; i+= tw)
	{
		for (int j = 0; j<width; j+= tw)
		{
			to.SetIndex(0,j);
			to.SetIndex(1,i);
			from.SetIndex(0,random_int(0,img_orig_.width()-tw));
			from.SetIndex(1,random_int(0,img_orig_.height()-tw));
			copy(from,to);
		}
	}
}

void Quilting::computeFittingTiles(int width, int height, int tw,  int ovl, int rc, int nbrt)
{
	// ensure image covered by complete tiles
	width = (width / ( tw - ovl) -1)* (tw - ovl ) + tw;
	height = (height / ( tw - ovl) -1)* (tw - ovl) + tw ;


	img_generated_.initItk(width,height);
	init(tw,ovl,rc,nbrt);

	Region from;
	from.SetSize(tile_sz_);
	Region to;
	to.SetSize(tile_sz_);

	int incr = tw - ovl;

	for (int i = 0; i<=height-tw; i+= incr)
	{
		for (int j = 0; j<=width-tw; j+= incr)
		{
			to.SetIndex(0,j);
			to.SetIndex(1,i);
			Index best = bestFittingTile(to.GetIndex());
			from.SetIndex(best);
			copy(from,to);
		}
	}

}


void Quilting::verticalPathCut(const Index& or_pos, const Index& gen_pos, std::vector<int>& minPos)
{
	// compute initial error
	ErrorCutPathBuffer e(overlap_width_,tile_width_);

	Region or_reg(or_pos,overlap_v_sz_);
	Region gen_reg(gen_pos,overlap_v_sz_);

	ImageType::IteratorIndexed it2 = img_generated_.beginIteratorIndexed(gen_reg);
	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(or_reg); !it.IsAtEnd(); ++it, ++it2)
	{
		double dR = it.Get()[0] - it2.Get()[0];
		double dG = it.Get()[1] - it2.Get()[1];
		double dB = it.Get()[2] - it2.Get()[2];
		double dP2 = dR*dR + dG*dG + dB*dB;
		e( it.GetIndex() - or_reg.GetIndex() ) = dP2;
	}

	// compute cumulative error
	ErrorCutPathBuffer E(overlap_width_,tile_width_);

	// first row
	for(int i=0;i<overlap_width_;++i)
		E(i,0) = e(i,0);

	// others rows
	for(int j=1;j<tile_width_;++j)
	{
		E(0,j) = e(0,j) + std::min(E(0,j-1),E(1,j-1));
		int i=1;
		while (i<overlap_width_-1)
		{
			E(i,j) = e(i,j) + std::min( std::min(E(i-1,j-1),E(i,j-1)), E(i+1,j-1)) ;
			++i;
		}
		E(i,j) = e(i,j) + std::min(E(i-1,j-1),E(i,j-1));
	}

	// up to store local position of min error cut
	int i=tile_width_-1;
	minPos[i] = E.minOfRow(i);
	while (i>0)
	{
		minPos[i-1] = E.minOfRow3(minPos[i],i-1);
		--i;
	}

//	std::cout << "POSCUT:" << std::endl;
//	for (int i=0;i<tile_width_;++i)
//	{
//		std::cout << minPos[i] << std::endl;
//	}

}


void Quilting::horizontalPathCut(const Index& or_pos, const Index& gen_pos, std::vector<int>& minPos)
{
	// compute initial error
	ErrorCutPathBuffer e(tile_width_,overlap_width_);

	Region or_reg(or_pos,overlap_h_sz_);
	Region gen_reg(gen_pos,overlap_h_sz_);

	ImageType::IteratorIndexed it2 = img_generated_.beginIteratorIndexed(gen_reg);
	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(or_reg); !it.IsAtEnd(); ++it, ++it2)
	{
		double dR = it.Get()[0] - it2.Get()[0];
		double dG = it.Get()[1] - it2.Get()[1];
		double dB = it.Get()[2] - it2.Get()[2];
		double dP2 = dR*dR + dG*dG + dB*dB;
		e( it.GetIndex() - or_reg.GetIndex() ) = dP2;
	}

	// compute cumulative error
	ErrorCutPathBuffer E(tile_width_,overlap_width_);
	// first row
	for(int i=0;i<overlap_width_;++i)
		E(0,i) = e(0,i);
	// others rows
	for(int j=1;j<tile_width_;++j)
	{
		E(j,0) = e(j,0) + std::min(E(j-1,0),E(j-1,1));
		int i=1;
		while (i<overlap_width_-1)
		{
			E(j,i) = e(j,i) + std::min( std::min(E(j-1,i-1),E(j-1,i)), E(j-1,i+1)) ;
			++i;
		}
		E(j,i) = e(j,i) + std::min(E(j-1,i-1),E(j-1,i));
	}

	// up pass to store local position of min cut
	int i=tile_width_-1;
	minPos[i] = E.minOfColumn(i);
	while (i>0)
	{
		minPos[i-1] = E.minOfColumn3(i-1,minPos[i]);
		--i;
	}
}


void Quilting::pathSlope(const std::vector<int>& pos, std::vector<float> slopes)
{
	int nb = pos.size();
	slopes.resize(nb);

	slopes[0] = pos[1] - pos[0];
	slopes[1] = 0.5f*(pos[2] - pos[0]);

	slopes[nb-1] = pos[nb-1] - pos[nb-2];
	slopes[nb-2] = 0.5f*(pos[nb-1] - pos[nb-3]);

	for (int i = 2; i<=nb-3; ++i)
	{
		slopes[i] = 1.0f/12.0f*(-pos[i+2]+8*pos[i+1]-8*pos[i-1]+pos[i-2]);
	}
}



void Quilting::copyPathCut(const Region& from, const Region& to, const std::vector<int>& leftMinCut,const std::vector<int>& upMinCut)
{
	const Index& origin = from.GetIndex();

	ImageType::IteratorIndexed it2 = img_generated_.beginIteratorIndexed(to);
	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(from); !it.IsAtEnd(); ++it, ++it2)
	{
		Offset off= it.GetIndex() - origin;


		if ( (off[0] == leftMinCut[off[1]]) || (off[1] == upMinCut[off[0]]))
		{
			itk::RGBPixel<double> pix;
			pix[1]=1.0; pix[1]=0.0; pix[2]=0;
			it2.Set(pix);
		}
		else
		if ( (off[0] >= leftMinCut[off[1]]) && (off[1] >= upMinCut[off[0]]))
		{
			it2.Set(it.Get());
		}
	}
}



void Quilting::copyPathCut_BlendDisk(const Region& from, const Region& to, const std::vector<int>& leftMinCut,const std::vector<int>& upMinCut, int rad)
{
	const Index& origin = from.GetIndex();

	ImageType::IteratorIndexed it2 = img_generated_.beginIteratorIndexed(to);
	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(from); !it.IsAtEnd(); ++it, ++it2)
	{
		Offset P= it.GetIndex() - origin;

		// simple distance computation (manhattan)
		int d1 = P[0] - leftMinCut[P[1]];
		int d2 = P[1] - upMinCut[P[0]];
		int d;
		if (d1<0 && d2<0)
			d = std::max(std::abs(d1),std::abs(d2));
		else
			d = std::min(std::abs(d1),std::abs(d2));


		if ( (P[0] > leftMinCut[P[1]]) && (P[1] > upMinCut[P[0]]))
		{	// we are inside cut path
			if ((d>rad) || (P[0]>=overlap_width_) || (P[1]>=overlap_width_))
				it2.Set(it.Get());
			else
			{
				// d=0 -> alpha=0.5
				// d=rad -> alpha=0.0
				double alpha = (1.0 - double(d)/rad) / 2.0;
				it2.Set(it2.Get()*alpha + it.Get()*(1.0-alpha));
			}

		}
		else if ( (P[0] < leftMinCut[P[1]]) || (P[1] < upMinCut[P[0]]))
		{// we are outside
			if (d<=rad)
			{
				// d=0 -> alpha=0.5
				// d=rad -> alpha=0.0
				double alpha = (1.0 -double(d)/rad) / 2.0;
				it2.Set(it.Get()*alpha + it2.Get()*(1.0-alpha));
			}
		}
		else
		{
//			itk::RGBPixel<double> pix;
//			pix[1]=1.0; pix[1]=0.0; pix[2]=0;
//			it2.Set(pix);
			it2.Set(it.Get()*0.5 + it2.Get()*0.5);
		}
	}
}

/*
void Quilting::copyPathCut_BlendDisk(const Region& from, const Region& to, const std::vector<int>& leftMinCut,const std::vector<int>& upMinCut, int rad)
{
	const Index& origin = from.GetIndex();

	ImageType::IteratorIndexed it2 = img_generated_.beginIteratorIndexed(to);
	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(from); !it.IsAtEnd(); ++it, ++it2)
	{
		itk::Offset<2> P= it.GetIndex() - origin;

		// simple distance computation (manhattan)
//		int d = std::min(std::abs(P[0] - leftMinCut[P[1]]),std::abs(P[1] - upMinCut[P[0]]));

		int d1 = P[0] - leftMinCut[P[1]];
		int d2 = P[1] - upMinCut[P[0]];
		int d;
		if (d1<0 && d2<0)
			d = std::max(-d1,-d2);
		else
			d = std::min(std::abs(d1),std::abs(d2));




		if ( (P[0] >= leftMinCut[P[1]]) && (P[1] >= upMinCut[P[0]]))
		{	// we are inside cut path
			it2.Set(it.Get());
		}
		else
		{// we are outside
			if (d<=rad)
			{
				// d=0 -> alpha=0 : cut
				// d=rad -> alpha=1.0
				double alpha = double(d)/rad;
				it2.Set(it.Get()*(1.0-alpha) + it2.Get()*alpha);
			}
		}
	}
}
*/



void Quilting::computeFittingTilesPathCut(int width, int height, int tw, int ovl, int rc, int nbrt)
{
	// ensure image covered by complete tiles
	width = (width / ( tw - ovl) -1)* (tw - ovl ) + tw;
	height = (height / ( tw - ovl) -1)* (tw - ovl) + tw ;

	img_generated_.initItk(width,height);
	init(tw,ovl,rc,nbrt);

	Region from;
	from.SetSize(tile_sz_);
	Region to;
	to.SetSize(tile_sz_);

	std::vector<int> minPosH;
	minPosH.resize(tile_width_);
	std::vector<int> minPosV;
	minPosV.resize(tile_width_);


	int incr = tw - ovl;

	for (int i = 0; i<=height-tw; i+= incr)
	{
		for (int j = 0; j<=width-tw; j+= incr)
		{
			to.SetIndex(0,j);
			to.SetIndex(1,i);
			Index best = bestFittingTile(to.GetIndex());
			from.SetIndex(best);

			if (j==0)
				minPosV.assign(tile_width_,0);
			else
				verticalPathCut(from.GetIndex(), to.GetIndex(), minPosV);

			if (i==0)
				minPosH.assign(tile_width_,0);
			else
				horizontalPathCut(from.GetIndex(), to.GetIndex(), minPosH);

//			copyPathCut(from,to,minPosV,minPosH);
			copyPathCut_BlendDisk(from,to,minPosV,minPosH,3);
		}
	}
}



