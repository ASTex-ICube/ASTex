/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/



#include "quilting_interp.h"
#include <ASTex/colorspace_filters.h>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iomanip>

using namespace ASTex;



int QuiltingInterp::random_int(int min, int max)
{
	double d = double(std::rand())/RAND_MAX;
	return  int(d*(max-min) + min);
}



void QuiltingInterp::loadInput(const std::string& imgName)
{
	ImageRGBu8 img;
	img.load(imgName);
	ColorSpace::FilterRGB255To01<ImageRGBu8::ItkImg,QuiltingInterp::ImageType::ItkImg>::Pointer filter =
			ColorSpace::FilterRGB255To01<ImageRGBu8::ItkImg,QuiltingInterp::ImageType::ItkImg>::New();
	filter->SetInput(img.itk());
	filter->Update();
	img_orig_.itk() = filter->GetOutput();
}

void QuiltingInterp::saveOutput(const std::string& imgName)
{
	ColorSpace::FilterRGB01To255<QuiltingInterp::ImageType::ItkImg,ImageRGBu8::ItkImg>::Pointer filter =
			ColorSpace::FilterRGB01To255<QuiltingInterp::ImageType::ItkImg,ImageRGBu8::ItkImg>::New();
	filter->SetInput(img_generated_.itk());
	ImageRGBu8 img(filter->GetOutput());
	img.save(imgName);

	filter->SetInput(img_dbg_interp.itk());
	ImageRGBu8 img2(filter->GetOutput());
	img2.save("/tmp/interp.png");
}




void QuiltingInterp::init(int tw, int ovl, int rc, int nbrt)
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


double QuiltingInterp::computeErrorOverlap(const Region& or_reg, const Region& gen_reg)
{
	double total=0.0;

	ImageType::IteratorIndexed it2 = img_generated_.beginIteratorIndexed(gen_reg);
	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(or_reg); !it.IsAtEnd(); ++it, ++it2)
	{
		double dR = it.Get()[0] - it2.Get()[0];
		double dG = it.Get()[1] - it2.Get()[1];
		double dB = it.Get()[2] - it2.Get()[2];
		double dP2 = dR*dR + dG*dG + dB*dB; // sqrt ??
		total += dP2;
	}

	return total/(or_reg.GetSize()[0]*or_reg.GetSize()[1]);
}



double QuiltingInterp::computeErrorTile(const Index& or_pos, const Index& gen_pos)
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


QuiltingInterp::Index QuiltingInterp::bestFittingTile(const Index& gen_pos, double error_max)
{
	if (gen_pos[0]+gen_pos[1] == 0)
	{
		Index best;
		best[0]=random_int(blend_rad_,img_orig_.width()-tile_width_-blend_rad_);
		best[1]=random_int(blend_rad_,img_orig_.height()-tile_width_-blend_rad_);

		return best;
	}

	multi_mins_.reset();

	std::size_t nbt = bufferPositions_.size();

	for (std::size_t i=0; i<nbt;++i)
	{
		bufferPositions_[i][0] = random_int(blend_rad_,img_orig_.width()-tile_width_-blend_rad_);
		bufferPositions_[i][1] = random_int(blend_rad_,img_orig_.height()-tile_width_-blend_rad_);
	}


	#pragma omp parallel for
	for (long i=0; i<nbt;++i)
		bufferErrors_[i] = computeErrorTile(bufferPositions_[i], gen_pos);

	for (std::size_t i=0; i<nbt;++i)
	{
		if (bufferErrors_[i] < error_max)
		{
			multi_mins_.try_add(bufferPositions_[i],bufferErrors_[i]);
		}
	}

	// take one of the best
	int r = random_int(0,multi_mins_.nbMaxElts_);


	std::multimap<double,Index>::iterator it = multi_mins_.mpos_.begin();
	for (int i=0;i<r;++i)
		++it;

	return it->second;
}





void QuiltingInterp::verticalPathCut(const Index& or_pos, const Index& gen_pos, std::vector<int>& minPos)
{
	// compute initial error
	ErrorCutPathBuffer e(overlap_width_,tile_width_+blend_rad_);

	Size sz = overlap_v_sz_;
	sz[1] += blend_rad_;
	Region or_reg(or_pos,sz);
	Region gen_reg(gen_pos,sz);

//	Region or_reg(or_pos,overlap_v_sz_);
//	Region gen_reg(gen_pos,overlap_v_sz_);

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
	ErrorCutPathBuffer E(overlap_width_,tile_width_ + blend_rad_);

	// first row
	for(int i=0;i<overlap_width_;++i)
		E(i,0) = e(i,0);

	// others rows
	for(int j=1;j<tile_width_ + blend_rad_;++j)
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
	int i=tile_width_ + blend_rad_ -1;
	minPos[i] = E.minOfRow(i);
	while (i>0)
	{
		minPos[i-1] = E.minOfRow3(minPos[i],i-1);
		--i;
	}
}


void QuiltingInterp::horizontalPathCut(const Index& or_pos, const Index& gen_pos, std::vector<int>& minPos)
{
	// compute initial error
	ErrorCutPathBuffer e(tile_width_+ blend_rad_,overlap_width_);

	Size sz = overlap_h_sz_;
	sz[0] += blend_rad_;
	Region or_reg(or_pos,sz);
	Region gen_reg(gen_pos,sz);

//	Region or_reg(or_pos,overlap_h_sz_);
//	Region gen_reg(gen_pos,overlap_h_sz_);

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
	ErrorCutPathBuffer E(tile_width_ + blend_rad_,overlap_width_);
	// first row
	for(int i=0;i<overlap_width_;++i)
		E(0,i) = e(0,i);
	// others rows
	for(int j=1;j<tile_width_ + blend_rad_;++j)
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
	int i=tile_width_+ blend_rad_-1;
	minPos[i] = E.minOfColumn(i);
	while (i>0)
	{
		minPos[i-1] = E.minOfColumn3(i-1,minPos[i]);
		--i;
	}
}



void QuiltingInterp::copyPathCut_BlendDisk(const Region& from_p, const Region& to_p, const std::vector<int>& leftMinCut,const std::vector<int>& upMinCut)
{
	Index origin = from_p.GetIndex();

	Index pt = to_p.GetIndex();
	Index pf = from_p.GetIndex();
	Size s = from_p.GetSize();

	int bl_rd = blend_rad_;

	if (pt[0]>=bl_rd)
	{
		pt[0] -= bl_rd;
		pf[0] -= bl_rd;
		s[0] += 2*bl_rd;
	}
	else
		s[0] += bl_rd;

	if (pt[1]>=bl_rd)
	{
		pt[1] -= bl_rd;
		pf[1] -= bl_rd;
		s[1] += 2*bl_rd;
	}
	else
		s[1] += bl_rd;

	Region from;
	from.SetIndex(pf);
	from.SetSize(s);

	Region to;
	to.SetIndex(pt);
	to.SetSize(s);

//	std::cout << "From: "<<from.GetIndex()[0] << ","<< from.GetIndex()[1] << ":"<< from.GetSize()[0] << ","<< from.GetSize()[1] ;
//	std::cout << "  To: "<<to.GetIndex()[0] << ","<< to.GetIndex()[1] << ":"<< to.GetSize()[0] << ","<< to.GetSize()[1] ;
//	std::cout << " Ori: "<<origin[0] << ","<< origin[1] << std::endl; ;


	ImageType::IteratorIndexed it3 = img_dbg_interp.beginIteratorIndexed(to);

	ImageType::IteratorIndexed it2 = img_generated_.beginIteratorIndexed(to);
	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(from); !it.IsAtEnd(); ++it, ++it2, ++it3)
	{
		itk::Offset<2> P= it.GetIndex() - origin;

		itk::Offset<2> Q= from.GetIndex() + from.GetSize() - it.GetIndex();

		// simple distance computation (manhattan)
		int d1 = P[0] - leftMinCut[P[1]];
		int d2 = P[1] - upMinCut[P[0]];

		if (P[0]<0)
			d2 = P[1] - upMinCut[0];
		if (P[1]<0)
			d1 = P[0] - leftMinCut[0];

		int d;
		if (d1<0 && d2<0)
			d = std::max(-d1,-d2);
		else if (d1>=0 && d2>=0)
			d = std::min(d1,d2);
		else if (d1<0)
			d = -d1;
		else
			d = -d2;


		double alpha = (1.0 - double(d)/blend_rad_) / 2.0;
		if (alpha < 0.0)
			alpha = 0.0;

		if ((P[1]<overlap_width_) && (P[1] > upMinCut[P[0]]) && (Q[0]<blend_rad_))
			alpha = (double(blend_rad_ - Q[0])/blend_rad_);

		if ((P[0]<overlap_width_) && (P[0] > upMinCut[P[1]]) && (Q[1]<blend_rad_))
			alpha = (double(blend_rad_ - Q[1])/blend_rad_);


		if ( (P[0] > leftMinCut[P[1]]) && (P[1] > upMinCut[P[0]]))
		{	// we are inside cut path
			it2.Set(it2.Get()*alpha + it.Get()*(1.0-alpha));

		}
		else
		{	// we are outside
			it2.Set(it.Get()*alpha + it2.Get()*(1.0-alpha));
		}

		itk::RGBPixel<double> pix;
		double beta = alpha + 0.5;

		if (d1<0 && d2<0)
		{
			pix[0]=1.0; pix[1]=0.0; pix[2]=0.0;
			it3.Set(pix*beta);
		}
		else if (d1<0 && d2>0)
		{
			pix[0]=0.0; pix[1]=1.0; pix[2]=0.0;
			it3.Set(pix*beta);
		}
		else if (d1>0 && d2<0)
		{
			pix[0]=0.0; pix[1]=0.0; pix[2]=1.0;
			it3.Set(pix*beta);
		}
		else if (d1>0 && d2>0)
		{
			pix[0]=1.0; pix[1]=0.0; pix[2]=1.0;
			it3.Set(pix*beta);
		}
		else if (d1==0 || d2==0)
			{
				pix[0]=1.0; pix[1]=1.0; pix[2]=1.0;
				it3.Set(pix*beta);
			}


	}
}




void QuiltingInterp::computeFittingTilesPathCutInterp(int width, int height, int tw, int ovl, int bl_rad ,int rc, int nbrt)
{
	// ensure image covered by complete tiles

	img_dbg_interp.initItk(width,height);

	img_generated_.initItk(width,height);
	init(tw,ovl,rc,nbrt);

	blend_rad_ = bl_rad;

	Region from;
	from.SetSize(tile_sz_);
	Region to;
	to.SetSize(tile_sz_);

	std::vector<int> minPosH;
	minPosH.resize(tile_width_+blend_rad_);
	std::vector<int> minPosV;
	minPosV.resize(tile_width_+blend_rad_);


	int K=0;

	int incr = tw - ovl;

	for (int i = 0; i<=height-tw-blend_rad_; i+= incr)
	{
//		if (K==1)
//			break;
//		else
//			K++;
		for (int j = 0; j<=width-tw-blend_rad_; j+= incr)
		{
			to.SetIndex(0,j);
			to.SetIndex(1,i);
			Index best = bestFittingTile(to.GetIndex());

			saveDebugOverlap(to.GetIndex(),best,K++);

			from.SetIndex(best);

			if (j==0)
				minPosV.assign(tile_width_,0);
			else
				verticalPathCut(from.GetIndex(), to.GetIndex(), minPosV);

			if (i==0)
				minPosH.assign(tile_width_,0);
			else
				horizontalPathCut(from.GetIndex(), to.GetIndex(), minPosH);

			copyPathCut_BlendDisk(from,to,minPosV,minPosH);
		}
	}
}

void QuiltingInterp::saveDebugOverlap(const Index& P, const Index& Q, int k)
{
	ImageRGBu8 img_out;

	img_out.initItk(2*overlap_width_ + tile_width_ + 16, tile_width_ + 8);

	for (ImageRGBu8::IteratorIndexed it = img_out.beginIteratorIndexed(); !it.IsAtEnd(); ++it)
		it.Set(RGBu8(0,0,0));

	Region rP;
	rP.SetIndex(P);
	rP.SetSize(0,overlap_width_);
	rP.SetSize(1,tile_width_);

	Region rQ;
	rQ.SetIndex(Q);
	rQ.SetSize(0,overlap_width_);
	rQ.SetSize(1,tile_width_);

	Region r1 = rP;
	r1.SetIndex(0,4);
	r1.SetIndex(1,4);

	Region r2 = rQ;
	r2.SetIndex(0,8+overlap_width_);
	r2.SetIndex(1,4);

	ColorSpace::fonctorRGB01To255<ImageType::PixelType,RGBu8> f2u;

	ImageRGBu8::IteratorIndexed it1 = img_out.beginIteratorIndexed(r1);
	for (ImageType::IteratorIndexed itP = img_generated_.beginIteratorIndexed(rP); !itP.IsAtEnd(); ++itP, ++it1)
		it1.Set(f2u(itP.Get()));

	ImageRGBu8::IteratorIndexed it2 = img_out.beginIteratorIndexed(r2);
	for (ImageType::IteratorIndexed itQ = img_orig_.beginIteratorIndexed(rQ); !itQ.IsAtEnd(); ++itQ, ++it2)
		it2.Set(f2u(itQ.Get()));

	rP.SetSize(1,overlap_width_);
	rP.SetSize(0,tile_width_);
	rQ.SetSize(1,overlap_width_);
	rQ.SetSize(0,tile_width_);

	r1 = rP;
	r1.SetIndex(0,12+2*overlap_width_);
	r1.SetIndex(1,4);

	r2 = rQ;
	r2.SetIndex(0,12+2*overlap_width_);
	r2.SetIndex(1,8+overlap_width_);

	it1 = img_out.beginIteratorIndexed(r1);
	for (ImageType::IteratorIndexed itP = img_generated_.beginIteratorIndexed(rP); !itP.IsAtEnd(); ++itP, ++it1)
		it1.Set(f2u(itP.Get()));

	it2 = img_out.beginIteratorIndexed(r2);
	for (ImageType::IteratorIndexed itQ = img_orig_.beginIteratorIndexed(rQ); !itQ.IsAtEnd(); ++itQ, ++it2)
		it2.Set(f2u(itQ.Get()));


	std::stringstream ss;
	ss << "/tmp/overlap_"<< k << ".png";

	img_out.save(ss.str());
}
