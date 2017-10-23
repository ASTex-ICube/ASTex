#include "quilting_all.h"
#include <ASTex/colorspace_filters.h>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iomanip>

using namespace ASTex;


QuiltingAll::QuiltingAll():
	tile_width_(64),
	overlap_width_(16),
	blend_rad_(7),
	level_comp_(-1)
{
	pos_first_tile[0]=-1;
	pos_first_tile[1]=-1;
}

int QuiltingAll::random_int(int min, int max)
{
	double d = double(std::rand()) / RAND_MAX;
	return  int(d*(max-min) + min);
}



void QuiltingAll::loadInput(const std::string& imgName)
{
	ImageRGBu8 img;
	img.load(imgName);
	ColorSpace::FilterRGB255To01<ImageRGBu8::ItkImg,QuiltingAll::ImageType::ItkImg>::Pointer filter =
			ColorSpace::FilterRGB255To01<ImageRGBu8::ItkImg,QuiltingAll::ImageType::ItkImg>::New();
	filter->SetInput(img.itk());
	filter->Update();
	img_orig_.itk() = filter->GetOutput();

	std::string basename = imgName.substr(0,imgName.size()-4);
	content_.load(basename);
}

void QuiltingAll::saveOutput(const std::string& imgNameOut)
{
//	{
//		ColorSpace::FilterGray01To255<ImageGrayd::ItkImg,ImageGrayu8::ItkImg>::Pointer filter =
//				ColorSpace::FilterGray01To255<ImageGrayd::ItkImg,ImageGrayu8::ItkImg>::New();
//		filter->SetInput(img_blend_coef_.itk());
//		ImageGrayu8 img(filter->GetOutput());
//		img.save("/tmp/blend_coefs.png");
//	}

	std::size_t slash = imgNameOut.find_last_of('/')+1;
	std::string imgName = "/tmp/quilting"+imgNameOut.substr(slash,imgNameOut.size()-slash-4);

	std::cout << "imgName:"<<imgName << std::endl;
	ColorSpace::FilterRGB01To255<QuiltingAll::ImageType::ItkImg,ImageRGBu8::ItkImg>::Pointer filter =
			ColorSpace::FilterRGB01To255<QuiltingAll::ImageType::ItkImg,ImageRGBu8::ItkImg>::New();
	filter->SetInput(img_generated_.itk());
	ImageRGBu8 img(filter->GetOutput());
	img.save(imgName+".png");


	filter->SetInput(img_generated_b_.itk());
	ImageRGBu8 img2(filter->GetOutput());
	std::string name2 = imgName+"_blend.png";
	img2.save(name2);

	filter->SetInput(img_generated_bl_.itk());
	ImageRGBu8 img3(filter->GetOutput());
	std::string name3 = imgName+"_blend_lab.png";
	img3.save(name3);

	filter->SetInput(img_generated2_.itk());
	ImageRGBu8 img4(filter->GetOutput());
	std::string name4 = imgName+"_cut.png";
	img4.save(name4);

	filter->SetInput(img_dbg_interp.itk());
	ImageRGBu8 img5(filter->GetOutput());
	std::string name5 = imgName+"_interp.png";
	img5.save(name5);
}




void QuiltingAll::init(int tw, int ovl, int rc, int nbrt, int blr)
{
	std::srand(std::time(0));

//	std::srand(0);

	blend_rad_ = blr;

	tile_width_ = tw;
	overlap_width_ = ovl;
	multi_mins_.nbMaxElts_ = rc;

	overlap_h_sz_[0] = tw+blend_rad_;
	overlap_h_sz_[1] = ovl;

	overlap_v_sz_[0] = ovl;
	overlap_v_sz_[1] = tw+blend_rad_;

	tile_sz_[0] = tw;
	tile_sz_[1] = tw;

	bufferPositions_.resize(nbrt);
	bufferErrors_.resize(nbrt); // ????

	label_used_for_choice_ = false;

}


double QuiltingAll::computeErrorOverlap(const Region& or_reg, const Region& gen_reg)
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

void QuiltingAll::computeErrorForCutting(const Region& or_reg, const Region& gen_reg, ErrorCutPathBuffer& e)
{
	ImageType::IteratorIndexed it2 = img_generated_.beginIteratorIndexed(gen_reg);
	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(or_reg); !it.IsAtEnd(); ++it, ++it2)
	{
		double dR = it.Get()[0] - it2.Get()[0];
		double dG = it.Get()[1] - it2.Get()[1];
		double dB = it.Get()[2] - it2.Get()[2];
		double dP2 = dR*dR + dG*dG + dB*dB;
		e( it.GetIndex() - or_reg.GetIndex() ) = dP2;
	}
}


double QuiltingAll::computeErrorOverlapLabels(const Region& or_reg, const Region& gen_reg)
{
	double total=0.0;

	ImageGrayu32::IteratorIndexed it_idx = img_gen_idx_.beginIteratorIndexed(gen_reg);
	ImageType::IteratorIndexed it2 = img_generated_.beginIteratorIndexed(gen_reg);
	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(or_reg); !it.IsAtEnd(); ++it, ++it2, ++it_idx)
	{

		Index p = idx_to_pos(it_idx.Get());
		int l = content_.levels_equality(it.GetIndex(),p);
		double local = 1.0 - double(l)/content_.get_number_of_levels();
		total += local;
	}
	return total/(or_reg.GetSize()[0]*or_reg.GetSize()[1]);
}


void QuiltingAll::computeErrorForCuttingLabels(const Region& or_reg, const Region& gen_reg, ErrorCutPathBuffer& e)
{
	ImageGrayu32::IteratorIndexed it_idx = img_gen_idx_.beginIteratorIndexed(gen_reg);
	ImageType::IteratorIndexed it2 = img_generated_.beginIteratorIndexed(gen_reg);
	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(or_reg); !it.IsAtEnd(); ++it, ++it2)
	{
		Index p = idx_to_pos(it_idx.Get());
		int l = content_.levels_equality(it.GetIndex(),p);
		double err = 1.0 - double(l)/content_.get_number_of_levels();
		e( it.GetIndex() - or_reg.GetIndex() ) = err;
	}
}



double QuiltingAll::computeErrorTile(const Index& or_pos, const Index& gen_pos)
{
	double errV = 0.0;
	// Vertical overlapping (except on left)
	if (gen_pos[0]>=overlap_width_)
	{
		Region r1;
		r1.SetIndex(or_pos);
		r1.SetSize(overlap_v_sz_);

		Region r2;
		r2.SetIndex(gen_pos);
		r2.SetSize(overlap_v_sz_);

		if (label_used_for_choice_)
			errV = computeErrorOverlapLabels(r1,r2);
		else
			errV = computeErrorOverlap(r1,r2);
	}

	double errH = 0.0;
	// horizontal overlapping (except on top)
	if (gen_pos[1]>=overlap_width_)
	{
		Region r3;
		r3.SetIndex(or_pos);
		r3.SetSize(overlap_h_sz_);

		Region r4;
		r4.SetIndex(gen_pos);
		r4.SetSize(overlap_h_sz_);

		if (label_used_for_choice_)
			errH = computeErrorOverlapLabels(r3,r4);
		else
			errH = computeErrorOverlap(r3,r4);
	}

	return std::max(errV,errH);
}


QuiltingAll::Index QuiltingAll::bestFittingTile(const Index& gen_pos, double error_max)
{
	if (gen_pos[0]+gen_pos[1] == 0)
	{
		if (pos_first_tile[0] <0 || pos_first_tile[1] <0 )
		{
			Index best;
			best[0]=random_int(blend_rad_,img_orig_.width()-tile_width_-blend_rad_);
			best[1]=random_int(blend_rad_,img_orig_.height()-tile_width_-blend_rad_);

			std::cout << "FIRST TILE = "<<best << std::endl;
			return best;
		}
		std::cout << "FIRST TILE = "<<pos_first_tile << std::endl;
		return pos_first_tile;
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




void QuiltingAll::verticalPathCut(const Index& or_pos, const Index& gen_pos, std::vector<int>& minPos)
{
	// compute initial error
	ErrorCutPathBuffer e(overlap_width_,tile_width_+blend_rad_);

	Size sz = overlap_v_sz_;
	Region or_reg(or_pos,sz);
	Region gen_reg(gen_pos,sz);

	if (label_used_for_choice_)
		computeErrorForCuttingLabels(or_reg,gen_reg,e);
	else
		computeErrorForCutting(or_reg,gen_reg,e);


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


//	i=tile_width_+ blend_rad_-1;
//	for (int j=0; j< overlap_width_; ++j,--i)
//	{
//		if (minPos[i]< (overlap_width_-1-j))
//			minPos[i] = (overlap_width_-1-j);
//	}
}


void QuiltingAll::horizontalPathCut(const Index& or_pos, const Index& gen_pos, std::vector<int>& minPos)
{
	// compute initial error
	ErrorCutPathBuffer e(tile_width_+ blend_rad_,overlap_width_);

	Size sz = overlap_h_sz_;
	Region or_reg(or_pos,sz);
	Region gen_reg(gen_pos,sz);

	if (label_used_for_choice_)
		computeErrorForCuttingLabels(or_reg,gen_reg,e);
	else
		computeErrorForCutting(or_reg,gen_reg,e);

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

//	i=tile_width_+ blend_rad_-1;
//	for (int j=0; j< overlap_width_; ++j,--i)
//	{
//		if (minPos[i]< (overlap_width_-1-j))
//			minPos[i] = (overlap_width_-1-j);
//	}
}



void QuiltingAll::copyPathCut_BlendDisk(const Region& from_p, const Region& to_p, const std::vector<int>& leftMinCut,const std::vector<int>& upMinCut)
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

	ImageType::IteratorIndexed it3 = img_dbg_interp.beginIteratorIndexed(to);
	ImageGrayu32::IteratorIndexed it_idx_to   = img_gen_idx_.beginIteratorIndexed(to);

	ImageType::IteratorIndexed it_blend   = img_generated_b_.beginIteratorIndexed(to);
	ImageType::IteratorIndexed it_b_lab = img_generated_bl_.beginIteratorIndexed(to);
	ImageType::IteratorIndexed it_gen2 = img_generated2_.beginIteratorIndexed(to);
	ImageType::IteratorIndexed it_gen = img_generated_.beginIteratorIndexed(to);
	for (ImageType::IteratorIndexed it = img_orig_.beginIteratorIndexed(from); !it.IsAtEnd(); ++it, ++it_gen, ++it_idx_to, ++it_blend, ++it_b_lab, ++it_gen2, ++it3)
	{
		const Index& or_pos = it.GetIndex();

		itk::Offset<2> P= it.GetIndex() - origin;

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

		// if necessary blend at end of overlapping zone
		itk::Offset<2> Q= from.GetIndex() + from.GetSize() - it.GetIndex();
		if ((P[1]<overlap_width_+blend_rad_) && (d2>0) && (Q[0]<blend_rad_))
		{
			alpha = std::max(alpha,(double(blend_rad_ - Q[0])/blend_rad_));
		}
		if ((P[0]<overlap_width_+blend_rad_) && (d1>0) && (Q[1]<blend_rad_))
		{
			alpha = std::max(alpha,(double(blend_rad_ - Q[1])/blend_rad_));
		}

		if (alpha < 0.0)
			alpha = 0.0;

		float beta=alpha;
		if (level_comp_ > -8)
		{
			if (is_an_index(it_idx_to.Get()))
			{
				Index p_l = idx_to_pos(it_idx_to.Get());
				if (!content_.has_same_label(it.GetIndex(),p_l,level_comp_))
					beta = 0.0;
			}
		}

		if ( (P[0] >= leftMinCut[P[1]]) && (P[1] >= upMinCut[P[0]]))
		{	// we are inside cut path
			it_blend.Set(it_blend.Get()*alpha + it.Get()*(1.0-alpha));
			it_b_lab.Set(it_b_lab.Get()*beta + it.Get()*(1.0-beta));

			it_gen.Set(it.Get());
			it_gen2.Set(it.Get());
			// set the idx
			it_idx_to.Set(pos_to_idx(or_pos));
		}
		else
		{	// we are outside
			it_blend.Set(it.Get()*alpha + it_blend.Get()*(1.0-alpha));
			it_b_lab.Set(it.Get()*beta + it_b_lab.Get()*(1.0-beta));
		}


		// debug image output
		itk::RGBPixel<double> pix;
		double gamma = beta + 0.5;
		if (gamma > 1.0)
			gamma = 1.0;

		if ((d1==0 && d2 >=0) || (d2==0 && d1>=0))
			it_gen2.Set(RGBd(1,0,0));


		if (d1<0 && d2<0)
		{
			pix[0]=1.0; pix[1]=0.0; pix[2]=0.0;
			it3.Set(pix*gamma);
		}
		else if (d1<0 && d2>0)
		{
			pix[0]=0.0; pix[1]=1.0; pix[2]=0.0;
			it3.Set(pix*gamma);
		}
		else if (d1>0 && d2<0)
		{
			pix[0]=0.0; pix[1]=0.0; pix[2]=1.0;
			it3.Set(pix*gamma);
		}
		else if (d1>0 && d2>0)
		{
			pix[0]=1.0; pix[1]=0.0; pix[2]=1.0;
			it3.Set(pix*gamma);
		}
		else if (d1==0 || d2==0)
			{
				pix[0]=1.0; pix[1]=1.0; pix[2]=1.0;
				it3.Set(pix*gamma);
			}
	}
}




void QuiltingAll::computeFittingTilesPathCutInterp(int width, int height, int tw, int ovl, int bl_rad ,int rc, int nbrt)
{

	img_dbg_interp.initItk(width,height);

	img_gen_idx_.initItk(width,height);
	img_generated_.initItk(width,height);
	img_generated2_.initItk(width,height);
	img_generated_b_.initItk(width,height);
	img_generated_bl_.initItk(width,height);


//	img_blend_coef_.initItk(width,height);


	init(tw,ovl,rc,nbrt,bl_rad);

	// fill idx with 0xffffffff
	for (ImageGrayu32::IteratorIndexed it = img_gen_idx_.beginIteratorIndexed(); !it.IsAtEnd(); ++it)
		it.Set(0xffffffff);

	Region from;
	from.SetSize(tile_sz_);
	Region to;
	to.SetSize(tile_sz_);

	std::vector<int> minPosH;
	minPosH.resize(tile_width_+blend_rad_);
	std::vector<int> minPosV;
	minPosV.resize(tile_width_+blend_rad_);

	int incr = tw - ovl;

	int k=0;
	for (int i = 0; i<=height-tw-blend_rad_; i+= incr)
	{
		int kk=0;

		for (int j = 0; j<=width-tw-blend_rad_; j+= incr)
		{
//			if ((k==2)&&(kk==1))
//				return;

			to.SetIndex(0,j);
			to.SetIndex(1,i);
			Index best = bestFittingTile(to.GetIndex());

//			saveDebugOverlap(to.GetIndex(), best, k++);

			from.SetIndex(best);

			if (j==0)
				minPosV.assign(minPosV.size(),-1);
			else
				verticalPathCut(from.GetIndex(), to.GetIndex(), minPosV);

			if (i==0)
				minPosH.assign(minPosH.size(),-1);
			else
				horizontalPathCut(from.GetIndex(), to.GetIndex(), minPosH);

			copyPathCut_BlendDisk(from,to,minPosV,minPosH);
			kk++;
		}

		k++;
	}
}



void QuiltingAll::saveDebugOverlap(const Index& P, const Index& Q, int k)
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
	if (label_used_for_choice_)
		ss << "/tmp/overlap_lab"<< k << ".png";
	else
		ss << "/tmp/overlap_"<< k << ".png";

	img_out.save(ss.str());
}

