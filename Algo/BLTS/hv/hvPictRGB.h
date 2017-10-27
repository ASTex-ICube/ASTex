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



// hvPictRGB.h: interface for the picture class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PICTRGB_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_)
#define AFX_PICTRGB_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>

#include <string.h>
#include "hvLinearTransform3.h"
#include "hvBitmap.h"
#include "hvField2.h"
#include "hvColor.h"
#include "hvPict.h"


namespace hview {

	// from Ward's RGBE images
const int MAX_SYNTHESIS_PIXEL_LIST = 5;
typedef struct {
  int valid;            /* indicate which fields are valid */
  char programtype[16]; /* listed at beginning of file to identify it 
                         * after "#?".  defaults to "RGBE" */ 
  float gamma;          /* image has already been gamma corrected with 
                         * given gamma.  defaults to 1.0 (no correction) */
  float exposure;       /* a value of 1.0 in an image corresponds to
			 * <exposure> watts/steradian/m^2. 
			 * defaults to 1.0 */
} rgbe_header_info;

/* flags indicating which fields in an rgbe_header_info are valid */
#define RGBE_VALID_PROGRAMTYPE 0x01
#define RGBE_VALID_GAMMA       0x02
#define RGBE_VALID_EXPOSURE    0x04

template <class T> class hvPictRGBA;
template <class R, unsigned int n> class hvQPictRGB;

////////////////////////////////////////////////////////////
template <class T> class hvPictRGB : public hvField2< hvColRGB<T> >  
////////////////////////////////////////////////////////////
{
public:
	hvPictRGB<T>() : hvField2< hvColRGB<T> >(),hvArray2< hvColRGB<T> >() { }
	hvPictRGB<T>(int sx, int sy, const hvColRGB<T> &nil) : hvField2< hvColRGB<T> >(sx, sy, nil),hvArray2< hvColRGB<T> >(sx, sy, nil) { }
	void reset(int sx, int sy,const hvColRGB<T> &nil)
	{
		hvField2< hvColRGB<T> >::reset(sx,sy, nil);
	}
	void reset()
	{
		hvField2< hvColRGB<T> >::reset();
	}
	void clone(const hvPictRGB<T> &pict,int x, int y, int sx, int sy)
	{
		hvField2< hvColRGB<T> >::reset(sx-x+1, sy-y+1, hvColRGB<T>(0));
		int i,j;
		for (i=x; i<=sx; i++) for (j=y; j<=sy; j++)
		{
			this->update(i-x,j-y,pict.get(i,j));
		}
	}
	template <class U> void clone(const hvPictRGB<U> &pict, U scalu, U gamma, U scalt, int x, int y, int sx, int sy)
	{
		hvField2< hvColRGB<T> >::reset(sx-x+1, sy-y+1, hvColRGB<T>(0));
		int i,j;
		for (i=x; i<= sx; i++) for (j=y; j<= sy; j++)
		{
			hvColRGB<U> col = pict.get(i,j);
			col.scale(1.0/scalu);
			col.gamma(1.0, gamma);
			this->update(i-x,j-y,hvColRGB<T>( (T)(scalt*col.RED()),(T)(scalt*col.GREEN()),(T)(scalt*col.BLUE()) ));
		}
	}
	hvPictRGB<T>(const hvBitmap &pict, const hvColRGB<T> &va, const hvColRGB<T> &vb): hvField2< hvColRGB<T> >(pict.sizeX(), pict.sizeY(), hvColRGB<T>(0)), hvArray2< hvColRGB<T> >(pict.sizeX(), pict.sizeY(), hvColRGB<T>(0))
	{
		int i,j;
		for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++)
		{
			if (pict.get(i,j)) this->update(i,j,va); else this->update(i,j,vb);
		}
	}
	void convert(const hvBitmap &pict, const hvColRGB<T> &va, const hvColRGB<T> &vb)
	{
		hvField2< hvColRGB<T> >::reset(pict.sizeX(), pict.sizeY(), hvColRGB<T>(0));
		//hvArray2< hvColRGB<T> >::reset(pict.sizeX(), pict.sizeY(), hvColRGB<T>(0));
		int i,j;
		for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++)
		{
			if (pict.get(i,j)) this->update(i,j,va); else this->update(i,j,vb);
		}
	}
	hvPictRGB<T>(const hvPictRGBA<T> &pict, hvPictComponent cc): hvField2< hvColRGB<T> >(pict.sizeX(), pict.sizeY(), hvColRGB<T>(0)), hvArray2< hvColRGB<T> >(pict.sizeX(), pict.sizeY(), hvColRGB<T>(0))
	{
		int i,j;
		for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++)
		{
			hvColRGB<T> v;
			hvColRGBA<T> co = pict.get(i,j);
			switch(cc) 
			{
				case HV_RED: v=hvColRGB<T>(co.RED()); break;
				case HV_GREEN: v=hvColRGB<T>(co.GREEN()); break;
				case HV_BLUE: v=hvColRGB<T>(co.BLUE()); break;
				case HV_ALPHA: v=hvColRGB<T>(co.ALPHA()); break;
				case HV_RGB: v=hvColRGB<T>(co); break;
				default: v = hvColRGB<T>(co.luminance());
			}
			this->update(i,j,v); 
		}
	}
	hvPictRGB<T>(const hvPictRGB<T> &pict, hvPictComponent cc): hvField2< hvColRGB<T> >(pict.sizeX(), pict.sizeY(), hvColRGB<T>(0)), hvArray2< hvColRGB<T> >(pict.sizeX(), pict.sizeY(), hvColRGB<T>(0))
	{
		int i,j;
		for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++)
		{
			hvColRGB<T> v;
			hvColRGB<T> co = pict.get(i,j);
			switch(cc) 
			{
				case HV_RED: v=hvColRGB<T>(co.RED()); break;
				case HV_GREEN: v=hvColRGB<T>(co.GREEN()); break;
				case HV_BLUE: v=hvColRGB<T>(co.BLUE()); break;
				case HV_RGB: v=hvColRGB<T>(co); break;
				default: v = hvColRGB<T>(co.luminance());
			}
			this->update(i,j,v); 
		}
	}
	template <class U, class V> hvPictRGB(const hvPict<U> &p, V scal);
	template <class U, class V> hvPictRGB(const hvPict<U> &p, V scal, V shift);
	template <class U, class V> hvPictRGB(const hvPict<U> &p, V scal, int x, int y, int sx, int sy);
	template <class U, class V> void convert(const hvPict<U> &p, V scal, int x, int y, int sx, int sy);
	template <class U, class V> void convertloga(const hvPict<U> &p, V loga, V max, V scal, int x, int y, int sx, int sy);
	template <class U, class V> hvPictRGB(const hvPict<U> &p, V scal, V min, V max);
	template <class U> hvPictRGB(const hvPict<U> &p, const std::vector<hvColRGB<unsigned char> > &lcol);


	// choose RGB components from two input
	void merge(const hvPictRGB<T> &pa, const hvPictRGB<T> &pb, bool rr, bool gg, bool bb)
	{
		int i,j;
		hvField2< hvColRGB<T> >::reset(pa.sizeX(), pa.sizeY(), hvColRGB<T>(0));
		for (i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> ca = pa.get(i,j);
			hvColRGB<T> cb = pb.get(i,j);
			hvColRGB<T> v=hvColRGB<T>(rr?ca.RED():cb.RED(), gg?ca.GREEN():cb.GREEN(), bb?ca.BLUE():cb.BLUE()); 
			this->update(i,j,v);
		}
	}

	// keep only some RGB components
	void extract(bool rr, bool gg, bool bb, T cr, T cg, T cb)
	{
		int i,j;
		for (i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> co = this->get(i,j);
			hvColRGB<T> v=hvColRGB<T>(rr?co.RED():cr,gg?co.GREEN():cg,bb?co.BLUE():cb); 
			this->update(i,j,v);
		}
	}
	// Convert into HSV color space
	void toHSV(const hvPictRGB<T> &p, T scal, hvColRGB<double> weights=hvColRGB<double>(1.0,1.0,1.0))
	{
		int i,j;
		hvField2< hvColRGB<T> >::reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		for (i=0; i<p.sizeX(); i++) for (j=0; j<p.sizeY(); j++)
		{
			hvColRGB<T> val = p.get(i,j);
			hvColRGB<T> valhsv; valhsv.tohsv(val, scal);
			this->update(i,j,hvColRGB<T>(T((double)valhsv.RED()*weights.RED()), T((double)valhsv.GREEN()*weights.GREEN()), T((double)valhsv.BLUE()*weights.BLUE())));
		}
	}
	// Convert into XYZ color space
	void toXYZ(const hvPictRGB<T> &p, T scal, hvColRGB<double> weights=hvColRGB<double>(1.0,1.0,1.0))
	{
		int i,j;
		hvField2< hvColRGB<T> >::reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		for (i=0; i<p.sizeX(); i++) for (j=0; j<p.sizeY(); j++)
		{
			hvColRGB<T> val = p.get(i,j);
			hvColRGB<T> valhsv; valhsv.toxyz(val, scal);
			this->update(i,j,hvColRGB<T>(T((double)valhsv.RED()*weights.RED()), T((double)valhsv.GREEN()*weights.GREEN()), T((double)valhsv.BLUE()*weights.BLUE())));
		}
	}
	// Convert into LUV color space
	void toLUV(const hvPictRGB<T> &p, T scal, hvColRGB<double> weights=hvColRGB<double>(1.0,1.0,1.0))
	{
		int i,j;
		hvField2< hvColRGB<T> >::reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		for (i=0; i<p.sizeX(); i++) for (j=0; j<p.sizeY(); j++)
		{
			hvColRGB<T> val = p.get(i,j);
			hvColRGB<T> valhsv; valhsv.toLuv(val, scal);
			this->update(i,j,hvColRGB<T>(T((double)valhsv.RED()*weights.RED()), T((double)valhsv.GREEN()*weights.GREEN()), T((double)valhsv.BLUE()*weights.BLUE())));
		}
	}
	// Convert from XYZ color space back to RGB
	void fromXYZ(const hvPictRGB<T> &p, T scal, hvColRGB<double> weights=hvColRGB<double>(1.0,1.0,1.0))
	{
		int i,j;
		hvField2< hvColRGB<T> >::reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		for (i=0; i<p.sizeX(); i++) for (j=0; j<p.sizeY(); j++)
		{
			hvColRGB<T> val = p.get(i,j);
			hvColRGB<T> valhsv; valhsv.fromxyz(val, scal);
			this->update(i,j,hvColRGB<T>(T((double)valhsv.RED()*weights.RED()), T((double)valhsv.GREEN()*weights.GREEN()), T((double)valhsv.BLUE()*weights.BLUE())));
		}
	}
	// Convert from LUV color space back to RGB
	void fromLUV(const hvPictRGB<T> &p, T scal, hvColRGB<double> weights=hvColRGB<double>(1.0,1.0,1.0))
	{
		int i,j;
		hvField2< hvColRGB<T> >::reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		for (i=0; i<p.sizeX(); i++) for (j=0; j<p.sizeY(); j++)
		{
			hvColRGB<T> val = p.get(i,j);
			hvColRGB<T> valhsv; valhsv.fromLuv(val, scal);
			this->update(i,j,hvColRGB<T>(T((double)valhsv.RED()*weights.RED()), T((double)valhsv.GREEN()*weights.GREEN()), T((double)valhsv.BLUE()*weights.BLUE())));
		}
	}
	// difference between two images of same resolution
	void difference(const hvPictRGB<T> &pia,const hvPictRGB<T> &pib)
	{
		hvField2< hvColRGB<T> >::reset(pia.sizeX(), pia.sizeY(), hvColRGB<T>(0));
		int i,j;
		for (i=0; i<pia.sizeX(); i++) for (j=0; j<pia.sizeY(); j++)
		{
			hvColRGB<T> cola = pia.get(i,j);
			hvColRGB<T> colb = pib.get(i,j);
			hvColRGB<T> cold; cold.difference(cola,colb);
			this->update(i,j,cold);
		}
	}
	// difference between two images of same resolution
	void difference(const hvPictRGB<T> &pia,const hvPictRGB<T> &pib, double scale, double offset)
	{
		hvField2< hvColRGB<T> >::reset(pia.sizeX(), pia.sizeY(), hvColRGB<T>(0));
		int i,j;
		for (i=0; i<pia.sizeX(); i++) for (j=0; j<pia.sizeY(); j++)
		{
			hvColRGB<T> cola = pia.get(i,j);
			hvColRGB<T> colb = pib.get(i,j);
			hvColRGB<T> cold; cold.difference(cola,colb, scale, offset);
			this->update(i,j,cold);
		}
	}
	// difference between two images of same resolution
	void differenceMask(double sqdthresh, const hvPictRGB<T> &pi, hvBitmap &diff)
	{
		diff.reset(pi.sizeX(), pi.sizeY(), true);
		int i,j;
		for (i=0; i<pi.sizeX(); i++) for (j=0; j<pi.sizeY(); j++)
		{
			hvColRGB<T> cola = pi.get(i,j);
			hvColRGB<T> colb = this->get(i,j);
			diff.set(i,j,(cola.squaredDifference(colb)>sqdthresh));
		}
	}

	void squaredDifference(int px, int py, int dx, int dy, const hvPictRGB<unsigned char> &pia, int ix, int iy, const hvPictRGB<unsigned char> &pib);

	double meansquareDifference(const hvPictRGB<T> &pi, hvBitmap &mask)
	{
		int i,j;
		int count=0;
		double sum=0.0;
		for (i=0; i<pi.sizeX(); i++) for (j=0; j<pi.sizeY(); j++)
		{
			if (mask.get(i,j))
			{
			hvColRGB<T> cola = pi.get(i,j);
			hvColRGB<T> colb = this->get(i,j);
			sum+=cola.squaredDifference(colb);
			count++;
			}
		}
		return sum/(double)count;
	}

	double meansquareDifference(int u, int v, const hvPictRGB<T> &pi, int x, int y, int sx, int sy)
	{
		int i,j;
		int count=0;
		double sum=0.0;
		for (i=-sx; i<=sx; i++) for (j=-sy; j<=sy; j++)
		{
			int xx = x+i; if (xx<0) xx+= pi.sizeX(); else if (xx>=pi.sizeX()) xx-=pi.sizeX();
			int yy = y+j; if (yy<0) yy+= pi.sizeY(); else if (yy>=pi.sizeY()) yy-=pi.sizeY();
			int uu = u+i; if (uu<0) uu+= this->sizeX(); else if (uu>=this->sizeX()) uu-=this->sizeX();
			int vv = v+j; if (vv<0) vv+= this->sizeY(); else if (vv>=this->sizeY()) vv-=this->sizeY();
			hvColRGB<T> cola = pi.get(xx,yy);
			hvColRGB<T> colb = this->get(uu,vv);
			sum+=cola.squaredDifference(colb)/3.0;
			count++;
		}
		return sum/(double)count;
	}

	// applies a blurring mask of size sx,sy
	// scal is the normalization factor, usually sx*sy
	void blur(T scal, int sx, int sy)
	{
		hvField2<hvColRGB<T> > source; source.clone(this);
		hvField2<hvColRGB<T> >:: template blur<hvColRGB<double> >(&source, sx, sy, hvColRGB<double>(1.0/(double)scal));
	}
	void gaussianBlur(int sx, int sy)
	{
		hvField2<hvColRGB<T> > source; source.clone(this);
		hvField2<hvColRGB<T> >:: template gaussianBlur<hvColRGB<double> >(&source, sx, sy);
	}
	// applies a deblurring mask of size sx,sy
	void deblur(int niter, int sx, int sy, double scal, double min, double max)
	{
		hvField2<hvColRGB<T> > source; source.clone(this);
		hvField2<hvColRGB<T> >:: template deblur<hvColRGB<double> >(&source, niter, sx, sy, scal, min, max);
	}
	// Difference of Gaussians
	void DoG(const hvPictRGB<T> &pia, int sx, int sy, int nrec=-1)
	{
		hvPictRGB<T> pi, pib;
		pi.clone(pia, 0, 0, pia.sizeX()-1, pia.sizeY()-1);
		if (nrec>=0) { pib.shrink(&pi,nrec); pi.clone(pib,0,0,pib.sizeX()-1,pib.sizeY()-1); }
		else pib.clone(pia, 0, 0, pia.sizeX()-1, pia.sizeY()-1);
		pib.gaussianBlur(sx,sy);
		this->difference(pi,pib);
	}
	// Difference of Gaussians
	void DoG(const hvPictRGB<T> &pia, int sx, int sy, double scale, double offset)
	{
		hvPictRGB<T> pib;
		pib.clone(pia, 0, 0, pia.sizeX()-1, pia.sizeY()-1);
		pib.gaussianBlur(sx,sy);
		this->difference(pia,pib, scale, offset);
	}
	// Difference of Gaussians
	void DoG(hvPict<T> &pires, int sx, int sy)
	{
		int i,j;
		hvPictRGB<T> pi;
		pi.clone(*this, 0, 0, this->sizeX()-1, this->sizeY()-1);
		pi.gaussianBlur(sx,sy);
		pires.reset(this->sizeX(), this->sizeY(),0);
		for (i=0; i<pi.sizeX(); i++) for (j=0; j<pi.sizeY(); j++)
		{
			pires.update(i,j,(T)(sqrt(pi.get(i,j).squaredDifference(this->get(i,j))/3.0)));
		}
	}

	// bilateral filter
	void bilateral(const hvPictRGB<T> &pia, const hvPictRGB<T> &pib, double sigma, T scal, int sx, int sy)
	{
		hvField2< hvColRGB<T> >::reset(pia.sizeX(), pia.sizeY(), hvColRGB<T>(0));
		int i,j, ii, jj;
		for (i=0; i<pia.sizeX(); i++) for (j=0; j<pia.sizeY(); j++)
		{
			hvColRGB<T> cola = pia.get(i,j);
			double ra=(double)cola.RED()/(double)scal;
			double rg=(double)cola.GREEN()/(double)scal;
			double rb=(double)cola.BLUE()/(double)scal;
			double ww=0.0;
			double resa=0.0, resg=0.0, resb=0.0;
			for (ii=-sx/2; ii<=sx/2; ii++) for (jj=-sy/2; jj<=sy/2; jj++)
			{
				if (i+ii>=0 && i+ii<pia.sizeX() && j+jj>=0 && j+jj<pia.sizeY())
				{
					cola = pia.get(i+ii,j+jj);
					double rra=(double)cola.RED()/(double)scal;
					double rrg=(double)cola.GREEN()/(double)scal;
					double rrb=(double)cola.BLUE()/(double)scal;
					cola = pib.get(i+ii,j+jj);
					double va=(double)cola.RED()/(double)scal;
					double vg=(double)cola.GREEN()/(double)scal;
					double vb=(double)cola.BLUE()/(double)scal;
					double dist = (va-ra)*(va-ra)+(vg-rg)*(vg-rg)+(vb-rb)*(vb-rb);
					dist = exp(-dist*sigma);
					ww+= dist;
					resa+=dist*rra; resg+=dist*rrg; resb+=dist*rrb;
				}
			}
			resa/=ww; resg/=ww; resb/=ww;
			this->update(i,j,hvColRGB<T>(T(resa*(double)scal), T(resg*(double)scal), T(resb*(double)scal)));
		}
	}
	// bilateral filter
	template <class U> void bilateral(const hvPictRGB<T> &pia, const hvPict<U> &pib, double sigma, T scal, int sx, int sy)
	{
		int i,j, k, ii, jj;
		/*
		int nb = pib.maxValue()+1;
		if (nb>200) hvFatal("too many classes for bilateral filtering");
		int nx = 2*pia.sizeX()/sx+1, ny=2*pia.sizeY()/sy+1;
		int tx = pia.sizeX()/nx+1, ty = pia.sizeY()/ny+1;
		printf("bilateral: filt size %d,%d blocs %d,%d (size: %d,%d)\n", sx,sy,nx,ny,tx,ty);
		hvArray2<hvList<hvVec2<int> > *> *blocs[200];
		for (k=0; k<nb; k++) blocs[k] = new hvArray2<hvList<hvVec2<int> > *>(nx,ny,0);
		for (k=0; k<nb; k++) for (i=0; i<nx; i++) for (j=0; j<ny; j++)
		{
			blocs[k]->update(i,j, new hvList<hvVec2<int> >(tx*ty/2));
		}
		for (i=0; i<pia.sizeX(); i++) for (j=0; j<pia.sizeY(); j++)
		{
			int q = pib.get(i,j);
			if (q<0 || q>=nb) hvFatal("index out of range in bilateral filter");
			hvList<hvVec2<int> > *ll = blocs[q]->get(i/tx,j/ty);
			if (ll->isFull()) { ll->resize(ll->getMaxElements()+100); blocs[q]->update(i/tx,j/ty,ll); }
			ll->pushBack(hvVec2<int>(i,j));
		}
		hvField2< hvColRGB<T> >::reset(pia.sizeX(), pia.sizeY(), hvColRGB<T>(0));
		for (i=0; i<pia.sizeX(); i++) for (j=0; j<pia.sizeY(); j++)
		{
			hvColRGB<T> cola = pia.get(i,j);
			double ra=(double)cola.RED()/(double)scal;
			double rg=(double)cola.GREEN()/(double)scal;
			double rb=(double)cola.BLUE()/(double)scal;
			double ww=0.0;
			double resa=0.0, resg=0.0, resb=0.0;
			int ida = pib.get(i,j);

			for (ii=-1; ii<=1; ii++) for (jj=-1; jj<=1; jj++)
			{
				int indx = i/tx+ii, indy = j/ty+jj;
				if (indx>=0 && indx<nx && indy>=0 && indy<ny)
				{
					hvList<hvVec2<int> > *ll = blocs[ida]->get(indx,indy);
					for (k=0; k<ll->length(); k++)
					{
						hvVec2<int> pix = ll->get(k);
						if (pix.X()>=i-sx/2 && pix.X()<=i+sx/2 && pix.Y()>=j-sy/2 && pix.Y()<=j+sy/2)
						{
							cola = pia.get(pix.X(),pix.Y());
							double rra=(double)cola.RED()/(double)scal;
							double rrg=(double)cola.GREEN()/(double)scal;
							double rrb=(double)cola.BLUE()/(double)scal;
							double q=exp(-2.0*(double)((pix.X()-i)*(pix.X()-i)+(pix.Y()-j)*(pix.Y()-j))/(double)(sx*sx+sy*sy));
							double dist = (rra-ra)*(rra-ra)+(rrg-rg)*(rrg-rg)+(rrb-rb)*(rrb-rb);
							dist = exp(-dist*sigma);
							ww+= dist*q;
							resa+=dist*rra*q; resg+=dist*rrg*q; resb+=dist*rrb*q;
						}
					}
				}
			}
			resa/=ww; resg/=ww; resb/=ww;
			update(i,j,hvColRGB<T>(T(resa*(double)scal), T(resg*(double)scal), T(resb*(double)scal)));
		}
		 for (k=0; k<nb; k++) for (i=0; i<nx; i++) for (j=0; j<ny; j++) delete blocs[k]->get(i,j);
		 for (k=0; k<nb; k++) delete blocs[k];
		*/
		hvField2< hvColRGB<T> >::reset(pia.sizeX(), pia.sizeY(), hvColRGB<T>(0));
		for (i=0; i<pia.sizeX(); i++) for (j=0; j<pia.sizeY(); j++)
		{
			hvColRGB<T> cola = pia.get(i,j);
			double ra=(double)cola.RED()/(double)scal;
			double rg=(double)cola.GREEN()/(double)scal;
			double rb=(double)cola.BLUE()/(double)scal;
			double ww=0.0;
			double resa=0.0, resg=0.0, resb=0.0;
			int ida = (int)pib.get(i,j);
			for (ii=-sx/2; ii<=sx/2; ii++) for (jj=-sy/2; jj<=sy/2; jj++)
			{
				if (i+ii>=0 && i+ii<pia.sizeX() && j+jj>=0 && j+jj<pia.sizeY())
				{
					cola = pia.get(i+ii,j+jj);
					double rra=(double)cola.RED()/(double)scal;
					double rrg=(double)cola.GREEN()/(double)scal;
					double rrb=(double)cola.BLUE()/(double)scal;
					if ((int)pib.get(i+ii,j+jj)==ida) 
					{ 
						double q=exp(-2.0*(double)(ii*ii+jj*jj)/(double)(sx*sx+sy*sy));
						double dist = (rra-ra)*(rra-ra)+(rrg-rg)*(rrg-rg)+(rrb-rb)*(rrb-rb);
						dist = exp(-dist*sigma);
						ww+= dist*q;
						resa+=dist*rra*q; resg+=dist*rrg*q; resb+=dist*rrb*q;
						//resa+=rra*q; resg+=rrg*q; resb+=rrb*q; ww+=q; 
					}
				}
			}	
			resa/=ww; resg/=ww; resb/=ww;
			this->update(i,j,hvColRGB<T>(T(resa*(double)scal), T(resg*(double)scal), T(resb*(double)scal)));
		}
	}

	template <class U> void bilateralApprox(const hvPictRGB<T> &pia, const hvPict<U> &pib, T scal, int sx, int sy)
	{
		int nb = pib.maxValue()+1;
		int i,j, k, ii, jj;
		int TT = 16;
		while (sx<TT*8 || sy<TT*8) TT/=2;
		//printf("bilateral: classes: %d, filt size %d,%d blocs %d,%d\n", nb, sx,sy,TT,TT);
		if (nb>200 || TT<=2) { this->bilateral(pia,pib,0.0,scal,sx,sy); return; }
		int tx = pia.sizeX()/TT+1, ty=pia.sizeY()/TT+1;
		//if (nb>200) hvFatal("too many classes for bilateral filtering");
		if (nb*tx*ty>20*1024*1024) { this->bilateral(pia,pib,0.0,scal,sx,sy); return; }
		hvPict<int> *count[200];
		hvPictRGB<double> *avgcol[200];
		for (k=0; k<nb; k++)
		{
			count[k] = new hvPict<int>(tx,ty,0);
			avgcol[k]= new hvPictRGB<double>(tx,ty,hvColRGB<double>(0.0));
		}
		for (i=0; i<pia.sizeX(); i+=TT) for (j=0; j<pia.sizeY(); j+=TT)
		{
			hvColRGB<double> col[200];
			int nn[200];
			for (k=0; k<nb; k++) { col[k]=hvColRGB<double>(0.0); nn[k]=0; }
			for (ii=0; ii<TT; ii++) for (jj=0; jj<TT;jj++)
			{
				if (i+ii<pia.sizeX() && j+jj<pia.sizeY())
				{
					int inda = (int)pib.get(i+ii,j+jj);
					nn[inda]++;
					hvColRGB<T> cola = pia.get(i+ii,j+jj);
					double ra=(double)cola.RED()/(double)scal;
					double rg=(double)cola.GREEN()/(double)scal;
					double rb=(double)cola.BLUE()/(double)scal;
					col[inda] = hvColRGB<double>(ra+col[inda].RED(),rg+col[inda].GREEN(),rb+col[inda].BLUE()); 
				}
			}
			for (k=0; k<nb; k++)
			{
				count[k]->update(i/TT,j/TT,nn[k]);
				avgcol[k]->update(i/TT,j/TT,col[k]);
			}
		}
		hvField2< hvColRGB<T> >::reset(pia.sizeX(), pia.sizeY(), hvColRGB<T>(0));
		for (i=0; i<pia.sizeX(); i++) for (j=0; j<pia.sizeY(); j++)
		{
			hvColRGB<T> cola = pia.get(i,j);
			double ra=(double)cola.RED()/(double)scal;
			double rg=(double)cola.GREEN()/(double)scal;
			double rb=(double)cola.BLUE()/(double)scal;
			double ww=0.0;
			double resa=0.0, resg=0.0, resb=0.0;
			int ida = pib.get(i,j);
			for (ii=-sx/2/TT; ii<=sx/2/TT; ii++) for (jj=-sy/2/TT; jj<=sy/2/TT; jj++)
			{
				if (i/TT+ii>=0 && i/TT+ii<tx && j/TT+jj>=0 && j/TT+jj<ty)
				{
					hvColRGB<double> cola = avgcol[ida]->get(i/TT+ii,j/TT+jj);
					double q=exp(-2.0*(double)(ii*ii*TT*TT+jj*jj*TT*TT)/(double)(sx*sx+sy*sy));
					ww+= (q*(double)count[ida]->get(i/TT+ii,j/TT+jj));
					resa+=cola.RED()*q; resg+=cola.GREEN()*q; resb+=cola.BLUE()*q;
				}
			}	
			resa/=ww; resg/=ww; resb/=ww;
			this->update(i,j,hvColRGB<T>(T(resa*(double)scal), T(resg*(double)scal), T(resb*(double)scal)));
		}
		for (k=0; k<nb; k++)
		{
			delete count[k] ;
			delete avgcol[k] ;
		}
	}

	void histogramm(double *histo, int NN, double norm, int x, int y, int fx, int fy, hvPictComponent comp)
	{
		int i,j;
		for (i=0; i<NN; i++) histo[i]=0.0;
		for (i=x; i<=fx; i++) for (j=y; j<=fy; j++)
		{
			hvColRGB<T> col = this->get(i,j);
			double v;
			switch(comp)
			{
				case HV_RED: v=col.RED(); break;
				case HV_GREEN: v=col.GREEN(); break;
				case HV_BLUE: v=col.BLUE(); break;
				default: v = col.luminance();
			}
			v /= norm;
			int bin = v*(double)NN; if (bin<0) bin=0; else if (bin>=NN-1) bin=NN-1;
			histo[bin]+=1.0;
		}
		for (i=0; i<NN; i++) histo[i]/=(double)((fx-x+1)*(fy-y+1));
	}

	// checks if R,G,B values are local max in mask size [-1,1] , returns value or 0 if not max
	hvColRGB<T> isLocalMaximum(const hvPictRGB<T> &dog2, int x, int y) const
	{
		hvColRGB<T> vxy = this->get(x,y);
		int i,j;
		bool rr=true, gg=true, bb=true;
		for (i=-1; i<=1; i++) for (j=-1; j<=1; j++) 
		{
			if (x+i>=0 && x+i<this->sizeX() && y+j>=0 && y+j<this->sizeY())
			{
				hvColRGB<T> v = this->get(x+i,y+j);
				if (vxy.RED()<v.RED()) rr=false;
				if (vxy.GREEN()<v.GREEN()) gg=false;
				if (vxy.BLUE()<v.BLUE()) bb=false;
				v = dog2.get((x+i)/2,(y+j)/2);
				if (vxy.RED()<v.RED()) rr=false;
				if (vxy.GREEN()<v.GREEN()) gg=false;
				if (vxy.BLUE()<v.BLUE()) bb=false;
			}
		}
		return hvColRGB<T>(rr?vxy.RED():T(0),gg?vxy.GREEN():T(0),bb?vxy.BLUE():T(0));
	}

	// computes local gradient using central differences along X and Y axes
	void gradient(int x, int y, T scale, hvColRGB<double> *dx, hvColRGB<double> *dy) const
	{
		if (x-1>=0 && x+1<this->sizeX())
		{
			hvColRGB<T> cola= this->get(x-1,y);
			hvColRGB<T> colb= this->get(x+1,y);
			*dx=hvColRGB<double>((double)colb.X()/(double)scale-(double)cola.X()/(double)scale,
				(double)colb.Y()/(double)scale-(double)cola.Y()/(double)scale,
				(double)colb.Z()/(double)scale-(double)cola.Z()/(double)scale);
		}
		if (y-1>=0 && y+1<this->sizeY())
		{
			hvColRGB<T> cola= this->get(x,y-1);
			hvColRGB<T> colb= this->get(x,y+1);
			*dy=hvColRGB<double>((double)colb.X()/(double)scale-(double)cola.X()/(double)scale,
				(double)colb.Y()/(double)scale-(double)cola.Y()/(double)scale,
				(double)colb.Z()/(double)scale-(double)cola.Z()/(double)scale);
		}
	}

	void mySIFTDescriptor(float angle, int x, int y, T scale, hvColRGB<double> descrx[8], hvColRGB<double> descry[8])
	{
		int i;
		int dx[8]={ 1,1,0,-1,-1,-1,0,1};
		int dy[8]={ 0,1,1,1,0,-1,-1,-1};
		int decal = (int)(angle/M_PI*4.0);
		//hvColRGB<T> colxy = get(x,y);
		for (i=0; i<8; i++)
		{
				hvColRGB<T> col(T(0));
				int ddx = dx[(i+decal)%8];
				int ddy = dy[(i+decal)%8];
				if (x+ddx>=0 && x+ddx<this->sizeX() && y+ddy>=0 && y+ddy<this->sizeY())
				{
					this->gradient(x+ddx,y+ddy, scale, descrx+i, descry+i);
				}
				else { descrx[i]=hvColRGB<double>(0.0); descry[i]=hvColRGB<double>(0.0); }
		}
	}
	void allmySIFTDescriptors(float angle, T scale, hvPictRGB<double> imgdx[8], hvPictRGB<double> imgdy[8])
	{
		int i,j,k;
		hvColRGB<double> descrx[8], descry[8];
		for (k=0; k<8; k++) imgdx[k].reset(this->sizeX(), this->sizeY(), hvColRGB<double>(0.0));
		for (k=0; k<8; k++) imgdy[k].reset(this->sizeX(), this->sizeY(), hvColRGB<double>(0.0));
		for (i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			this->mySIFTDescriptor(angle, i,j,scale,descrx, descry);
			for (k=0; k<8; k++) imgdx[k].update(i,j,descrx[k]);
			for (k=0; k<8; k++) imgdy[k].update(i,j,descry[k]);
		}

	}
	static double minSIFTDescriptorError(int decal, hvColRGB<double> descr1x[8], hvColRGB<double> descr1y[8], 
		hvColRGB<double> descr2x[8], hvColRGB<double> descr2y[8], hvColRGB<double> &minerr)
	{
		int i,j;
		double minv=10.0;
		//int decal = (int)(delta/M_PI*4.0);
		for (i=-(decal+1)+1; i<decal+1; i++)
		{
			hvColRGB<double> err(0.0,0.0,0.0);
			for (j=0; j<8; j++)
			{
					hvColRGB<double> diffx, diffy; 
					diffx.sub(descr1x[j],descr2x[(j+(i+8))%8]);
					diffx.square();
					diffy.sub(descr1y[j],descr2y[(j+(i+8))%8]);
					diffy.square();
					err += diffx; err += diffy;
			}
			double vv = err.RED()+err.GREEN()+err.BLUE();
			if (vv<minv) { minv=vv; minerr=err; }
			/*
			err=hvColRGB<double>(0.0,0.0,0.0);
			for (j=0; j<8; j++)
			{
					hvColRGB<double> diff; 
					diff.sub(descr1[j],descr2[(((8-j)%8)+(i+8))%8]);
					diff.square();
					err += diff;
			}
			vv = err.RED()+err.GREEN()+err.BLUE();
			if (vv<minv) { minv=vv; minerr=err; }
			*/
		}
		return minv;
	}
	static double minSIFTDescriptorError2(int decal, hvColRGB<double> descrx1[8], hvColRGB<double> descrx1b[8],
		hvColRGB<double> descry1[8], hvColRGB<double> descry1b[8],
		hvColRGB<double> descrx2[8], hvColRGB<double> descrx2b[8], 
		hvColRGB<double> descry2[8], hvColRGB<double> descry2b[8], hvColRGB<double> &minerr)
	{
		int i,j;
		double minv=10.0;
		//int decal = (int)(delta/M_PI*4.0);
		for (i=-(decal+1)+1; i<decal+1; i++)
		{
			hvColRGB<double> err(0.0,0.0,0.0);
			for (j=0; j<8; j++)
			{
					hvColRGB<double> diffx, diffy,diffxb, diffyb; 
					diffx.sub(descrx1[j],descrx2[(j+(i+8))%8]);
					diffx.square();
					diffy.sub(descry1[j],descry2[(j+(i+8))%8]);
					diffy.square();
					diffxb.sub(descrx1b[j],descrx2b[(j+(i+8))%8]);
					diffxb.square();
					diffyb.sub(descry1b[j],descry2b[(j+(i+8))%8]);
					diffyb.square();
					err += diffx; err += diffxb; err += diffy; err += diffyb;
			}
			double vv = err.RED()+err.GREEN()+err.BLUE();
			if (vv<minv) { minv=vv; minerr=err; }
			/*
			err=hvColRGB<double>(0.0,0.0,0.0);
			for (j=0; j<8; j++)
			{
					hvColRGB<double> diff, diffb; 
					diff.sub(descr1[j],descr2[(((8-j)%8)+(i+8))%8]);
					diff.square();
					diffb.sub(descr1b[j],descr2b[(((8-j)%8)+(i+8))%8]);
					diffb.square();
					err += diff;  err += diffb;
			}
			vv = err.RED()+err.GREEN()+err.BLUE();
			if (vv<minv) { minv=vv; minerr=err; }
			*/
		}
		return minv;
	}
	void mySIFTImage(int decal, float alpha, T scale, hvPictRGB<double> &errpi, 
		hvColRGB<double> descrx[8], hvColRGB<double> descry[8], double &minv, double &maxv)
	{
		errpi.reset(this->sizeX(), this->sizeY(), hvColRGB<double>(0.0));
		int i,j;
		hvColRGB<double> ddx[8], ddy[8], minerr;
		double vv;
		minv=1.0;maxv=0.0;
		for (i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			mySIFTDescriptor(alpha,i,j,scale,ddx, ddy);
			vv = hvPictRGB<T>::minSIFTDescriptorError(decal, descrx, descry,ddx, ddy, minerr);
			if(vv<minv && vv!=0.0) minv=vv;
			if (vv>maxv) maxv=vv;
			errpi.update(i,j,minerr);
		}
	}
	void mySIFTImage(int decal, T scale, hvPictRGB<double> &errpi, hvColRGB<double> descrx[8], hvColRGB<double> descry[8], 
		hvPictRGB<double> siftdescrx[8], hvPictRGB<double> siftdescry[8], double &minv, double &maxv)
	{
		errpi.reset(this->sizeX(), this->sizeY(), hvColRGB<double>(0.0));
		int i,j,k;
		hvColRGB<double> ddx[8], ddy[8], minerr;
		double vv;
		minv=1.0;maxv=0.0;
		for (i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			for (k=0; k<8; k++) ddx[k] = siftdescrx[k].get(i,j);
			for (k=0; k<8; k++) ddy[k] = siftdescry[k].get(i,j);
			vv = hvPictRGB<T>::minSIFTDescriptorError(decal,descrx,descry,ddx, ddy, minerr);
			if(vv<minv && vv!=0.0) minv=vv;
			if (vv>maxv) maxv=vv;
			errpi.update(i,j,minerr);
		}
	}

	// FFT
	void statfft(const hvPictRGB<T> &pict, hvPictComponent cc, double ss, double loga, bool centred, double scal, double offset, int pow_2, int niter, bool amplitude = false)
	{
		hvArray1<double> *resfft;
		hvField2<double> func(pict.sizeX(), pict.sizeY(), 0.0);
		int i, j;
		for (i = 0; i<pict.sizeX(); i++) for (j = 0; j<pict.sizeY(); j++)
		{
			double v;
			hvColRGB<T> co = pict.get(i, j);
			switch (cc)
			{
			case HV_RED: v = (double)co.RED(); break;
			case HV_GREEN: v = (double)co.GREEN(); break;
			case HV_BLUE: v = (double)co.BLUE(); break;
			default: v = (double)co.luminance();
			}
			func.update(i, j, v);
		}
		int nn = 1 << pow_2;
		resfft = func.statfft(centred, scal, offset, pow_2, niter, amplitude);
		this->hvField2< hvColRGB<T> >::reset(nn, nn, hvColRGB<T>(0));
		double max = resfft->get(0);
		for (i = 1; i<nn*nn; i++) if (resfft->get(i)>max) max = resfft->get(i);
		for (i = 0; i<nn; i++) for (j = 0; j<nn; j++)
		{
			double val = resfft->get(i + j*nn);
			double vv = log(1.0 + loga*(double)val / max) / log(loga + 1.0);
			if (vv>1.0) vv = 1.0;
			this->update(i, j, hvColRGB<T>(T(vv*ss)));
		}

	}



	// devides size by factor of 2^(nrec+1), nrec=0 means factor 2, nrec=1 means factor 4, etc.
	void shrink(hvPictRGB<T> *p, int nrec=0)
	{
		if (nrec>0) { hvPictRGB<T> source; source.shrink(p,0); this->shrink(&source, nrec-1); }
		else 
		{
			this->reset(p->sizeX()/2, p->sizeY()/2, hvColRGB<T>());
			int ii,jj;
			for (ii=0; ii<this->sizeX(); ii++) for (jj=0; jj<this->sizeY(); jj++)
			{
				hvColRGB<T> cc[4];
				cc[0]=p->get(2*ii,2*jj); cc[1]=p->get(2*ii+1,2*jj); cc[2]=p->get(2*ii+1,2*jj+1); cc[3]=p->get(2*ii,2*jj+1);
				hvColRGB<T> col; col.mean(4,cc);
				this->update(ii,jj,col);
			}
		}
	}
	void shrinkx(hvPictRGB<T> *p, int nrec=0)
	{
		if (nrec>0) { hvPictRGB<T> source; source.shrinkx(p,0); this->shrinkx(&source, nrec-1); }
		else 
		{
			this->reset(p->sizeX()/2, p->sizeY(), hvColRGB<T>());
			int ii,jj;
			for (ii=0; ii<this->sizeX(); ii++) for (jj=0; jj<this->sizeY(); jj++)
			{
				hvColRGB<T> cc[2];
				cc[0]=p->get(2*ii,jj); cc[1]=p->get(2*ii+1,jj);
				hvColRGB<T> col; col.mean(2,cc);
				this->update(ii,jj,col);
			}
		}
	}
	void shrinky(hvPictRGB<T> *p, int nrec=0)
	{
		if (nrec>0) { hvPictRGB<T> source; source.shrinky(p,0); this->shrinky(&source, nrec-1); }
		else 
		{
			this->reset(p->sizeX(), p->sizeY()/2, hvColRGB<T>());
			int ii,jj;
			for (ii=0; ii<this->sizeX(); ii++) for (jj=0; jj<this->sizeY(); jj++)
			{
				hvColRGB<T> cc[2];
				cc[0]=p->get(ii,2*jj); cc[1]=p->get(ii,2*jj+1); 
				hvColRGB<T> col; col.mean(2,cc);
				this->update(ii,jj,col);
			}
		}
	}
	// multiplies size by factor of 2^(nrec+1), nrec=0 means factor 2, nrec=1 means factor 4, etc.
	void enlarge(hvPictRGB<T> *p, int nrec=0)
	{
		if (nrec>0) { hvPictRGB<T> source; source.hvField2<hvColRGB<T> >::enlarge(p); this->enlarge(&source, nrec-1); }
		else hvField2<hvColRGB<T> >::enlarge(p);
	}
	// applies a local operator mm on a mask of size [-size,+size]
	hvVec3<double> maskPixel(T scal, int size, int x, int y, hvPictMask mm) const
	{
		int ii,jj,i,j;
		hvVec3<double> res(0.0);
		for (ii=-size; ii<=size; ii++) for (jj=-size; jj<=size; jj++)
		{
			int px = x+ii; if (px<0) px=0; else if (px>= this->sizeX()) px= this->sizeX()-1;
			int py = y+jj; if (py<0) py=0; else if (py>= this->sizeY()) py= this->sizeY()-1;
			hvColRGB<T> v = this->get(px,py);
			hvVec3<double> vd((double)v.RED()/(double)scal,(double)v.GREEN()/(double)scal,(double)v.BLUE()/(double)scal);
			switch(mm)
			{
			case HV_EDGE_VERT: if (ii<0) res-=vd; else if (ii>0) res+=vd; break;
			case HV_EDGE_HORIZ: if (jj<0) res-=vd; else if (jj>0) res+=vd; break;
			case HV_EDGE_DIAG1: if (jj>ii) res-=vd; else if (jj<ii) res+=vd; break;
			case HV_EDGE_DIAG2: if (jj>-ii) res-=vd; else if (jj<-ii) res+=vd; break;
			case HV_DESCR_VERT: if (ii<0) res-=vd; else res+=vd; break;
			case HV_DESCR_HORIZ: if (jj<0) res-=vd; else res+=vd; break;
			case HV_DESCR_DIAG1: if ((ii<0 && jj<0)||(ii>=0 && jj>=0)) res-=vd; else res+=vd; break;
			case HV_DESCR_DIAG2: if ((ii<0 && jj>0)||(ii>=0 && jj<=0)) res-=vd; else res+=vd; break;
			//case HV_DESCR_DOTS: i=ii<0?-ii:ii; j=jj<0?-jj:jj; i/=2; j/=2; if ((i%2==0 && j%2!=0)||(i%2!=0 && j%2==0)) res-=vd; else res+=vd; break;
			case HV_DESCR_DOTS: i=ii<0?-ii:ii; j=jj<0?-jj:jj; if ((i+j)%2==0) res-=vd; else res+=vd; break;
			default: res+=vd; break;
			}
		}
		return res;
	}
	void mask(const hvPictRGB<T> &p, T scal, int size, hvPictMask mm, double norm)
	{
		this->reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvVec3<double> res=p.maskPixel(scal,size,i,j,mm);
			res.scale(1.0/norm);
			res.abs();
			this->update(i,j,hvColRGB<T>((T)(res.X()*(double)scal),(T)(res.Y()*(double)scal),(T)(res.Z()*(double)scal))); 
		}
	}
	void erase(const hvBitmap &mm, const hvColRGB<T> &vv)
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			if (!mm.get(i,j)) this->update(i,j,vv);
		}
	}
	void edge(const hvPictRGB<T> &p, T scal, int size)
	{
		this->reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvVec3<double> vert=p.maskPixel(scal,size,i,j,HV_EDGE_VERT);
			hvVec3<double> horiz=p.maskPixel(scal,size,i,j,HV_EDGE_HORIZ);
			hvVec3<double> diag1=p.maskPixel(scal,size,i,j,HV_EDGE_DIAG1);
			hvVec3<double> diag2=p.maskPixel(scal,size,i,j,HV_EDGE_DIAG2);
			vert.scale(1.0/(double)(size*(size*2+1)));
			vert.abs();
			horiz.scale(1.0/(double)(size*(size*2+1)));
			horiz.abs();
			diag1.scale(1.0/(double)(size*(size*2+1)));
			diag1.abs();
			diag2.scale(1.0/(double)(size*(size*2+1)));
			diag2.abs();
			vert.keepMax(vert, horiz);
			vert.keepMax(vert, diag1);
			vert.keepMax(vert, diag2);
			this->update(i,j,hvColRGB<T>((T)(vert.X()*(double)scal),(T)(vert.Y()*(double)scal),(T)(vert.Z()*(double)scal)));
		}
	}
	void discont(const hvPictRGB<T> &p, const hvColRGB<T> c1, const hvColRGB<T> c2)
	{
		this->reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		int i,j,ii,jj;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> cc = p.get(i,j);
			bool eq=true;
			for (ii=-1; ii<=1;ii++) for (jj=-1; jj<=1; jj++)
			{
				int x = i+ii, y = j+jj;
				if (x<0) x+= this->sizeX(); x %= this->sizeX();
				if (y<0) y+= this->sizeY(); y %= this->sizeY();
				if (!cc.equals(p.get(x,y))) eq=false;
			}
			this->update(i,j,eq?c1:c2);
		}
	}

	hvColRGB<T> avg() const 
	{
		int i,j;
		double vr=0.0, vg=0.0, vb=0.0;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			vr+=(double)v.RED(); vg+=(double)v.GREEN(); vb+=(double)v.BLUE();
		}
		return hvColRGB<T>((T)(vr/(double)(this->sizeX()*this->sizeY())),(T)(vg/(double)(this->sizeX()*this->sizeY())),(T)(vb/(double)(this->sizeX()*this->sizeY())) ) ;
	}
	hvColRGB<T> avg(const hvBitmap &bm) const 
	{
		int i,j,n=0;
		double vr=0.0, vg=0.0, vb=0.0;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			if (bm.get(i,j))
			{
				hvColRGB<T> v = this->get(i,j);
				vr+=(double)v.RED(); vg+=(double)v.GREEN(); vb+=(double)v.BLUE();
				n++;
			}
		}
		return hvColRGB<T>((T)(vr/(double)(n)),(T)(vg/(double)(n)),(T)(vb/(double)(n)) ) ;
	}

	void gammaNormalizedMax(T scal, double power)
	{
		int i,j;
		hvColRGB<T> min,max;
		this->minmax(min,max);

		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			v.gammaNormalizedMax(max,scal, power);
			this->update(i,j,v);
		}
	}
	void minmax(hvColRGB<T> &min, hvColRGB<T> &max)
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			if (i==0 && j==0) { min=v; max=v; }
			min.keepMin(min,v);
			max.keepMax(max,v);
		}
	}
	void normalize(const hvColRGB<T> &min, const hvColRGB<T> &max, double scal)
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			v.normalize(min,max, scal);
			this->update(i,j,v);
		}
	}
	void luminance(T scal, double factor)
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			double rr=(double)v.RED()*factor; if (rr>(double)scal) rr=(double)scal;
			double gg=(double)v.GREEN()*factor; if (gg>(double)scal) gg=(double)scal;
			double bb=(double)v.BLUE()*factor; if (bb>(double)scal) bb=(double)scal;
			this->update(i,j,hvColRGB<T>(T((rr+gg+bb)/3.0)));
		}
	}
	void scale(T scal, double factor)
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			double rr=(double)v.RED()*factor; if (rr>(double)scal) rr=(double)scal;
			double gg=(double)v.GREEN()*factor; if (gg>(double)scal) gg=(double)scal;
			double bb=(double)v.BLUE()*factor; if (bb>(double)scal) bb=(double)scal;
			this->update(i,j,hvColRGB<T>(T(rr),T(gg),T(bb)));
		}
	}
	template <class X> void scale(T scal, const hvPict<X> &sc, double norm)
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			double factor = (double)sc.get(i,j)/norm;
			double rr=(double)v.RED()*factor; if (rr>(double)scal) rr=(double)scal;
			double gg=(double)v.GREEN()*factor; if (gg>(double)scal) gg=(double)scal;
			double bb=(double)v.BLUE()*factor; if (bb>(double)scal) bb=(double)scal;
			this->update(i,j,hvColRGB<T>(T(rr),T(gg),T(bb)));
		}
	}
	void step(T scal, T ss)
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			this->update(i,j,hvColRGB<T>(v.RED()<ss?0:scal, v.GREEN()<ss?0:scal, v.BLUE()<ss?0:scal));
		}
	}
	void drawRect(int px, int py, int sx, int sy, hvColRGB<T> v)
	{
		int i,j;
		for (i=px; i<px+sx; i++) for (j=py; j<py+sy; j++)
		{
			if (i>=0 && i<this->sizeX() && j>=0 && j<this->sizeY()) this->update(i,j,v);
		}
	}
	void drawRectBlend(int px, int py, int sx, int sy, hvColRGB<T> v, double alpha)
	{
		int i,j;
		for (i=px; i<px+sx; i++) for (j=py; j<py+sy; j++)
		{
			hvColRGB<unsigned char> cc; cc.blend(v,this->get(i,j),alpha);
			if (i>=0 && i<this->sizeX() && j>=0 && j<this->sizeY()) this->update(i,j,cc);
		}
	}

	void copy(int x, int y, const hvPictRGB<T> &pict)
	{
		int i,j;
		for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++)
		{
			this->update(x+i,y+j,pict.get(i,j));
		}
	}
	void copyRect(int px, int py, int x, int y, int sx, int sy, const hvPictRGB<T> &pict, const hvBitmap &mask)
	{
		int i,j;
		for (i=0; i<sx; i++) for (j=0; j<sy; j++)
		{
			if (mask.get(x+i, y+j))
			{
				if (px+i>=0 && px+i<this->sizeX() && py+j>=0 && py+j<this->sizeY()) this->update(px+i,py+j,pict.get(x+i,y+j));
			}
		}
	}
	void copyRect(int px, int py, int x, int y, int sx, int sy, const hvPictRGB<T> &pict)
	{
		int i,j;
		for (i=0; i<sx; i++) for (j=0; j<sy; j++)
		{
				if (px+i>=0 && px+i<this->sizeX() && py+j>=0 && py+j<this->sizeY()) this->update(px+i,py+j,pict.get(x+i,y+j));
		}
	}
	void copyRect(int px, int py, int x, int y, int sx, int sy, const hvPictRGB<T> &pict, hvColRGB<unsigned char> col, int dd=2, int bscale=1)
	{
		int i,j;
		for (i=0; i<sx; i++) for (j=0; j<sy; j++)
		{
				if (px+i>=0 && px+i<this->sizeX() && py+j>=0 && py+j<this->sizeY())
				{
					if (i<=dd || j<=dd || i>=sx-1-dd || j>=sy-1-dd) this->update(px+i,py+j,col);
					else this->update(px+i,py+j,pict.get(x+i*bscale,y+j*bscale));
				}
		}
	}
	void copyRectShadow(int px, int py, int x, int y, int sx, int sy, const hvPictRGB<T> &pict, hvColRGB<T> col, int hh)
	{
		int i,j;
		for (i=0; i<sx; i++) for (j=0; j<sy; j++)
		{
				if (px+i>=0 && px+i<this->sizeX() && py+j>=0 && py+j<this->sizeY())
				{
					if (i==0 || j==0 || i==sx-1 || j==sy-1) this->update(px+i,py+j,col);
					else this->update(px+i,py+j,pict.get(x+i,y+j));
				}
		}
		for (i=0; i<sx; i++) for (j=1; j<=hh; j++)
		{
				if (px+i+j>=0 && px+i+j<this->sizeX() && py-j>=0 && py-j<this->sizeY())
				{
					hvColRGB<T> cc = this->get(px+i+j,py-j);
					cc.scale(1.0-(double)j/(double)(hh+1));
					this->update(px+i+j,py-j,cc);
				}
		}
	}

	double minShift(double scale, const hvPictRGB<T> &pi, int posx, int posy, int deltax, int deltay, int &minx, int &miny)
	{
		minx=0; miny=0;
		double minv=1000000.0;
		int i,j;
		for (i=-deltax; i<=deltax; i++) for (j=-deltay; j<=deltay; j++)
		{
			hvColRGB<double> rr;
			this->squareDifference(scale,this->sizeX()/2-posx, this->sizeY()/2-posy,pi.sizeX()/2-posx+i,pi.sizeY()/2-posy+j,10,10,pi,rr);
			double vv = rr.RED()+rr.GREEN()+rr.BLUE();
			if (vv<minv) { minx=i, miny=j; minv=vv; }
		}
		return minv;
	}

	void copyWangRect(int px, int py, int x, int y, int sx, int sy, const hvPictRGB<T> &pict, const hvBitmap &mask)
	{
		int i,j;
		for (i=0; i<sx; i++) for (j=0; j<sy; j++)
		{
			if (mask.get(x+i, y+j))
			{
				if (px+i>=0 && px+i<this->sizeX() && py+j>=0 && py+j<this->sizeY())
				{
					hvColRGB<unsigned char> cc = pict.get(x+i,y+j);
					int tx=(px+i)/pict.sizeX(), ty=(py+j)/pict.sizeY();
					hvColRGB<unsigned char> col(cc.RED(), cc.GREEN(), (unsigned char)((tx * 4 + ty * 16 + 128) % 256));
					this->update(px+i,py+j,col);
				}
			}
		}
	}
	void copyWangRect(int px, int py, int x, int y, int sx, int sy, int id, const hvBitmap &mask)
	{
		int i,j;
		for (i=0; i<sx; i++) for (j=0; j<sy; j++)
		{
			if (mask.get(x+i, y+j))
			{
				if (px+i>=0 && px+i<this->sizeX() && py+j>=0 && py+j<this->sizeY())
				{
					int tx=(px+i)/mask.sizeX(), ty=(py+j)/mask.sizeY();
					hvColRGB<unsigned char> col((unsigned char)((id * 7) % 128 + 128), (unsigned char)(((id * 11) / 256) % 128 + 128), (unsigned char)((tx * 4 + ty * 16 + 128) % 256));
					this->update(px+i,py+j,col);
				}
			}
		}
	}

	void apply(T scal, const hvLinearTransform3<double> &t)
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			hvVec3<double> col((double)v.RED()/(double)scal, (double)v.GREEN()/(double)scal,(double)v.BLUE()/(double)scal);
			col = t.apply(col);
			v=hvColRGB<T>((T)(col.X()*(double)scal),(T)(col.Y()*(double)scal),(T)(col.Z()*(double)scal));
			this->update(i,j,v);
		}
	}
	void applyinverse(T scal, const hvFrame3<double> &fr, double offset, double rescal)
	{
		hvLinearTransform3<double> t(fr);
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			hvVec3<double> col(((double)v.RED()/(double)scal-offset)/rescal, ((double)v.GREEN()/(double)scal-offset)/rescal,((double)v.BLUE()/(double)scal-offset)/rescal);
			col = t.apply(col);
			v=hvColRGB<T>((T)(col.X()*(double)scal),(T)(col.Y()*(double)scal),(T)(col.Z()*(double)scal));
			this->update(i,j,v);
		}
	}
	void apply(T scal, const hvFrame3<double> &fr, double offset, double rescal)
	{
		hvLinearTransform3<double> t; t.inverseFrame3(fr);
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			hvVec3<double> col((double)v.RED()/(double)scal, (double)v.GREEN()/(double)scal,(double)v.BLUE()/(double)scal);
			col = t.apply(col);
			double rr = (col.X()*rescal+offset);
			if (rr<0.0) rr=0.0; else if (rr>1.0) rr=1.0;
			double gg = (col.Y()*rescal+offset);
			if (gg<0.0) gg=0.0; else if (gg>1.0) gg=1.0;
			double bb = (col.Z()*rescal+offset);
			if (bb<0.0) bb=0.0; else if (bb>1.0) bb=1.0;
			v=hvColRGB<T>((T)(rr*(double)scal),(T)(gg*(double)scal),(T)(bb*(double)scal));
			this->update(i,j,v);
		}
	}

	hvFrame3<double> pca(T scal, hvBitmap *mask=0) const
	{
		int i,j;
		hvVec3<double> sum;
		int npix=0;
		/*** computing mean value ***/
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			hvVec3<double> col;
			bool ok=true;
			if (mask!=0) ok = mask->get(i,j);
			if (ok)
			{ 
				npix++; 
				col = hvVec3<double>((double)v.RED()/(double)scal, (double)v.GREEN()/(double)scal,(double)v.BLUE()/(double)scal);
				sum += col;
			}
		}
		sum /= (double)(npix);
		
		/*** computing Covariance matrix covar ***/
		hvMat3<double> covar;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			bool ok=true;
			if (mask!=0) ok = mask->get(i,j);
			if (ok)
			{ 
				hvColRGB<T> v = this->get(i,j);
				hvVec3<double> col((double)v.RED()/(double)scal, (double)v.GREEN()/(double)scal,(double)v.BLUE()/(double)scal);
				col -= sum;
				covar += hvMat3<double>( 
					hvVec3<double>(col.X()*col.X(), col.X()*col.Y(), col.X()*col.Z()), 
					hvVec3<double>(col.Y()*col.X(), col.Y()*col.Y(), col.Y()*col.Z()),
					hvVec3<double>(col.Z()*col.X(), col.Z()*col.Y(), col.Z()*col.Z())   );
			}
		}
		hvMat3<double> rr = covar.eigen();
		return hvFrame3<double>(sum, rr);
	}


	static void readPPMLine(FILE *fd, char buff[256]) {
	do {
        fgets(buff,255,fd);
        } while(buff[0]=='#');
	}

	void loadPPM(FILE *fd, T norm)
	{
		int  sx, sy;
		int i,j, type;
		char buff[256];
		hvColRGB<T> co;

		readPPMLine(fd,buff);
		if (strcmp(buff,"P6\n")==0) type = 0;
		else if (strcmp(buff,"P3\n")==0) type = 1;
		else { type = 2; printf("unknown picture PPM type=%d (%s)\n", type,buff); } 
		readPPMLine(fd,buff);
		sscanf(buff,"%d %d",&sx,&sy);
		readPPMLine(fd,buff);
		if (strcmp(buff,"255\n")!=0){ printf("type=%d\n", type); reset(); }
		reset(sx, sy, hvColRGB<T>()); 
		for (i=0; i<sy; i++)
		for (j=0; j<sx; j++)
			{
				unsigned char r,g,b;
				if (type==0)
				{
					fread(&r,1,sizeof(unsigned char),fd);
					fread(&g,1,sizeof(unsigned char),fd);
					fread(&b,1,sizeof(unsigned char),fd);
				}
				else if (type==1)
				{
					int rr, gg, bb;
					fscanf(fd, "%d %d %d", &rr, &gg, &bb);
					r= (unsigned char)rr;
					g= (unsigned char)gg;
					b= (unsigned char)bb;
				}
				else { r=0; g=0; b=0; }
				hvArray2< hvColRGB<T> >::update(j,sy-i-1,hvColRGB<T>((T)r/norm, (T)g/norm, (T)b/norm));
			}
	}

	void savePPM(FILE *fd, T norm)
	{
		int i,j;
		hvColRGB<T> co;
		unsigned char v;

		fprintf(fd,"P6\n");
		fprintf(fd,"%d %d\n",this->sizeX(),this->sizeY());
		fprintf(fd,"255\n");
		for (i=0; i<this->sizeY(); i++)
		for (j=0; j<this->sizeX(); j++)
			{
			co = hvArray2< hvColRGB<T> >::get(j, this->sizeY()-i-1);
			v = (unsigned char)((T)co.RED()*norm);
			fwrite(&v,1,sizeof(unsigned char),fd);
			v = (unsigned char)((T)co.GREEN()*norm);
			fwrite(&v,1,sizeof(unsigned char),fd);
			v = (unsigned char)((T)co.BLUE()*norm);
			fwrite(&v,1,sizeof(unsigned char),fd);
			}
	}

	static void RGBE_WriteHeader(FILE *fp, int width, int height, rgbe_header_info *info)
	{
	 const char *programtype = "RGBE";

	  if (info && (info->valid & RGBE_VALID_PROGRAMTYPE)) programtype = info->programtype;
	  fprintf(fp,"#?%s\n",programtype);
	  if (info && (info->valid & RGBE_VALID_GAMMA)) { fprintf(fp,"GAMMA=%g\n",info->gamma); }
	  if (info && (info->valid & RGBE_VALID_EXPOSURE)) { fprintf(fp,"EXPOSURE=%g\n",info->exposure); }
	  fprintf(fp,"FORMAT=32-bit_rle_rgbe\n\n");
	  fprintf(fp, "-Y %d +X %d\n", height, width);
	}

	bool RGBE_ReadHeader(FILE *fp, int *width, int *height, rgbe_header_info *info)
	{
	  char buf[256];
	  int found_format;
	  float tempf;
	  int i;

	  found_format = 0;
	  if (info) 
	  {
		info->valid = 0;
		info->programtype[0] = 0;
		info->gamma = 1.0;
		info->exposure = 1.0;
	  }
	  else hvFatal("no rgbe_header_info pointer given");
	  if (fgets(buf,256,fp) == 0) return false;
	  if ((buf[0] != '#')||(buf[1] != '?')) return false;
	  info->valid |= RGBE_VALID_PROGRAMTYPE;
	  for(i=0;i<sizeof(info->programtype)-1;i++) 
		{
		  if ((buf[i+2] == 0) || buf[i+2]==' ') break;
		  info->programtype[i] = buf[i+2];
		}
	  info->programtype[i] = 0;
	  bool cont=true;
	  while(cont) 
	  {
		if (fgets(buf,256,fp) == 0) return false;
		if (strcmp(buf,"FORMAT=32-bit_rle_rgbe\n") == 0) { }
		else if (sscanf(buf,"GAMMA=%g",&tempf) == 1) 
		{
		  info->gamma = tempf;
		  info->valid |= RGBE_VALID_GAMMA;
		}
		else if (sscanf(buf,"EXPOSURE=%g",&tempf) == 1) 
		{
		  info->exposure = tempf;
		  info->valid |= RGBE_VALID_EXPOSURE;
		}
		else if (sscanf(buf,"-Y %d +X %d",height,width) == 2) cont=false;
		//if (cont) if (fgets(buf,256,fp) == 0) return false;
	  }
	  return true;
	}

	void RGBE_WritePixels(FILE *fp, T norm)
	{
		unsigned char rgbe[4];
		int i,j;
		hvColRGB<T> co;
		unsigned char v;

		for (i=0; i<this->sizeY(); i++)
		for (j=0; j<this->sizeX(); j++)
			{
			co = hvArray2< hvColRGB<T> >::get(j, this->sizeY()-i-1);
			co.torgbe(norm, rgbe);
			fwrite(rgbe, 4, 1, fp);
			}
	}

	bool RGBE_ReadPixels_RLE(FILE *fp, T norm, int scanline_width, int num_scanlines)
	{
	  unsigned char rgbe[4], *scanline_buffer, *ptr, *ptr_end;
	  int i, count;
	  unsigned char buf[2];

	  if ((scanline_width < 8)||(scanline_width > 0x7fff)) { printf("not RLE encoded\n"); return false; } //RGBE_ReadPixels(fp,data,scanline_width*num_scanlines);
	  reset(scanline_width, num_scanlines, hvColRGB<T>());
	  scanline_buffer = 0;
	  /* read in each successive scanline */
	  while(num_scanlines > 0) 
	  {
		    if (fread(rgbe,sizeof(rgbe),1,fp) < 1) { printf("file corrupt in line %d\n",num_scanlines);  if (scanline_buffer!=0) free(scanline_buffer); return false; }
			if ((rgbe[0] != 2)||(rgbe[1] != 2)||(rgbe[2] & 0x80)) 
			{
				printf("file is not run length encoded in line %d\n",num_scanlines);
				/* this file is not run length encoded */
				//rgbe2float(&data[0],&data[1],&data[2],rgbe);
				//data += RGBE_DATA_SIZE;
				if (scanline_buffer!=0) free(scanline_buffer);
				return false; //RGBE_ReadPixels(fp,data,scanline_width*num_scanlines-1);
			}
			if ((((int)rgbe[2])<<8 | rgbe[3]) != scanline_width) 
			{
				if (scanline_buffer!=0)  free(scanline_buffer);
				printf("wrong scanline width in line %d\n",num_scanlines);
				return false;
			}
			if (scanline_buffer == 0) scanline_buffer = (unsigned char *)malloc(sizeof(unsigned char)*4*scanline_width);	    
			ptr = &scanline_buffer[0];
		/* read each of the four channels for the scanline into the buffer */
			for(i=0;i<4;i++) 
			{
				ptr_end = &scanline_buffer[(i+1)*scanline_width];
				while(ptr < ptr_end) 
				{
					if (fread(buf,sizeof(buf[0])*2,1,fp) < 1) { free(scanline_buffer); printf("file corrupt\n"); return false; }
					if (buf[0] > 128) 
					{
						/* a run of the same value */
						count = buf[0]-128;
						if ((count == 0)||(count > ptr_end - ptr)) { free(scanline_buffer); printf("bad scanline data"); return false; }
						while(count-- > 0) *ptr++ = buf[1];
					}
					else 
					{
						/* a non-run */
						count = buf[0];
						if ((count == 0)||(count > ptr_end - ptr)) { free(scanline_buffer); printf("bad scanline data"); return false; }
						*ptr++ = buf[1];
						if (--count > 0) 
						{
							if (fread(ptr,sizeof(*ptr)*count,1,fp) < 1) { free(scanline_buffer); printf("file corrupt\n"); return false; }
							ptr += count;
						}
					}
				}
			}
		/* now convert data from buffer into floats */
		for(i=0;i<scanline_width;i++) 
		{
		  rgbe[0] = scanline_buffer[i];
		  rgbe[1] = scanline_buffer[i+scanline_width];
		  rgbe[2] = scanline_buffer[i+2*scanline_width];
		  rgbe[3] = scanline_buffer[i+3*scanline_width];
		  hvColRGB<T> col; col.fromrgbe(norm, rgbe);
		  hvArray2< hvColRGB<T> >::update(i, num_scanlines-1,  col);
		}
		num_scanlines--;
	  }
	  free(scanline_buffer);
	  return true;
	}

	// applies a quantification into nr colors, but on the PCA of the image
	void segmentPCA(hvPict<int> &pi, std::vector<hvColRGB<unsigned char> > &lcol, int sx, int sy, int nr) const;


	// Find the best neigbor list match on (x,y) using a simple square difference
	// bm is the set of pixels that can be used to compute the match
	// neighborhood mask arround (x,y) is given by bseed 
	// use bcheck to check only among these pixels (bcheck is included into bm)
	int bestNeighborMatch(const hvBitmap &bcheck, const hvBitmap &bseed, const hvBitmap &bm, int neighbor, int x, int y, int nb, int rx[], int ry[])
	{
		int i,j,k, nx, ny;
		std::vector<double> lerr(MAX_SYNTHESIS_PIXEL_LIST);
		std::vector<int> lx(MAX_SYNTHESIS_PIXEL_LIST), ly(MAX_SYNTHESIS_PIXEL_LIST);

		lerr.clear(); lx.clear(); ly.clear();
		if (nb>MAX_SYNTHESIS_PIXEL_LIST) nb=MAX_SYNTHESIS_PIXEL_LIST;

		for (i=0; i<bcheck.sizeX(); i++)
		for (j=0; j<bcheck.sizeY(); j++)
		{
			if (bcheck.get(i,j)==true)
			{
				hvColRGB<double> res;
				int count=squareDifference(x,y,i,j,neighbor, bseed, bm, res);
				//if (count>=neighbor*neighbor/4)
				if (count>=neighbor*neighbor)
				{
					res.scale(1.0/(double)count);
					double err = sqrt(res.luminance())*255.0;
					for (k=0; k<lerr.size(); k++) if (err<lerr[k]) break;
					if (k<nb)
					{
						if (lerr.size()==nb) { lx.pop_back(); ly.pop_back(); lerr.pop_back();  }
						lx.insert(lx.begin()+k,i); ly.insert(ly.begin() + k,j); lerr.insert(lerr.begin() + k, err);
					}
				}
			}
		}
		if (lx.size()==0) { hvFatal("cannot select a best match pixel in hvPictRGB"); }
		//printf("%d,%d -> errors:",x,y);
		//for (k=0; k<lx.size(); k++) printf("%d(%d,%d) ",(unsigned int)lerr[k],lx[k],sizeY()-ly[k]-1);
		//printf("\n");
		for(k=0; k<lx.size(); k++) { rx[k]=lx[k]; ry[k]=ly[k]; }

		for (nx=-neighbor;nx<=neighbor; nx++)
		{
		for (ny=-neighbor;ny<=neighbor; ny++)
		{
			int px, py, epx, epy;
			px = x+nx; if (px<0) px+=this->sizeX(); else if (px>= this->sizeX()) px-= this->sizeX();
			py = y+ny; if (py<0) py+= this->sizeY(); else if (py>= this->sizeY()) py-= this->sizeY();
			epx = rx[0]+nx; if (epx<0) epx+= this->sizeX(); else if (epx>= this->sizeX()) epx-= this->sizeX();
			epy = ry[0]+ny; if (epy<0) epy+= this->sizeY(); else if (epy>= this->sizeY()) epy-= this->sizeY();
			if (bseed.get(px,py)==true && bm.get(epx, epy)==true)
			{
				double r,g,b;
				hvColRGB<unsigned char> col = this->get(px,py);
				r = (double)col.RED()/255.0;
				g = (double)col.GREEN()/255.0;
				b = (double)col.BLUE()/255.0;
				col = this->get(epx,epy);
				r = r - (double)col.RED()/255.0;
				g = g - (double)col.GREEN()/255.0;
				b = b - (double)col.BLUE()/255.0;
				//printf("%03d ", (unsigned char)(sqrt((r*r+g*g+b*b)/3.0)*255.0));
			}
			//else if (bseed.get(px,py)==true) printf("/// ");
			//else if (bm.get(px,py)==true) printf("### ");
			//else printf("--- ");
		}
		//printf("\n");
		}

		return lx.size();
	}

	// Find the best neigbor list match on (x,y) within a Gaussian Pyramid using a simple square difference
	// bm is the set of pixels that can be used to compute the match
	// neighborhood mask arround (x,y) is given by bseed 
	// use bcheck to check only among these pixels (bcheck is included into bm)
	int bestNeighborMatchPyramid(hvPictRGB<unsigned char> exlod[], int nlevels, int curlevel, const hvBitmap &bcheck, const hvBitmap &bseed, const hvBitmap &bm, int neighbor, int x, int y, int nb, int rx[], int ry[])
	{
		int i,j,k, nx, ny;
		std::vector<double> lerr(MAX_SYNTHESIS_PIXEL_LIST);
		std::vector<int> lx(MAX_SYNTHESIS_PIXEL_LIST), ly(MAX_SYNTHESIS_PIXEL_LIST);

		lerr.clear(); lx.clear(); ly.clear();

		if (nb>MAX_SYNTHESIS_PIXEL_LIST) nb=MAX_SYNTHESIS_PIXEL_LIST;

		for (i=neighbor; i<bcheck.sizeX()-neighbor; i++)
		for (j=neighbor; j<bcheck.sizeY()-neighbor; j++)
		{
			if (bcheck.get(i,j)==true)
			{
				hvColRGB<double> res;
				int count=squareDifference(255,x,y,i,j,neighbor, bseed, bm, res);
				//if (count>=neighbor*neighbor/4)
				if (count>=neighbor*neighbor)
				{
					// add square difference of next higher level
					int plx=x/2, ply=y/2;
					int pli=i/2, plj=j/2;
					int ll=curlevel+1;
					if (ll<nlevels)

					//for (int ll=curlevel+1; ll<nlevels; ll++)
					{
						int neigh = neighbor - ll;
						if (neigh<2) neigh=2;
						hvColRGB<double> pres;
						int nr = exlod[ll].squareDifference(255,plx, ply, pli, plj, neigh, pres);
						count+=nr;
						res += pres;
						//plx/=2; ply/=2; pli/=2; plj/=2;
					}
					// stack final error to keep lowest
					res.scale(1.0/(double)count);
					double err = sqrt(res.luminance())*255.0;
					for (k=0; k<lerr.size(); k++) if (err<lerr[k]) break;
					if (k<nb)
					{
						if (lerr.size()==nb) { lx.pop_back(); ly.pop_back(); lerr.pop_back();  }
						lx.insert(lx.begin() + k, i); ly.insert(ly.begin() + k, j); lerr.insert(lerr.begin() + k, err);
					}
				}
			}
		}
		if (lx.size()==0) { hvFatal("cannot select a best match pixel in hvPictRGB"); }
		//printf("%d,%d -> errors:",x,y);
		//for (k=0; k<lx.length(); k++) printf("%d(%d,%d) ",(unsigned int)lerr[k],lx[k],sizeY()-ly[k]-1);
		//printf("\n");
		for(k=0; k<lx.size(); k++) { rx[k]=lx[k]; ry[k]=ly[k]; }
		return lx.size();
	}

	// Use a random walk like strategy to find best match
	void bestNeighborMatch(const hvPictRGB<unsigned char> &ex, int x, int y, int neighbor, int &rx, int &ry, int SX = 0, int SY = 0, int DX = 0, int DY = 0)
	{
		const int NSAMPLES = 5;
		const int DEPTH = 2;
		int NITER = NSAMPLES*(DX*DY+1);
		if (SX == 0) { SX = ex.sizeX(); }
		if (SY == 0) { SY = ex.sizeY(); }
		int RADIUS = SX < SY ? SX / 4 : SY / 4;
		if (RADIUS < 1) RADIUS = 1;
		int i, j, k;
		double minerr, searcherr, besterr;
		int spx, spy, bestx, besty;
		for (i = 0; i < NITER; i++)
		{
			for (j = 0; j < NSAMPLES; j++) // first level is random
			{
				int px = (int)((double)rand() / (double)RAND_MAX*(double)SX);
				if (px < neighbor) px = neighbor;
				if (px >= SX-neighbor) px = SX - neighbor-1;
				int py = (int)((double)rand() / (double)RAND_MAX*(double)SY);
				if (py < neighbor) py = neighbor;
				if (py >= SY-neighbor) py = SY - neighbor-1;
				int deltax = (int)((double)rand() / (double)RAND_MAX*(double)DX);
				if (deltax >= DX) deltax = DX - 1;
				int deltay = (int)((double)rand() / (double)RAND_MAX*(double)DY);
				if (deltay >= DY) deltay = DY - 1;
				px += SX*deltax; py += SY*deltay;
				double err = this->meanSquareDifference(ex, 255.0, x, y, px, py, neighbor);
				if (j == 0) { searcherr = err; spx = px; spy = py; }
				else if (err < searcherr) { searcherr = err; spx = px; spy = py; }
			}
			minerr = searcherr; bestx = spx; besty = spy;
			for (k = 0; k < DEPTH && RADIUS>1; k++)
			{
				int deltax = bestx / SX, deltay = besty / SY;
				for (j = 0; j < NSAMPLES; j++) // next levels are closer and closer to previous best
				{
					int px = bestx - deltax*SX - RADIUS + (int)(2.0*(double)rand() / (double)RAND_MAX*(double)RADIUS);
					if (px < neighbor) px = neighbor;
					if (px >= SX-neighbor) px = SX -neighbor - 1;
					int py = besty - deltay*SY - RADIUS + (int)(2.0*(double)rand() / (double)RAND_MAX*(double)RADIUS);
					if (py < neighbor) py = neighbor;
					if (py >= SY-neighbor) py = SY - neighbor - 1;
					px += SX*deltax; py += SY*deltay;
					double err = this->meanSquareDifference(ex, 255.0, x, y, px, py, neighbor);
					if (err < searcherr) { searcherr = err; spx = px; spy = py; }
				}
				if (searcherr < minerr) { minerr = searcherr; bestx = spx; besty = spy; }
				RADIUS /= 2; if (RADIUS < 1) RADIUS = 1;
			}
			if (i == 0) { besterr = minerr; rx = bestx; ry = besty; }
			else if (besterr>minerr) { besterr = minerr; rx = bestx; ry = besty; }
		}
	}
	// Use a random walk like strategy to find best match
	void refineBestNeighborMatch(const hvPictRGB<unsigned char> &ex, int x, int y, int neighbor, int ix[], int iy[], int nn, int &rx, int &ry, int SX = 0, int SY = 0, int DX = 0, int DY = 0)
	{
		const int NSAMPLES = 5;
		const int DEPTH = 1;
		if (SX == 0) { SX = ex.sizeX(); }
		if (SY == 0) { SY = ex.sizeY(); }
		int RADIUS = neighbor*neighbor;
		if (RADIUS < 1) RADIUS = 1;
		int i, j, k;
		double minerr, searcherr, besterr;
		int spx, spy, bestx, besty;
		for (i = 0; i < nn; i++)
		{
			int deltax = ix[i] / SX, deltay = iy[i] / SY;
			for (j = 0; j < NSAMPLES; j++) // first level is random
			{
				int px = ix[i] - deltax*SX - RADIUS + (int)(2.0*(double)rand() / (double)RAND_MAX*(double)RADIUS);
				if (px < neighbor) px = neighbor;
				if (px >= SX - neighbor) px = SX - neighbor - 1;
				int py = iy[i] - deltay*SY - RADIUS + (int)(2.0*(double)rand() / (double)RAND_MAX*(double)RADIUS);
				if (py < neighbor) py = neighbor;
				if (py >= SY - neighbor) py = SY - neighbor - 1;
				px += SX*deltax; py += SY*deltay;
				double err = this->meanSquareDifference(ex, 255.0, x, y, px, py, neighbor);
				if (j == 0) { searcherr = err; spx = px; spy = py; }
				else if (err < searcherr) { searcherr = err; spx = px; spy = py; }
			}
			minerr = searcherr; bestx = spx; besty = spy;
			RADIUS /= 2; if (RADIUS < 1) RADIUS = 1;
			for (k = 0; k < DEPTH && RADIUS>1; k++)
			{
				int deltax = bestx / SX, deltay = besty / SY;
				for (j = 0; j < NSAMPLES; j++) // next levels are closer and closer to previous best
				{
					int px = bestx - deltax*SX - RADIUS + (int)(2.0*(double)rand() / (double)RAND_MAX*(double)RADIUS);
					if (px < neighbor) px = neighbor;
					if (px >= SX - neighbor) px = SX - neighbor - 1;
					int py = besty - deltay*SY - RADIUS + (int)(2.0*(double)rand() / (double)RAND_MAX*(double)RADIUS);
					if (py < neighbor) py = neighbor;
					if (py >= SY - neighbor) py = SY - neighbor - 1;
					px += SX*deltax; py += SY*deltay;
					double err = this->meanSquareDifference(ex, 255.0, x, y, px, py, neighbor);
					if (err < searcherr) { searcherr = err; spx = px; spy = py; }
				}
				if (searcherr < minerr) { minerr = searcherr; bestx = spx; besty = spy; }
				RADIUS /= 2; if (RADIUS < 1) RADIUS = 1;
			}
			if (i == 0) { besterr = minerr; rx = bestx; ry = besty; }
			else if (besterr > minerr) { besterr = minerr; rx = bestx; ry = besty; }
		}
	}
	////////////////////////////////////////////////////////////
	// compute square differences
	////////////////////////////////////////////////////////////
	int squareDifference(double scale, int x, int y, int i, int j, int neighbor, const hvBitmap &bseed, const hvBitmap &bm, hvColRGB<double> &res)
	{
		int count=0, nx, ny;
		double errr=0.0, errg=0.0, errb=0.0;
		for (nx=-neighbor;nx<=neighbor; nx++)
		for (ny=-neighbor;ny<=neighbor; ny++)
		{
			int px, py, epx, epy;
			px = x+nx; if (px<0) px+=this->sizeX(); else if (px>= this->sizeX()) px-= this->sizeX();
			py = y+ny; if (py<0) py+= this->sizeY(); else if (py>= this->sizeY()) py-= this->sizeY();
			epx = i+nx; if (epx<0) epx+= this->sizeX(); else if (epx>= this->sizeX()) epx-= this->sizeX();
			epy = j+ny; if (epy<0) epy+= this->sizeY(); else if (epy>= this->sizeY()) epy-= this->sizeY();
			if (bseed.get(px,py)==true && bm.get(epx, epy)==true)
			{
				double r,g,b;
				hvColRGB<unsigned char> col = this->get(px,py);
				r = (double)col.RED()/scale;
				g = (double)col.GREEN()/scale;
				b = (double)col.BLUE()/scale;
				col = this->get(epx,epy);
				r = r - (double)col.RED()/scale;
				g = g - (double)col.GREEN()/scale;
				b = b - (double)col.BLUE()/scale;
				errr += r*r; errg += g*g; errb += b*b;
				count++;
			}
		}
		res = hvColRGB<double>(errr, errg, errb);
		return count;
	}
	
	double meanSquareDifference(const hvPictRGB<unsigned char> &ex, double scale, int x, int y, int i, int j, int neighbor)
	{
		int count = 0, nx, ny;
		double errr = 0.0, errg = 0.0, errb = 0.0;
		for (nx = -neighbor; nx <= neighbor; nx++)
			for (ny = -neighbor; ny <= neighbor; ny++)
			{
				int px, py, epx, epy;
				px = x + nx; if (px<0) px += this->sizeX(); else if (px >= this->sizeX()) px -= this->sizeX();
				py = y + ny; if (py<0) py += this->sizeY(); else if (py >= this->sizeY()) py -= this->sizeY();
				epx = i + nx; if (epx<0) epx += this->sizeX(); else if (epx >= this->sizeX()) epx -= this->sizeX();
				epy = j + ny; if (epy<0) epy += this->sizeY(); else if (epy >= this->sizeY()) epy -= this->sizeY();
				if (px>=0 && px<this->sizeX() && py >= 0 && py<this->sizeY() &&
					epx>=0 && epx<ex.sizeX() && epy >= 0 && epy<ex.sizeY() )
				{
					double r, g, b;
					hvColRGB<unsigned char> col = this->get(px, py);
					r = (double)col.RED() / scale;
					g = (double)col.GREEN() / scale;
					b = (double)col.BLUE() / scale;
					col = ex.get(epx, epy);
					r = r - (double)col.RED() / scale;
					g = g - (double)col.GREEN() / scale;
					b = b - (double)col.BLUE() / scale;
					errr += r*r; errg += g*g; errb += b*b;
					count++;
				}
			}
		if (count == 0) hvFatal("could not find neighborhood");
		return (errr+errg+errb)/3.0/(double)count;
	}

	int squareDifference(double scale, int x, int y, int i, int j, int neighbor, hvColRGB<double> &res)
	{
		int count=0, nx, ny;
		double errr=0.0, errg=0.0, errb=0.0;
		for (nx=-neighbor;nx<=neighbor; nx++)
		for (ny=-neighbor;ny<=neighbor; ny++)
		{
			int px, py, epx, epy;
			px = x+nx; if (px<0) px+= this->sizeX(); else if (px>= this->sizeX()) px-= this->sizeX();
			py = y+ny; if (py<0) py+= this->sizeY(); else if (py>= this->sizeY()) py-= this->sizeY();
			epx = i+nx; if (epx<0) epx+= this->sizeX(); else if (epx>= this->sizeX()) epx-= this->sizeX();
			epy = j+ny; if (epy<0) epy+= this->sizeY(); else if (epy>= this->sizeY()) epy-= this->sizeY();
			double r,g,b;
			hvColRGB<unsigned char> col = this->get(px,py);
			r = (double)col.RED()/scale;
			g = (double)col.GREEN()/scale;
			b = (double)col.BLUE()/scale;
			col = this->get(epx,epy);
			r = r - (double)col.RED()/scale;
			g = g - (double)col.GREEN()/scale;
			b = b - (double)col.BLUE()/scale;
			errr += r*r; errg += g*g; errb += b*b;
			count++;
		}
		res = hvColRGB<double>(errr, errg, errb);
		return count;
	}
	void squareDifference(double scale, int px, int py, int x, int y, int sx, int sy, const hvPictRGB<T> &pi, const hvBitmap &mask, hvColRGB<double> &res)
	{
		int count=0, i,j;
		double errr=0.0, errg=0.0, errb=0.0;
		for (i=0; i<sx;i++)
		for (j=0; j<sy;j++)
		{
			if (mask.get(x+i,y+j))
			{
				double r,g,b;
				//if (px+i<0 || px+i>=sizeX() || py+j<0 || py+j>=sizeY()) { printf("out of this picture range: %d,%d\n", px+i,py+j); }
				int ppx=px+i, ppy=py+j;
				if (ppx<0) ppx+= this->sizeX(); else if (ppx>= this->sizeX()) ppx-= this->sizeX();
				if (ppy<0) ppy+= this->sizeY(); else if (ppy>= this->sizeY()) ppy-= this->sizeY();
				hvColRGB<unsigned char> col = this->get(ppx,ppy);
				r = (double)col.RED()/scale;
				g = (double)col.GREEN()/scale;
				b = (double)col.BLUE()/scale;
				//if (x+i<0 || x+i>=pi.sizeX() || y+j<0 || y+j>=pi.sizeY()) { printf("out of pi picture range: %d,%d\n", x+i,y+j); }
				col = pi.get(x+i,y+j);
				r = r - (double)col.RED()/scale;
				g = g - (double)col.GREEN()/scale;
				b = b - (double)col.BLUE()/scale;
				errr += r*r; errg += g*g; errb += b*b;
			}
		}
		res = hvColRGB<double>(errr, errg, errb);
	}
	void squareDifference(double scale, int px, int py, int x, int y, int sx, int sy, const hvPictRGB<T> &pi, hvColRGB<double> &res, int step=1)
	{
		int count=0, i,j;
		double errr=0.0, errg=0.0, errb=0.0;
		for (i=0; i<sx;i+=step)
		for (j=0; j<sy;j+=step)
		{
				double r,g,b;
				//if (px+i<0 || px+i>=sizeX() || py+j<0 || py+j>=sizeY()) { printf("out of this picture range: %d,%d\n", px+i,py+j); }
				int ppx=px+i, ppy=py+j;
				if (ppx<0) ppx+= this->sizeX(); else if (ppx>= this->sizeX()) ppx-= this->sizeX();
				if (ppy<0) ppy+= this->sizeY(); else if (ppy>= this->sizeY()) ppy-= this->sizeY();
				hvColRGB<unsigned char> col = this->get(ppx,ppy);
				r = (double)col.RED()/scale;
				g = (double)col.GREEN()/scale;
				b = (double)col.BLUE()/scale;
				//if (x+i<0 || x+i>=pi.sizeX() || y+j<0 || y+j>=pi.sizeY()) { printf("out of pi picture range: %d,%d\n", x+i,y+j); }
				col = pi.get(x+i,y+j);
				r = r - (double)col.RED()/scale;
				g = g - (double)col.GREEN()/scale;
				b = b - (double)col.BLUE()/scale;
				errr += r*r; errg += g*g; errb += b*b;
		}
		res = hvColRGB<double>(errr, errg, errb);
	}
	void squareDifferenceBorder(double scale, int px, int py, int x, int y, int sx, int sy, const hvPictRGB<T> &pi, hvColRGB<double> &res, int depth = 1)
	{
		int count = 0, i, j;
		double errr = 0.0, errg = 0.0, errb = 0.0;
		for (i = 0; i<sx; i += 1)
			for (j = 0; j<sy; j += 1)
			{
				if (i < depth || i >= sx - depth || j < sy || j >= sy - depth)
				{
					double r, g, b;
					//if (px+i<0 || px+i>=sizeX() || py+j<0 || py+j>=sizeY()) { printf("out of this picture range: %d,%d\n", px+i,py+j); }
					int ppx = px + i, ppy = py + j;
					if (ppx < 0) ppx += this->sizeX(); else if (ppx >= this->sizeX()) ppx -= this->sizeX();
					if (ppy < 0) ppy += this->sizeY(); else if (ppy >= this->sizeY()) ppy -= this->sizeY();
					hvColRGB<unsigned char> col = this->get(ppx, ppy);
					r = (double)col.RED() / scale;
					g = (double)col.GREEN() / scale;
					b = (double)col.BLUE() / scale;
					//if (x+i<0 || x+i>=pi.sizeX() || y+j<0 || y+j>=pi.sizeY()) { printf("out of pi picture range: %d,%d\n", x+i,y+j); }
					col = pi.get(x + i, y + j);
					r = r - (double)col.RED() / scale;
					g = g - (double)col.GREEN() / scale;
					b = b - (double)col.BLUE() / scale;
					errr += r*r; errg += g*g; errb += b*b;
				}
			}
		res = hvColRGB<double>(errr, errg, errb);
	}

	void weightedSquareDifference(double scale, int px, int py, int x, int y, int sx, int sy, const hvPictRGB<T> &pi, const hvPict<double> &weight, hvColRGB<double> &res)
	{
		int count=0, i,j;
		double errr=0.0, errg=0.0, errb=0.0;
		for (i=0; i<sx;i++)
		for (j=0; j<sy;j++)
		{
			if (weight.get(x+i,y+j)>0.0)
			{
				double r,g,b;
				double ww = weight.get(x+i,y+j);
				//if (px+i<0 || px+i>=sizeX() || py+j<0 || py+j>=sizeY()) { printf("out of this picture range: %d,%d\n", px+i,py+j); }
				int ppx=px+i, ppy=py+j;
				if (ppx<0) ppx+= this->sizeX(); else if (ppx>= this->sizeX()) ppx-= this->sizeX();
				if (ppy<0) ppy+= this->sizeY(); else if (ppy>= this->sizeY()) ppy-= this->sizeY();
				hvColRGB<unsigned char> col = this->get(ppx,ppy);
				r = (double)col.RED()/scale;
				g = (double)col.GREEN()/scale;
				b = (double)col.BLUE()/scale;
				//if (x+i<0 || x+i>=pi.sizeX() || y+j<0 || y+j>=pi.sizeY()) { printf("out of pi picture range: %d,%d\n", x+i,y+j); }
				col = pi.get(x+i,y+j);
				r = r - (double)col.RED()/scale;
				g = g - (double)col.GREEN()/scale;
				b = b - (double)col.BLUE()/scale;
				errr += r*r*ww; errg += g*g*ww; errb += b*b*ww;
			}
		}
		res = hvColRGB<double>(errr, errg, errb);
	}
	/*
	// synthesis of textures: using the lapped technique
	void lapped(const hvPictRGB<T> &example, const hvBitmap &pict, double lpower)
	{
		const int MAX_LAPPED=50;
		hvColRGB<T> med; med=example.avg(pict);
		reset(example.sizeX(), example.sizeY(), med);
		hvBitmap mask; mask.erosion(pict,5,5); 
		if (mask.count()<20) return;
		hvBitmap alldots; alldots.erosion(mask,3,3);
		hvBitmap alldotsinv; alldotsinv = alldots; ~alldotsinv;
		hvPict<unsigned char> alpha(alldots, 5, 255);
		hvBitmap dots[MAX_LAPPED];
		int baryx[MAX_LAPPED], baryy[MAX_LAPPED];
		int i,j,k, nlapped=0;
		bool cont = true;
		int avgsize=0;
		printf("extracting features for lapped technique...\n");
		for (k=0; k<MAX_LAPPED && cont; k++)
		{
			dots[k].reset(example.sizeX(), example.sizeY(), false);
			cont=false;
			for (i=0; i<alldots.sizeX() && !cont; i++) 
			for (j=0; j<alldots.sizeY() && !cont; j++) 
				if (alldotsinv.get(i,j)==false) cont=true;
			if (cont)
			{
				printf("feature level %d at %d,%d\n", k, i-1,j-1);
				alldotsinv.seedfill(i-1,j-1,dots[k]);
				if (dots[k].count()<10) { k--; }
				else { dots[k].bary(baryx[k], baryy[k]); nlapped++; avgsize+=dots[k].count(); }
			}
		}
		if (nlapped==0) return;
		avgsize /= nlapped;
		printf("found %d features (avg size %d)\n", nlapped, avgsize);
		printf("doing lapped...\n");
		//blend(example, alpha, (unsigned char)255, 0.25, alldots);
		hvBitmap already(alldots.sizeX(), alldots.sizeY(), false);
		alldotsinv.clear(false);
		cont = true;
		int count=2*alldots.sizeX()*alldots.sizeY()/avgsize;
		printf("setting randomly %d features now\n", count);
		// put randomly features
		while (count>0)
		{
			int px,py, cc=0;
			do { 
				px = (int)((double)rand()/(double)RAND_MAX*(double)alldots.sizeX()); if (px>=alldots.sizeX()) px = alldots.sizeX()-1;
				py = (int)((double)rand()/(double)RAND_MAX*(double)alldots.sizeY()); if (py>=alldots.sizeY()) py = alldots.sizeY()-1;
				cc++;
			} while (alldotsinv.get(px,py) && cc<50);
			already.clear(false);
			k = (int)((double)rand()/(double)RAND_MAX*(double)nlapped);
			if (k<0) k=0; else if (k>=nlapped) k=nlapped-1;
			//printf("set on %d,%d, feature %d\n", px,py,k);
			shiftedblend(px-baryx[k], py-baryy[k], example, alpha, (unsigned char)255, lpower, dots[k], already);
			alldotsinv |= already;
			count--;
			if (count%10==0) { printf("."); fflush(stdout); }
		}
		// fill the rest
		printf("\ncompleting the remaining pixels...\n");
		while (cont)
		{
			cont=false;
			for (i=0; i<alldots.sizeX() && !cont; i++) for (j=0; j<alldots.sizeY() && !cont; j++) if (alldotsinv.get(i,j)==false) cont=true;
			if (cont)
			{
				if (count%10==0) { printf("."); fflush(stdout); }
				already.clear(false);
				k = (int)((double)rand()/(double)RAND_MAX*(double)nlapped);
				//printf("set on %d,%d, feature %d\n", i-1,j-1,k);
				shiftedblend(i-1-baryx[k], j-1-baryy[k], example, alpha, (unsigned char)255, lpower, dots[k], already);
				alldotsinv |= already;
				count++;
			}
		}
		printf("finish\n");
		for (k=0; k<MAX_LAPPED; k++) dots[k].reset();
	}
	*/
	/////////////////////////////////////////////////////////
	// Texture synthesis algorithms
	/////////////////////////////////////////////////////////
		// synthesis of textures: using the lapped technique
	void lapped(const hvPictRGB<T> &example, const hvBitmap &pict, double lpower, hvBitmap *subp=0)
	{
//		std::cout << "doing lapped inpainting...\n"; // (%d, %d), ex(%d, %d), map %d, %d...\n", sizeX(),sizeY(),example.sizeX(), example.sizeY(), pict.sizeX(), pict.sizeY());
		int i,j,k, nlapped=0;
		const int MAX_LAPPED=100;
		//hvColRGB<T> med; med=example.avg(pict);
		//this->hvArray2<hvColRGB<T> >::clear(med);
		//reset(example.sizeX(), example.sizeY(), med);
		hvBitmap mask; mask=pict; //mask.erosion(pict,3,3); 
		if (mask.count()<5) { std::cout<<"only "<< mask.count() <<"pixels in lapped\n"; return; }
		hvBitmap alldots; alldots.erosion(mask,3,3);
		hvBitmap alldotsinv; alldotsinv = alldots; ~alldotsinv;
		hvPict<unsigned char> alpha(alldots, 5, 255);
		hvBitmap dots[MAX_LAPPED];
		int baryx[MAX_LAPPED], baryy[MAX_LAPPED];
		bool cont = true;
		int avgsize=0;
//		std::cout<<"extracting features for lapped technique...\n";
		for (k=0; k<MAX_LAPPED && cont; k++)
		{
			dots[k].reset(example.sizeX(), example.sizeY(), false);
			cont=false;
			for (i=0; i<alldots.sizeX() && !cont; i++) 
			for (j=0; j<alldots.sizeY() && !cont; j++) 
				if (alldotsinv.get(i,j)==false) cont=true;
			if (cont)
			{
				//printf("feature level %d at %d,%d\n", k, i-1,j-1);
//				char buff[256]; fgets(buff, 10, stdin);
				alldotsinv.seedfill(i-1,j-1,dots[k]);
				if (dots[k].count()<4) { k--; }
				else 
				{ 
					dots[k].bary(baryx[k], baryy[k]); 
					int cc = 500;
					int bx = baryx[k], by = baryy[k];
					while (!dots[k].get(baryx[k], baryy[k]) && cc>0)
					{
						int px = bx - 5 + (int)((double)rand() / (double)RAND_MAX*10.0) / 2; if (px < 0) px = 0;  if (px >= dots[k].sizeX()) px = dots[k].sizeX() - 1;
						int py = by - 5 + (int)((double)rand() / (double)RAND_MAX*10.0) / 2; if (py < 0) py = 0;  if (py >= dots[k].sizeY()) py = dots[k].sizeY() - 1;
						baryx[k] = px;  baryy[k] = py;
						cc--;
					}
					if (cc == 0)
					{
						cc = 500;
						while (!dots[k].get(baryx[k], baryy[k]) && cc > 0)
						{
							int px = i - 6 + (int)((double)rand() / (double)RAND_MAX*10.0) / 2; if (px < 0) px = 0;  if (px >= dots[k].sizeX()) px = dots[k].sizeX() - 1;
							int py = j - 6 + (int)((double)rand() / (double)RAND_MAX*10.0) / 2; if (py < 0) py = 0; if (py >= dots[k].sizeY()) py = dots[k].sizeY() - 1;
							baryx[k] = px;  baryy[k] = py;
							cc--;
						}
					}
					if (cc == 0) { baryx[k] = i - 1;  baryy[k] = j - 1; if (dots[k].get(baryx[k], baryy[k])) cc = 1; }
					if (cc == 0) hvFatal("cannot find a point for inpainting");
					nlapped++; 
					avgsize+=dots[k].count(); 
				}
			}
		}
		if (nlapped==0) return;
		avgsize /= nlapped;
		//std::cout<<"found "<<nlapped<<" features (avg size "<<avgsize<<")\n";
		//std::cout<<"now doing lapped...\n";
		//blend(example, alpha, (unsigned char)255, 0.25, alldots);
		hvBitmap already(this->sizeX(), this->sizeY(), false);
		if (subp!=0) alldotsinv=*subp;
		else alldotsinv.reset(this->sizeX(), this->sizeY(), false);
		cont = true;
		int count=2*alldotsinv.sizeX()*alldotsinv.sizeY()/avgsize;
//		std::cout<<"setting randomly "<<count<<" features\n";
		// put randomly features
		while (count>0)
		{
			int px,py, cc=0;
			do { 
				px = (int)((double)rand()/(double)RAND_MAX*(double)alldotsinv.sizeX()); if (px>=alldotsinv.sizeX()) px = alldotsinv.sizeX()-1;
				py = (int)((double)rand()/(double)RAND_MAX*(double)alldotsinv.sizeY()); if (py>=alldotsinv.sizeY()) py = alldotsinv.sizeY()-1;
				cc++;
			} while (alldotsinv.get(px,py) && cc<50);
			already.clear(false);
			k = (int)((double)rand()/(double)RAND_MAX*(double)nlapped);
			if (k<0) k=0; else if (k>=nlapped) k=nlapped-1;
			//printf("set on %d,%d, feature %d\n", px,py,k);
			shiftedblend(px-baryx[k], py-baryy[k], example, alpha, (unsigned char)255, lpower, dots[k], already);
			alldotsinv |= already;
			count--;
			//if (count%10==0) { printf("."); fflush(stdout); }
		}
		// fill the rest
//		std::cout<<"\ncompleting the remaining pixels...\n";
		while (cont)
		{
			cont = false;
			for (i = 0; i<alldotsinv.sizeX() && !cont; i++) for (j = 0; j<alldotsinv.sizeY() && !cont; j++) if (alldotsinv.get(i, j) == false) cont = true;
			if (cont)
			{
				//if (count % 10 == 0) { printf("."); fflush(stdout); }
				already.clear(false);
				k = (int)((double)rand() / (double)RAND_MAX*(double)nlapped);
				if (k >= nlapped) k = nlapped - 1;
				int siftx = (int)((double)rand() / (double)RAND_MAX*(double)5);
				int sifty = (int)((double)rand() / (double)RAND_MAX*(double)5);
				//printf("set on %d,%d, feature %d\n", i-1,j-1,k);
				shiftedblend(i - siftx - baryx[k], j - sifty - baryy[k], example, alpha, (unsigned char)255, lpower, dots[k], already);
				alldotsinv |= already;
				count++;
			}
		}
//		std::cout<<"finish\n";
		for (k=0; k<MAX_LAPPED; k++) dots[k].reset();
	}

	void searchBestBlocPos(const hvPictRGB<T> &example, int bsize,int px, int py, std::vector<hvVec2<int> > &ll, int &posx, int &posy)
	{
		int i,j,k;
		double mindiff=0.0;
		const double PROBA = 0.9;

		posx=0; posy=0;
		for (i=0; i<example.sizeX()-bsize; i++) for (j=0; j<example.sizeY()-bsize;j++)
		{
			if (i==0 && j==0)
			{
				double diff=0.0;
				int nn=0;
				for (k=0; k<ll.size(); k++)
				{
					hvVec2<int> pp = ll.at(k);
					if (i+pp.X()<example.sizeX() && j+pp.Y()<example.sizeY())
					{
						nn++;
						if (px+pp.X()<this->sizeX() && py+pp.Y()<this->sizeY())
						{
							diff+=example.get(i+pp.X(),j+pp.Y()).squaredDifference(this->get(px+pp.X(),py+pp.Y()));
						}
					}
					else hvFatal("example to small for bloc size");
				}
				posx=0; posy=0; mindiff=diff;
			}
			else
			{
				double diff=0.0;
				int nna=0, nnb=0;
				for (k=0; k<ll.size(); k++)
				{
					hvVec2<int> pp = ll.at(k);
					if (i+pp.X()<example.sizeX() && j+pp.Y()<example.sizeY())
					{
						nnb++;
						if (px+pp.X()<this->sizeX() && py+pp.Y()<this->sizeY()) 
						{
							nna++;
							diff+=example.get(i+pp.X(),j+pp.Y()).squaredDifference(this->get(px+pp.X(),py+pp.Y()));
						}
					}
				}
				if (nnb==ll.size() && diff<mindiff && (double)rand() / (double)RAND_MAX<PROBA)
				{
					posx=i; posy=j; mindiff=diff;
				}
			}
		}
	}
	// synthesis of textures: using the chaos blocs technique / quilting
	void chaosblocs(const hvPictRGB<T> &example, int bsize, int blends, bool wblending=true)
	{
		int i,j,bx,by,nx,ny;
		nx = this->sizeX()/(bsize-blends);
		ny = this->sizeY()/(bsize-blends);
		if (nx<1 || ny<1) hvFatal("bloc size too large or example too small in chaosblocs");
		std::vector<hvVec2<int> > lx(bsize*blends); lx.clear();
		for (i=0; i<bsize; i++) for (j=0; j<blends; j++) lx.push_back(hvVec2<int>(i,j));
		std::vector<hvVec2<int> > ly(bsize*blends); ly.clear();
		for (i=0; i<bsize; i++) for (j=0; j<blends; j++) ly.push_back(hvVec2<int>(j,i));
		std::vector<hvVec2<int> > lxy(bsize*blends + (bsize - blends)*blends); lxy.clear();
		for (i=0; i<bsize; i++) for (j=0; j<blends; j++) lxy.push_back(hvVec2<int>(i,j));
		for (i=0; i<blends; i++) for (j=blends; j<bsize; j++) lxy.push_back(hvVec2<int>(i,j));
		for (bx=0; bx<=nx; bx++) for (by=0; by<=ny; by++)
		{
			if (bx==0 && by==0)
			{
				for (i=0; i<bsize; i++) for (j=0; j<bsize; j++)
				{
					if (i<this->sizeX() && j<this->sizeY()) this->update(i,j, example.get(example.sizeX()/2-bsize/2+i,example.sizeY()/2-bsize/2+j));
				}
			}
			else if (bx==0)
			{
				int posx,posy;
				this->searchBestBlocPos(example, bsize, 0, by*(bsize-blends), lx, posx, posy);
				for (i=0; i<bsize; i++) for (j=0; j<bsize; j++)
				{
					if (i<this->sizeX() && j+by*(bsize-blends)<this->sizeY()) 
					{
						hvColRGB<unsigned char> col; 
						double alpha = 0.0;
						if (j<blends) alpha = 1.0-(double)(j)/(double)(blends);
						if (!wblending) alpha = alpha > 0.5 ? 1.0 : 0.0;
						col.blend(this->get(i,j+by*(bsize-blends)),example.get(posx+i,posy+j),alpha);
						this->update(i,j+by*(bsize-blends), col);
					}
				}
			}
			else if (by==0)
			{
				int posx,posy;
				this->searchBestBlocPos(example, bsize,bx*(bsize-blends), 0, ly, posx, posy);
				for (i=0; i<bsize; i++) for (j=0; j<bsize; j++)
				{
					if (i+bx*(bsize-blends)<this->sizeX() && j<this->sizeY()) 
					{
						hvColRGB<unsigned char> col; 
						double alpha = 0.0;
						if (i<blends) alpha = 1.0-(double)(i)/(double)(blends);
						if (!wblending) alpha = alpha > 0.5 ? 1.0 : 0.0;
						col.blend(this->get(i+bx*(bsize-blends),j),example.get(posx+i,posy+j),alpha);
						this->update(i+bx*(bsize-blends),j, col);
					}
				}
			}
			
			else
			{
				int posx,posy;
				this->searchBestBlocPos(example, bsize,bx*(bsize-blends), by*(bsize-blends), lxy, posx, posy);
				for (i=0; i<bsize; i++) for (j=0; j<bsize; j++)
				{
					if (i+bx*(bsize-blends)<this->sizeX() && j+by*(bsize-blends)<this->sizeY()) 
					{
						hvColRGB<unsigned char> col; 
						double alphax = 0.0;
						if (j<blends) alphax = 1.0-(double)(j)/(double)(blends);
						if (!wblending) alphax = alphax > 0.5 ? 1.0 : 0.0;
						double alphay = 0.0;
						if (i<blends) alphay = 1.0-(double)(i)/(double)(blends);
						if (!wblending) alphay = alphay > 0.5 ? 1.0 : 0.0;
						col.blend(this->get(i+bx*(bsize-blends),j+by*(bsize-blends)),example.get(posx+i,posy+j),alphax>alphay?alphax:alphay);
						this->update(i+bx*(bsize-blends),j+by*(bsize-blends), col);
					}
				}
			}
			
		}
	}

	void optimizeNeighborhoods(hvPictRGB<T> &example, int neighbor, int SX=0, int SY=0, int DX=0, int DY=0)
	{
		int i, j, ii, jj;
		hvPictRGB<T> input;
		input.clone(*this, 0, 0, this->sizeX() - 1, this->sizeY() - 1);
		hvArray2<hvVec2<int> > pos(this->sizeX() / neighbor + 1, this->sizeY() / neighbor + 1, hvVec2<int>(0));
		for (i = 0; j < neighbor; i++) for (j = 0; j < neighbor; j++)
		{
			if (i==0 && j==0) for (ii = i; ii < this->sizeX(); ii += neighbor) for (jj = j; jj < this->sizeY(); jj += neighbor)
			{
				int px, py;
				input.bestNeighborMatch(example, ii, jj, neighbor, px, py, SX,SY,DX,DY);
				this->update(ii, jj, example.get(px, py));
				pos.update(ii / neighbor, jj / neighbor, hvVec2<int>(px, py));
			}
			else for (ii = i; ii < this->sizeX(); ii += neighbor) for (jj = j; jj < this->sizeY(); jj += neighbor)
			{
				int px, py;
				int ix[4], iy[4];
				ix[0] = pos.get(ii / neighbor, jj / neighbor).X(); iy[0] = pos.get(ii / neighbor, jj / neighbor).Y();
				ix[1] = pos.get(ii / neighbor+1, jj / neighbor).X(); iy[1] = pos.get(ii / neighbor + 1, jj / neighbor).Y();
				ix[2] = pos.get(ii / neighbor, jj / neighbor+1).X(); iy[2] = pos.get(ii / neighbor, jj / neighbor + 1).Y();
				ix[3] = pos.get(ii / neighbor+1, jj / neighbor+1).X(); iy[3] = pos.get(ii / neighbor + 1, jj / neighbor + 1).Y();
				input.refineBestNeighborMatch(example, ii, jj, neighbor, ix, iy, 4, px, py, SX, SY, DX, DY);
				this->update(ii, jj, example.get(px, py));
			}
		}
	}

	// synthesis of textures: using the chaos blocs technique / quilting
	void chaosblocs2(hvPictRGB<T> &res2, const hvPictRGB<T> &example, const hvPictRGB<T> &example2, int bsize, int blends, bool wblending=true)
	{
		int i,j,bx,by,nx,ny;
		nx = this->sizeX()/(bsize-blends);
		ny = this->sizeY()/(bsize-blends);
		if (nx<1 || ny<1) hvFatal("bloc size too large or example too small in chaosblocs");
		if (blends == 0)
		{
			for (bx = 0; bx <= nx; bx++) for (by = 0; by <= ny; by++)
			{
				int ppx = (int)((double)rand() / (double)RAND_MAX*(double)(example.sizeX()-bsize));
				if (ppx >= example.sizeX() - bsize) ppx = example.sizeX() - bsize - 1;
				int ppy = (int)((double)rand() / (double)RAND_MAX*(double)(example.sizeY() - bsize));
				if (ppy >= example.sizeY() - bsize) ppy = example.sizeY() - bsize - 1;
				for (i = 0; i < bsize; i++) for (j = 0; j < bsize; j++)
				{
					if (bx*bsize+i < this->sizeX() && by*bsize+j < this->sizeY())
					{
						this->update(bx*bsize + i, by*bsize + j, example.get(ppx+i,ppy+j));
						res2.update(bx*bsize + i, by*bsize + j, example2.get(ppx + i, ppy + j));
					}
				}
			}
			return;
		}
		std::vector<hvVec2<int> > lx(bsize*blends); lx.clear();
		for (i=0; i<bsize; i++) for (j=0; j<blends; j++) lx.push_back(hvVec2<int>(i,j));
		std::vector<hvVec2<int> > ly(bsize*blends); ly.clear();
		for (i=0; i<bsize; i++) for (j=0; j<blends; j++) ly.push_back(hvVec2<int>(j,i));
		std::vector<hvVec2<int> > lxy(bsize*blends + (bsize - blends)*blends); lxy.clear();
		for (i=0; i<bsize; i++) for (j=0; j<blends; j++) lxy.push_back(hvVec2<int>(i,j));
		for (i=0; i<blends; i++) for (j=blends; j<bsize; j++) lxy.push_back(hvVec2<int>(i,j));
		for (bx=0; bx<=nx; bx++) for (by=0; by<=ny; by++)
		{
			if (bx==0 && by==0)
			{
				for (i=0; i<bsize; i++) for (j=0; j<bsize; j++)
				{
					if (i<this->sizeX() && j<this->sizeY()) 
					{
							this->update(i,j, example.get(example.sizeX()/2-bsize/2+i,example.sizeY()/2-bsize/2+j));
							res2.update(i,j, example2.get(example2.sizeX()/2-bsize/2+i,example2.sizeY()/2-bsize/2+j));
					}
				}
			}
			else if (bx==0)
			{
				int posx,posy;
				this->searchBestBlocPos(example, bsize, 0, by*(bsize-blends), lx, posx, posy);
				for (i=0; i<bsize; i++) for (j=0; j<bsize; j++)
				{
					if (i<this->sizeX() && j+by*(bsize-blends)<this->sizeY()) 
					{
						hvColRGB<unsigned char> col; 
						double alpha = 0.0;
						if (j<blends) alpha = 1.0-(double)(j)/(double)(blends);
						if (!wblending) alpha = alpha > 0.5 ? 1.0 : 0.0;
						col.blend(this->get(i,j+by*(bsize-blends)),example.get(posx+i,posy+j),alpha);
						this->update(i,j+by*(bsize-blends), col);
						col.blend(res2.get(i,j+by*(bsize-blends)),example2.get(posx+i,posy+j),alpha);
						res2.update(i,j+by*(bsize-blends), col);
					}
				}
			}
			else if (by==0)
			{
				int posx,posy;
				this->searchBestBlocPos(example, bsize,bx*(bsize-blends), 0, ly, posx, posy);
				for (i=0; i<bsize; i++) for (j=0; j<bsize; j++)
				{
					if (i+bx*(bsize-blends)<this->sizeX() && j<this->sizeY()) 
					{
						hvColRGB<unsigned char> col; 
						double alpha = 0.0;
						if (i<blends) alpha = 1.0-(double)(i)/(double)(blends);
						if (!wblending) alpha = alpha > 0.5 ? 1.0 : 0.0;
						col.blend(this->get(i+bx*(bsize-blends),j),example.get(posx+i,posy+j),alpha);
						this->update(i+bx*(bsize-blends),j, col);
						col.blend(res2.get(i+bx*(bsize-blends),j),example2.get(posx+i,posy+j),alpha);
						res2.update(i+bx*(bsize-blends),j, col);
					}
				}
			}
			
			else
			{
				int posx,posy;
				this->searchBestBlocPos(example, bsize,bx*(bsize-blends), by*(bsize-blends), lxy, posx, posy);
				for (i=0; i<bsize; i++) for (j=0; j<bsize; j++)
				{
					if (i+bx*(bsize-blends)<this->sizeX() && j+by*(bsize-blends)<this->sizeY()) 
					{
						hvColRGB<unsigned char> col; 
						double alphax = 0.0;
						if (j<blends) alphax = 1.0-(double)(j)/(double)(blends);
						if (!wblending) alphax = alphax > 0.5 ? 1.0 : 0.0;
						double alphay = 0.0;
						if (i<blends) alphay = 1.0-(double)(i)/(double)(blends);
						if (!wblending) alphay = alphay > 0.5 ? 1.0 : 0.0;
						col.blend(this->get(i+bx*(bsize-blends),j+by*(bsize-blends)),example.get(posx+i,posy+j),alphax>alphay?alphax:alphay);
						this->update(i+bx*(bsize-blends),j+by*(bsize-blends), col);
						col.blend(res2.get(i+bx*(bsize-blends),j+by*(bsize-blends)),example2.get(posx+i,posy+j),alphax>alphay?alphax:alphay);
						res2.update(i+bx*(bsize-blends),j+by*(bsize-blends), col);
					}
				}
			}
			
		}

	}

	void optimizeNeighborhoods2(hvPictRGB<T> &res2, const hvPictRGB<T> &example, const hvPictRGB<T> &example2, int neighbor, int SX = 0, int SY = 0, int DX = 0, int DY = 0)
	{
		int i, j, ii, jj;
		hvPictRGB<T> input;
		input.clone(*this, 0, 0, this->sizeX() - 1, this->sizeY() - 1);
		hvArray2<hvVec2<int> > pos(this->sizeX() / neighbor + 1, this->sizeY() / neighbor + 1, hvVec2<int>(0));
		for (i = 0; i < neighbor; i++) for (j = 0; j < neighbor; j++)
		{
			if (i == 0 && j == 0) for (ii = i; ii < this->sizeX(); ii += neighbor) for (jj = j; jj < this->sizeY(); jj += neighbor)
			{
				int px, py;
				input.bestNeighborMatch(example, ii, jj, neighbor, px, py, SX, SY, DX, DY);
				this->update(ii, jj, example.get(px, py));
				res2.update(ii, jj, example2.get(px, py));
				pos.update(ii / neighbor, jj / neighbor, hvVec2<int>(px, py));
			}
			else for (ii = i; ii < this->sizeX(); ii += neighbor) for (jj = j; jj < this->sizeY(); jj += neighbor)
			{
				int px, py;
				if (ii / neighbor + 1 >= pos.sizeX() || jj / neighbor + 1 >= pos.sizeY())
				{
					input.bestNeighborMatch(example, ii, jj, neighbor, px, py, SX, SY, DX, DY);
				}
				else
				{
					int ix[4], iy[4];
					ix[0] = pos.get(ii / neighbor, jj / neighbor).X(); iy[0] = pos.get(ii / neighbor, jj / neighbor).Y();
					ix[1] = pos.get(ii / neighbor + 1, jj / neighbor).X(); iy[1] = pos.get(ii / neighbor + 1, jj / neighbor).Y();
					ix[2] = pos.get(ii / neighbor, jj / neighbor + 1).X(); iy[2] = pos.get(ii / neighbor, jj / neighbor + 1).Y();
					ix[3] = pos.get(ii / neighbor + 1, jj / neighbor + 1).X(); iy[3] = pos.get(ii / neighbor + 1, jj / neighbor + 1).Y();
					input.refineBestNeighborMatch(example, ii, jj, neighbor, ix, iy, 4, px, py, SX, SY, DX, DY);
				}
				this->update(ii, jj, example.get(px, py));
				res2.update(ii, jj, example2.get(px, py));
			}
		}
	}

	template <class U> void blend(const hvPictRGB<T> &example, const hvPict<U> &alpha, U scal, double power, const hvBitmap &mask)
	{
		int i,j;
		for (i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			if (mask.get(i,j))
			{
				double coeff = pow((double)alpha.get(i,j)/(double)scal, power);
				hvColRGB<T> col = this->get(i,j);
				hvColRGB<T> colex = example.get(i,j);
				hvColRGB<T> colres; colres.blend(colex,col, coeff);
				this->update(i,j,colres);
			}
		}
	}
	template <class U> void blendRect(int px, int py, int x, int y, int sx, int sy, const hvPictRGB<T> &example, const hvPict<U> &alpha, U scal, double power, const hvBitmap &mask, bool mshift=true)
	{
		int i,j;
		for (j=0; j<sy; j++) for (i=0; i<sx; i++) 
		{
			if (mask.get((mshift?x:0)+i,(mshift?y:0)+j))
			{
				if (px+i>=0 && px+i<this->sizeX() && py+j>=0 && py+j<this->sizeY())
				{
					double coeff = pow((double)alpha.get((mshift?x:0)+i,(mshift?y:0)+j)/(double)scal, power);
					hvColRGB<T> col = this->get(px+i,py+j);
					hvColRGB<T> colex = example.get(x+i,y+j);
					hvColRGB<T> colres; colres.blend(colex,col, coeff);
					//if (coeff<1.0) colres=hvColRGB<T>(T(255),T(255),T(0));
					//colres = hvColRGB<T>(0);
					this->update(px+i,py+j,colres);
				}
			}
		}
	}
	template <class U> void shiftedblend(int dx, int dy, const hvPictRGB<T> &example, const hvPict<U> &alpha, U scal, double power, const hvBitmap &mask, hvBitmap &affected)
	{
		int i,j;
		for (i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			int x = i-dx; if (x<0) x += mask.sizeX(); else if (x>=mask.sizeX()) x -= mask.sizeX();
			int y = j-dy; if (y<0) y += mask.sizeY(); else if (y>=mask.sizeY()) y -= mask.sizeY();
			if (mask.get(x,y))
			{
				double coeff = pow((double)alpha.get(x,y)/(double)scal, power);
				hvColRGB<T> col = this->get(i,j);
				hvColRGB<T> colex = example.get(x,y);
				hvColRGB<T> colres; colres.blend(colex,col, coeff);
				this->update(i,j,colres);
				affected.set(i,j,true);
			}
		}
	}
	
	void seedfill(const hvColRGB<T> &col, int x, int y, hvBitmap &bm, hvVec2<int> &min, hvVec2<int> &max) const
	{
		if (x<0 || y<0 || x>= this->sizeX() || y>= this->sizeY()) return;
		if (!(this->get(x,y).equals(col))) return;
		if (bm.get(x,y)) return;
		bm.set(x,y,true);
		min.keepMin(min,hvVec2<int>(x,y));
		max.keepMax(max,hvVec2<int>(x,y));
		int i,a,b;
		for (i=x+1; i<this->sizeX() && (!bm.get(i,y)) && this->get(i,y).equals(col); i++) { bm.set(i,y, true); max.keepMax(max,hvVec2<int>(i,y)); }
		b=i-1;
		for (i=x-1; i>=0 && (!bm.get(i,y)) && this->get(i,y).equals(col); i--) { bm.set(i,y, true); min.keepMin(min,hvVec2<int>(i,y)); }
		a = i+1;
		for (i=a; i<=b; i++) 
		{ 
			seedfill(col,i,y-1,bm,min,max); 
			seedfill(col,i,y+1,bm,min,max); 
		} 
	}
	void seedfill(const hvColRGB<T> &col, int x, int y, hvBitmap &bm, hvVec2<int> &min, hvVec2<int> &max, std::vector<hvVec2<unsigned short> > &lpts) const
	{
		if (x<0 || y<0 || x>= this->sizeX() || y>= this->sizeY()) return;
		if (!(this->get(x,y).equals(col))) return;
		if (bm.get(x,y)) return;
		bm.set(x,y,true);
		lpts.push_back(hvVec2<unsigned short>((unsigned short)x,(unsigned short)y));
		min.keepMin(min,hvVec2<int>(x,y));
		max.keepMax(max,hvVec2<int>(x,y));
		int i,a,b;
		for (i=x+1; i<this->sizeX() && (!bm.get(i,y)) && this->get(i,y).equals(col); i++) { bm.set(i,y, true); max.keepMax(max,hvVec2<int>(i,y)); }
		b=i-1;
		for (i=x-1; i>=0 && (!bm.get(i,y)) && this->get(i,y).equals(col); i--) { bm.set(i,y, true); min.keepMin(min,hvVec2<int>(i,y)); }
		a = i+1;
		for (i=a; i<=b; i++) 
		{ 
			seedfill(col,i,y-1,bm,min,max); 
			seedfill(col,i,y+1,bm,min,max); 
		} 
	}
	void seedfilltorus(const hvColRGB<T> &col, int xx, int yy, hvBitmap &bm, hvVec2<int> &min, hvVec2<int> &max) const
	{
		//if (x<0 || y<0 || x>=sizeX() || y>=sizeY()) return;
		int x = xx; while (x<0) x+= this->sizeX(); while (x>= this->sizeX()) x-= this->sizeX();
		int y = yy; while (y<0) y+= this->sizeY(); while (y>= this->sizeY()) y-= this->sizeY();
		if (!(this->get(x,y).equals(col))) return;
		if (bm.get(x,y)) return;
		bm.set(x,y,true);
		min.keepMin(min,hvVec2<int>(xx,yy));
		max.keepMax(max,hvVec2<int>(xx,yy));
		int i,a,b;
		for (i=x+1; i<this->sizeX() && (!bm.get(i,y)) && this->get(i,y).equals(col); i++) { bm.set(i,y, true); max.keepMax(max,hvVec2<int>(i-x+xx,yy)); }
		b=i-1;
		for (i=x-1; i>=0 && (!bm.get(i,y)) && this->get(i,y).equals(col); i--) { bm.set(i,y, true); min.keepMin(min,hvVec2<int>(i-x+xx,yy)); }
		a = i+1;
		for (i=a; i<=b; i++) 
		{ 
			seedfill(col,i-x+xx,yy-1,bm,min,max); 
			seedfill(col,i-x+xx,yy+1,bm,min,max); 
		} 
	}	
	// creates an index picture "pf" by enumerating all connected components (fragments) in input picture
	// the index corresponds to a fragment number
	// returns the max index (total number of fragments found)
	// the hvVec3<U> contains: frag index, displement (x,y) to frag center
	template <class U> U fragment(hvPict<hvVec3<U> > &pf)
	{
		std::cout<<"fragmenting picture...\n";
		pf.reset(this->sizeX(), this->sizeY(), U(-1));
		hvBitmap yet(this->sizeX(), this->sizeY(), false);
		hvBitmap clas(this->sizeX(), this->sizeY(), false);
		U count = U(0);
		int i,j, x, y;
		bool cont;
		do
		{
			//if (count%10==0) { printf("."); fflush(stdout); }
			cont=true;
			for (i=0; i<this->sizeX() && cont; i++) for (j=0; j<this->sizeY() && cont; j++) if (!yet.get(i,j)) { cont=false; x=i; y=j; }
			if (!cont)
			{
				hvColRGB<T> col = this->get(x,y);
				clas.clear(false);
				hvVec2<int> min(x,y), max(x,y);
				this->seedfill(col, x,y,clas,min,max);
				double meanx=0.0, meany=0.0;
				int c=0;
				for (i=min.X(); i<=max.X(); i++) for (j=min.Y(); j<=max.Y(); j++) if (clas.get(i,j)) 
				{ 
					yet.set(i,j,true); 
					meanx += (double)i; meany+= (double)j;
					c++;
				}
				meanx /= (double)c; meany /= (double)c;
				for (i=min.X(); i<=max.X(); i++) for (j=min.Y(); j<=max.Y(); j++) if (clas.get(i,j)) 
				{
					pf.update(i,j,hvVec3<U>(U(count), U(meanx), U(meany))); 
				}
				count++;
			}
		} while (!cont);
		return count;
	}

	// creates an index picture "pf" by enumerating all connected components (fragments) in input picture
	// the index corresponds to a fragment number
	// returns the max index (total number of fragments found)
	template <class U> U fragment(hvPict<U> &pf, bool periodic=true)
	{
		std::cout<<"fragmenting picture...\n";
		pf.reset(this->sizeX(), this->sizeY(), U(-1));
		hvBitmap yet(this->sizeX(), this->sizeY(), false);
		U count = U(0);
		int i,j, x, y;
		bool cont;
		do
		{
			//if (count%10==0) { printf("."); fflush(stdout); }
			cont=true;
			for (i=0; i<this->sizeX() && cont; i++) for (j=0; j<this->sizeY() && cont; j++) if (!yet.get(i,j)) { cont=false; x=i; y=j; }
			if (!cont)
			{
				hvColRGB<T> col = this->get(x,y);
				hvBitmap clas(this->sizeX(), this->sizeY(), false);
				hvVec2<int> min(x,y), max(x,y);
				this->seedfill(col, x,y,clas,min,max);
				for (i=min.X(); i<=max.X(); i++) for (j=min.Y(); j<=max.Y(); j++) if (clas.get(i,j)) { yet.set(i,j,true); pf.update(i,j,U(count)); }
				if (periodic)
				{
					if (x==0 && y==0) 
					{
							x = this->sizeX()-1; y=0;
							col = this->get(x,y);
							clas.hvBoolArray2::clear(false);
							min=hvVec2<int>(x,y); max=hvVec2<int>(x,y);
							this->seedfill(col, x,y,clas,min,max);
							for (i=min.X(); i<=max.X(); i++) for (j=min.Y(); j<=max.Y(); j++) if (clas.get(i,j)) { yet.set(i,j,true); pf.update(i,j,U(count)); }
							x = 0; y= this->sizeY()-1;
							col = this->get(x,y);
							clas.hvBoolArray2::clear(false);
							min=hvVec2<int>(x,y); max=hvVec2<int>(x,y);
							this->seedfill(col, x,y,clas,min,max);
							for (i=min.X(); i<=max.X(); i++) for (j=min.Y(); j<=max.Y(); j++) if (clas.get(i,j)) { yet.set(i,j,true); pf.update(i,j,U(count)); }
							x = this->sizeX()-1; y= this->sizeY()-1;
							col = this->get(x,y);
							clas.hvBoolArray2::clear(false);
							min=hvVec2<int>(x,y); max=hvVec2<int>(x,y);
							this->seedfill(col, x,y,clas,min,max);
							for (i=min.X(); i<=max.X(); i++) for (j=min.Y(); j<=max.Y(); j++) if (clas.get(i,j)) { yet.set(i,j,true); pf.update(i,j,U(count)); }
					}
					else if (x==0)
					{
							x = this->sizeX()-1; y=0;
							col = this->get(x,y);
							clas.hvBoolArray2::clear(false);
							min=hvVec2<int>(x,y); max=hvVec2<int>(x,y);
							this->seedfill(col, x,y,clas,min,max);
							for (i=min.X(); i<=max.X(); i++) for (j=min.Y(); j<=max.Y(); j++) if (clas.get(i,j)) { yet.set(i,j,true); pf.update(i,j,U(count)); }
					}
					else if (y==0)
					{
							x = 0; y= this->sizeY()-1;
							col = this->get(x,y);
							clas.hvBoolArray2::clear(false);
							min=hvVec2<int>(x,y); max=hvVec2<int>(x,y);
							this->seedfill(col, x,y,clas,min,max);
							for (i=min.X(); i<=max.X(); i++) for (j=min.Y(); j<=max.Y(); j++) if (clas.get(i,j)) { yet.set(i,j,true); pf.update(i,j,U(count)); }
					}
				}
				count++;
			}
		} while (!cont);
		return count;
	}

	// creates a hvVec4 picture by browsing connected components (fragments) in input picture
	// the hvVec4<U> contains: frag index, displement (x,y) to frag center, mean variational value
	template <class U> void fragmentRegular(hvPict<hvVec4<U> > &pf, const hvPictRGB<unsigned char> &pictinput, const hvQPictRGB<unsigned char, 64> &pquant) const;

	// creates a hvVec3 picture by browsing connected components (fragments) in input picture
	// the hvVec3<U> contains: frag index, displement (x,y) to frag center
	template <class U> void fragmentRegular(hvPict<hvVec3<U> > &pf, hvBitmap &border) const
	{
		pf.reset(this->sizeX(), this->sizeY(), hvVec3<U>(U(-1)));
		border.clear(true);
		hvBitmap yet(this->sizeX(), this->sizeY(), false);
		hvBitmap clas(this->sizeX(), this->sizeY(), false);
		U count = U(0);
		int i,j, x, y, ii, jj;
		hvColRGB<T> col;
		hvVec2<int> gmin, gmax, min, max;
		float meanx=0.0, meany=0.0;
		int num=0;
		// features
		bool cont=false;
		do		
		{
			cont=false;
			for (x=0; x<this->sizeX();x++) for (y=0; y<this->sizeY(); y++)
			{
				if (!yet.get(x,y))
				{
					bool onborder=false;
					num=0;
					col = this->get(x,y);
					clas.hvBoolArray2::clear(false);
					min=hvVec2<int>(x,y); max=hvVec2<int>(x,y);
					this->seedfill(col, x,y,clas,min,max);
					int bcount=0;
					for (i=min.X(); i<=max.X(); i++) for (j=min.Y(); j<=max.Y(); j++) if (clas.get(i,j)) 
					{ 
						if (i==0 || i== this->sizeX()-1 || j==0 || j== this->sizeY()-1) bcount++;
					}
					if (bcount>5) onborder=true;
					for (i=min.X(); i<=max.X(); i++) for (j=min.Y(); j<=max.Y(); j++) if (clas.get(i,j)) 
					{ 
						yet.set(i,j,true); 
						if (onborder) border.set(i,j,false);
						pf.update(i,j,hvVec3<U>(U(count))); 
						num++;
					}
					meanx=(float)(max.X()+min.X())*0.5; meany=(float)(max.Y()+min.Y())*0.5;
					for (i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
					{ 
						if (pf.get(i,j).X()==U(count))  
						{ 
							float vx=meanx, vy=meany;
							pf.update(i,j,hvVec3<U>(U(count),U(vx/(float)this->sizeX()*127.0+127.0),U(vy/(float)this->sizeY()*127.0+127.0)));
						}
					}
					count++;
					cont=true;
				}
			}
		} while(cont);
	}

	template <class U> U fragmentWang(hvPict<hvVec4<U> > &pf, const hvPictRGB<unsigned char> &pictinput, const hvPictRGB<unsigned char> &psynth, const hvQPictRGB<unsigned char, 64> &pquant);
	
	// searches for matching picture parts in list of rectangles
	void rematchParticles(const hvPictRGB<unsigned char> &part, const hvPictRGB<unsigned char> &nouv, hvPictRGB<unsigned char> &nouvpart)
	{
		bool cont=true;
		int px=0,py=0,sx,sy, maxy;
		int lx,ly;
		int i,j;
		do {
			cont=true;
			maxy=0;
			// search next particle
			bool found=false;
			sx=4; sy=4; 
			while (!found)
			{
				while(!part.get(px+sx,py).isNull() && px+sx<part.sizeX()) sx++;
				while(!part.get(px,py+sy).isNull() && py+sy<part.sizeY()) sy++;
				found = true;
				for (i=0; i<sx && found; i++) if (!part.get(px+i,py+sy).isNull()) found=false;
				for (i=0; i<sy && found; i++) if (!part.get(px+sx,py+i).isNull()) found=false;
			}
			lx=0; ly=0;
			found=false;
			int startx=0, starty=0;
			while (!found)
			{
				for (i=startx; i<this->sizeX() && !found; i++) for (j=(i==startx?starty:0); j<this->sizeY() && !found; j++)
				{
					if (part.get(px,py).equals(this->get(i,j)) ) found=true;
				}
				if (!found) { printf("inconsistency in particles and example image!\n"); hvFatal("must stop here"); }
				lx=i-1; ly=j-1;
				for (i=0; i<sx && found; i+=2) for (j=0; j<sy && found; j+=2) 
				{
					if (lx+i<this->sizeX() && ly+j<this->sizeY())
						if (!part.get(px+i,py+j).equals(this->get(lx+i,ly+j))) found=false;
				}
				if (!found) { startx=lx; starty=ly+1; if (starty== this->sizeY()) { startx=lx+1; starty=0; } }
			}
			printf("found particle at: %d,%d-%d,%d, matching %d,%d\n", px,py,px+sx-1,py+sy-1,lx,ly);
			for (i=0; i<sx; i++) for (j=0; j<sy; j++)
			{
				nouvpart.update(px+i,py+j,nouv.get((lx+i)% this->sizeX(),(ly+j)% this->sizeY()));
			}
			if (maxy<sy) maxy=sy;
			found=true;
			if (px+sx+1<part.sizeX())
				{
					for (i=0; i<sx && found; i++) if (!part.get(px+sx+1,py).isNull()) found=false;
				}
			if (!found) // start search on same line
			{
				px = px+sx+1;
			}
			else
			{
				cont=false;
				for (i=0; i<part.sizeX() && !cont; i++) if (!part.get(i,py+maxy+1).isNull()) cont=true;
				px=0; py=py+maxy+1;
			}
		} while (cont);
	}
	/*
	hvPictRGB<T>(const hvPictRGB<T> &example, const hvBitmap &pict, int sizex, int sizey, T scal)
	{
		reset(sizex, sizey, hvColRGB<T>(0));
		hvAList<hvVec2<int>,500> blist;
		int iter=2;
		int sx, sy;
		do {
			blist.popAll();
			sx = pict.sizeX()/iter;
			sy = pict.sizeY()/iter;
			pict.extractBoxes(blist,sx,sy);
			printf("%d BOXES of size %d,%d possible\n", blist.length(), sx, sy);
			if (blist.length()<50) { iter = iter+2; }
		} while(blist.length()<50);
		lapped(example, blist, sx, sy, scal);
	}

	void blend(int px, int py, const hvPictRGB<T> &example, int bbx, int bby, int bbsx, int bbsy, const hvPict<T> &mask, hvBitmap &done, T scal)
	{
		int i,j;
		for (i=0; i<bbsx; i++) for (j=0; j<bbsy; j++)
		{
			int mx = (int)((double)i/(double)bbsx*(double)mask.sizeX());
			int my = (int)((double)j/(double)bbsy*(double)mask.sizeY());
			int rx = px-bbsx/2+i; if (rx<0) rx+=sizeX(); else if (rx>=sizeX()) rx-=sizeX();
			int ry = py-bbsy/2+j; if (ry<0) ry+=sizeY(); else if (ry>=sizeY()) ry-=sizeY();
			if (mask.get(mx, my)!=T(0))
			{
				hvColRGB<T> colex = example.get(bbx+i, bby+j);
				hvColRGB<T> colres = get(rx,ry); 
				if (done.get(rx,ry))
				{
					double alpha = (double)mask.get(mx, my)/(double)scal;
					hvColRGB<T> col; col.blend(colex, colres, scal, alpha);
					update(rx,ry,col);
				}
				else 
				{
					update(rx,ry,colex);
				}
				done.set(rx,ry,true);
			}
		}
	}
	template <int nn> void lapped(const hvPictRGB<T> &example, const hvAList<hvVec2<int>,nn> &blist, int sx, int sy, T scal)
	{
		// make a mask for lapped
		hvBitmap bm(128,128, false);
		bm.drawEllipse(64,64, 40, 40, 0.0);
		hvBitmap dist; dist.noiseDistortion(bm, 0.1, 0.1, 20, 20);
		dist.median(4,4);
		hvPict<unsigned char> pi(dist,6, 255);
		pi.gamma(255,0.9);
		// create a bitmap for already visited
		hvBitmap done(sizeX(), sizeY(), false);
		int i,j;
		do
		{
			hvVec2<int> pos;
			done.randomPosition(false, pos);
			int quel=(int)((double)rand()/(double)RAND_MAX*(double)blist.length());
			hvVec2<int> bb = blist.get(quel);
			//printf("particle %d, pos %d,%d, box %d\n", i, pos.X(), pos.Y(), quel);
			blend(pos.X(), pos.Y(), example,  bb.X(), bb.Y(), sx, sy, pi, done, scal);
			printf("done count = %g\n", (double)done.count()/(double)(sizeX()*sizeY()));
		} while(done.count()!=sizeX()*sizeY());
	}
*/

	void makeIcon(int isize, int nback, hvColRGB<T> colback[], FILE *fd, T scal)
	{
		hvBitmap mask(this->sizeX(), this->sizeY(), false);
		int i,j,k, nn=0;
		for (i=0; i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> cc = this->get(i,j);
			for (k=0; k<nback; k++) if (cc.equals(colback[k])) break;
			if (k==nback) {nn++; mask.set(i,j,k==nback); }
		}
		mask.fillholes();
		hvPict<float> pi(mask, 3, 1.0f);
		hvPict<float> picl;
		hvPictRGB<T> pict,pcl;
		pict.clone(*this, 0,0, this->sizeX()-1, this->sizeY()-1);
		while(pi.sizeX()>isize || pi.sizeY()>isize)
		{
			picl.clone(pi); pi.shrink(&picl); 
			pcl.clone(pict, 0,0,pict.sizeX()-1, pict.sizeY()-1); pict.shrink(&pcl);
		}
		for (j=0; j<isize; j++) for (i=0; i<isize; i++) 
		{
			if (i==0) fprintf(fd,"\n");
			hvColRGB<T> cc;
			if (i<pict.sizeX() && j<pict.sizeY()) cc = pict.get(i,j); else cc = hvColRGB<T>();
			float alpha;
			if (i<pict.sizeX() && j<pict.sizeY()) alpha = pi.get(i,j); else alpha = 0.0f;
			fprintf(fd,"%d,%d,%d,%d,",(unsigned char)(cc.RED()*scal), (unsigned char)(cc.GREEN()*scal), (unsigned char)(cc.BLUE()*scal), (unsigned char)(alpha*scal));
		}
	}

	void makeFragFeature(const hvPict<hvVec3<int> > &fragcenters, double radius, int x, int y, hvBitmap &feature) const
	{
		int i;
		feature.clear(false);
		for (i=0; i<this->sizeX(); i++)
		{
			for (int j=0; j<this->sizeY(); j++)
			{
				hvVec3<int> col = fragcenters.get(i, j);
				int px = (int)((double)(col.Y()-127)/127.0*(double)this->sizeX());
				int py = (int)((double)(col.Z()-127)/127.0*(double)this->sizeY());
				double dx = (double)px/(double)this->sizeX()-(double)x/(double)this->sizeX();
				double dy = (double)py/(double)this->sizeY()-(double)y/(double)this->sizeY();
				double dist = sqrt(dx*dx+dy*dy);
				feature.set(i,j,dist<radius?true:false);
			}
		}
	}

	hvVec2<int> chooseMinSquareDiff(int bx, int by, int bsx, int bsy, int x,int y,int sx,int sy, const hvPictRGB<T> &inpict, hvColRGB<double> &minerr) 
	{
		int i,j;
		hvColRGB<double> res;
		double minv=0.0;
		hvVec2<int> minpos(bx,by);
		bool first=true;

		for (i=bx; i<bx+bsx; i+=2) for (j=by; j<by+bsy; j+=2)
		{
				//printf("mask at %d,%d, %s\n", cx-i+x,cy-j+y, mask.get(cx-i+x,cy-j+y)?"true":"flase");
				this->squareDifference(255.0,i,j,x,y,sx,sy,inpict,res,4);
				double vv = res.norm();
				if (first) { first=false; minv=vv; minpos=hvVec2<int>(i,j); minerr=res; }
				else if (vv<minv) { minv=vv; minpos=hvVec2<int>(i,j); minerr=res; }
		}
		return minpos;
	}
	hvVec2<int> chooseMinSquareDiffBorder(int bx, int by, int bsx, int bsy, int x, int y, int sx, int sy, const hvPictRGB<T> &inpict, hvColRGB<double> &minerr)
	{
		int i, j;
		hvColRGB<double> res;
		double minv = 0.0;
		hvVec2<int> minpos(bx, by);
		bool first = true;

		for (i = bx; i<bx + bsx; i += 2) for (j = by; j<by + bsy; j += 2)
		{
			//printf("mask at %d,%d, %s\n", cx-i+x,cy-j+y, mask.get(cx-i+x,cy-j+y)?"true":"flase");
			this->squareDifferenceBorder(255.0, i, j, x, y, sx, sy, inpict, res, 4);
			double vv = res.norm();
			if (first) { first = false; minv = vv; minpos = hvVec2<int>(i, j); minerr = res; }
			else if (vv<minv) { minv = vv; minpos = hvVec2<int>(i, j); minerr = res; }
		}
		return minpos;
	}

	hvVec2<int> chooseMinSquareDiff(int bx, int by, int bsx, int bsy, int x,int y,int sx,int sy, const hvPictRGB<T> &inpict, const hvBitmap &mask, int cx, int cy, hvColRGB<double> &minerr) 
	{
		int i,j;
		hvColRGB<double> res;
		bool first=true;
		double minv=0.0;
		hvVec2<int> minpos(bx,by);

		for (i=bx; i<bx+bsx; i+=2) for (j=by; j<by+bsy; j+=2)
		{
			//if (mask.get(cx-i+x,cy-j+y) && mask.get(cx-i+x-1,cy-j+y) && mask.get(cx-i+x+1,cy-j+y) && mask.get(cx-i+x,cy-j+y-1) && mask.get(cx-i+x,cy-j+y+1))
			if (mask.get(cx-i+x,cy-j+y)) 
			{
				//printf("mask at %d,%d, %s\n", cx-i+x,cy-j+y, mask.get(cx-i+x,cy-j+y)?"true":"flase");
				this->squareDifference(255.0,i,j,x,y,sx,sy,inpict,mask,res);
				double vv = res.norm();
				if (first) { first=false; minv=vv; minpos=hvVec2<int>(i,j); minerr=res; }
				else if (vv<minv) { minv=vv; minpos=hvVec2<int>(i,j); minerr=res; }
			}
		}
		if (first) { printf("warning cannot find best minsquare diff on %d,%d\n", cx,cy); }
		return minpos;
	}
	hvVec2<int> chooseWeightedMinSquareDiff(int bx, int by, int bsx, int bsy, int x,int y,int sx,int sy, const hvPictRGB<T> &inpict,  const hvBitmap &mask, const hvPict<double> &weight, int cx, int cy, hvColRGB<double> &minerr) 
	{
		int i,j;
		hvColRGB<double> res;
		bool first=true;
		double minv=0.0;
		hvVec2<int> minpos(bx,by);

		for (i=bx; i<bx+bsx; i+=2) for (j=by; j<by+bsy; j+=2)
		{
			//if (mask.get(cx-i+x,cy-j+y) && mask.get(cx-i+x-1,cy-j+y) && mask.get(cx-i+x+1,cy-j+y) && mask.get(cx-i+x,cy-j+y-1) && mask.get(cx-i+x,cy-j+y+1))
			if (mask.get(cx-i+x,cy-j+y))
			{
				//printf("mask at %d,%d, %s\n", cx-i+x,cy-j+y, mask.get(cx-i+x,cy-j+y)?"true":"flase");
				this->weightedSquareDifference(255.0,i,j,x,y,sx,sy,inpict,weight,res);
				double vv = res.norm();
				if (first) { first=false; minv=vv; minpos=hvVec2<int>(i,j); minerr=res; }
				else if (vv<minv) { minv=vv; minpos=hvVec2<int>(i,j); minerr=res; }
			}
		}
		if (first) { printf("warning cannot find best minsquare diff on %d,%d\n", cx,cy); }
		return minpos;
	}
	// Making Wang tiles using fragmentation
	void makeWangtiles16(double radius, double pow, const hvPictRGB<T> &inpict,  const hvPictRGB<T> &frag, const hvPict<unsigned char> &heightfield, hvPictRGB<T> &wangfrag, hvPict<unsigned char> &wangheight)
	{
		int i,j,k;
		int pixrad = (int)((1.0-radius)*(double)(inpict.sizeX()<inpict.sizeY()?inpict.sizeX():inpict.sizeY()));
		reset(inpict.sizeX()*4, inpict.sizeY()*4, hvColRGB<T>(0));
		wangfrag.reset(inpict.sizeX()*4, inpict.sizeY()*4, hvColRGB<T>(0));
		wangheight.reset(inpict.sizeX()*4, inpict.sizeY()*4, 0);
		hvBitmap wangyet(inpict.sizeX()*4, inpict.sizeY()*4, false);
		for (i=0; i<inpict.sizeX()*4;i++) for (j=0; j<inpict.sizeY()*4; j++) 
		{
			update(i,j,inpict.get(i%inpict.sizeX(),j%inpict.sizeY()));
			wangfrag.update(i,j,frag.get(i%inpict.sizeX(),j%inpict.sizeY()));
			wangheight.update(i,j,heightfield.get(i%inpict.sizeX(),j%inpict.sizeY()));
		}
		hvPict<hvVec3<int> > fragcenters;
		hvBitmap featborder(frag.sizeX(), frag.sizeY(), true);
		hvBitmap feature(frag.sizeX(), frag.sizeY(), false);
		hvBitmap fragfeature(frag.sizeX(), frag.sizeY(), false);
		frag.fragmentRegular(fragcenters, featborder);
		hvVec2<int> fmin, fmax;
		
		// fill from corner
		bool search=true, bordersdone=false;
		int tilex=0, tiley=0, tilecounter=0;
		int counter=0;
		do {
		//for (k=0; k<3; k++)
		//{
			int modulox=0; int moduloy=0;
			int locx=inpict.sizeX()-1,locy=inpict.sizeY()-1;
			search=true;

			// first manage to fill the corner and the borders
			if (!bordersdone)
			{
				
				locx=inpict.sizeX()-1; locy=inpict.sizeY()-1;
				while (wangyet.get(locx,locy) && locx>0) locx--;
				if (locx>0) search=false;
				
				if (search) 
				{ 
					locx=inpict.sizeX()-1; locy=inpict.sizeY()-1; 
					while (wangyet.get(locx,locy) && locy>0) locy--;
					if (locy>0) search=false;
				}
				
				if (search) 
				{ 
					locx=inpict.sizeX()-1; locy=2*inpict.sizeY()-1;
					while (wangyet.get(locx,locy) && locx>0) locx--;
					if (locx>0) { search=false; moduloy=1; }
				}
				
				if (search) 
				{ 
					locx=2*inpict.sizeX()-1; locy=inpict.sizeY()/2;
					while (wangyet.get(locx,locy) && locy>0) locy--;
					if (locy>0) { search=false; modulox=1; }
				}
				
				if (search) { bordersdone=true; printf("borders are done...\n"); }
			}
			// manage the interior of the tiles
			
			if (bordersdone)
			{
				if (tilex<4 && tiley<4)
				{
					int cc=0; 
					do {
						cc++;
						locx=tilex*inpict.sizeX()+inpict.sizeX()/2+int((2.0*(double)rand()/(double)RAND_MAX-1.0)*(double)(inpict.sizeX())/4.0);
						locy=tiley*inpict.sizeY()+inpict.sizeY()/2+int((2.0*(double)rand()/(double)RAND_MAX-1.0)*(double)(inpict.sizeY())/4.0);
						modulox=tilex; moduloy=tiley;
					} while(cc<20 && wangyet.get(locx,locy));
					search=false;
				}
				tilecounter++; if (tilecounter==2) { tilecounter=0; tilex++; if (tilex>=4) { tilex=0; tiley++; } }
			}
			
			if (!search)
			{
				printf("search to fill position: %d,%d\n", locx,locy);
				int kk;
				hvVec2<int> minpos;
				for (kk=0; kk<1; kk++)
				{
					int dimx,dimy,posx,posy;
					frag.makeFragFeature(fragcenters, radius, inpict.sizeX()/2+int((2.0*(double)rand()/(double)RAND_MAX-1.0)*(double)(pixrad/2)), inpict.sizeY()/2+int((2.0*(double)rand()/(double)RAND_MAX-1.0)*(double)(pixrad/2)), feature);
					feature &= featborder;
					fragfeature = feature;
					feature.dilatation(3,3);
					feature.fillholes();
					feature.box(fmin,fmax);
					dimx = (fmax.X()-fmin.X());
					dimy = (fmax.Y()-fmin.Y());
					posx=-dimx/2;
					posy=-dimy/2;
					// choose position to minimize squared difference
					minpos = chooseMinSquareDiff(locx-dimx, locy-dimy, dimx, dimy, fmin.X(),fmin.Y(),dimx,dimy, inpict, feature, locx, locy); 
				}				
				//printf("Feature: %d,%d - %d,%d (%d,%d)\n", fmin.X(), fmin.Y(), fmax.X(), fmax.Y(),dimx, dimy);
				counter++;

				// blend the patch
				hvBitmap ff; ff=feature;
				//ff.dilatation(6,6);
				hvPict<unsigned char> pfeat(ff, 5, 255);
				bool oncorner=false;
				bool onxa=false, onxb=false, onya=false, onyb=false;
				for (i=0; i<fmax.X()-fmin.X(); i++) for (j=0; j<fmax.Y()-fmin.Y(); j++)
				{
					if (ff.get(fmin.X()+i, fmin.Y()+j))
					{
						int vx=minpos.X()+i;  while (vx<0) vx+=inpict.sizeX();  vx%=inpict.sizeX();
						int vy=minpos.Y()+j;  while (vy<0) vy+=inpict.sizeY();  vy%=inpict.sizeY();
						if ( (vx==0 || vx==inpict.sizeX()-1) && (vy==0 || vy==inpict.sizeX()-1) ) oncorner=true;
						vx=minpos.X()+i; if (vx<0) vx+=4*inpict.sizeX();
						if (vx==0 || vx==inpict.sizeX()-1 || vx==inpict.sizeX() || vx==4*inpict.sizeX()-1) onxa=true;
						if (vx==2*inpict.sizeX()-1 || vx==2*inpict.sizeX() || vx==3*inpict.sizeX()-1 || vx==3*inpict.sizeX()) onxb=true;
						vy=minpos.Y()+j; if (vy<0) vy+=4*inpict.sizeY();
						if (vy==0 || vy==inpict.sizeY()-1 || vy==inpict.sizeY() || vy==4*inpict.sizeY()-1) onya=true;
						if (vy==2*inpict.sizeY()-1 || vy==2*inpict.sizeY() || vy==3*inpict.sizeY()-1 || vy==3*inpict.sizeY()) onyb=true;
					}
				}
				printf("feature pos %d,%d  (%d,%d) ", minpos.X(), minpos.Y(), modulox, moduloy);
				if (oncorner) printf("on corner ");
				else
				{
					if (onxa) printf("on XA "); if (onxb) printf("on XB ");
					if (onya) printf("on YA "); if (onyb) printf("on YB ");
				}
				printf("\n");
				if (!bordersdone && (oncorner || onxa || onxb || onya || onyb) )
				{
					minpos=hvVec2<int>(minpos.X()-modulox*inpict.sizeX(),minpos.Y()-moduloy*inpict.sizeY());
					for (i=0; i<=4; i++) for (j=0; j<=4; j++)
					{
						bool doblend=false;
						if (oncorner) doblend=true;
						else
						{
							if (onxa && !onya && !onyb && (i==0 || i==1 || i==4)) doblend=true;
							if (onxb && !onyb && !onya && (i==2 || i==3)) doblend=true;
							if (onya && !onxa && !onxb && (j==0 || j==1 || j==4)) doblend=true;
							if (onyb && !onxb && !onxa && (j==2 || j==3)) doblend=true;
						}
						if (doblend)
						{
							blendRect(inpict.sizeX()*(i-1)+minpos.X(),inpict.sizeY()*(j-1)+minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),inpict,pfeat, (unsigned char)(255), pow, ff);
							wangfrag.copyWangRect(inpict.sizeX()*(i-1)+minpos.X(),inpict.sizeY()*(j-1)+minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),frag,fragfeature);
							wangheight.blendRect(inpict.sizeX()*(i-1)+minpos.X(),inpict.sizeY()*(j-1)+minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),heightfield,pfeat, (unsigned char)(255), pow, ff);
							wangyet.operatorOr(inpict.sizeX()*(i-1)+minpos.X(),inpict.sizeY()*(j-1)+minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(), fragfeature);
						}
					}
				}
				else if (bordersdone && !(oncorner || onxa || onxb || onya || onyb) )
				{
					
							blendRect(minpos.X(),minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),inpict,pfeat, (unsigned char)(255), pow, ff);
							wangfrag.copyWangRect(minpos.X(),minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),frag,fragfeature);
							wangheight.blendRect(minpos.X(),minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),heightfield,pfeat, (unsigned char)(255), pow, ff);
							wangyet.operatorOr(minpos.X(),minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(), fragfeature);							
				}
				else printf("feature eliminated.\n");
			}
		} while (!search);	

		// do ax edge
	}

};

	


////////////////////////////////////////////////////////////
template <class T> class hvPictRGBA : public hvField2< hvColRGBA<T> >  
////////////////////////////////////////////////////////////
{
public:
	hvPictRGBA<T>() : hvField2< hvColRGBA<T> >() { }
	hvPictRGBA<T>(int sx, int sy, const hvColRGBA<T> &nil) : hvField2< hvColRGBA<T> >(sx, sy, nil),hvArray2< hvColRGBA<T> >(sx, sy, nil) { }
	
	void clone(const hvPictRGBA<T> &pict,int x, int y, int sx, int sy)
	{
		hvField2< hvColRGBA<T> >::reset(sx-x+1, sy-y+1, hvColRGBA<T>(0));
		int i,j;
		for (i=x; i<=sx; i++) for (j=y; j<=sy; j++)
		{
			this->update(i-x,j-y,pict.get(i,j));
		}
	}
	void clone(const hvPictRGB<T> &pict,const hvPict<T> &pa, int x, int y, int sx, int sy)
	{
		hvField2< hvColRGBA<T> >::reset(sx-x+1, sy-y+1, hvColRGBA<T>(0));
		int i,j;
		for (i=x; i<=sx; i++) for (j=y; j<=sy; j++)
		{
			hvColRGB<T> col = pict.get(i,j);
			hvColRGBA<T> cc(col, pa.get(i,j));
			//printf("clone %d,%d -> %d,%d,%d,%d\n", i,j,(int)cc.RED(), (int)cc.GREEN(), (int)cc.BLUE(), (int)cc.ALPHA());
			this->update(i-x,j-y,cc);
		}
	}
	void copy(int x, int y, const hvPictRGBA<T> &pict)
	{
		int i,j;
		for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++)
		{
			this->update(x+i,y+j,pict.get(i,j));
		}
	}
	void copyRect(int px, int py, int x, int y, int sx, int sy, const hvPictRGBA<T> &pict, const hvBitmap &mask)
	{
		int i,j;
		for (i=0; i<sx; i++) for (j=0; j<sy; j++)
		{
			if (mask.get(x+i, y+j))
			{
				this->update(px+i,py+j,pict.get(x+i,y+j));
			}
		}
	}

	void gamma(T scal, double power)
	{
		int i,j;
		for (i=0;i<this->sizeX(); i++) for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<T> v = this->get(i,j);
			v.gamma(scal, power);
			this->update(i,j,v);
		}
	}
	void convert(const hvBitmap &pict, const hvColRGBA<T> &va, const hvColRGBA<T> &vb)
	{
		hvField2< hvColRGBA<T> >::reset(pict.sizeX(), pict.sizeY(), hvColRGBA<T>(0));
		//hvArray2< hvColRGB<T> >::reset(pict.sizeX(), pict.sizeY(), hvColRGB<T>(0));
		int i,j;
		for (i=0; i<pict.sizeX(); i++) for (j=0; j<pict.sizeY(); j++)
		{
			if (pict.get(i,j)) this->update(i,j,va); else this->update(i,j,vb);
		}
	}
	void loadPPM(FILE *fd, T norm, T alpha)
	{
		int  sx, sy;
		int i,j, type;
		char buff[256];
		hvColRGB<T> co;

		hvPictRGB<T>::readPPMLine(fd,buff);
		if (strcmp(buff,"P6\n")==0) type = 0;
		else if (strcmp(buff,"P3\n")==0) type = 1;
		else { type = 2; printf("unknown picture PPM type=%d (%s)\n", type,buff); } 
		hvPictRGB<T>::readPPMLine(fd,buff);
		sscanf(buff,"%d %d",&sx,&sy);
		hvPictRGB<T>::readPPMLine(fd,buff);
		if (strcmp(buff,"255\n")!=0){ printf("type=%d\n", type); hvFatal("Not the right PPM Format"); }
		this->reset(sx, sy, hvColRGBA<T>());
		for (i=0; i<sy; i++)
		for (j=0; j<sx; j++)
			{
				unsigned char r,g,b;
				if (type==0)
				{
					fread(&r,1,sizeof(unsigned char),fd);
					fread(&g,1,sizeof(unsigned char),fd);
					fread(&b,1,sizeof(unsigned char),fd);
				}
				else if (type==1)
				{
					int rr, gg, bb;
					fscanf(fd, "%d %d %d", &rr, &gg, &bb);
					r= (unsigned char)rr;
					g= (unsigned char)gg;
					b= (unsigned char)bb;
				}
				else { r=0; g=0; b=0; }
				hvArray2< hvColRGBA<T> >::update(j,sy-i-1,hvColRGBA<T>((T)r/norm, (T)g/norm, (T)b/norm, alpha));
			}
	}
	void loadPPMA(FILE *fd, T norm)
	{
		int  sx, sy;
		int i,j, type;
		char buff[256];
		hvColRGB<T> co;

		hvPictRGB<T>::readPPMLine(fd,buff);
		if (strcmp(buff,"P6\n")==0) type = 0;
		else if (strcmp(buff,"P3\n")==0) type = 1;
		else { type = 2; printf("unknown picture PPM type=%d (%s)\n", type,buff); } 
		hvPictRGB<T>::readPPMLine(fd,buff);
		sscanf(buff,"%d %d",&sx,&sy);
		hvPictRGB<T>::readPPMLine(fd,buff);
		if (strcmp(buff,"255\n")!=0){ printf("type=%d\n", type); hvFatal("Not the right PPM Format"); }
		this->reset(sx, sy, hvColRGBA<T>());
		for (i=0; i<sy; i++)
		for (j=0; j<sx; j++)
			{
				unsigned char r,g,b,aa;
				if (type==0)
				{
					fread(&r,1,sizeof(unsigned char),fd);
					fread(&g,1,sizeof(unsigned char),fd);
					fread(&b,1,sizeof(unsigned char),fd);
					fread(&aa,1,sizeof(unsigned char),fd);
				}
				else if (type==1)
				{
					int rr, gg, bb, alpha;
					fscanf(fd, "%d %d %d %d", &rr, &gg, &bb, &alpha);
					r= (unsigned char)rr;
					g= (unsigned char)gg;
					b= (unsigned char)bb;
					aa = (unsigned char)alpha;
				}
				else { r=0; g=0; b=0; aa=0; }
				hvArray2< hvColRGBA<T> >::update(j,sy-i-1,hvColRGBA<T>((T)r/norm, (T)g/norm, (T)b/norm, (T)aa/norm));
			}
	}
	void savePPM(FILE *fd, T norm, bool alpha=false)
	{
		int i,j;
		hvColRGBA<T> co;
		unsigned char v;

		fprintf(fd,"P6\n");
		fprintf(fd,"%d %d\n", this->sizeX(), this->sizeY());
		fprintf(fd,"255\n");
		for (i=0; i<this->sizeY(); i++)
		for (j=0; j<this->sizeX(); j++)
			{
			co = hvArray2< hvColRGBA<T> >::get(j, this->sizeY()-i-1);
			v = (unsigned char)((T)(alpha?co.ALPHA():co.RED())*norm);
			fwrite(&v,1,sizeof(unsigned char),fd);
			v = (unsigned char)((T)(alpha?co.ALPHA():co.GREEN())*norm);
			fwrite(&v,1,sizeof(unsigned char),fd);
			v = (unsigned char)((T)(alpha?co.ALPHA():co.BLUE())*norm);
			fwrite(&v,1,sizeof(unsigned char),fd);
			}
	}
	void savePPMA(FILE *fd, T norm)
	{
		int i,j;
		hvColRGBA<T> co;
		unsigned char v;

		fprintf(fd,"P6\n");
		fprintf(fd,"%d %d\n", this->sizeX(), this->sizeY());
		fprintf(fd,"255\n");
		for (i=0; i<this->sizeY(); i++)
		for (j=0; j<this->sizeX(); j++)
			{
			co = hvArray2< hvColRGBA<T> >::get(j, this->sizeY()-i-1);
			v = (unsigned char)((T)co.RED()*norm);
			fwrite(&v,1,sizeof(unsigned char),fd);
			v = (unsigned char)((T)co.GREEN()*norm);
			fwrite(&v,1,sizeof(unsigned char),fd);
			v = (unsigned char)((T)co.BLUE()*norm);
			fwrite(&v,1,sizeof(unsigned char),fd);
			v = (unsigned char)((T)co.ALPHA()*norm);
			fwrite(&v,1,sizeof(unsigned char),fd);			}
	}
};




}

#endif // !efined(AFX_PICTRGB_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_)
