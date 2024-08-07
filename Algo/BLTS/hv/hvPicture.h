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



// hvPicture.h: interface for the picture class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PICTURE_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_)
#define AFX_PICTURE_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "hvBitmap.h"
#include "hvPict.h"
#include "hvPictRGB.h"
#include "hvQPicture.h"

namespace hview {

////////////////////////////////////////////////////////////
// template <class T> class hvPict : public hvField2< T > 


	template <class T> template <class X, class Y>  hvPict<T>::hvPict(const hvPictRGB<X> &pict, Y scal, int x, int y, int sx, int sy) : hvField2< T >(sx - x + 1, sy - y + 1, T(0)), hvArray2< T >(sx - x + 1, sy - y + 1, T(0))
	{
		this->reset(pict.sizeX(), pict.sizeY(), T(0));
		int i, j;
		for (i = x; i <= sx; i++) for (j = y; j <= sy; j++)
		{
			this->update(i - x, j - y, T(Y((pict.get(i, j)).luminance())*scal));
		}
	}
	template <class T> template <class X, class Y> void hvPict<T>::convert(const hvPictRGB<X> &pict, hvPictComponent cc, Y scal)
	{
		this->reset(pict.sizeX(), pict.sizeY(), T(0));
		int i, j;
		for (i = 0; i<pict.sizeX(); i++) for (j = 0; j<pict.sizeY(); j++)
		{
			Y v;
			hvColRGB<X> co = pict.get(i, j);
			switch (cc)
			{
			case HV_RED: v = Y(co.RED()); break;
			case HV_GREEN: v = Y(co.GREEN()); break;
			case HV_BLUE: v = Y(co.BLUE()); break;
			default: v = Y(co.luminance());
			}
			this->update(i, j, T(v*scal));
		}
	}
	template <class T> template <class X, class Y> hvPict<T>::hvPict(const hvPictRGB<X> &pict, hvPictComponent cc, Y scal, int x, int y, int sx, int sy) : hvField2< T >(sx - x + 1, sy - y + 1, T(0)), hvArray2< T >(sx - x + 1, sy - y + 1, T(0))
	{
		int i, j;
		for (i = x; i <= this->sizeX(); i++) for (j = y; j <= this->sizeY(); j++)
		{
			Y v;
			hvColRGB<X> co = pict.get(i, j);
			switch (cc)
			{
			case HV_RED: v = Y(co.RED()); break;
			case HV_GREEN: v = Y(co.GREEN()); break;
			case HV_BLUE: v = Y(co.BLUE()); break;
			default: v = Y(co.luminance());
			}
			this->update(i - x, j - y, T(v*scal));
			//printf("%d,%d=%d\n",i,j,get(i-x,j-y));
		}
	}
	template <class T> template <class X, class Y> hvPict<T>::hvPict(const hvPictRGBA<X> &pict, Y scal, int x, int y, int sx, int sy) : hvField2< T >(sx - x + 1, sy - y + 1, T(0)), hvArray2< T >(sx - x + 1, sy - y + 1, T(0))
	{
		int i, j;
		for (i = x; i <= this->sizeX(); i++) for (j = y; j <= this->sizeY(); j++)
		{
			this->update(i - x, j - y, T(Y((pict.get(i, j)).luminance())*scal));
		}
	}
	template <class T> template <class X, class Y> hvPict<T>::hvPict(const hvPictRGBA<X> &pict, hvPictComponent cc, Y scal, int x, int y, int sx, int sy) : hvField2< T >(sx - x + 1, sy - y + 1, T(0)), hvArray2< T >(sx - x + 1, sy - y + 1, T(0))
	{
		int i, j;
		for (i = x; i <= this->sizeX(); i++) for (j = y; j <= this->sizeY(); j++)
		{
			Y v;
			hvColRGBA<X> co = pict.get(i, j);
			switch (cc)
			{
			case HV_RED: v = Y(co.RED()); break;
			case HV_GREEN: v = Y(co.GREEN()); break;
			case HV_BLUE: v = Y(co.BLUE()); break;
			case HV_ALPHA: v = Y(co.ALPHA()); break;
			default: v = Y(co.luminance());
			}
			this->update(i - x, j - y, T(v*scal));
		}
	}

	template <class T> void hvPict<T>::squaredDifference(int px, int py, int dx, int dy, const hvPictRGB<unsigned char> &pia, int ix, int iy, const hvPictRGB<unsigned char> &pib)
	{
		int i, j;
		this->reset(dx, dy, T(0));
		for (j = 0; j<dy; j++) for (i = 0; i<dx; i++) 
		{
			int kx = px + i; while (kx<0) kx += pia.sizeX(); while (kx >= pia.sizeX()) kx -= pia.sizeX();
			int ky = py + j; while (ky<0) ky += pia.sizeY(); while (ky >= pia.sizeY()) ky -= pia.sizeY();
			this->update(i, j, T(pia.get(kx, ky).squaredDifference(pib.get(ix + i, iy + j))));
		}
	}

////////////////////////////////////////////////////////////
//template <class T> class hvPictRGB : public hvField2< hvColRGB<T> >

	template <class T>  template <class U, class V> hvPictRGB<T>::hvPictRGB(const hvPict<U> &p, V scal)
	{
		int i, j;
		V val;
		hvArray2< hvColRGB<T> >::reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		for (i = 0; i<p.sizeX(); i++) for (j = 0; j<p.sizeY(); j++)
		{
			val = V(p.get(i, j));
			val *= scal;
			this->update(i, j, hvColRGB<T>(T(val), T(val), T(val)));
		}
	}
	template <class T> template <class U, class V> hvPictRGB<T>::hvPictRGB(const hvPict<U> &p, V scal, V shift)
	{
		int i, j;
		V val;
		hvArray2< hvColRGB<T> >::reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		for (i = 0; i<p.sizeX(); i++) for (j = 0; j<p.sizeY(); j++)
		{
			val = V(p.get(i, j));
			val *= scal;
			val += shift;
			update(i, j, hvColRGB<T>(T(val), T(val), T(val)));
		}
	}
	template <class T> template <class U, class V> hvPictRGB<T>::hvPictRGB(const hvPict<U> &p, V scal, int x, int y, int sx, int sy)
	{
		int i, j;
		V val;
		hvArray2< hvColRGB<T> >::reset(sx - x + 1, sy - y + 1, hvColRGB<T>(0));
		for (i = x; i <= sx; i++) for (j = y; j <= sy; j++)
		{
			val = V(p.get(i, j));
			val *= scal;
			this->update(i - x, j - y, hvColRGB<T>(T(val), T(val), T(val)));
		}
	}
	template <class T> template <class U, class V> void hvPictRGB<T>::convert(const hvPict<U> &p, V scal, int x, int y, int sx, int sy)
	{
		int i, j;
		V val;
		this->reset(sx - x + 1, sy - y + 1, hvColRGB<T>(T(0)));
		for (i = x; i <= sx; i++) for (j = y; j <= sy; j++)
		{
			val = V(p.get(i, j));
			val *= scal;
			this->update(i - x, j - y, hvColRGB<T>(T(val), T(val), T(val)));
		}
	}
	template <class T> template <class U, class V> void hvPictRGB<T>::convertloga(const hvPict<U> &p, V loga, V max, V scal, int x, int y, int sx, int sy)
	{
		int i, j;
		V val;
		this->reset(sx - x + 1, sy - y + 1, hvColRGB<T>(T(0)));
		for (i = x; i <= sx; i++) for (j = y; j <= sy; j++)
		{
			val = V(p.get(i, j));
			val = (V)(log(1.0 + (double)loga*(double)val / (double)max) / log((double)loga + 1.0));
			if (val>V(1)) val = V(1);
			val *= scal;
			this->update(i - x, j - y, hvColRGB<T>(T(val), T(val), T(val)));
		}
	}
	template <class T> template <class U, class V> hvPictRGB<T>::hvPictRGB(const hvPict<U> &p, V scal, V min, V max)
	{
		int i, j;
		V val;
		hvArray2< hvColRGB<T> >::reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		for (i = 0; i<p.sizeX(); i++) for (j = 0; j<p.sizeY(); j++)
		{
			val = V(p.get(i, j));
			val = scal*(val - min) / (max - min);
			this->update(i, j, hvColRGB<T>(T(val), T(val), T(val)));
		}
	}
	template <class T> template <class U> hvPictRGB<T>::hvPictRGB(const hvPict<U> &p, const std::vector<hvColRGB<unsigned char> > &lcol)
	{
		int i, j;
		hvArray2< hvColRGB<T> >::reset(p.sizeX(), p.sizeY(), hvColRGB<T>(0));
		for (i = 0; i<p.sizeX(); i++) for (j = 0; j<p.sizeY(); j++)
		{
			U val = p.get(i, j);
			this->update(i, j, lcol.at(val));
		}
	}

	template <class T> void hvPictRGB<T>::segmentPCA(hvPict<int> &pi, std::vector<hvColRGB<unsigned char> > &lcol, int sx, int sy, int nr) const
	{
		int i, j, k;
		hvFrame3<double> fr;
		hvLinearTransform3<double> t;
		int nnr, nng, nnb;
		std::cout << "segmenting with PCA...\n";

		fr = this->pca(255);
		t = hvLinearTransform3<double>(fr);
		hvPictRGB<unsigned char> pict; pict.clone(*this, 0, 0, this->sizeX() - 1, this->sizeY() - 1);
		pict.apply(255, fr, 0.5, 0.7);
		//printf("sizex=%d, sizey=%d\n", pict.sizeX(), pict.sizeY());

		pict.setLoop(0, hvField<hvColRGB<unsigned char> >::MIRROR);
		pict.setLoop(1, hvField<hvColRGB<unsigned char> >::MIRROR);
		pict.blur(9, 3, 3);
		//hvPictRGB<unsigned char> pb; pb.clone(*this, 0,0,sizeX()-1, sizeY()-1);
		//pb.blur<hvColRGB<double> >(&pict,sx,sy);

		hvQPictRGB<int, 64> pquant(pict, nr);
		pquant.median(3, 3);
		pquant.convert(pi, lcol);
		int l = lcol.size();
		lcol.clear();
		for (i = 0; i<l; i++)
		{
			hvBitmap bm(pi, hvBitmap::EQUAL, i);
			hvColRGB<unsigned char> v = avg(bm);
			lcol.push_back(v);
		}
	}

	// creates a hvVec4 picture by browsing connected components (fragments) in input picture
	// the hvVec4<U> contains: frag index, displement (x,y) to frag center, mean variational value
	template <class T> template <class U> void hvPictRGB<T>::fragmentRegular(hvPict<hvVec4<U> > &pf, const hvPictRGB<unsigned char> &pictinput, const hvQPictRGB<unsigned char, 64> &pquant) const
	{
		pf.reset(this->sizeX(), this->sizeY(), hvVec4<U>(U(-1)));
		hvBitmap yet(this->sizeX(), this->sizeY(), false);
		hvBitmap clas(this->sizeX(), this->sizeY(), false);
		U count = U(0);
		int i, j, x, y, ii, jj;
		hvColRGB<T> col;
		hvVec2<int> gmin, gmax, min, max;
		float meanx = 0.0, meany = 0.0, meanv = 0.0;
		int num = 0;

		std::vector<hvFrame3<double> > lfr(20); lfr.clear();
		pquant.plsr(255, pictinput, lfr);
		printf("PLSRinto %d colors\n", pquant.ncolors());
		for (i = 0; i<lfr.size(); i++)
		{
			hvFrame3<double> fr = lfr.at(i);
			printf("frame %d : %g,%g,%g\n", i, fr.origin().X()*255.0, fr.origin().Y()*255.0, fr.origin().Z()*255.0);
		}

		// features
		bool cont = false;
		do
		{
			cont = false;
			for (x = 0; x<this->sizeX(); x++) for (y = 0; y<this->sizeY(); y++)
			{
				if (!yet.get(x, y))
				{
					num = 0;
					meanv = 0.0;
					col = this->get(x, y);
					clas.hvBoolArray2::clear(false);
					min = hvVec2<int>(x, y); max = hvVec2<int>(x, y);
					this->seedfill(col, x, y, clas, min, max);
					for (i = min.X(); i <= max.X(); i++) for (j = min.Y(); j <= max.Y(); j++) if (clas.get(i, j))
					{
						yet.set(i, j, true);
						pf.update(i, j, hvVec4<U>(U(count)));
						int qid;
						hvColRGB<unsigned char> plsrcol = pquant.apply(255, pictinput.get(i, j), lfr, 0.5, 0.7, qid);
						num++;
						meanv += (float)plsrcol.RED();
					}
					meanv /= (float)num;
					meanx = (float)(max.X() + min.X())*0.5; meany = (float)(max.Y() + min.Y())*0.5;
					for (i = 0; i<this->sizeX(); i++) for (j = 0; j<this->sizeY(); j++)
					{
						if (pf.get(i, j).X() == U(count))
						{
							float vx = meanx, vy = meany;
							pf.update(i, j, hvVec4<U>(U(count), U(vx / (float)this->sizeX()*127.0 + 127.0), U(vy / (float)this->sizeY()*127.0 + 127.0), U(meanv)));
						}
					}
					printf("feature: %d->  %g,%g (%g)\n", num, meanx, meany, meanv);
					count++;
					cont = true;
				}
			}
		} while (cont);
	}

	template <class T> template <class U> U hvPictRGB<T>::fragmentWang(hvPict<hvVec4<U> > &pf, const hvPictRGB<unsigned char> &pictinput, const hvPictRGB<unsigned char> &psynth, const hvQPictRGB<unsigned char, 64> &pquant)
	{
		hvPictRGB<T> logpict;
		logpict.clone(*this, 0, 0, this->sizeX() - 1, this->sizeY() - 1);
		printf("fragmenting Wang tiles...\n");
		if (this->sizeX() % 4 != 0 || this->sizeY() % 4 != 0) { printf("picture resol must be multiple of four!\n"); return U(-1); }

		pf.reset(this->sizeX(), this->sizeY(), hvVec4<U>(U(-1)));
		hvBitmap yet(this->sizeX(), this->sizeY(), false);
		hvBitmap clas(this->sizeX(), this->sizeY(), false);
		U count = U(0);
		int rx = this->sizeX() / 4, ry = this->sizeY() / 4;
		int i, j, x, y, tx, ty;
		hvColRGB<T> col;
		hvVec2<int> gmin, gmax, min, max;
		float meanx = 0.0, meany = 0.0, meanv = 0.0;
		int num = 0;

		//pquant.quantizeHue(pictinput, ncolseg);
		std::vector<hvFrame3<double> > lfr(20); lfr.clear();
		pquant.plsr(255, pictinput, lfr);
		printf("PLSRinto %d colors\n", pquant.ncolors());
		for (i = 0; i<lfr.size(); i++)
		{
			hvFrame3<double> fr = lfr.at(i);
			//printf("frame %d : %g,%g,%g\n", i, fr.origin().X()*255.0,fr.origin().Y()*255.0,fr.origin().Z()*255.0);
		}

		// corners
		for (tx = 0; tx<4; tx++) for (ty = 0; ty<4; ty++)
		{
			//gmin=hvVec2<int>(tx*rx,ty*ry); gmax=hvVec2<int>(tx*rx,ty*ry);
			for (int delta = 0; delta<4; delta++)
			{
				x = tx*rx + (delta == 0 || delta == 2 ? 0 : -1); y = ty*ry + (delta == 0 || delta == 1 ? 0 : -1);
				if (x == -1) x += this->sizeX();
				if (y == -1) y += this->sizeY();
				if (!yet.get(x, y))
				{
					col = this->get(x, y);
					clas.hvBoolArray2::clear(false);
					min = hvVec2<int>(x, y); max = hvVec2<int>(x, y);
					this->seedfill(col, x, y, clas, min, max);
					for (i = min.X(); i <= max.X(); i++) for (j = min.Y(); j <= max.Y(); j++) if (clas.get(i, j)) { yet.set(i, j, true); pf.update(i, j, hvVec3<U>(U(count))); }
					min = hvVec2<int>(min.X() - tx*rx - (x == this->sizeX() - 1 ? this->sizeX() : 0), min.Y() - ty*ry - (y == this->sizeY() - 1 ? this->sizeY() : 0)); max = hvVec2<int>(max.X() - tx*rx - (x == this->sizeX() - 1 ? this->sizeX() : 0), max.Y() - ty*ry - (y == this->sizeY() - 1 ? this->sizeY() : 0));
					if (delta == 0) { gmin = min; gmax = max; }
					else { gmin.keepMin(gmin, min); gmax.keepMax(gmax, max); }
				}
			}
			//printf("seed %d,%d: %d,%d <-> %d,%d\n", tx,ty, gmin.X(),gmin.Y(),gmax.X(),gmax.Y());
			//printf("corner %d,%d: %g,%g\n", tx,ty, (float)(gmax.X()+gmin.X())*0.5, (float)(gmax.Y()+gmin.Y())*0.5);
			meanx += (float)(gmax.X() + gmin.X())*0.5; meany += (float)(gmax.Y() + gmin.Y())*0.5;
		}
		meanx /= 16.0; meany /= 16.0;
		num = 0;
		meanv = 0.0;
		for (i = 0; i<this->sizeX(); i++) for (j = 0; j<this->sizeY(); j++)
		{
			if (pf.get(i, j).X() == U(count))
			{
				int qid;
				hvColRGB<unsigned char> plsrcol = pquant.apply(255, psynth.get(i, j), lfr, 0.5, 0.7, qid);
				num++;
				meanv += (float)plsrcol.RED();
			}
		}
		meanv /= (float)num;
		for (i = 0; i<this->sizeX(); i++) for (j = 0; j<this->sizeY(); j++)
		{
			if (pf.get(i, j).X() == U(count))
			{
				float vx = meanx, vy = meany;
				if (i%rx>rx / 2) vx += rx;
				if (j%ry>ry / 2) vy += ry;
				pf.update(i, j, hvVec4<U>(U(count), U(vx / (float)rx*127.0*0.75 + 127.0), U(vy / (float)ry*127.0*0.75 + 127.0), U(meanv)));
				logpict.update(i, j, hvColRGB<T>(255, 255, 255));
			}
		}
		printf("corner feature: %d->  %g,%g (%g)\n", num, meanx, meany, meanv);
		count++;

		// top and bottom border features
		bool cont;
		for (int wline = 0; wline<2; wline++)
			do
			{
				printf("wline tb border=%d\n", wline);
				meanx = 0.0; meany = 0.0;
				//search for largest feature
				cont = true;
				hvSortedList< hvColRGB<unsigned char> > lf(100); lf.clear();
				hvArray1<int> ncol(100, 0);
				int iline = (wline == 0 ? 0 : this->sizeY() / 2);
				for (i = 0; i<rx; i++) if (!yet.get(i, iline))
				{
					col = this->get(i, iline);
					int lfp = lf.search(col);
					if (lfp == -1) { lf.push_back(col); lfp = lf.size() - 1; } //printf("push new color %d,%d,%d at %d,%d, %s\n", col.RED(), col.GREEN(), col.BLUE(), i, iline,yet.get(i,iline)?"yet":"! yet" ); }
					ncol[lfp]++;
				}
				//printf("colors: "); for (i=0; i<lf.length(); i++) printf("%d  ", ncol[i]); printf("\n");
				bool stillfeature = false;
				for (i = 0; i<lf.size(); i++) if (ncol[i]>3) stillfeature = true;
				if (stillfeature) // there is still a feature on border
				{
					//printf("found %d different features on border\n", lf.length()); 
					int nmax = 0, idmax = -1; for (i = 0; i<lf.size(); i++) if (ncol[i]>nmax) { nmax = ncol[i]; idmax = i; }
					col = lf.at(idmax);
					int bleft = rx, bright = 0;
					for (i = 0; i<rx; i++) if (!yet.get(i, iline) && this->get(i, iline) == col) { if (bleft>i) bleft = i; if (bright<i) bright = i; }
					if (bleft == 0 || bright == rx - 1) { printf("Warning: inconstitency in Wang tiles fragments! Corners already done!\n"); }
					if (bleft>bright) { printf("Warning: inconsistency in Wang tiles fragments! no element of color found!\n");  hvFatal("stop"); }
					//printf("largest feature is between %d,%d for color %d,%d,%d\n", bleft, bright, col.RED(), col.GREEN(), col.BLUE());
					//if (bright-bleft>=4) 
					{
						//largest feature is between [bleft,bright]
						for (tx = 0; tx<4; tx++) for (ty = 0; ty<2; ty++)
						{
							for (i = bleft; i <= bright; i++) logpict.update(i, iline, hvColRGB<T>(255, 255, 128));
							bool partaok = false;
							int startx = tx*rx; int starty = ty*ry + iline;
							lf.clear(); for (i = 0; i<100; i++) ncol[i] = 0;
							for (i = startx + bleft - 1; i<startx + bright + 1; i++) if (!yet.get(i, starty))
							{
								col = this->get(i, starty);
								int lfp = lf.search(col);
								if (lfp == -1) { lf.push_back(col); lfp = lf.size() - 1; }
								ncol[lfp]++;
							}
							nmax = 0; idmax = -1; for (i = 0; i<lf.size(); i++) if (ncol[i]>nmax) { nmax = ncol[i]; idmax = i; }
							//if (nmax==0) printf("warning no corresponding frag found in part a of tile %d,%d at pixels (%d-%d,%d)\n",tx,ty,startx+bleft-1,startx+bright+1,starty);
							if (nmax>0)
							{
								//printf("inconsistency in tile: %d,%d\n", tx,ty); hvFatal("stop"); }
								col = lf.at(idmax);
								//printf("tile(1) %d,%d color: %d,%d,%d\n", tx,ty, col.RED(), col.GREEN(), col.BLUE());
								bool doseed;
								int nloops = 0;
								do {
									doseed = false;
									//printf("round %d (%d-%d) color: %d,%d,%d\n", nloops,startx+bleft ,startx+bright-bleft+2, col.RED(), col.GREEN(), col.BLUE());
									for (i = startx + bleft - 1; i<startx + bright + 1; i++) if ((!yet.get(i, starty)) && this->get(i, starty) == col) { doseed = true; break; }
									if (doseed)
									{
										x = i; y = starty;
										clas.hvBoolArray2::clear(false);
										if (nloops == 0) { min = hvVec2<int>(x, y); max = hvVec2<int>(x, y); }
										//printf("max col %d,%d,%d ->starting seed loop %d at %d,%d for col %d,%d,%d\n", col.RED(), col.GREEN(), col.BLUE(),nloops,x,y,get(x,y).RED(),get(x,y).GREEN(),get(x,y).BLUE());
										this->seedfill(col, x, y, clas, min, max);
										//printf("seed box: %d,%d  -  %d,%d\n", min.X(),max.X(),min.Y(),max.Y());
										for (i = min.X(); i <= max.X(); i++) for (j = min.Y(); j <= max.Y(); j++)
											if (clas.get(i, j)) { yet.set(i, j, true); pf.update(i, j, hvVec3<U>(U(count))); } //logpict.update(i,j,hvColRGB<T>(255,255,255)); }
										nloops++;
									}
								} while (doseed);
								//printf("Seed provided %d pixels\n", clas.count());
								min = hvVec2<int>(min.X() - tx*rx, min.Y() - ty*ry - iline); max = hvVec2<int>(max.X() - tx*rx, max.Y() - ty*ry - iline);
								gmin = min; gmax = max;
								partaok = true;
							}
							//if (!partaok) { printf("no seed : lf=%d : ",lf.length());  for (i=0; i<lf.length(); i++)  printf("%d,",ncol[i]); printf("\n"); }

							startx = tx*rx; starty = ty*ry - 1 + iline; if (starty == -1) starty += this->sizeY();
							lf.clear(); for (i = 0; i<100; i++) ncol[i] = 0;
							bool need = true;
							if (partaok) for (i = startx + bleft - 1; i<startx + bright + 1; i++) if (this->get(i, starty) == col) { need = false; }
							if (need) // need because the feature does not share the same identifier on the other side
							{
								//printf("need because the feature does not share the same identifier\n");
								for (i = startx + bleft - 1; i<startx + bright + 1; i++) if (!yet.get(i, starty))
								{
									col = this->get(i, starty);
									int lfp = lf.search(col);
									if (lfp == -1) { lf.push_back(col); lfp = lf.size() - 1; }
									ncol[lfp]++;
								}
								nmax = 0; idmax = -1; for (i = 0; i<lf.size(); i++) if (ncol[i]>nmax) { nmax = ncol[i]; idmax = i; }
								//if (nmax==0) printf("warning no corresponding frag found in part b of tile %d,%d at pixels (%d-%d,%d)\n",tx,ty,startx+bleft-1,startx+bright+1,starty);
								if (nmax>0)
								{
									col = lf.at(idmax);
									//printf("tile(2) %d,%d color: %d,%d,%d\n", tx,ty, col.RED(), col.GREEN(), col.BLUE());
									bool doseed;
									int nloops = 0;
									do {
										doseed = false;
										for (i = startx + bleft; i<startx + bright + 1; i++) if (!yet.get(i, starty) && this->get(i, starty) == col) { doseed = true; break; }
										if (doseed)
										{
											x = i; y = starty;
											clas.hvBoolArray2::clear(false);
											if (nloops == 0) { min = hvVec2<int>(x, y); max = hvVec2<int>(x, y); }
											//printf("max col %d,%d,%d ->starting seed loop %d at %d,%d for col %d,%d,%d\n", col.RED(), col.GREEN(), col.BLUE(),nloops,x,y,get(x,y).RED(),get(x,y).GREEN(),get(x,y).BLUE());
											this->seedfill(col, x, y, clas, min, max);
											//printf("seed box: %d,%d  -  %d,%d\n", min.X(),max.X(),min.Y(),max.Y());
											for (i = min.X(); i <= max.X(); i++) for (j = min.Y(); j <= max.Y(); j++)
												if (clas.get(i, j)) { yet.set(i, j, true); pf.update(i, j, hvVec3<U>(U(count))); }//logpict.update(i,j,hvColRGB<T>(255,255,255)); }
											nloops++;
										}
									} while (doseed);
									min = hvVec2<int>(min.X() - tx*rx, min.Y() - ty*ry - iline - (starty == this->sizeY() - 1 ? this->sizeY() : 0)); max = hvVec2<int>(max.X() - tx*rx, max.Y() - ty*ry - iline - (starty == this->sizeY() - 1 ? this->sizeY() : 0));
									if (partaok)
									{
										gmin.keepMin(gmin, min); gmax.keepMax(gmax, max);
									}
									else
									{
										gmin = min; gmax = max;
										partaok = true;
									}
								}
							}
							//printf("seed %d,%d: %d,%d <-> %d,%d\n", tx,ty, gmin.X(),gmin.Y(),gmax.X(),gmax.Y());
							//printf("border %d,%d: %g,%g\n", tx,ty, (float)(gmax.X()+gmin.X())*0.5, (float)(gmax.Y()+gmin.Y())*0.5);
							meanx += (float)(gmax.X() + gmin.X())*0.5; meany += (float)(gmax.Y() + gmin.Y())*0.5;
							if (!partaok)
							{
								for (i = startx + bleft - 1; i<startx + bright + 1; i++) logpict.update(i, starty, hvColRGB<T>(128, 128, 200));
								FILE *fd = fopen("errorlog.ppm", "wb");
								logpict.savePPM(fd, 1);
								fclose(fd);
								printf("inconsistency in tile: %d,%d\n", tx, ty);
								hvFatal("see errorlog.ppm file. Stop");
							}
						}
						meanx /= 8.0; meany /= 8.0;
						num = 0;
						meanv = 0.0;
						for (i = 0; i<this->sizeX(); i++) for (j = 0; j<this->sizeY(); j++)
						{
							if (pf.get(i, j).X() == U(count))
							{
								int qid;
								hvColRGB<unsigned char> plsrcol = pquant.apply(255, psynth.get(i, j), lfr, 0.5, 0.7, qid);
								num++;
								meanv += (float)plsrcol.RED();
							}
						}
						meanv /= (float)num;
						for (i = 0; i<this->sizeX(); i++) for (j = 0; j<this->sizeY(); j++)
						{
							if (pf.get(i, j).X() == U(count))
							{
								float vx = meanx, vy = meany;
								//if (i%rx>rx/2) vx+=rx;
								if (j%ry>ry / 2) vy += ry;
								pf.update(i, j, hvVec4<U>(U(count), U(vx / (float)rx*127.0*0.75 + 127.0), U(vy / (float)ry*127.0*0.75 + 127.0), U(meanv)));
							}
						}
						//printf("new border feature: %d->  %g,%g\n", num, meanx, meany);
						count++;
					}
					//else cont=false;
				}
				else cont = false;
			} while (cont);

			// left and right border features
			for (int wline = 0; wline<2; wline++)
				do
				{
					//printf("wline lr border=%d\n\n",wline);
					meanx = 0.0; meany = 0.0;
					//search for largest feature
					cont = true;
					hvSortedList< hvColRGB<unsigned char> > lf(100); lf.clear();
					hvArray1<int> ncol(100, 0);
					int iline = (wline == 0 ? 0 : this->sizeX() / 2);
					for (j = 0; j<ry; j++) if (!yet.get(iline, j))
					{
						col = this->get(iline, j);
						int lfp = lf.search(col);
						if (lfp == -1) { lf.push_back(col); lfp = lf.size() - 1; } //printf("push new color %d,%d,%d at %d,%d, %s\n", col.RED(), col.GREEN(), col.BLUE(), iline,j, yet.get(iline,j)?"yet":"! yet" ); }
						ncol[lfp]++;
					}
					//printf("colors: "); for (i=0; i<lf.length(); i++) printf("%d  ", ncol[i]); printf("\n");
					bool stillfeature = false;
					for (i = 0; i<lf.size(); i++) if (ncol[i]>3) stillfeature = true;
					if (stillfeature) // there is still a feature on border
					{
						//printf("found %d different features on border\n", lf.length()); 
						int nmax = 0, idmax = -1; for (i = 0; i<lf.size(); i++) if (ncol[i]>nmax) { nmax = ncol[i]; idmax = i; }
						col = lf.at(idmax);
						int bleft = ry, bright = 0;
						for (j = 0; j<ry; j++) if (!yet.get(iline, j) && this->get(iline, j) == col) { if (bleft>j) bleft = j; if (bright<j) bright = j; }
						if (bleft == 0 || bright == ry - 1) { printf("Warning: inconstitency in Wang tiles fragments! Corners already done!\n"); }
						if (bleft>bright) { printf("Warning: inconsistency in Wang tiles fragments! no element of color found!\n");  hvFatal("stop"); }
						//printf("largest feature is between %d,%d for color %d,%d,%d\n", bleft, bright, col.RED(), col.GREEN(), col.BLUE());
						//if (bright-bleft>=4) 
						{
							//largest feature is between [bleft,bright]
							for (ty = 0; ty<4; ty++) for (tx = 0; tx<2; tx++)
							{
								for (j = bleft; j <= bright; j++) logpict.update(iline, j, hvColRGB<T>(255, 255, 128));
								bool partaok = false;
								int startx = tx*rx + iline; int starty = ty*ry;
								lf.clear(); for (i = 0; i<100; i++) ncol[i] = 0;
								for (j = starty + bleft - 1; j<starty + bright + 1; j++) if (!yet.get(startx, j))
								{
									col = this->get(startx, j);
									int lfp = lf.search(col);
									if (lfp == -1) { lf.push_back(col); lfp = lf.size() - 1; }
									ncol[lfp]++;
								}
								nmax = 0; idmax = -1; for (i = 0; i<lf.size(); i++) if (ncol[i] >= nmax) { nmax = ncol[i]; idmax = i; }
								//if (nmax==0) printf("warning no corresponding frag found in part a of tile %d,%d at pixels (%d,%d-%d)\n",tx,ty,startx,starty+bleft-1,starty+bright+1);
								if (nmax>0)
								{ //printf("inconsistency in tile: %d,%d\n", tx,ty); hvFatal("stop"); }
									col = lf.at(idmax);
									//printf("tile(1) %d,%d color: %d,%d,%d\n", tx,ty, col.RED(), col.GREEN(), col.BLUE());
									bool doseed;
									int nloops = 0;
									do {
										doseed = false;
										//printf("round %d (%d-%d) color: %d,%d,%d\n", nloops,startx+bleft ,startx+bright-bleft+2, col.RED(), col.GREEN(), col.BLUE());
										for (j = starty + bleft; j<starty + bright + 1; j++) if ((!yet.get(startx, j)) && this->get(startx, j) == col) { doseed = true; break; }
										if (doseed)
										{
											x = startx; y = j;
											clas.hvBoolArray2::clear(false);
											if (nloops == 0) { min = hvVec2<int>(x, y); max = hvVec2<int>(x, y); }
											//printf("max col %d,%d,%d ->starting seed loop %d at %d,%d for col %d,%d,%d\n", col.RED(), col.GREEN(), col.BLUE(),nloops,x,y,get(x,y).RED(),get(x,y).GREEN(),get(x,y).BLUE());
											this->seedfill(col, x, y, clas, min, max);
											//printf("seed box: %d,%d  -  %d,%d\n", min.X(),max.X(),min.Y(),max.Y());
											for (i = min.X(); i <= max.X(); i++) for (j = min.Y(); j <= max.Y(); j++)
												if (clas.get(i, j)) { yet.set(i, j, true); pf.update(i, j, hvVec3<U>(U(count))); }// logpict.update(i,j,hvColRGB<T>(255,255,255)); }
											nloops++;
										}
									} while (doseed);
									//printf("Seed provided %d pixels\n", clas.count());
									min = hvVec2<int>(min.X() - tx*rx - iline, min.Y() - ty*ry); max = hvVec2<int>(max.X() - tx*rx - iline, max.Y() - ty*ry);
									gmin = min; gmax = max;
									partaok = true;
								}
								//if (!partaok) { printf("tile %d,%d no seed : lf=%d : ",tx, ty, lf.length());  for (i=0; i<lf.length(); i++)  printf("%d,",ncol[i]); printf("\n"); }

								startx = tx*rx - 1 + iline; starty = ty*ry; if (startx == -1) startx += this->sizeX();
								lf.clear(); for (i = 0; i<100; i++) ncol[i] = 0;
								bool need = true;
								if (partaok) for (j = starty + bleft; j<starty + bright + 1; j++) if (this->get(startx, j) == col) { need = false; }
								if (need) // need because the feature does not share the same identifier on the other side
								{
									for (j = starty + bleft - 1; j<starty + bright + 1; j++) if (!yet.get(startx, j))
									{
										col = this->get(startx, j);
										int lfp = lf.search(col);
										if (lfp == -1) { lf.push_back(col); lfp = lf.size() - 1; }
										ncol[lfp]++;
									}
									nmax = 0; idmax = -1; for (i = 0; i<lf.size(); i++) if (ncol[i]>nmax) { nmax = ncol[i]; idmax = i; }
									//if (nmax==0) printf("warning no corresponding frag found in part b of tile %d,%d at pixels (%d,%d-%d)\n",tx,ty,startx,starty+bleft-1,starty+bright+1);
									if (nmax>0)
									{
										col = lf.at(idmax);
										//printf("tile(2) %d,%d color: %d,%d,%d\n", tx,ty, col.RED(), col.GREEN(), col.BLUE());
										bool doseed;
										int nloops = 0;
										do {
											doseed = false;
											for (j = starty + bleft; j<starty + bright + 1; j++) if (!yet.get(startx, j) && this->get(startx, j) == col) { doseed = true; break; }
											if (doseed)
											{
												x = startx; y = j;
												clas.hvBoolArray2::clear(false);
												if (nloops == 0) { min = hvVec2<int>(x, y); max = hvVec2<int>(x, y); }
												//printf("max col %d,%d,%d ->starting seed loop %d at %d,%d for col %d,%d,%d\n", col.RED(), col.GREEN(), col.BLUE(),nloops,x,y,get(x,y).RED(),get(x,y).GREEN(),get(x,y).BLUE());
												this->seedfill(col, x, y, clas, min, max);
												//printf("seed box: %d,%d  -  %d,%d\n", min.X(),max.X(),min.Y(),max.Y());
												for (i = min.X(); i <= max.X(); i++) for (j = min.Y(); j <= max.Y(); j++)
													if (clas.get(i, j)) { yet.set(i, j, true); pf.update(i, j, hvVec3<U>(U(count))); }// logpict.update(i,j,hvColRGB<T>(255,255,255)); }
												nloops++;
											}
										} while (doseed);
										min = hvVec2<int>(min.X() - tx*rx - iline - (startx == this->sizeX() - 1 ? this->sizeX() : 0), min.Y() - ty*ry); max = hvVec2<int>(max.X() - tx*rx - iline - (startx == this->sizeX() - 1 ? this->sizeX() : 0), max.Y() - ty*ry);
										if (partaok)
										{
											gmin.keepMin(gmin, min); gmax.keepMax(gmax, max);
										}
										else
										{
											gmin = min; gmax = max;
											partaok = true;
										}
									}
								}
								//printf("seed %d,%d: %d,%d <-> %d,%d\n", tx,ty, gmin.X(),gmin.Y(),gmax.X(),gmax.Y());
								if (!partaok)
								{
									for (j = starty + bleft - 1; j<starty + bright + 1; j++) logpict.update(startx, j, hvColRGB<T>(128, 128, 200));
									FILE *fd = fopen("errorlog.ppm", "wb");
									logpict.savePPM(fd, 1);
									fclose(fd);
									printf("inconsistency in tile: %d,%d\n", tx, ty);
									hvFatal("see errorlog.ppm file. Stop");
								}
								//printf("border %d,%d: %g,%g\n", tx,ty, (float)(gmax.X()+gmin.X())*0.5, (float)(gmax.Y()+gmin.Y())*0.5);
								meanx += (float)(gmax.X() + gmin.X())*0.5; meany += (float)(gmax.Y() + gmin.Y())*0.5;
							}
							meanx /= 8.0; meany /= 8.0;
							num = 0;
							meanv = 0.0;
							for (i = 0; i<this->sizeX(); i++) for (j = 0; j<this->sizeY(); j++)
							{
								if (pf.get(i, j).X() == U(count))
								{
									int qid;
									hvColRGB<unsigned char> plsrcol = pquant.apply(255, psynth.get(i, j), lfr, 0.5, 0.7, qid);
									num++;
									meanv += (float)plsrcol.RED();
								}
							}
							meanv /= (float)num;
							for (i = 0; i<this->sizeX(); i++) for (j = 0; j<this->sizeY(); j++)
							{
								if (pf.get(i, j).X() == U(count))
								{
									float vx = meanx, vy = meany;
									if (i%rx>rx / 2) vx += rx;
									//if (j%ry>ry/2) vy+=ry;
									pf.update(i, j, hvVec4<U>(U(count), U(vx / (float)rx*127.0*0.75 + 127.0), U(vy / (float)ry*127.0*0.75 + 127.0), U(meanv)));
								}
							}
							//printf("new border feature: %d->  %g,%g\n", num, meanx, meany);
							count++;
						}
						//else cont=false;
					}
					else cont = false;
				} while (cont);

				// all other features

				for (tx = 0; tx<4; tx++) for (ty = 0; ty<4; ty++)
				{
					printf("tile %d,%d...\n", tx, ty);
					do
					{
						cont = false;
						for (i = tx*rx; i<(tx + 1)*rx && !cont; i++) for (j = ty*ry; j<(ty + 1)*ry && !cont; j++)
						{
							if (!yet.get(i, j))
							{
								cont = true;
								x = i; y = j;
							}
						}
						if (cont)
						{
							col = this->get(x, y);
							clas.hvBoolArray2::clear(false);
							min = hvVec2<int>(x, y); max = hvVec2<int>(x, y);
							this->seedfill(col, x, y, clas, min, max);
							for (i = min.X(); i <= max.X(); i++) for (j = min.Y(); j <= max.Y(); j++) if (clas.get(i, j)) { yet.set(i, j, true); pf.update(i, j, hvVec3<U>(U(count))); }
							meanx = (float)(((max.X() + min.X()) / 2) % rx); meany = (float)(((max.Y() + min.Y()) / 2) % ry);
							meanv = 0.0;
							num = 0;
							for (i = tx*rx; i<(tx + 1)*rx; i++) for (j = ty*ry; j<(ty + 1)*ry; j++)
							{
								if (pf.get(i, j).X() == U(count))
								{
									int qid;
									hvColRGB<unsigned char> plsrcol = pquant.apply(255, psynth.get(i, j), lfr, 0.5, 0.7, qid);
									num++;
									meanv += (float)plsrcol.RED();
								}
							}
							meanv /= (float)num;
							for (i = tx*rx; i<(tx + 1)*rx; i++) for (j = ty*ry; j<(ty + 1)*ry; j++)
							{
								if (pf.get(i, j).X() == U(count))
								{
									pf.update(i, j, hvVec4<U>(U(count), U(meanx / (float)rx*127.0*0.75 + 127.0), U(meany / (float)ry*127.0*0.75 + 127.0), U(meanv)));
								}
							}
							//printf("new feature: ->  %g,%g\n", meanx, meany);
							count++;
						}
					} while (cont);
				}

				return count;
	}

////////////////////////////////////////////////////////////
//template <class T, unsigned int n> class hvQPict : public hvArray2< T >, public std::vector<unsigned char>  // T is an unsigned integer type (char, short or int)

template <class T, unsigned int n> void hvQPictRGB<T, n>::update(const hvPictRGB<unsigned char> &pi)
{
	int i, j;
	for (i = 0; i<pi.sizeX(); i++)
		for (j = 0; j<pi.sizeY(); j++)
		{
			hvColRGB<unsigned char> val = pi.get(i, j);
			hvArray2< T >::update(i, j, (T)closest(val));
		}
}
template <class T, unsigned int n> void hvQPictRGB<T, n>::updateTable(const hvPictRGB<unsigned char> &pi)
{
	int i, j;
	std::vector<hvColRGB<double> > coltab(n); coltab.clear();
	std::vector<int> count(n); count.clear();
	for (i = 0; i<std::vector<hvColRGB<unsigned char> >::size(); i++) { coltab.push_back(hvColRGB<double>(0.0)); count.push_back(0); }
	for (i = 0; i<pi.sizeX(); i++)
		for (j = 0; j<pi.sizeY(); j++)
		{
			hvColRGB<double> val = hvColRGB<double>(pi.get(i, j));
			val += coltab.at(hvArray2< T >::get(i, j));
			coltab[hvArray2< T >::get(i, j)] = val;
			count[hvArray2< T >::get(i, j)] = count.at(hvArray2< T >::get(i, j)) + 1;
		}
	for (i = 0; i<std::vector<hvColRGB<unsigned char> >::size(); i++)
	{
		hvColRGB<double> val = coltab.at(i);
		val /= (double)count.at(i);
		//printf("table %d: %g,%g,%g, count=%d\n", i, val.RED(), val.GREEN(), val.BLUE(), count.get(i));
		std::vector<hvColRGB<unsigned char> >::operator[](i) = hvColRGB<unsigned char>(val);
	}
}

template <class T, unsigned int n> hvQPictRGB<T, n>::hvQPictRGB(const hvPictRGB<unsigned char> &pi, int nbcol) : hvArray2< T >(pi.sizeX(), pi.sizeY(), T(0)), std::vector<hvColRGB<unsigned char> >(n)
{
	std::vector<hvColRGB<unsigned char> >::clear();
	this->quantize(pi, nbcol);
}
template <class T, unsigned int n> hvQPictRGB<T, n>::hvQPictRGB(const hvPictRGB<unsigned char> &pi, int nbcol, int level) : hvArray2< T >(pi.sizeX(), pi.sizeY(), T(0)), std::vector<hvColRGB<unsigned char> >(n)
{
	std::vector<hvColRGB<unsigned char> >::clear();
	this->quantize(pi, nbcol, level);
}
template <class T, unsigned int n> hvQPictRGB<T, n>::hvQPictRGB(const hvPictRGB<unsigned char> &pi, const hvPictRGB<unsigned char> &ps, int nbcol) : hvArray2< T >(pi.sizeX(), pi.sizeY(), T(0)), std::vector<hvColRGB<unsigned char> >(n)
{
	std::vector<hvColRGB<unsigned char> >::clear();
	quantize(pi, ps, nbcol);
}

template <class T, unsigned int n> void hvQPictRGB<T, n>::quantizeHue(const hvPictRGB<unsigned char> &pi, int nbcol, hvColRGB<double> ww, int level)
{
	hvPictRGB<unsigned char> pihsv;
	//pihsv.toHSV(pi, 255, ww);
	pihsv.toLUV(pi, 255, ww);
	quantize(pihsv, nbcol, level);
	updateTable(pi);
}
template <class T, unsigned int n> void hvQPictRGB<T, n>::quantizeLuv(const hvPictRGB<unsigned char> &pi, int nbcol, hvColRGB<double> ww, int level)
{
	hvPictRGB<unsigned char> pihsv;
	//pihsv.toHSV(pi, 255, ww);
	pihsv.toLUV(pi, 255, ww);
	quantize(pihsv, nbcol, level);
	updateTable(pi);
}


template <class T, unsigned int n> void hvQPictRGB<T, n>::quantize(const hvPictRGB<unsigned char> &pi, int nbcol, int level)
{
	//if (sizeX()!=pi.sizeX() || sizeY()!=pi.sizeY()) { hvFatal("hvQPictRGB<T,n>::quantize pictures must have same resolution"); return; }
	this->reset(pi.sizeX(), pi.sizeY(), T(0));
	if (nbcol <= 1) { std::vector<hvColRGB<unsigned char> >::push_back(pi.avg()); return; }
	int nc = (nbcol >= 8 ? nbcol : 8);
	hvOctree<hvQPictRGBVal> *troot = hvOctree<hvQPictRGBVal>::createOctree(hvQPictRGBVal());
	int nleaf = 0;
	hvPictRGB<unsigned char> ps;
	int i, j;
	//if (level>0) step=1<<level;
	if (level <= 1) ps.clone(pi, 0, 0, pi.sizeX() - 1, pi.sizeY() - 1);
	else {
		hvPictRGB<unsigned char> pbuf; pbuf.clone(pi, 0, 0, pi.sizeX() - 1, pi.sizeY() - 1);
		for (i = 1; i<level; i++) { ps.shrink(&pbuf); pbuf.clone(ps, 0, 0, ps.sizeX() - 1, ps.sizeY() - 1); }
	}
	for (i = 0; i<ps.sizeX(); i++)
		for (j = 0; j<ps.sizeY(); j++)
		{
			// insert the value
			hvColRGB<unsigned char> val = ps.get(i, j);
			hvQPictRGB<T, n>::insertValue(troot, val);
			// test the number of leafs
			hvQPictRGB<T, n>::updateTree(troot, nc);
		}
	// create the table
	std::vector<hvColRGB<unsigned char> >::clear();
	this->updateTable(troot);
	hvOctree<hvQPictRGBVal>::destroy(troot);
	this->update(pi);
	this->updateTable(pi);
	while (std::vector<hvColRGB<unsigned char> >::size()>nbcol)
	{
		this->reduce();
		this->update(pi);
		this->updateTable(pi);
	}
	//printf("Ncolors=%d\n", hvAList<hvColRGB<unsigned char>,n>::length());
}
template <class T, unsigned int n> void hvQPictRGB<T, n>::quantize(const hvPictRGB<unsigned char> &pi, const hvPictRGB<unsigned char> &ps, int nbcol)
{
	//if (sizeX()!=pi.sizeX() || sizeY()!=pi.sizeY()) { hvFatal("hvQPictRGB<T,n>::quantize pictures must have same resolution"); return; }
	this->reset(pi.sizeX(), pi.sizeY(), T(0));
	if (nbcol <= 1) { std::vector<hvColRGB<unsigned char> >::push_back(pi.avg()); return; }
	int nc = (nbcol >= 4 ? nbcol : 4);
	hvOctree<hvQPictRGBVal> *troot = hvOctree<hvQPictRGBVal>::createOctree(hvQPictRGBVal());
	int nleaf = 0;
	int i, j;
	for (i = 0; i<ps.sizeX(); i++)
		for (j = 0; j<ps.sizeY(); j++)
		{
			// insert the value
			hvColRGB<unsigned char> val = ps.get(i, j);
			hvQPictRGB<T, n>::insertValue(troot, val);
			// test the number of leafs
			hvQPictRGB<T, n>::updateTree(troot, nc);
		}
	// create the table
	std::vector<hvColRGB<unsigned char> >::clear();
	this->updateTable(troot);
	hvOctree<hvQPictRGBVal>::destroy(troot);
	this->update(pi);
	this->updateTable(pi);
	while (std::vector<hvColRGB<unsigned char> >::size()>nbcol)
	{
		this->reduce();
		this->update(pi);
		this->updateTable(pi);
	}
	//printf("Ncolors=%d\n", hvAList<hvColRGB<unsigned char>,n>::length());
}

template <class T, unsigned int n> void hvQPictRGB<T, n>::apply(T scal, hvPictRGB<unsigned char> &pi, std::vector<hvFrame3<double> > &lfr, double offset, double rescal)
{
	int i, j;
	for (i = 0; i<pi.sizeX(); i++) for (j = 0; j<pi.sizeY(); j++)
	{
		hvColRGB<T> v = pi.get(i, j);
		hvVec3<double> col((double)v.RED() / (double)scal, (double)v.GREEN() / (double)scal, (double)v.BLUE() / (double)scal);
		int q = this->hvArray2< T >::get(i, j);
		hvFrame3<double> fr = lfr.at(q);
		hvLinearTransform3<double> t; t.inverseFrame3(fr);
		col = t.apply(col);
		double rr = (col.X()*rescal + offset);
		if (rr<0.0) rr = 0.0; else if (rr>1.0) rr = 1.0;
		double gg = (col.Y()*rescal + offset);
		if (gg<0.0) gg = 0.0; else if (gg>1.0) gg = 1.0;
		double bb = (col.Z()*rescal + offset);
		if (bb<0.0) bb = 0.0; else if (bb>1.0) bb = 1.0;
		v = hvColRGB<T>((T)(rr*(double)scal), (T)(gg*(double)scal), (T)(bb*(double)scal));
		pi.update(i, j, v);
	}
}
template <class T, unsigned int n> void hvQPictRGB<T, n>::applyInverse(T scal, hvPictRGB<unsigned char> &pi, std::vector<hvFrame3<double> > &lfr, double offset, double rescal)
{
	int i, j;
	for (i = 0; i<pi.sizeX(); i++) for (j = 0; j<pi.sizeY(); j++)
	{
		int q = this->getIndex(i, j);
		hvFrame3<double> fr = lfr.at(q);
		hvLinearTransform3<double> t; t = hvLinearTransform3<double>(fr);
		hvColRGB<unsigned char> cc = pi.get(i, j);
		hvVec3<double> col(((double)cc.RED() / (double)scal - offset) / rescal, ((double)cc.GREEN() / (double)scal - offset) / rescal, ((double)cc.BLUE() / (double)scal - offset) / rescal);
		//hvVec3<double> col( ((double)cc.RED()/255.0-0.5)/0.7, 0.0,0.0);
		col = t.apply(col);
		double rr = col.X()*(double)scal; if (rr<0.0) rr = 0.0; else if (rr>(double)scal) rr = (double)scal;
		double gg = col.Y()*(double)scal; if (gg<0.0) gg = 0.0; else if (gg>(double)scal) gg = (double)scal;
		double bb = col.Z()*(double)scal; if (bb<0.0) bb = 0.0; else if (bb>(double)scal) bb = (double)scal;
		pi.update(i, j, hvColRGB<unsigned char>((unsigned char)(rr), (unsigned char)(gg), (unsigned char)(bb)));
	}
}

template <class T, unsigned int n> void hvQPictRGB<T, n>::update(double ww, T scal, hvPictRGB<unsigned char> &pi, std::vector<hvFrame3<double> > &lfr, double offset, double rescal)
{
	int i, j, q;
	T indmin;
	double errmin;
	for (i = 0; i<pi.sizeX(); i++)
		for (j = 0; j<pi.sizeY(); j++)
		{
			hvColRGB<T> v = pi.get(i, j);
			hvVec3<double> col((double)v.RED() / (double)scal, (double)v.GREEN() / (double)scal, (double)v.BLUE() / (double)scal);
			for (q = 0; q<this->ncolors(); q++)
			{
				hvColRGB<unsigned char> tabv = std::vector<hvColRGB<unsigned char> >::at(q);
				hvVec3<double> tabval((double)tabv.RED() / (double)scal, (double)tabv.GREEN() / (double)scal, (double)tabv.BLUE() / (double)scal);
				hvFrame3<double> fr = lfr.at(q);
				hvLinearTransform3<double> t;
				t.inverseFrame3(fr);
				hvVec3<double> collfr = col;
				collfr = t.apply(collfr);
				collfr = hvVec3<double>(collfr.X(), 0.0, 0.0);
				t = hvLinearTransform3<double>(fr);
				collfr = t.apply(collfr);
				hvVec3<double> errp; errp.PVec(col, collfr);
				hvVec3<double> errptab; errptab.PVec(tabval, col);
				if (q == 0) { errmin = ww*errp.norm() + (1.0 - ww)*errptab.norm(); indmin = 0; }
				else if (ww*errp.norm() + (1.0 - ww)*errptab.norm()<errmin) { errmin = ww*errp.norm() + (1.0 - ww)*errptab.norm(); indmin = q; }
			}
			//if (hvArray2< T >::get(i,i)!=indmin) printf("pixel %d,%d from %d to %d\n", i,j,hvArray2< T >::get(i,i),indmin);
			hvArray2< T >::update(i, j, (T)indmin);
		}
}

template <class T, unsigned int n> void hvQPictRGB<T, n>::plsr(T scal, const hvPictRGB<unsigned char> &pi, std::vector<hvFrame3<double> > &lfr) const
{
	lfr.clear();
	int i;
	for (i = 0; i<this->ncolors(); i++)
	{
		hvBitmap mm; this->convert(mm, i);
		hvFrame3<double> fr = pi.pca(255, &mm);
		lfr.push_back(fr);
	}
}
template <class T, unsigned int n> hvVec3<double> hvQPictRGB<T, n>::plsrmean(const hvPictRGB<unsigned char> &pi, const std::vector<hvFrame3<double> > &lfr, std::vector<double> &lvar) const
{
	lvar.clear();
	hvVec3<double> mean(0.0);
	for (int i = 0; i<this->ncolors(); i++)
	{
		hvBitmap mm; this->convert(mm, i);
		int count = mm.count();
		hvFrame3<double> fr = lfr.at(i);
		hvVec3<double> vv = fr.origin();
		vv.scale((double)count / (double)(this->sizeX()*this->sizeY()));
		mean += vv;
		int ii, jj;
		double var = 0.0;
		for (ii = 0; ii<pi.sizeX(); ii++) for (jj = 0; jj<pi.sizeY(); jj++)
		{
			if (mm.get(ii, jj))
			{
				hvColRGB<T> v = pi.get(ii, jj);
				hvVec3<double> col((double)v.RED() / (double)255, (double)v.GREEN() / (double)255, (double)v.BLUE() / (double)255);
				hvLinearTransform3<double> t; t.inverseFrame3(fr);
				col = t.apply(col);
				if (col.X()>var) var = col.X();
				//var+=abs(col.X());
			}
		}
		//var /= (double)count;
		lvar.push_back(var);
	}
	return mean;
}
template <class T, unsigned int n> void hvQPictRGB<T, n>::load(FILE *fd)
{
	int  sx, sy;
	int i, j, ncol;
	char buff[256];
	hvColRGB<T> co;

	hvPictRGB<T>::readPPMLine(fd, buff);
	if (strcmp(buff, "PQ\n") != 0) { hvFatal("cannot read PictQuant, not correct format"); return; }
	hvPictRGB<T>::readPPMLine(fd, buff);
	sscanf(buff, "%d %d %d", &ncol, &sx, &sy);
	if (ncol>n) { hvFatal("cannot read PictQuant, too many color indices"); return; }
	hvPictRGB<T>::readPPMLine(fd, buff);
	if (strcmp(buff, "255\n") != 0) { hvFatal("Not the right PictQuant Format"); }
	hvArray2< T >::reset(sx, sy);
	std::vector<hvColRGB<T> >::clear();
	for (i = 0; i<ncol; i++)
	{
		unsigned char r, g, b;
		fread(&r, 1, sizeof(unsigned char), fd);
		fread(&g, 1, sizeof(unsigned char), fd);
		fread(&b, 1, sizeof(unsigned char), fd);
		std::vector<hvColRGB<T> >::push_back(hvColRGB<T>((T)r, (T)g, (T)b));
	}
	for (i = 0; i<sy; i++)
		for (j = 0; j<sx; j++)
		{
			unsigned char v;
			fread(&v, 1, sizeof(unsigned char), fd);
			hvArray2<T>::update(j, sy - i - 1, (T)v);
		}
}

/////////////////////////////////////////////////////
// class hvBitmap
template <class T, class U> hvBitmap::hvBitmap(const hvPict<T> &p, hvBitmap::operation op, U value) : hvBoolArray2(p.sizeX(), p.sizeY(), false)
{
	int i, j;
	for (i = 0; i<p.sizeX(); i++) for (j = 0; j<p.sizeY(); j++)
	{
		switch (op)
		{
		case LESS: if (U(p.get(i, j))<value) set(i, j, true); break;
		case LEQUAL: if (U(p.get(i, j)) <= value) set(i, j, true); break;
		case EQUAL: if (U(p.get(i, j)) == value) set(i, j, true); break;
		case GEQUAL: if (U(p.get(i, j)) >= value) set(i, j, true); break;
		case GREATER: if (U(p.get(i, j))>value) set(i, j, true); break;
		case NOTEQUAL: if (U(p.get(i, j)) != value) set(i, j, true); break;
		default: break;
		}
	}
}
template <class T, class U> void hvBitmap::convert(const hvPict<T> &p, hvBitmap::operation op, U value)
{
	int i, j;
	hvBoolArray2::reset(p.sizeX(), p.sizeY(), false);
	for (i = 0; i<p.sizeX(); i++) for (j = 0; j<p.sizeY(); j++)
	{
		switch (op)
		{
		case LESS: if (U(p.get(i, j))<value) set(i, j, true); break;
		case LEQUAL: if (U(p.get(i, j)) <= value) set(i, j, true); break;
		case EQUAL: if (U(p.get(i, j)) == value) set(i, j, true); break;
		case GEQUAL: if (U(p.get(i, j)) >= value) set(i, j, true); break;
		case GREATER: if (U(p.get(i, j))>value) set(i, j, true); break;
		case NOTEQUAL: if (U(p.get(i, j)) != value) set(i, j, true); break;
		default: break;
		}
	}
}
template <class T, class U> void hvBitmap::convert(const hvPict<T> &p, hvBitmap::operation op, int nn, U value[])
{
	int i, j, k;
	hvBoolArray2::reset(p.sizeX(), p.sizeY(), false);
	for (i = 0; i<p.sizeX(); i++) for (j = 0; j<p.sizeY(); j++)
	{
		for (k = 0; k<nn && get(i, j) == false; k++)
		{
			switch (op)
			{
			case LESS: if (U(p.get(i, j))<value[k]) set(i, j, true); break;
			case LEQUAL: if (U(p.get(i, j)) <= value[k]) set(i, j, true); break;
			case EQUAL: if (U(p.get(i, j)) == value[k]) set(i, j, true);  break;
			case GEQUAL: if (U(p.get(i, j)) >= value[k]) set(i, j, true); break;
			case GREATER: if (U(p.get(i, j))>value[k]) set(i, j, true); break;
			case NOTEQUAL: if (U(p.get(i, j)) != value[k]) set(i, j, true); break;
			default: break;
			}
		}
	}
}

}

#endif // !efined(AFX_COLOR_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_)
