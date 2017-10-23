// hvQPicture.h: interface for the quantized picture class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_QPICTURE_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_)
#define AFX_QPICTURE_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000



#include "hvBinTree.h"
#include "hvOctree.h"

namespace hview {

template <class T> class hvPictRGB;

class hvQPictVal 
{
	bool e;
	int n;
	int v;
public:
	hvQPictVal() { e=false; n=0; v=0; }
	hvQPictVal(bool x) { e=x; n=0; v=0; }
	void incr() { n++; }
	void incr(int x) { n+=x; }
	void operator+=(int x) { v += x; }
	void set(bool x) { e=x; }
	int nleafs() const { return n; }
	int value() const { return v; }
	bool isEnd() const { return e; }
};

////////////////////////////////////////////////////////////
template <class T, unsigned int n> class hvQPict : public hvArray2< T >, public std::vector<unsigned char>  // T is an unsigned integer type (char, short or int)
////////////////////////////////////////////////////////////
{
	static int NCOLORS(hvBinTree<hvQPictVal> *o)
	{
		int i,nb;

		if (o==0) return(0);
		if (hvBinTree<hvQPictVal>::getValue(o).isEnd()) return(hvBinTree<hvQPictVal>::getValue(o).nleafs());
		nb=hvQPict<T,n>::NCOLORS(hvBinTree<hvQPictVal>::getLeft(o));
		nb+=hvQPict<T,n>::NCOLORS(hvBinTree<hvQPictVal>::getRight(o));
		return(nb);
	}

	static hvBinTree<hvQPictVal> *reduce(hvBinTree<hvQPictVal> *o, int start_depth, int *res_depth)
	{
		int dl, dr;
		if (o==0) { *res_depth=8; return 0; }
		bool red=false;
		if (! hvBinTree<hvQPictVal>::getValue(o).isEnd()) 
		{
			red=true;
			if (hvBinTree<hvQPictVal>::getLeft(o)!=0) { if (! hvBinTree<hvQPictVal>::getValue(hvBinTree<hvQPictVal>::getLeft(o)).isEnd()) red=false; }
			if (hvBinTree<hvQPictVal>::getRight(o)!=0) { if (! hvBinTree<hvQPictVal>::getValue(hvBinTree<hvQPictVal>::getRight(o)).isEnd()) red=false; }
		}
		if (red) { *res_depth=start_depth; return o; }
		hvBinTree<hvQPictVal> *ol = hvQPict<T,n>::reduce(hvBinTree<hvQPictVal>::getLeft(o),start_depth-1,&dl);
		hvBinTree<hvQPictVal> *ori = hvQPict<T,n>::reduce(hvBinTree<hvQPictVal>::getRight(o),start_depth-1,&dr);
		if (dr<dl) { *res_depth=dr; return ori; }
		else if (dl<dr) { *res_depth=dl; return ol; }
		if (hvQPict<T,n>::NCOLORS(ol)<hvQPict<T,n>::NCOLORS(ori)) { *res_depth=dl; return ol; }
		else { *res_depth=dr; return ori; }
	}
	static void updateTree(hvBinTree<hvQPictVal> *troot, int nbcol)
	{
		while (hvBinTree<hvQPictVal>::countAllLeafs(troot)>nbcol)
		{
			//printf("reduction:%d\n", hvBinTree<hvQPictVal>::countAllLeafs(troot));
			int dd;
			hvBinTree<hvQPictVal> *prev = hvQPict<T,n>::reduce(troot,7,&dd);
			if (prev==0) { hvFatal("null in reduce!"); }
			hvQPictVal vv(true);
			hvQPictVal vl;
			if (hvBinTree<hvQPictVal>::getLeft(prev)!=0)
			{
				vl = hvBinTree<hvQPictVal>::getValue(hvBinTree<hvQPictVal>::getLeft(prev));
				vv.incr(vl.nleafs()); vv+=vl.value();
			}
			if (hvBinTree<hvQPictVal>::getRight(prev)!=0)
			{
				vl = hvBinTree<hvQPictVal>::getValue(hvBinTree<hvQPictVal>::getRight(prev));
				vv.incr(vl.nleafs()); vv+=vl.value();
			}
			hvBinTree<hvQPictVal>::update(prev,vv); 
			hvBinTree<hvQPictVal>::destroy(hvBinTree<hvQPictVal>::getLeft(prev));
			hvBinTree<hvQPictVal>::destroy(hvBinTree<hvQPictVal>::getRight(prev));
			hvBinTree<hvQPictVal>::setLeft(prev,0);
			hvBinTree<hvQPictVal>::setRight(prev,0);
		}
	}
	static void insertValue(hvBinTree<hvQPictVal> *troot, unsigned char val)
	{
		int depth = 7;
		bool cont = true;
		hvBinTree<hvQPictVal> *tt, *prev;
		prev = troot;
		int quel=0; 
		if ((val & (unsigned char)(1<<depth))==0) tt=hvBinTree<hvQPictVal>::getLeft(troot); else { quel=1; tt=hvBinTree<hvQPictVal>::getRight(troot); }
		while (depth>=0 && cont)
		{
			if (tt==0) 
			{
				tt = hvBinTree<hvQPictVal>::createBinTree(hvQPictVal(depth==0));
				if (quel==0) hvBinTree<hvQPictVal>::setLeft(prev, tt); else hvBinTree<hvQPictVal>::setRight(prev, tt);
			}
			hvQPictVal vv = hvBinTree<hvQPictVal>::getValue(tt);
			if (vv.isEnd())
			{
				vv.incr(); vv+=(int)val; 
				hvBinTree<hvQPictVal>::update(tt, vv);
				cont = false;
			}
			else { depth--; prev=tt; if ((val & (unsigned char)(1<<depth))==0) { quel=0; tt=hvBinTree<hvQPictVal>::getLeft(tt); } else { quel=1; tt=hvBinTree<hvQPictVal>::getRight(tt); } }
			//printf("i=%d,j=%d,val=%d, depth=%d, count=%d\n", i,j,(int)val,depth, hvBinTree<hvQPictVal>::countAllLeafs(troot));
		}
	}
	void updateTable(hvBinTree<hvQPictVal> *troot)
	{
		if (troot==0) return;
		hvQPictVal vv = hvBinTree<hvQPictVal>::getValue(troot);
		if (vv.isEnd()) { std::vector<unsigned char>::push_back((unsigned char)(vv.value()/vv.nleafs())); return; }
		this->updateTable(hvBinTree<hvQPictVal>::getLeft(troot));
		this->updateTable(hvBinTree<hvQPictVal>::getRight(troot));
	}
	void update(const hvPict<unsigned char> &pi)
	{
		int i,j;
		for (i=0; i<pi.sizeX(); i++)
		for (j=0; j<pi.sizeY(); j++)
		{
			unsigned char val = pi.get(i,j);
			hvArray2< T >::update(i,j, (T)closest(val));
		}
	}

public:
	hvQPict<T, n>() : hvArray2< T >(), std::vector<unsigned char>(n) { std::vector<unsigned char>::clear();  }
	hvQPict<T,n>(int sx, int sy, T nil) : hvArray2< T >(sx, sy, nil), std::vector<unsigned char>(n) { std::vector<unsigned char>::clear(); }
	hvQPict<T,n>(const hvPict<unsigned char> &pi, int nbcol, bool regular=false) : hvArray2< T >(pi.sizeX(), pi.sizeY(), T(0)), std::vector<unsigned char>(n)
	{ 
		std::vector<unsigned char>::clear();
		if (nbcol<=1) { std::vector<unsigned char>::push_back((unsigned char)pi.avg()); return; }
		if (regular)
		{
			unsigned char min, max; 
			pi.minmax(min, max);
			printf("quant min=%d, max=%d\n",(int)min, (int)max);
			for (int i=0; i<nbcol; i++) std::vector<unsigned char>::push_back((unsigned char)((double)min+(double)(max-min)/(double)nbcol*((double)i+0.5)));
		}
		else
		{
			hvBinTree<hvQPictVal> *troot = hvBinTree<hvQPictVal>::createBinTree(hvQPictVal());
			int nleaf = 0;
			int i,j;
			for (i=0; i<pi.sizeX(); i++)
			for (j=0; j<pi.sizeY(); j++)
			{
				// insert the value
				unsigned char val = pi.get(i,j);
				hvQPict<T,n>::insertValue(troot, val);
				// test the number of leafs
				hvQPict<T,n>::updateTree(troot, nbcol);
			}
			// create the table
			std::vector<unsigned char>::clear();
			this->updateTable(troot);
			hvBinTree<hvQPictVal>::destroy(troot);
		}
		//printf("Ncolors=%d\n", hvAList<unsigned char,n>::length());
		this->update(pi);
	}

	int closest(unsigned char col) const
	{
		int idmin=0;
		unsigned char diffmin = (col<std::vector<unsigned char>::at(0)?std::vector<unsigned char>::at(0)-col:col- std::vector<unsigned char>::at(0));
		int i;
		for (i=1; i<std::vector<unsigned char>::size(); i++)
		{
			unsigned char diff = (col<std::vector<unsigned char>::at(i)? std::vector<unsigned char>::at(i)-col:col- std::vector<unsigned char>::at(i));
			if (diff<diffmin) { diffmin=diff; idmin=i; }
		}
		return idmin;
	}

	template <class X, class Y> void convert(hvPict<X> &pi, Y scal)
	{
		pi.reset(this->sizeX(), this->sizeY(), X(0));
		int i,j;
		for (i=0; i<this->sizeX(); i++)
		for (j=0; j<this->sizeY(); j++)
		{
			pi.update(i,j, X(Y(std::vector<unsigned char>::at(hvArray2< T >::get(i,j)))*scal));
		}
	}

	void convert(hvBitmap &bm, unsigned char value, hvBitmap::operation op)
	{
		int i,j;
		bm.reset(this->sizeX(), this->sizeY(), false);
		for (i=0; i<this->sizeX(); i++)
		for (j=0; j<this->sizeY(); j++)
		{
			unsigned char v= std::vector<unsigned char>::at(hvArray2< T >::get(i,j));
			switch(op)
			{
			case hvBitmap::LESS: if (v<value) bm.set(i,j,true); break;
			case hvBitmap::LEQUAL: if (v<=value) bm.set(i,j,true);  break;
			case hvBitmap::EQUAL: if (v==value) bm.set(i,j,true);  break; 
			case hvBitmap::GEQUAL: if (v>=value) bm.set(i,j,true);  break; 
			case hvBitmap::GREATER: if (v>value) bm.set(i,j,true);  break; 
			case hvBitmap::NOTEQUAL: if (v!=value) bm.set(i,j,true);  break;
			default: break;
			}
		}
	}

};

class hvQPictRGBVal 
{
	bool e;
	int n;
	hvVec3<int> v;
public:
	hvQPictRGBVal() { e=false; n=0; v=hvVec3<int>(0); }
	hvQPictRGBVal(bool x) { e=x; n=0; v=hvVec3<int>(0); }
	void incr() { n++; }
	void incr(int x) { n+=x; }
	void operator+=(const hvVec3<int> &x) { v += x; }
	void set(bool x) { e=x; }
	int nleafs() const { return n; }
	hvVec3<int> value() const { return v; }
	bool isEnd() const { return e; }
};

////////////////////////////////////////////////////////////
template <class T, unsigned int n> class hvQPictRGB : public virtual hvArray2< T >, public virtual std::vector<hvColRGB<unsigned char> >  // T is an unsigned integer type (char, short or int)
////////////////////////////////////////////////////////////
{
	static int NCOLORS(hvOctree<hvQPictRGBVal> *o)
	{
		int i,nb=0;

		if (o==0) return(0);
		if (hvOctree<hvQPictRGBVal>::getValue(o).isEnd()) return(hvOctree<hvQPictRGBVal>::getValue(o).nleafs());
		for (i=0; i<8; i++) nb+=hvQPictRGB<T,n>::NCOLORS(hvOctree<hvQPictRGBVal>::getChildNode(o,i));
		return(nb);
	}

	static hvOctree<hvQPictRGBVal> *reduce(hvOctree<hvQPictRGBVal> *o, int start_depth, int *res_depth)
	{
		int i, d[8],best;
		hvOctree<hvQPictRGBVal> *oc[8];

		if (o==0) { *res_depth=8; return 0; }
		bool red=false;
		if (! hvOctree<hvQPictRGBVal>::getValue(o).isEnd()) 
		{
			red=true;
			for (i=0; i<8; i++) if (hvOctree<hvQPictRGBVal>::getChildNode(o,i)!=0) { if (! hvOctree<hvQPictRGBVal>::getValue(hvOctree<hvQPictRGBVal>::getChildNode(o,i)).isEnd()) red=false; }
		}
		if (red) { *res_depth=start_depth; return o; }
		for (i=0; i<8; i++)
		{
			oc[i] = hvQPictRGB<T,n>::reduce(hvOctree<hvQPictRGBVal>::getChildNode(o,i),start_depth-1,d+i);
		}
		best=0;
		for (i=1; i<8; i++)
		{
			if (d[i]<d[best]) best=i;
			else if (d[i]==d[best] && hvQPictRGB<T,n>::NCOLORS(oc[i])<hvQPictRGB<T,n>::NCOLORS(oc[best])) best = i;
		}
		*res_depth = d[best];
		return(oc[best]);
	}

	static void updateTree(hvOctree<hvQPictRGBVal> *troot, int nbcol)
	{
		while (hvOctree<hvQPictRGBVal>::countAllLeafs(troot)>nbcol)
		{
			//printf("reduction:%d\n", hvBinTree<hvQPictVal>::countAllLeafs(troot));
			int dd;
			hvOctree<hvQPictRGBVal> *prev = hvQPictRGB<T,n>::reduce(troot,7,&dd);
			if (prev==0) { hvFatal("null in reduce!"); }
			hvQPictRGBVal vv(true);
			hvQPictRGBVal vl;
			int i;
			for (i=0; i<8; i++)
			{
				if (hvOctree<hvQPictRGBVal>::getChildNode(prev,i)!=0)
				{
					vl = hvOctree<hvQPictRGBVal>::getValue(hvOctree<hvQPictRGBVal>::getChildNode(prev,i));
					vv.incr(vl.nleafs()); vv+=vl.value();
				}
			}
			hvOctree<hvQPictRGBVal>::update(prev,vv); 
			for (i=0; i<8; i++) 
			{
				hvOctree<hvQPictRGBVal>::destroy(hvOctree<hvQPictRGBVal>::getChildNode(prev,i));
				hvOctree<hvQPictRGBVal>::setChildNode(prev,i,0);
			}
		}
	}

	static int BRANCH(const hvColRGB<unsigned char> &col, int depth)
	{
	int s,b;

	b = 1<<depth;
	s = 0;
	if (col.RED()&b) s |= 1;
	if (col.GREEN()&b) s |= 2;
	if (col.BLUE()&b) s |= 4;
	return(s);
	}

	static void insertValue(hvOctree<hvQPictRGBVal> *troot, const hvColRGB<unsigned char> &val)
	{
		int depth = 7;
		bool cont = true;
		hvOctree<hvQPictRGBVal> *tt, *prev;
		prev = troot;
		int quel=BRANCH(val,depth); 
		tt=hvOctree<hvQPictRGBVal>::getChildNode(troot,quel);
		while (depth>=0 && cont)
		{
			if (tt==0) 
			{
				tt = hvOctree<hvQPictRGBVal>::createOctree(hvQPictRGBVal(depth==0));
				hvOctree<hvQPictRGBVal>::setChildNode(prev, quel, tt);
			}
			hvQPictRGBVal vv = hvOctree<hvQPictRGBVal>::getValue(tt);
			if (vv.isEnd())
			{
				vv.incr(); vv+=hvVec3<int>((int)val.RED(),(int)val.GREEN(), (int)val.BLUE()); 
				hvOctree<hvQPictRGBVal>::update(tt, vv);
				cont = false;
			}
			else { depth--; prev=tt; quel=BRANCH(val,depth); tt=hvOctree<hvQPictRGBVal>::getChildNode(tt,quel); } 
			//printf("i=%d,j=%d,val=%d, depth=%d, count=%d\n", i,j,(int)val,depth, hvBinTree<hvQPictVal>::countAllLeafs(troot));
		}
	}
	public:
	void updateTable(hvOctree<hvQPictRGBVal> *troot)
	{
		if (troot==0) return;
		hvQPictRGBVal vv = hvOctree<hvQPictRGBVal>::getValue(troot);
		if (vv.isEnd()) { std::vector<hvColRGB<unsigned char> >::push_back(hvColRGB<unsigned char>((unsigned char)(vv.value().X()/vv.nleafs()),(unsigned char)(vv.value().Y()/vv.nleafs()),(unsigned char)(vv.value().Z()/vv.nleafs()))); return; }
		int i; for (i=0; i<8; i++) updateTable(hvOctree<hvQPictRGBVal>::getChildNode(troot,i));
	}
	void update(const hvPictRGB<unsigned char> &pi);
	void updateTable(const hvPictRGB<unsigned char> &pi);

	void newTable(const std::vector<hvColRGB<unsigned char> > &lcol)
	{
		std::vector<hvColRGB<unsigned char> >::clear();
		int i;
		for (i=0; i<lcol.size(); i++) std::vector<hvColRGB<unsigned char> >::push_back(lcol.at(i));
	}
public:
	hvQPictRGB<T, n>() : hvArray2< T >(), std::vector<hvColRGB<unsigned char> >(n) { std::vector<hvColRGB<unsigned char> >::clear();  }
	hvQPictRGB<T,n>(int sx, int sy, T nil) : hvArray2< T >(sx, sy, nil), std::vector<hvColRGB<unsigned char> >(n) { std::vector<hvColRGB<unsigned char> >::clear(); }

	hvQPictRGB<T, n>(const hvPictRGB<unsigned char> &pi, int nbcol);
	hvQPictRGB<T, n>(const hvPictRGB<unsigned char> &pi, int nbcol, int level);
	hvQPictRGB<T, n>(const hvPictRGB<unsigned char> &pi, const hvPictRGB<unsigned char> &ps, int nbcol);

	hvQPictRGB<T,n>(hvQPictRGB<T,n> &q) { hvFatal("no temporary creation of hvQPictRGB!"); }
	
	hvQPictRGB<T,n> &operator=(const hvQPictRGB<T,n> &q)
	{
		((hvArray2< T > *)this)->clone(q);
		std::vector<hvColRGB<unsigned char> >::operator=((std::vector<hvColRGB<unsigned char> >)q);
		return *this;
	}
	void reset(int sx, int sy, T nil) { hvArray2< T >::reset(sx, sy, nil); std::vector<hvColRGB<unsigned char> >::clear(); }
	int ncolors() const { return std::vector<hvColRGB<unsigned char> >::size(); }
	T getIndex(int x, int y) const { return hvArray2< T >::get(x,y); }
	hvColRGB<unsigned char> getColor(int ind) const { return std::vector<hvColRGB<unsigned char> >::at(ind); }

	void quantizeHue(const hvPictRGB<unsigned char> &pi, int nbcol, hvColRGB<double> ww = hvColRGB<double>(1.0, 1.0, 1.0), int level = 0 );
	void quantizeLuv(const hvPictRGB<unsigned char> &pi, int nbcol, hvColRGB<double> ww = hvColRGB<double>(1.0, 1.0, 1.0), int level = 0 );
	void quantize(const hvPictRGB<unsigned char> &pi, int nbcol, int level = 0 );
	void quantize(const hvPictRGB<unsigned char> &pi, const hvPictRGB<unsigned char> &ps, int nbcol);

	void median(int sx, int sy)
	{
		int i,j,k,ii,jj;
		int rx = sx/2;
		int ry = sy/2;
		int count[n];
		hvArray2< T > npict(this->sizeX(), this->sizeY(), T(0));
		for (i=0; i<this->sizeX(); i++)
		for (j=0; j<this->sizeY(); j++)
		{
			for (k=0; k< std::vector<hvColRGB<unsigned char> >::size(); k++) count[k]=0;
			for (ii=-rx; ii<=rx; ii++) for (jj=-ry; jj<=ry; jj++)
			{
				if (i+ii>=0 && i+ii<this->sizeX() && j+jj>=0 && j+jj<this->sizeY())
				{
					T id = hvArray2< T >::get(i+ii,j+jj);
					count[id]++;
				}
			}
			int imax = 0; int max=count[0];
			for (k=1; k< std::vector<hvColRGB<unsigned char> >::size(); k++) if (count[k]>max) {max=count[k]; imax=k; }
			npict.update(i,j,imax);
		}
		for (i=0; i<this->sizeX(); i++)
		for (j=0; j<this->sizeY(); j++)
		{
			hvArray2< T >::update(i,j,npict.get(i,j));
		}
	}

	int closest(const hvColRGB<unsigned char> &col) const
	{
		int idmin=0;
		hvColRGB<unsigned char> val = std::vector<hvColRGB<unsigned char> >::at(0);
		hvColRGB<unsigned char> diffval; diffval.subabs(val,col);
		unsigned char diffmin = diffval.luminance();
		int i;
		for (i=1; i<std::vector<hvColRGB<unsigned char> >::size(); i++)
		{
			val = std::vector<hvColRGB<unsigned char> >::at(i);
			diffval.subabs(val,col);
			unsigned char diff = diffval.luminance();
			if (diff<diffmin) { diffmin=diff; idmin=i; }
		}
		return idmin;
	}

	void reduce()
	{
		int i,j, inda, indb;
		unsigned char err = 255, diffmin=0;
		for(i=0; i<std::vector<hvColRGB<unsigned char> >::size(); i++)
		{
			hvColRGB<unsigned char> val = std::vector<hvColRGB<unsigned char> >::at(i);
			for (j=i+1; j<std::vector<hvColRGB<unsigned char> >::size(); j++)
			{
				hvColRGB<unsigned char> col = std::vector<hvColRGB<unsigned char> >::at(j);
				hvColRGB<unsigned char> diffval; diffval.subabs(val,col);
				diffmin = diffval.luminance();
				if (diffmin<err) { inda = i; indb=j; err=diffmin; }
			}
		}
		std::vector<hvColRGB<unsigned char> >::erase(std::vector<hvColRGB<unsigned char> >::begin()+indb);
	}

	template <class X, class Y> void convert(hvPictRGB<X> &pi, Y scal)
	{
		pi.reset(this->sizeX(), this->sizeY(), hvColRGB<X>(0));
		int i,j;
		for (i=0; i<this->sizeX(); i++)
		for (j=0; j<this->sizeY(); j++)
		{
			hvColRGB<unsigned char> col = std::vector<hvColRGB<unsigned char> >::at(hvArray2< T >::get(i,j));
			hvColRGB<Y> cc = hvColRGB<Y>(col);
			cc.scale(scal);
			pi.update(i,j, hvColRGB<X>(cc));
		}
	}
	void convert(hvBitmap &pi, T q) const
	{
		pi.reset(this->sizeX(), this->sizeY(), false);
		int i,j;
		for (i=0; i<this->sizeX(); i++)
		for (j=0; j<this->sizeY(); j++)
		{
			T v = hvArray2< T >::get(i,j);
			if (v==q) pi.set(i,j, true);
		}
	}
	template <class X> void convert(hvPict<X> &pi)
	{
		pi.reset(this->sizeX(), this->sizeY(), X(0));
		int i,j;
		for (i=0; i<this->sizeX(); i++)
		for (j=0; j<this->sizeY(); j++)
		{
			pi.update(i,j,X(hvArray2< T >::get(i,j)));
		}
	}

	template <class X> void convert(hvPict<X> &pi, std::vector<hvColRGB<unsigned char> > &lcol)
	{
		this->convert(pi);
		lcol.clear();
		for (int i=0; i<std::vector<hvColRGB<unsigned char> >::size(); i++) lcol.push_back(std::vector<hvColRGB<unsigned char> >::at(i));
	}

	void gamma(double power) 
	{
		int i,j;
		for (i=0;i<std::vector<hvColRGB<unsigned char> >::size(); i++)
		{
			hvColRGB<unsigned char> v = std::vector<hvColRGB<unsigned char> >::at(i);
			v.gamma(255,power);
			std::vector<hvColRGB<unsigned char> >::operator[](i)=v;
		}
	}

	hvColRGB<T> apply(T scal, const hvColRGB<T> &v, std::vector<hvFrame3<double> > &lfr, double offset, double rescal, int &q) const
	{
			hvVec3<double> col((double)v.RED(), (double)v.GREEN(),(double)v.BLUE());
			q = (int)this->closest(v);
			hvFrame3<double> fr = lfr.at(q);
			hvLinearTransform3<double> t; t.inverseFrame3(fr);
			col = hvVec3<double>((double)v.RED()/(double)scal, (double)v.GREEN()/(double)scal,(double)v.BLUE()/(double)scal);
			col = t.apply(col);
			double rr = (col.X()*rescal+offset);
			if (rr<0.0) rr=0.0; else if (rr>1.0) rr=1.0;
			double gg = (col.Y()*rescal+offset);
			if (gg<0.0) gg=0.0; else if (gg>1.0) gg=1.0;
			double bb = (col.Z()*rescal+offset);
			if (bb<0.0) bb=0.0; else if (bb>1.0) bb=1.0;
			return hvColRGB<T>((T)(rr*(double)scal),(T)(gg*(double)scal),(T)(bb*(double)scal));
	}

	void apply(T scal, hvPictRGB<unsigned char> &pi, std::vector<hvFrame3<double> > &lfr, double offset, double rescal);
	void applyInverse(T scal, hvPictRGB<unsigned char> &pi, std::vector<hvFrame3<double> > &lfr, double offset, double rescal);
	void update(double ww, T scal, hvPictRGB<unsigned char> &pi, std::vector<hvFrame3<double> > &lfr, double offset, double rescal);
	void plsr(T scal, const hvPictRGB<unsigned char> &pi, std::vector<hvFrame3<double> > &lfr) const;
	hvVec3<double> plsrmean(const hvPictRGB<unsigned char> &pi, const std::vector<hvFrame3<double> > &lfr, std::vector<double> &lvar) const;

	void savePPM(FILE *fd)
	{
		int i,j;
		hvColRGB<unsigned char> co;
		unsigned char v;

		fprintf(fd,"P6\n");
		fprintf(fd,"%d %d\n",this->sizeX(),this->sizeY());
		fprintf(fd,"255\n");
		for (i=0; i<this->sizeY(); i++)
		for (j=0; j<this->sizeX(); j++)
			{
			co = std::vector<hvColRGB<unsigned char> >::at(hvArray2< T >::get(j,this->sizeY()-i-1));
			v = (unsigned char)(co.RED());
			fwrite(&v,1,sizeof(unsigned char),fd);
			v = (unsigned char)(co.GREEN());
			fwrite(&v,1,sizeof(unsigned char),fd);
			v = (unsigned char)(co.BLUE());
			fwrite(&v,1,sizeof(unsigned char),fd);
			}
	}

	void save(FILE *fd)
	{
		int i,j;
		unsigned char v;

		fprintf(fd,"PQ\n");
		fprintf(fd,"%d %d %d\n", std::vector<hvColRGB<T> >::size(),this->sizeX(), this->sizeY());
		fprintf(fd,"255\n");
		for (i=0; i<std::vector<hvColRGB<T> >::size(); i++)
		{
			hvColRGB<T> co;
			co = std::vector<hvColRGB<T> >::at(i);
			v = (unsigned char)(co.RED());
			fwrite(&v,1,sizeof(unsigned char),fd);
			v = (unsigned char)(co.GREEN());
			fwrite(&v,1,sizeof(unsigned char),fd);
			v = (unsigned char)(co.BLUE());
			fwrite(&v,1,sizeof(unsigned char),fd);
		}
		for (i=0; i<this->sizeY(); i++)
		for (j=0; j<this->sizeX(); j++)
			{
			v = (unsigned char)hvArray2< T >::get(j,this->sizeY()-i-1);
			fwrite(&v,1,sizeof(unsigned char),fd);
			}
	}

	void load(FILE *fd);

};


}

#endif // !efined(AFX_QPICTURE_H__098453F0_1C38_49E9_A6F4_AABF90AA55E8__INCLUDED_)
