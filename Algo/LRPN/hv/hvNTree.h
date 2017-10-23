// N-tuple Tree.h: interface for the tree class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_NTREE_H__07428FA9_06A1_45AC_A60F_C632153045C6__INCLUDED_)
#define AFX_NTREE_H__07428FA9_06A1_45AC_A60F_C632153045C6__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "hvError.h"

namespace hview {

template <class T, int N> class hvNTree  
{
protected:
	T				data;
	hvNTree<T,N>	*node[N];
	
	hvNTree(T d)
	{
		for (int i=0; i<N; i++) node[i]=0;
		data = d;
	}

public:
	static hvNTree<T,N> *createNTree() { return 0; }
	static hvNTree<T,N> *createNTree(T d) { return new hvNTree<T,N>(d); }

	static bool isVoid(hvNTree<T,N> *t) { return t==0; }
	static bool isLeaf(hvNTree<T,N> *t)
	{
		if (t==0) { hvFatal("tree is void!"); }
		for (int i=0; i<N; i++) if (t->node[i]!=0) return false;
		return true;
	}

	static T getValue(hvNTree<T,N> *t)
	{
		if (t==0) { hvFatal("tree is void!"); }
		return t->data;
	}
	static void update(hvNTree<T,N> *t, const T &d)
	{
		if (t==0) { hvFatal("hvNTree is void!"); }
		t->data = d;
	}

	static hvNTree<T,N> *getChildNode(hvNTree<T,N> *t, int i)
	{
		if (t==0) { hvFatal("tree is void!"); }
		if (i<0 || i>=N) { printf("node index %d out of range!",i); hvFatal("out of range!"); }
		return t->node[i];
	}

	static void setChildNode(hvNTree<T,N> *t, int i, hvNTree<T,N> *n)
	{
		if (t==0) { hvFatal("tree is void!"); }
		if (i<0 || i>=N) { printf("node index %d out of range!",i); hvFatal("out of range!"); }
		t->node[i]=n;
	}

	static void destroy(hvNTree<T,N> *t)
	{
		if (t==0) return;
		for (int i=0; i<N; i++) if (t->node[i]!=0) destroy(t->node[i]);
		delete t;
	}

	static int countAllLeafs(hvNTree<T,N> *t)
	{
		if (t==0) return 0;
		if (hvNTree<T,N>::isLeaf(t)) return 1;
		int nb=0;
		for (int i=0; i<N; i++) if (t->node[i]!=0) nb += hvNTree<T,N>::countAllLeafs(t->node[i]);
		return nb;
	}

	static int height(hvNTree<T,N> *t)
	{
		if (t==0) return 0;
		int nb=0, nn;
		for (int i=0; i<N; i++) 
			if (t->node[i]!=0) 
			{ 
				nn = hvNTree<T,N>::height(t->node[i]);
				if (nn>nb) nb=nn;
			}
		return nb+1;
	}

};

}

#endif // !defined(AFX_NTREE_H__07428FA9_06A1_45AC_A60F_C632153045C6__INCLUDED_)
