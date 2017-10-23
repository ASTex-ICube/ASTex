// tree.h: interface for the tree class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_OCTREE_H__07428FA9_06A1_45AC_A60F_C632153045C6__INCLUDED_)
#define AFX_OCTREE_H__07428FA9_06A1_45AC_A60F_C632153045C6__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "hvNTree.h"

namespace hview {

template <class T> class hvOctree : public hvNTree<T,8>  
{
protected:
	hvOctree(T d):hvNTree<T,8>(d) { };
public:

	static hvOctree<T> *createOctree() { return 0; }
	static hvOctree<T> *createOctree(T d) { return new hvOctree<T>(d); }
	static hvOctree<T> *getChildNode(hvOctree<T> *t, int i)
	{
		if (t==0) { hvFatal("octree is void!"); }
		if (i<0 || i>=8) { hvFatal("octree node index out of range!"); }
		return (hvOctree<T> *)t->node[i];
	}

	static void setChildNode(hvOctree<T> *t, int i, hvOctree<T> *n)
	{
		if (t==0) { hvFatal("tree is void!"); }
		if (i<0 || i>=8) { hvFatal("octree node index out of range!"); }
		t->node[i]=n;
	}
};

}

#endif // !defined(AFX_OCTREE_H__07428FA9_06A1_45AC_A60F_C632153045C6__INCLUDED_)
