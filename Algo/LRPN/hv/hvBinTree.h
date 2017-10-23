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



// bintree.h: interface for the bintree class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BINTREE_H__1B175543_6EB8_40BB_A497_0932B70C00BE__INCLUDED_)
#define AFX_BINTREE_H__1B175543_6EB8_40BB_A497_0932B70C00BE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "hvError.h"
#include "hvNTree.h"

namespace hview {

template <class T> class hvBinTree : public hvNTree<T, 2> 
{

protected:

	hvBinTree<T>(T d): hvNTree<T, 2>(2) { hvNTree<T, 2>::data = d; }
	hvBinTree<T>(T d, hvBinTree<T> *l, hvBinTree<T> *r) { hvNTree<T, 2>::data = d; hvNTree<T, 2>::node[0]=l; hvNTree<T, 2>::node[1]=r; }

public:

	// no copy, no affectation
	hvBinTree<T>(const hvBinTree<T> &a)
	{
		hvFatal("No temporary creation of hvBinTree!");
	}
	hvBinTree<T> &operator=(hvBinTree<T> &a)
	{
		hvFatal("No affectation of hvBinTree!");
	}

	static hvBinTree<T> *createBinTree() { return 0; }
	static bool isVoid(hvBinTree<T> *t) { return t==0; }
	static bool isLeaf(hvBinTree<T> *t) { return t->node[0]==0 && t->node[1]==0; }

	static hvBinTree<T> *createBinTree(T d, hvBinTree<T> *l, hvBinTree<T> *r)
	{
		hvBinTree<T> *n=new hvBinTree<T>(d, l, r);
		return n;
	}
	static hvBinTree<T> *createBinTree(T d)
	{
		hvBinTree<T> *n=new hvBinTree<T>(d);
		return n;
	}

	static hvBinTree<T> *getLeft(hvBinTree<T> *t)
	{
		if (t==0) { hvFatal("hvBinTree is void!"); }
		return (hvBinTree<T> *)t->node[0];
	}
	static hvBinTree<T> *getRight(hvBinTree<T> *t)
	{
		if (t==0) { hvFatal("hvBinTree is void!"); }
		return (hvBinTree<T> *)t->node[1];
	}
	static void setLeft(hvBinTree<T> *t, hvBinTree<T> *l)
	{
		if (t==0) { hvFatal("hvBinTree is void!"); }
		t->node[0] = l;
	}
	static void setRight(hvBinTree<T> *t, hvBinTree<T> *r)
	{
		if (t==0) { hvFatal("hvBinTree is void!"); }
		t->node[1] = r;
	}

	static int height(hvBinTree<T> *t)
	{
		if (t==0) return 0;
		int nl = hvBinTree<T>::height((hvBinTree<T> *)t->node[0]);
		int nr = hvBinTree<T>::height((hvBinTree<T> *)t->node[1]);
		if (nl>nr) return nl+1; else return nr+1;
	}
	
	static void destroy(hvBinTree<T> *t)
	{
		if (t==0) return;
		if (t->node[0]!=0) hvBinTree<T>::destroy((hvBinTree<T> *)t->node[0]);
		if (t->node[1]!=0) hvBinTree<T>::destroy((hvBinTree<T> *)t->node[1]);
		delete t;
	}

};

}
#endif // !defined(AFX_BINTREE_H__1B175543_6EB8_40BB_A497_0932B70C00BE__INCLUDED_)
