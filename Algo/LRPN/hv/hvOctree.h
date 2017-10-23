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
