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



// hviTransform.h: interface (abstract object) for a Transformation.
//
// JMD 1/01/2007
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TRANSFORM_H__B6AC0A32_75EF_428E_BC10_6219F619FA29__INCLUDED_)
#define AFX_TRANSFORM_H__B6AC0A32_75EF_428E_BC10_6219F619FA29__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

////////////////////////////////////////////
/*
A Transformation T consists of a conversion of an object U to another object of the same data type U
The interface implements following issues: 
the application of the transform to an object U
and the composition of two transforms
Neutral transform is called Identity.
*/ 
////////////////////////////////////////////
////////////////////////////////////////////

namespace hview {

template <class T, class U> class hviTransform  
{
public:
	// apply a transform to an object of type T, returns another object of type T
	virtual U apply(const U  &a) const =0;

	// compose two transforms
	virtual void compose(const T &a, const T &b) =0;

	// set identity
	virtual void setIdentity() =0;
};

}
#endif // !defined(AFX_ALGEBRA_H__B6AC0A32_75EF_428E_BC10_6219F619FA29__INCLUDED_)
