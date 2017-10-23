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



// LinearTransform2.h: interface for the LinearTransform2 class.
//
// Define a linear transform in 2D space (scaling, translation, rotation, etc.) 
//
// LinearTransform2 derives from matrix 3x3 (homogenous coordinates)
// A main operation consists in applying the transform to 2D points
//
// By JMD 18/1/08
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LINEARTRANSFORM2_H__33F894B9_F2FB_4902_97D2_D69639938055__INCLUDED_)
#define AFX_LINEARTRANSFORM2_H__33F894B9_F2FB_4902_97D2_D69639938055__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "hvMat3.h"
#include "hviTransform.h"

namespace hview {

template <class T> class hvLinearTransform2 : public hvMat3<T>, public hviTransform<hvLinearTransform2<T>, hvVec2<T> > 
{
public:
	// constructor: defines identity transform
	hvLinearTransform2() { }

	template <class X> void cast(const hvLinearTransform2<X> &v) { hvMat3<T>::i.cast<X>(v.I()); hvMat3<T>::j.cast<X>(v.J()); hvMat3<T>::k.cast<X>(v.K());  }

	void setIdentity() { hvMat3<T>::i=hvVec3<T>(1,0,0); hvMat3<T>::j=hvVec3<T>(0,1,0); hvMat3<T>::k=hvVec3<T>(0,0,1); }

	// List of transforms: each is composed with the current one

	// shift of vector v
	void translation(const hvVec2<T> &v)
	{
		hvMat3<T> m(hvVec3<T>(1,0,0), hvVec3<T>(0,1,0), hvVec3<T>(v.X(),v.Y(),1));
		this->mult(*this, m);
	}

	// scale the three components by factors given by v
	void scale(const hvVec2<T> &v)
	{
		hvMat3<T> m(hvVec3<T>(v.X(),0,0), hvVec3<T>(0,v.Y(),0), hvVec3<T>(0,0,1));
		this->mult(*this, m);
	}


	// rotation by angle a in radiant
	void rot(T a)
	{
		hvMat3<T> m(	hvVec3<T>(cos(a), sin(a), 0.0),
						hvVec3<T>(-sin(a), cos(a), 0.0),
						hvVec3<T>(0.0, 0.0, 1.0) );
		this->mult(*this, m);
	}


	// compose two transforms with each other
	void compose(const hvLinearTransform2<T> &a, const hvLinearTransform2<T> &b) { this->mult(a,b); }

	// apply a transform to a 2D point, result is a 2D point
	hvVec2<T> apply(const hvVec2<T> &v) const { hvVec3<T> u(v.X(), v.Y(), 1); hvVec3<T> r=this->mult(u); return hvVec2<T>(r.X()/r.Z(),r.Y()/r.Z()); }
	hvVec3<T> apply(const hvVec3<T> &v) const { return this->mult(v); }
	
	// apply transform to a vector, result is a new vector
	// this is the same as for point, but no translation is applied
	// NOTE: the norm of resulting vector is not preserved
	hvVec2<T> applyVec2(const hvVec2<T> &v) const
	{
		hvMat3<T> m = *this;
		hvMat2<T> mm= (hvMat2<T>)m;
		return mm.mult(v);
	}
	

};

}

#endif // !defined(AFX_TRANSFORM2_H__33F894B9_F2FB_4902_97D2_D69639938055__INCLUDED_)
