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



// line3.h: interface for the line3 class.
//
// line3 defines a line in a 3D space: L= O+Dl, where O is origin and D direction
// 
// By JMD 10/8/06
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LINE3_H__EC87CD56_6E08_4390_B876_CEC6D44EEA86__INCLUDED_)
#define AFX_LINE3_H__EC87CD56_6E08_4390_B876_CEC6D44EEA86__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "hvVec3.h"

namespace hview {

	template <class T> class line3
	{
	protected:
		hvVec3<T>		o;	// origin
		hvVec3<T>		d;	// direction
	public:
		// constructors
		line3<T>() { }

		// defines a line by a point pt (origin) and a vector dir (direction)
		// dir must be normalized
		line3<T>(const hvVec3<T> &pt, const hvVec3<T> &dir) { o = pt; d = dir; }

		// defines a line from origin (0,0,0) in direction dir
		// dir must be normalized
		line3<T>(const hvVec3<T> &dir) : o(T(0)) { d = dir; }

		// selectors
		hvVec3<T> origin() const { return o; }
		hvVec3<T> direction() const { return d; }

		// Compute a point (result) on the line at distance l from the origin
		hvVec3<T> pointOfLine(T l) const
		{
			hvVec3<T> p; p.scale(d, l);
			p += o;
			return p;
		}

		void reverse() { d.reverse(); }

	};

}
#endif // !defined(AFX_LINE3_H__EC87CD56_6E08_4390_B876_CEC6D44EEA86__INCLUDED_)
