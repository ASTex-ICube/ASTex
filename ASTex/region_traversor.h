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



#ifndef __ASTEX_WALKINTHEIMAGE__
#define __ASTEX_WALKINTHEIMAGE__

#include <ASTex/image_gray.h>

/** \file
 * \brief Kind of an iterator which returns square regions in the image.
*/

/**
 * \brief namespace ASTex
 */
namespace ASTex
{

/**
 * \class RegionTraversor
 * \brief Kind of an iterator which returns square regions in the image.
 * The region is [x;x+T-1]*[y;y+T-1] where
 * T is the largest power of two < input width and height.
 * The top left corner (x,y) takes all valid (i*step, j*step) values.
 */
class RegionTraversor
{
public :
	/**
	 * \brief constructor
	 * \param [in] imgSz an image size
	 * \param [in] step by which the region is shifted at each iteration
	 */
	RegionTraversor(const itk::Size<2>& imgSz, int step)
	{
		m_W = imgSz[0];
		m_H = imgSz[1];

		uint64_t M = std::min(m_W,m_H);
		m_T = 1;
		while ((2*m_T) < M)
		{
			m_T*=2;
		}

		m_region.SetIndex({{0l,0l}});
		m_region.SetSize({{m_T,m_T}});
		m_step=step;
	}

	/** \brief tests if the abscissa is valid */
	inline bool isXValid () const { return ( (getX()+m_T) <= m_W ); }

	/** \brief tests if the ordinate is valid */
	inline bool isYValid () const { return ( (getY()+m_T) <= m_H ); }

	/** \brief tests if the position is valid */
	inline bool isValid() { return isXValid() && isYValid(); }

	/** \brief moves the region of interest one step forward */
	inline RegionTraversor& operator ++ ()
	{
//		long& Y = m_region.GetModifiableIndex()[1];
		m_region.GetModifiableIndex()[1] += m_step;
		if (!isYValid())
		{
			m_region.GetModifiableIndex()[1] = 0;
			m_region.GetModifiableIndex()[0] += m_step;
		}
		return *this;
	}

	inline void init()
	{
		m_region.SetIndex({{0l,0l}});
	}

	/** \brief gets current abscissa */
	inline int getX() const { return m_region.GetIndex()[0]; }

	/** \brief gets current ordinate */
	inline int getY() const { return m_region.GetIndex()[1];}

	/** \brief gets region size */
	inline int getSize() const {return m_T;}

	/** \brief gets current region */
	inline const Region& operator * () const { return m_region; }

private :
	/** \brief gets current abscissa */
//	inline long& X() {return m_region.GetModifiableIndex()[0];}

	/** \brief gets current ordinate */
//	inline long& Y() {return m_region.GetModifiableIndex()[1];}

	uint64_t m_W;
	uint64_t m_H;
	uint64_t m_T;
	uint64_t m_step;
	Region m_region;
};

} // namespace ASTex

#endif
