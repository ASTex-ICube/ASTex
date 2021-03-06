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



#include "error_buffer.h"

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iomanip>

using namespace ASTex;



ErrorCutPathBuffer::ErrorCutPathBuffer(int w, int h):
		width_(w),height_(h)
{
	data_ = new double[w*h];
}

ErrorCutPathBuffer::~ErrorCutPathBuffer()
{
	delete[] data_;
}


int ErrorCutPathBuffer::minOfRow(int r)
{
	double *ptr = data_+(r*width_);
	double m = *ptr++;
	int mi=0;
	for (int i=1;i<width_;++i)
	{
		if (*ptr < m)
		{
			m = *ptr;
			mi = i;
		}
		++ptr;
	}
	return mi;
}

int ErrorCutPathBuffer::minOfColumn(int c)
{
	double *ptr = data_+c;
	double m = *ptr;
	int mi=0;
	for (int i=1;i<height_;++i)
	{
		ptr += width_;
		if (*ptr < m)
		{
			m = *ptr;
			mi = i;
		}
	}
	return mi;
}

int ErrorCutPathBuffer::minOfRow3(int i, int j)
{
	double m = (*this)(i,j);
	int mi = i;
	if (i>0)
	{
		double v = (*this)(i-1,j);
		if (v<m)
		{
			m = v;
			mi = i-1;
		}
	}
	i++;
	if (i<width_)
	{
		double v = (*this)(i,j);
		if (v<m)
		{
			mi = i;
		}
	}
	return mi;
}

int ErrorCutPathBuffer::minOfColumn3(int i, int j)
{
	double m = (*this)(i,j);
	int mj = j;
	if (j>0)
	{
		double v = (*this)(i,j-1);
		if (v<m)
		{
			m = v;
			mj = j-1;
		}
	}
	j++;
	if (j<height_)
	{
		double v = (*this)(i,j);
		if (v<m)
		{
			mj = j;
		}
	}
	return mj;
}



void ContentLabel::load(const std::string& bname)
{
	label_images_.reserve(32);
	std::ifstream* f = nullptr;
	int i=1;
	do
	{
		std::stringstream ss;
		ss << "_cl_"<<i<<"_uhd.png";
		std::string name = bname+ss.str();
		if (f)
		{
			f->close();
			delete f;
		}
		f = new std::ifstream(name.c_str());
		std::cout << "trying "<< name << std::endl;
		if (f->good())
		{
			std::cout << "loading "<< name << std::endl;
			label_images_.resize(i);
			label_images_.back().load(name);
		}
		++i;
	} while (f->good());

	f->close();
	delete f;
}


int ContentLabel::levels_equality(const Index& p, const Index& q)
{
	int nb = label_images_.size();

	for (int i=0; i<nb;++i)
	{
		if (! has_same_label(p,q,i))
			return i;
	}

	return nb;
}

