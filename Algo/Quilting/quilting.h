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



#ifndef __ASTEX_QUILTING__
#define __ASTEX_QUILTING__

#include <ASTex/image_rgb.h>
#include <iostream>
#include <map>
#include <vector>


#include <ASTex/image_indexed.h>
#include "error_buffer.h"

using namespace ASTex;

class Quilting
{
public:
	using ImageType = ImageRGBd;

protected:

	int tile_width_;
	int overlap_width_;

	Size overlap_h_sz_;
	Size overlap_v_sz_;
	Size tile_sz_;

	ImageType img_orig_;

	ImageType img_generated_;

	std::vector<Index> bufferPositions_;

	std::vector<float> bufferErrors_;


	class MultiMin
	{
	public:
		std::multimap<double,Index> mpos_;
		std::size_t nbMaxElts_;
		double max_of_mins_;

		inline MultiMin():
			nbMaxElts_(8), max_of_mins_(0.0) {}

		inline MultiMin(int nbElts):
			nbMaxElts_(nbElts), max_of_mins_(0.0) {}

		inline void reset()
		{
			max_of_mins_ = 0.0;
			mpos_.clear();
		}

		inline void try_add(const Index& pos, double v)
		{
			// insert if possible
			if ((v<max_of_mins_) || (mpos_.size() < nbMaxElts_))
				mpos_.insert(std::pair<double,Index>(v,pos));

			// remove last if too many
			if ((mpos_.size() > nbMaxElts_))
			{
				std::multimap<double,Index>::iterator it = mpos_.end();
				--it;
				mpos_.erase(it);
			}

			// update max
			max_of_mins_ = mpos_.rbegin()->first;
		}
	};

	MultiMin multi_mins_;


	static int random_int(int min, int max);


	/**
	 * @brief init
	 * @param tw Tile width
	 * @param ovl overlap width
	 * @param rc number of fitting tile kept before choosing randomly one of them
	 * @param nbrt number of tiles randomly choosen on original image tested for boundary fitting (0: all possiblee tile tested)
	 */
	void init(int tw, int ovl, int rc, int nbrt);

	/**
	 * @brief compute error between 2 regions of the 2 images (origin & generated)
	 * @param or_reg region in original image
	 * @param gen_reg region in generated image
	 * @return error in [0,1]
	 */
	double computeErrorOverlap(const Region& or_reg, const Region& gen_reg);

	/**
	 * @brief compute the error of boundary overlap matching between a tile in original image and the generated image
	 * @param or_pos position of tile in original image
	 * @param gen_pos position of insertion in generated image
	 * @return error in [0,1]
	 */
	double computeErrorTile(const Index& or_pos, const Index& gen_pos);

	/**
	 * @brief find the best tiles (minimal overlap error)
	 * @param gen_pos position of tile to insert in generated image
	 * @param error_max
	 * @return position of choosen tile in original image
	 */
	Index bestFittingTile(const Index& gen_pos, double error_max = 2.0 );

	/**
	 * @brief simple copy
	 * @param from region in original image
	 * @param to region in generated image
	 */
	void copy(const Region& from, const Region& to);


	/**
	 * @brief compute minimal errot cutting path for vertical overlapping
	 * @param or_pos position in original image
	 * @param gen_pos position in generated image
	 * @param minPos vector of position of minimal cut path
	 */
	void verticalPathCut(const Index& or_pos, const Index& gen_pos, std::vector<int>& minPos);

	/**
	 * @brief compute minimal errot cutting path for horizontal overlapping
	 * @param or_pos position in original image
	 * @param gen_pos position in generated image
	 * @param minPos vector of position of minimal cut path
	 */
	void horizontalPathCut(const Index& or_pos, const Index& gen_pos, std::vector<int>& minPos);

	/**
	 * @brief copy region of image with respect of cutting paths
	 * @param from region in original image
	 * @param to region in generated image
	 * @param leftMinCut cutting path of left overlapping
	 * @param upMinCut cutting path of up overlapping
	 */
	void copyPathCut(const Region& from, const Region& to, const std::vector<int>& leftMinCut,const std::vector<int>& upMinCut);


	void copyPathCut_BlendDisk(const Region& from, const Region& to, const std::vector<int>& leftMinCut,const std::vector<int>& upMinCut, int rad);


	/**
	 * @brief compute path slope
	 * @warning need to be inverted for vertical path
	 * @param pos position (y = f(x))
	 * @param slopes computed slopes
	 */
	void pathSlope(const std::vector<int>& pos, std::vector<float> slopes);


public:
	/**
	 * @brief load input image
	 * @param imgName
	 */
	void loadInput(const std::string& imgName);

	/**
	 * @brief save output image
	 * @param imgName
	 */
	void saveOutput(const std::string& imgName);

	/**
	 * @brief compute image by random choosen tile
	 * @param width width of generated image
	 * @param height height of generated image
	 * @param tw tile width (tiles are squares)
	 */
	void computeRandom(int width, int height, int tw);

	/**
	 * @brief compute image by choosing the tiles that fit boundaries
	 * @param width width of generated image
	 * @param height height of generated image
	 * @param tw tile width (tiles are squares)
	 * @param ovl overlap width
	 * @param rc number of best fitting tiles to keep
	 * @param nbrt number of tiles to test
	 */
	void computeFittingTiles(int width, int height, int tw,  int ovl = 6, int rc = 5, int nbrt = 4096);


	/**
	 * @brief compute image by choosing the tiles that fit boundaries with error cutting paths
	 * @param width width of generated image
	 * @param height height of generated image
	 * @param tw tile width (tiles are squares)
	 * @param ovl overlap width
	 * @param rc number of best fitting tiles to keep
	 * @param nbrt number of tiles to test
	 */
	void computeFittingTilesPathCut(int width, int height, int tw, int ovl = 6, int rc = 5, int nbrt = 4096);


};



#endif

