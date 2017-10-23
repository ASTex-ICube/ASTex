#ifndef __ASTEX_QUILTING_INTERP_CONT__
#define __ASTEX_QUILTING_INTERP_CONT__


#include <ASTex/image_rgb.h>
#include <iostream>
#include <map>
#include <vector>

#include <ASTex/image_indexed.h>
#include "error_buffer.h"

using namespace ASTex;


class QuiltingInterpCont
{
public:
	typedef ImageRGBd ImageType;
	typedef itk::ImageRegion<2> Region;
	typedef itk::Index<2> Index;
	typedef itk::Size<2> Size;

protected:

	ContentLabel content_;

	int tile_width_;
	int overlap_width_;
	int blend_rad_;

	itk::Size<2> overlap_h_sz_;
	itk::Size<2> overlap_v_sz_;
	itk::Size<2> tile_sz_;

	ImageType img_orig_;

	ImageType img_generated_;

	ImageGrayu32 img_gen_idx_;

	ImageType img_dbg_interp;

	Index pos_first_tile;

	bool label_used_for_choice_;

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

	inline uint32_t pos_to_idx(const Index& p)
	{
		 return p[0] | p[1]<<16;
	}

	inline Index idx_to_pos(uint32_t idx)
	{
		Index P;
		P[0] = idx & 0xffff;
		P[1] = idx >> 16;
		return P;
	}

	inline bool is_an_index(uint32_t idx)
	{
		return idx != 0xffffffff;
	}

//	inline label_origin(const Index& P, uint32_t lev)
//	{
//		return content_.ee
//	}


	/**
	 * @brief random_int
	 * @param min
	 * @param max
	 * @return
	 */
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


	double computeErrorOverlapLabels(const Region& or_reg, const Region& gen_reg);

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



	void computeErrorForCutting(const Region& or_reg, const Region& gen_reg, ErrorCutPathBuffer& e);

	void computeErrorForCuttingLabels(const Region& or_reg, const Region& gen_reg, ErrorCutPathBuffer& e);

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

//	/**
//	 * @brief copy region of image with respect of cutting paths
//	 * @param from region in original image
//	 * @param to region in generated image
//	 * @param leftMinCut cutting path of left overlapping
//	 * @param upMinCut cutting path of up overlapping
//	 */
//	void copyPathCut(const Region& from, const Region& to, const std::vector<int>& leftMinCut,const std::vector<int>& upMinCut);


	void copyPathCut_BlendDisk(const Region& from, const Region& to, const std::vector<int>& leftMinCut,const std::vector<int>& upMinCut);

	void saveDebugOverlap(const Index& P, const Index& Q, int k);

public:
	QuiltingInterpCont();

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

//	/**
//	 * @brief compute image by random choosen tile
//	 * @param width width of generated image
//	 * @param height height of generated image
//	 * @param tw tile width (tiles are squares)
//	 */
//	void computeRandom(int width, int height, int tw);

//	/**
//	 * @brief compute image by choosing the tiles that fit boundaries
//	 * @param width width of generated image
//	 * @param height height of generated image
//	 * @param tw tile width (tiles are squares)
//	 * @param ovl overlap width
//	 * @param rc number of best fitting tiles to keep
//	 * @param nbrt number of tiles to test
//	 */
//	void computeFittingTiles(int width, int height, int tw,  int ovl = 6, int rc = 5, int nbrt = 4096);

	inline void setFirstTile(int x, int y)
	{
		pos_first_tile[0] = x;
		pos_first_tile[1] = y;
	}

	/**
	 * @brief compute image by choosing the tiles that fit boundaries with error cutting paths
	 * @param width width of generated image
	 * @param height height of generated image
	 * @param tw tile width (tiles are squares)
	 * @param ovl overlap width
	 * @param bl_rad blending radius
	 * @param rc number of best fitting tiles to keep
	 * @param nbrt number of tiles to test
	 */
	void computeFittingTilesPathCutInterp(int width, int height, int tw, int ovl = 6, int bl_rad = 5, int rc = 5, int nbrt = 4096);


};



#endif

