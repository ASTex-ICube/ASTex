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



#ifndef __ASTEX_GRAPHCUT__
#define __ASTEX_GRAPHCUT__

#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include <iostream>
#include <vector>
#include <ctime>

#include "maxflow/graph.h"

using namespace ASTex;

// After:
//
// Vivek Kwatra, Arno Schödl, Irfan Essa, Greg Turk and Aaron Bobick.
// Graphcut Textures: Image and Video Synthesis Using Graph Cuts.
// ACM Transactions on Graphics (Proc. SIGGRAPH), 22(3):227-286, 2003.
//
// Xuejie Qin and Yee-Hong Yang. Theoretical Analysis of Graphcut Textures.
// Technical Report 05-26, Department of Computer Science, University of Alberta, 2005.
//
// Kwatra et al . do not prove that solving the min-cut problem is equivalent to obtaining
// optimal labeling of the graph including their definition of seam nodes. This gap is filled by
// X. Qin and Y.-H. Yang, "Theoretical Analysis of Graphcut Textures", Technical Report 05-26,
// Department of Computer Science, University of Alberta, 2005, who gave the mathematical proof of
// the optimality of label assignments using the graph-cut formulation described by Kwatra et al.
//

// Data-structure for seam management in the output
// texture, to be used in conjunction with min s-t cut
//
//   X--right--o
//   |
// bottom
//   |
//   o
//
class GlobalNode
{
public:
	GlobalNode() {
		empty_ = true; rightCost = 0.0; bottomCost = 0.0; seamRight = false; seamBottom = false; maxFlow = 0.0;
		newSeam = false;
		csbc = 0;
		csrc = 0;
	}

	void setColor(const itk::RGBPixel<unsigned char>& c) { color = c; empty_ = false; }
	itk::RGBPixel<unsigned char> getColor() const { return color; }

	bool empty() const { return empty_; }
	bool notEmpty() const { return !empty_; }

	void setRightCost(double rc) { rightCost = rc; csrc++;  }
	int getCountSetRightCost() { return csrc; }
	double getRightCost() const { return rightCost; }
	void setBottomCost(double bc) { bottomCost = bc; csbc++; }
	int getCountSetBottomCost() { return csbc;  }
	double getBottomCost() const { return bottomCost; }

	void unsetSeams() { seamRight = false; seamBottom = false; }
	void unsetSeamRight() { seamRight = false; }
	void setSeamRight(double mf) { seamRight = true; maxFlow = mf; }
	bool onSeamRight() const { return seamRight; }
	void unsetSeamBottom() { seamBottom = false; }
	void setSeamBottom(double mf) { seamBottom = true; maxFlow = mf; }
	bool onSeamBottom() const { return seamBottom; }

	void setColorOtherPatch(const itk::RGBPixel<unsigned char>& c) { colorOtherPatch = c; }
	itk::RGBPixel<unsigned char> getColorOtherPatch() const { return colorOtherPatch; }

	void setGradXOtherPatch(uint8_t g) { gradXOtherPatch = g; }
	uint8_t getGradXOtherPatch() const { return gradXOtherPatch; }

	void setGradYOtherPatch(uint8_t g) { gradYOtherPatch = g; }
	uint8_t getGradYOtherPatch() const { return gradYOtherPatch; }

	bool onNewSeam() const { return newSeam; }
	void setNewSeam() { newSeam = true; }
	void unsetNewSeam() { newSeam = false; }

private:
	itk::RGBPixel<unsigned char> color;
	bool empty_;
	double rightCost;
	double bottomCost;
	bool seamRight;
	bool seamBottom;
	double maxFlow;
	itk::RGBPixel<unsigned char> colorOtherPatch;
	uint8_t gradXOtherPatch;
	uint8_t gradYOtherPatch;

	bool newSeam;

	int csbc;
	int csrc;
};

// Data-structure for the management of old seam nodes
class SeamNode
{
public:
	SeamNode(int s, int e, double c1, double c2, double c3, int o) {
		start = s;
		end = e;
		capacity1 = c1;
		capacity2 = c2;
		capacity3 = c3;
		orientation = o;
	}

	void setSeam(int s) { seam = s; }
	int getStart() const { return start; }
	int getEnd() const { return end; }
	int getSeam() const { return seam; }
	double getCapacity1() const { return capacity1; }
	double getCapacity2() const { return capacity2; }
	double getCapacity3() const { return capacity3; }
	int getOrientation() const { return orientation; }

private:
	int start;
	int end;
	int seam;
	double capacity1;
	double capacity2;
	double capacity3;
	int orientation; // right: 0, bottom: 1
};

// Main class for graph-cut texture synthesis
class GCTexture
{
public:
	GCTexture();
	virtual ~GCTexture();
	void initialize(const ImageRGBu8& iImage, const std::string& tName, int o_w, int o_h);
	void topDownRandomFilling(int bandX, int bandY, bool blending);
	void topDownEntireMatchingMinErrorFilling(int bandX, int bandY, bool blending);
	void topDownEntireMatchingRandomnessFilling(int bandX, int bandY, double k, bool blending);
	void randomRefinement(int bandX, int bandY, int iter, bool blending);
	void seamsErrorCenterRefinement(int maxIter, int radius, bool blending);
	void seamsErrorRefinement(int maxIter, int radius, bool blending);
	void seamsErrorFullPatchRefinement(int maxIter, int radius, bool blending);
	void seamsErrorSubPatchRefinementA(int maxIter, int radius, bool blending);
	void seamsErrorSubPatchRefinement(int maxIter, int radius, bool blending);
	void revealSeams();
	void revealSeamsError();
	void revealSeamsMaxError(int radius);
	void writeOutputImage(const std::string& fileName, const std::string& fileNameS);

	inline ImageRGBu8 getOutputImage() const { return outputImage; }

	inline void setSeamSize(int n) { seamSize = n; }

	inline int getSeamSize() const { return seamSize; }

	inline void setBlendingParams(int r, int n) { blradius = r; bln = n; }

	inline int getBlendingRadius() const { return blradius; }

	inline int getBlendingN() const { return bln; }

	inline void setSSDStep(int s) { ssdStep = s; }

	inline int getSSDStep() const { return ssdStep; }

	inline void setBorderSize(int s) { borderSize = s; }

	inline int getBorderSize() const { return borderSize; }

private:
	std::string textureName;

	ImageRGBu8 inputImage;
	ImageRGBu8 croppedInputImage;
	ImageGrayu8 inputImageGX;
	ImageGrayu8 croppedInputImageGX;
	ImageGrayu8 inputImageGY;
	ImageGrayu8 croppedInputImageGY;
	ImageRGBu8 outputImage;

	typedef Graph<double, double, double> GraphType;
	GraphType *g;
	double infiniteCap;
	double minCap;

	int input_w;
	int input_h;
	int croppedInput_w;
	int croppedInput_h;
	int output_w;
	int output_h;
	GlobalNode* globalNode;
	std::vector<SeamNode> seamNode;

	int patchCount;
	int seamSize;

	bool isInitialized;
	bool isFilled;

	int ssdStep;
	int maxErrNodeNbGlobal;
	int blradius;
	int bln;

	int borderSize;

private:
	void fillOutputImage();
	void writeOutputImagePatch();

	int getNodeNbGlobal(int x, int y, int i, int j);
	int getNodeNbLocal(int i, int j);
	void getPixel(int nodeNb, int& i, int& j);
	void getPixelGlobal(int nodeNbGlobal, int& x, int& y);

	int insertPatch(int x, int y, int px, int py, int sx, int sy, bool filling, bool blending, int radius, int inputX, int inputY, int tX, int tY);
	int insertPatch(int x, int y, bool filling, bool blending, int radius = 0, int inputX = 0, int inputY = 0, int tX = 0, int tY = 0);

	double computeSeamsAverageError();
	double updateSeamsMaxError();
	double computeSSD(int x, int y);
	double computeLocalSSD(int x, int y, int inputX, int inputY, int sizeX, int sizeY);
};

#endif
