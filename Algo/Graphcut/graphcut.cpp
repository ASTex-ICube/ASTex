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



#include "graphcut.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <queue>
#include <limits>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "itkCastImageFilter.h"

//#define _DEBUG

using namespace ASTex;

typedef itk::RGBPixel<uint8_t > itkRGBu8;

// ---- Utility functions ----

// 2*sqrt(3)
const double two_sqrt3 = 3.4641016151377545870548926830117;

// Distance between two colors in RGB space
double distRGB(const itkRGBu8& c1, const itkRGBu8& c2)
{
	double r1 = c1.GetRed() / 255.0;
	double g1 = c1.GetGreen() / 255.0;
	double b1 = c1.GetBlue() / 255.0;

	double r2 = c2.GetRed() / 255.0;
	double g2 = c2.GetGreen() / 255.0;
	double b2 = c2.GetBlue() / 255.0;

	return sqrt((r2 - r1)*(r2 - r1) + (g2 - g1)*(g2 - g1) + (b2 - b1)*(b2 - b1));
}

// Squared distance between two colors in RGB space
double sqDistRGB(const itkRGBu8& c1, const itkRGBu8& c2)
{
	double r1 = c1.GetRed() / 255.0;
	double g1 = c1.GetGreen() / 255.0;
	double b1 = c1.GetBlue() / 255.0;

	double r2 = c2.GetRed() / 255.0;
	double g2 = c2.GetGreen() / 255.0;
	double b2 = c2.GetBlue() / 255.0;

	return ((r2 - r1)*(r2 - r1) + (g2 - g1)*(g2 - g1) + (b2 - b1)*(b2 - b1));
}

// Distance between two colors based on their luminance (result in [0, 1])
double distRGBLum(const itkRGBu8& c1, const itkRGBu8& c2)
{
	double r1 = c1.GetRed() / 255.0;
	double g1 = c1.GetGreen() / 255.0;
	double b1 = c1.GetBlue() / 255.0;

	double l1 = (0.2126 * r1 + 0.7152 * g1 + 0.0722 * b1);
	if (l1 > 1.0) l1 = 1.0;

	double r2 = c2.GetRed() / 255.0;
	double g2 = c2.GetGreen() / 255.0;
	double b2 = c2.GetBlue() / 255.0;

	double l2 = (0.2126 * r2 + 0.7152 * g2 + 0.0722 * b2);
	if (l2 > 1.0) l2 = 1.0;

	return std::abs(l2 - l1);
}

// Distance between two colors (result in [0, 1])
// - lum: if lum is false, distance is computed in RGB space,
// luminance is used otherwise
double distRGB(const itkRGBu8& c1, const itkRGBu8& c2, bool lum)
{
	if (lum)
		return distRGBLum(c1, c2);
	else
		return distRGB(c1, c2);
}

double distColors(const itkRGBu8& c1, const itkRGBu8& c2)
{
	return distRGB(c1, c2, false);
	//return distRGB(c1, c2, true);
}

// Conversion from RGB to luminance
int rgbToLum(const itkRGBu8& c)
{
	double r = c.GetRed() / 255.0;
	double g = c.GetGreen() / 255.0;
	double b = c.GetBlue() / 255.0;

	double l = (0.2126 * r + 0.7152 * g + 0.0722 * b);
	if (l > 1.0) l = 1.0;

	return (int)(255.0 * l);
}

// Compute gradient images for a given input image
// - image: image to be processed
// - imageGX: horizontal gradients
// - imageGY: vertical gradients
void computeGradientImage(const ImageRGBu8& image, ImageGrayu8& imageGX, ImageGrayu8& imageGY)
{
	typedef itk::Image< float, 2 > FloatImageType;
	typedef itk::Image< uint8_t , 2 > UnsignedCharImageType;
	typedef itk::RGBToLuminanceImageFilter< ImageRGBu8::ItkImg, UnsignedCharImageType > ToGrayFilterType;
	typedef itk::DerivativeImageFilter< UnsignedCharImageType, FloatImageType > DerivativeFilterType;	
	typedef itk::RescaleIntensityImageFilter< FloatImageType, UnsignedCharImageType > NormalizeFilterType;

	// Convert image to gray
	ToGrayFilterType::Pointer toGrayFilter = ToGrayFilterType::New();
	toGrayFilter->SetInput(image.itk());

	// Create and setup derivative filter along X
	DerivativeFilterType::Pointer derivativeFilterX = DerivativeFilterType::New();
	derivativeFilterX->SetDirection(0);
	derivativeFilterX->SetInput(toGrayFilter->GetOutput());
	NormalizeFilterType::Pointer normalizerX = NormalizeFilterType::New();
	normalizerX->SetOutputMinimum(0);
	normalizerX->SetOutputMaximum(255);
	normalizerX->SetInput(derivativeFilterX->GetOutput());
	normalizerX->Update();
	//imageGX.itk()->CopyInformation(normalizerX->GetOutput());
	imageGX.itk() = ImageGrayu8(normalizerX->GetOutput()).itk();


	// Create and setup derivative filter along Y
	DerivativeFilterType::Pointer derivativeFilterY = DerivativeFilterType::New();
	derivativeFilterY->SetDirection(1);
	derivativeFilterY->SetInput(toGrayFilter->GetOutput());
	NormalizeFilterType::Pointer normalizerY = NormalizeFilterType::New();
	normalizerY->SetOutputMinimum(0);
	normalizerY->SetOutputMaximum(255);
	normalizerY->SetInput(derivativeFilterY->GetOutput());
	normalizerY->Update();
	//imageGY.itk()->CopyInformation(normalizerY->GetOutput());
	imageGY.itk() = ImageGrayu8(normalizerY->GetOutput()).itk();

}

// ---- Private member functions ----

// Fill the output texture with globalNode
void GCTexture::fillOutputImage()
{
	for (int j = 0; j < output_h; j++)
	{
		for (int i = 0; i < output_w; i++)
		{
			if (globalNode[getNodeNbGlobal(0, 0, i, j)].notEmpty())
			{
				outputImage.pixelAbsolute(i, j) = globalNode[getNodeNbGlobal(0, 0, i, j)].getColor();
			}
			else
			{
				outputImage.pixelAbsolute(i, j) = RGB<uint8_t >(0, 0, 0);
			}
		}
	}
}

// Write the output texture, with patch number
void GCTexture::writeOutputImagePatch()
{
	std::ostringstream resultName, resultNameSeams;
	resultName << textureName.c_str() << "_" << std::setfill('0') << std::setw(3) << patchCount << ".png";
	resultNameSeams << textureName.c_str() << "_s_" << std::setfill('0') << std::setw(3) << patchCount << ".png";
	outputImage.save(resultName.str());
	revealSeamsError();
	revealSeamsMaxError(8);
	outputImage.save(resultNameSeams.str());
}

// Locates the seam pixel that has the maximum error value
// and returns its index maxErrNodeNbGlobal in output texture
// space, as well as the error value
double GCTexture::updateSeamsMaxError()
{
	int nodeNbGlobal;
	maxErrNodeNbGlobal = -1;
	double maxErr = -1.0;


	// Prevents the max error point to be located near the
	// boundaries of the synthesized image (some issues
	// appear when max error is located near borders,
	// check this).
	int size_2 = borderSize;

	for (int j = size_2; j < output_h - size_2; j++)
	{
		for (int i = size_2; i < output_w - size_2; i++)
		{
			nodeNbGlobal = getNodeNbGlobal(0, 0, i, j);
			if (globalNode[nodeNbGlobal].notEmpty())
			{
				if (globalNode[nodeNbGlobal].onSeamRight())
				{
					if (globalNode[nodeNbGlobal].getRightCost() > maxErr)
					{
						maxErr = globalNode[nodeNbGlobal].getRightCost();
						maxErrNodeNbGlobal = nodeNbGlobal;
					}
				}

				if (globalNode[nodeNbGlobal].onSeamBottom())
				{
					if (globalNode[nodeNbGlobal].getBottomCost() > maxErr)
					{
						maxErr = globalNode[nodeNbGlobal].getBottomCost();
						maxErrNodeNbGlobal = nodeNbGlobal;
					}
				}
			}
		}
	}

#ifdef _DEBUG
	if (maxErrNodeNbGlobal != -1)
	{
		int maxErrX, maxErrY;
		getPixelGlobal(maxErrNodeNbGlobal, maxErrX, maxErrY);
		std::cout << "Max seam error: " << maxErr << " at (" << maxErrX << ", " << maxErrY << ")\n";
		computeSeamsAverageError();
	}
	else
	{
		std::cout << "No seam error\n";
	}
#endif

	return maxErr;
}

// Average seam error
double GCTexture::computeSeamsAverageError()
{
	int nodeNbGlobal;
	int avNodes = 0;
	double avErr = 0.0;

	for (int j = 0; j < output_h; j++)
	{
		for (int i = 0; i < output_w; i++)
		{
			nodeNbGlobal = getNodeNbGlobal(0, 0, i, j);
			if (globalNode[nodeNbGlobal].notEmpty())
			{
				if (globalNode[nodeNbGlobal].onSeamRight())
				{
					avErr += globalNode[nodeNbGlobal].getRightCost();
					avNodes++;
				}

				if (globalNode[nodeNbGlobal].onSeamBottom())
				{
					avErr += globalNode[nodeNbGlobal].getBottomCost();
					avNodes++;
				}
			}
		}
	}

#ifdef _DEBUG
	std::cout << "Average seam error: " << (avErr / avNodes) << std::endl;
#endif
	
	return (avErr / avNodes);
}

// Assumes that updateSeamsMaxError has been called before.
void GCTexture::revealSeamsMaxError(int radius)
{
	// No seam error to draw
	if (maxErrNodeNbGlobal == -1)
		return;

	int maxErrX, maxErrY;
	getPixelGlobal(maxErrNodeNbGlobal, maxErrX, maxErrY);

	for (int sj = -radius; sj < radius; sj++)
	{
		for (int si = -radius; si < radius; si++)
		{
			int tj = maxErrY + sj;
			int ti = maxErrX + si;
			if ((ti >= 0 && ti < output_w) && (tj >= 0 && tj < output_h))
			{
				outputImage.pixelAbsolute(ti, tj) = RGB<uint8_t >(255, 0, 255);
			}
		}
	}
}

// Given a pixel (i, j) in the cropped intput patch and an insertion position
// (x, y) in the output texture, returns the index of the corresponding node
// in the global graph, -1 otherwise
int GCTexture::getNodeNbGlobal(int x, int y, int i, int j)
{
	if ((x + i < output_w) && (x + i >= 0) && (y + j < output_h) && (y + j >= 0))
	{
		return (x + i) + output_w * (y + j);
	}
	else
	{
		return -1;
	}
}

// Given a node index in globalNode, returns the corresponding
// pixel position (x, y) in the output texture
void GCTexture::getPixelGlobal(int nodeNbGlobal, int& x, int& y)
{
	x = nodeNbGlobal % output_w;
	y = nodeNbGlobal / output_h;
}

// Given a pixel position (i, j) in the cropped input patch,
// returns the corresponding node index in the local graph
int GCTexture::getNodeNbLocal(int i, int j)
{
	return i + croppedInput_w * j;
}

// Given a node index nodeNbLocal in the local graph, returns
// the corresponding pixel position (i, j) in the cropped
// input patch
void GCTexture::getPixel(int nodeNbLocal, int& i, int& j)
{
	i = nodeNbLocal % croppedInput_w;
	j = nodeNbLocal / croppedInput_h;
}

class Point
{
public:
	Point(int x_, int y_) { px = x_; py = y_; }
	int x() const { return px; }
	int y() const { return py; }

private:
	int px;
	int py;
};

float distPoint(const Point& p1, const Point& p2)
{
	return std::sqrt((p2.x() - p1.x())*(p2.x() - p1.x()) + (p2.y() - p1.y())*(p2.y() - p1.y()));
}

float alphaFunction(float dist, float radius, int n)
{
	if (dist > radius)
		return 0.0f;
	else
		return std::pow(1 - ((dist*dist) / (radius*radius)), n);
}

// Insert a new patch in the output texture at position (x, y).
// The origin of a patch is a upper left corner.
// filling: true, initial texture filling
// blending: if true, blending is enabled
// radius, intputX, inputY, tX, tY: to be used in conjunction
// with refinement, define the region tha is forced to be inserted in the
// output texture (tx, ty: translation that gives the insertion point
// in the output texture)
int GCTexture::insertPatch(int x, int y, bool filling, bool blending, int radius, int inputX, int inputY, int tX, int tY)
{
	return insertPatch(x, y, 0, 0, input_w, input_h, filling, blending, radius, inputX, inputY, tX, tY);
}

int GCTexture::insertPatch(int x, int y, int px, int py, int sx, int sy, bool filling, bool blending, int radius, int inputX, int inputY, int tX, int tY)
{
#ifdef _DEBUG
	std::cout << "Insert patch\n";
#endif

	int nodeNbGlobal;

	// ---- Crop patch to match output texture ----

	// Bounding box of the region of the new patch included
	// in the output texture (coordinates in input image
	// space, i.e. with origin at input image top left corner)
	int v1X = std::max(0, x) - x;
	int v1Y = std::max(0, y) - y;
	int v2X = std::min(output_w - 1, x + sx - 1) - x;
	int v2Y = std::min(output_h - 1, y + sy - 1) - y;

	croppedInput_w = v2X - v1X + 1;
	croppedInput_h = v2Y - v1Y + 1;

#ifdef _DEBUG
	std::cout << "Cropped patch size: " << croppedInput_w << "*" << croppedInput_h << std::endl;
#endif

	if (croppedInput_w == 0 || croppedInput_h == 0)
	{
#ifdef _DEBUG
		std::cout << "Patch lies outside output texture\n";
#endif
		return -1;
	}

	croppedInputImage.initItk(croppedInput_w, croppedInput_h);
	croppedInputImageGX.initItk(croppedInput_w, croppedInput_h);
	croppedInputImageGY.initItk(croppedInput_w, croppedInput_h);

	for (int j = v1Y, cj = 0; j <= v2Y; j++, cj++)
	{
		for (int i = v1X, ci = 0; i <= v2X; i++, ci++)
		{
			croppedInputImage.pixelAbsolute(ci, cj) = inputImage.pixelAbsolute(i + px, j + py);
			croppedInputImageGX.pixelAbsolute(ci, cj) = inputImageGX.pixelAbsolute(i + px, j + py);
			croppedInputImageGY.pixelAbsolute(ci, cj) = inputImageGY.pixelAbsolute(i + px, j + py);
		}
	}

	// Update origin coordinates
	x = std::max(0, x);
	y = std::max(0, y);

	// ---- Fill globalNode in non-overlapping region ---

	int noOverlap = 0;
	for (int j = 0; j < croppedInput_h; j++)
	{
		for (int i = 0; i < croppedInput_w; i++)
		{
			nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
			if (globalNode[nodeNbGlobal].empty())
				noOverlap++;
		}
	}

	if (noOverlap == croppedInput_w * croppedInput_h)
	{
#ifdef _DEBUG
		std::cout << "No overlap detected\n";
#endif

		for (int j = 0; j < croppedInput_h; j++)
		{
			for (int i = 0; i < croppedInput_w; i++)
			{
				nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
				if (globalNode[nodeNbGlobal].empty())
					globalNode[nodeNbGlobal].setColor(croppedInputImage.pixelAbsolute(i, j));
			}
		}

		return 0;
	}
	else if ((noOverlap == 0) && filling)
	{
#ifdef _DEBUG
		std::cout << "Patch does not contribute\n";
#endif
		return -1;
	}

#ifdef _DEBUG
	std::cout << "Overlap detected\n";
#endif

	// ---- Computation of the gradients of the current output texture ----
	// for old pixels overlapping the new patch
	// Note: regions of the texture overlapped by the new patch
	// that do not contain old pixels are filled in black

	ImageRGBu8 image; image.initItk(croppedInput_w, croppedInput_h);
	ImageGrayu8 imageGX; imageGX.initItk(croppedInput_w, croppedInput_h);
	ImageGrayu8 imageGY; imageGY.initItk(croppedInput_w, croppedInput_h);

	for (int j = 0; j < croppedInput_h; j++)
	{
		for (int i = 0; i < croppedInput_w; i++)
		{
			nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
			if (globalNode[nodeNbGlobal].notEmpty())
			{
				image.pixelAbsolute(i, j) = globalNode[nodeNbGlobal].getColor();
			}
			else
			{
				image.pixelAbsolute(i, j) = RGBu8(0, 0, 0);
			}
		}
	}

	computeGradientImage(image, imageGX, imageGY);

	// ---- Graph construction ----

	// Initial number of nodes
	int nbNodes = croppedInput_w * croppedInput_h;
	// Initial number of graph edges (this number does not need to include edges to source and sink)
	int nbEdges = (croppedInput_w - 1)*(croppedInput_h - 1) * 2 + (croppedInput_w - 1) + (croppedInput_h - 1);
	// Number of seam nodes
	int nbSeamNodes = nbEdges;
	// Number of edges (add edges from seam nodes to artificial new nodes)
	nbEdges += nbSeamNodes;

	// Prepare graph data
	g = new GraphType(nbNodes + 2 * nbSeamNodes, nbEdges);
#ifdef _DEBUG
	std::cout << "Graph allocated\n";
#endif

	// Create nodes for texture pixels
	for (int i = 0; i < nbNodes; i++)
		g->add_node();

	int nodeNbGlobalRight, nodeNbGlobalBottom;
	int nodeNbLocal, nodeNbLocalRight, nodeNbLocalBottom;
	double grad;
	double d1, d2, d3, d4;
	double capRight, capBottom;
	double capacity1, capacity2, capacity3;

	int overlap = 0;

	// Create edges
	for (int j = 0; j < croppedInput_h; j++)
	{
		for (int i = 0; i < croppedInput_w; i++)
		{
			nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
			nodeNbLocal = getNodeNbLocal(i, j);

			if (i < croppedInput_w - 1)
			{
				nodeNbGlobalRight = getNodeNbGlobal(x, y, i + 1, j);
				nodeNbLocalRight = getNodeNbLocal(i + 1, j);

				if (globalNode[nodeNbGlobal].notEmpty())
				{
					if (globalNode[nodeNbGlobalRight].notEmpty())
					{
						// Overlap
						d1 = distColors(globalNode[nodeNbGlobal].getColor(), croppedInputImage.pixelAbsolute(i, j));
						d2 = distColors(globalNode[nodeNbGlobalRight].getColor(), croppedInputImage.pixelAbsolute(i + 1, j));
						if (globalNode[nodeNbGlobal].onSeamRight())
						{
							// Old seam: a seam node will created
							//
							//               Sink
							//                 o
							//                 |
							//             capacity1
							//                 |
							//   X--capacity2--+--capacity3--o
							//
							capacity1 = globalNode[nodeNbGlobal].getRightCost();
							d3 = distColors(globalNode[nodeNbGlobalRight].getColorOtherPatch(), croppedInputImage.pixelAbsolute(i + 1, j));
							d4 = distColors(croppedInputImage.pixelAbsolute(i, j), globalNode[nodeNbGlobal].getColorOtherPatch());
							grad = ((croppedInputImageGX.pixelAbsolute(i, j) / 255.0)
								+ (croppedInputImageGX.pixelAbsolute(i + 1, j) / 255.0)
								+ (imageGX.pixelAbsolute(i, j) / 255.0)
								+ (globalNode[nodeNbGlobalRight].getGradXOtherPatch() / 255.0));
							grad += 1.0;
							capacity2 = (d1 + d3) / grad;
							//capacity2 = d1 + d3;
							grad = (croppedInputImageGX.pixelAbsolute(i, j) / 255.0)
								+ (croppedInputImageGX.pixelAbsolute(i + 1, j) / 255.0)
								+ (imageGX.pixelAbsolute(i + 1, j) / 255.0)
								+ (globalNode[nodeNbGlobal].getGradXOtherPatch() / 255.0);
							grad += 1.0;
							capacity3 = (d4 + d2) / grad;
							//capacity3 = d4 + d2;
							capacity2 += minCap;
							capacity3 += minCap;
							seamNode.push_back(SeamNode(nodeNbLocal, nodeNbLocalRight, capacity1, capacity2, capacity3, 0));
							
							//globalNode[nodeNbGlobal].unsetSeamRight();
						}
						else
						{
							// No old seam
							grad = (croppedInputImageGX.pixelAbsolute(i, j) / 255.0)
								+ (croppedInputImageGX.pixelAbsolute(i + 1, j) / 255.0)
								+ (imageGX.pixelAbsolute(i, j) / 255.0)
								+ (imageGX.pixelAbsolute(i + 1, j) / 255.0);
							grad += 1.0;
							capRight = (d1 + d2) / grad;
							capRight += minCap;
							//capRight = d1 + d2;
							g->add_edge(nodeNbLocal, nodeNbLocalRight, capRight, capRight);
							globalNode[nodeNbGlobal].setRightCost(capRight);
						}
						overlap++;
					}
					else
					{
						// No overlap
						g->add_edge(nodeNbLocal, nodeNbLocalRight, 0.0, 0.0);
						globalNode[nodeNbGlobal].setRightCost(0.0);
					}
				}
				else
				{
					// No overlap
					g->add_edge(nodeNbLocal, nodeNbLocalRight, 0.0, 0.0);
					globalNode[nodeNbGlobal].setRightCost(0.0);
				}
			}

			if (j < croppedInput_h - 1)
			{
				nodeNbGlobalBottom = getNodeNbGlobal(x, y, i, j + 1);
				nodeNbLocalBottom = getNodeNbLocal(i, j + 1);

				if (globalNode[nodeNbGlobal].notEmpty())
				{
					if (globalNode[nodeNbGlobalBottom].notEmpty())
					{
						// Overlap
						d1 = distColors(globalNode[nodeNbGlobal].getColor(), croppedInputImage.pixelAbsolute(i, j));
						d2 = distColors(globalNode[nodeNbGlobalBottom].getColor(), croppedInputImage.pixelAbsolute(i, j + 1));
						if (globalNode[nodeNbGlobal].onSeamBottom())
						{
							// Old seam: a seam node will created
							//
							//        X
							//        |
							//    capacity2
							//        |
							//        +--capacity1--o Sink
							//        |
							//    capacity3
							//        |
							//        o
							//
							capacity1 = globalNode[nodeNbGlobal].getBottomCost();
							d3 = distColors(globalNode[nodeNbGlobalBottom].getColorOtherPatch(), croppedInputImage.pixelAbsolute(i, j + 1));
							d4 = distColors(croppedInputImage.pixelAbsolute(i, j), globalNode[nodeNbGlobal].getColorOtherPatch());
							grad = (croppedInputImageGY.pixelAbsolute(i, j) / 255.0)
								+ (croppedInputImageGY.pixelAbsolute(i, j + 1) / 255.0)
								+ (imageGY.pixelAbsolute(i, j) / 255.0)
								+ (globalNode[nodeNbGlobalBottom].getGradYOtherPatch() / 255.0);
							grad += 1.0;
							capacity2 = (d1 + d3) / grad;
							//capacity2 = d1 + d3;
							grad = (croppedInputImageGY.pixelAbsolute(i, j) / 255.0)
								+ (croppedInputImageGY.pixelAbsolute(i, j + 1) / 255.0)
								+ (imageGY.pixelAbsolute(i, j + 1) / 255.0)
								+ (globalNode[nodeNbGlobal].getGradYOtherPatch() / 255.0);
							grad += 1.0;
							capacity3 = (d4 + d2) / grad;
							//capacity3 = d4 + d2;
							capacity2 += minCap;
							capacity3 += minCap;
							seamNode.push_back(SeamNode(nodeNbLocal, nodeNbLocalBottom, capacity1, capacity2, capacity3, 1));
							
							//globalNode[nodeNbGlobal].unsetSeamBottom();
						}
						else
						{
							// No old seam
							grad = (croppedInputImageGY.pixelAbsolute(i, j) / 255.0)
								+ (croppedInputImageGY.pixelAbsolute(i, j + 1) / 255.0)
								+ (imageGY.pixelAbsolute(i, j) / 255.0)
								+ (imageGY.pixelAbsolute(i, j + 1) / 255.0);
							grad += 1.0;
							capBottom = (d1 + d2) / grad;
							capBottom += minCap;
							//capBottom = d1 + d2;
							g->add_edge(nodeNbLocal, nodeNbLocalBottom, capBottom, capBottom);
							globalNode[nodeNbGlobal].setBottomCost(capBottom);
						}
						overlap++;
					}
					else
					{
						// No overlap
						g->add_edge(nodeNbLocal, nodeNbLocalBottom, 0.0, 0.0);
						globalNode[nodeNbGlobal].setBottomCost(0.0);
					}
				}
				else
				{
					// No overlap
					g->add_edge(nodeNbLocal, nodeNbLocalBottom, 0.0, 0.0);
					globalNode[nodeNbGlobal].setBottomCost(0.0);
				}
			}
		}
	}

#ifdef _DEBUG
	std::cout << "Number of seam nodes: " << seamNode.size() << std::endl;
#endif

	// Add seam nodes
	// Note: Seam nodes are inserted at the end so that the indices
	// of globalNodes still match the indices of graph nodes.
	int nodeOldSeam;
	for (unsigned int i = 0; i < seamNode.size(); i++)
	{
		nodeOldSeam = g->add_node();
		seamNode[i].setSeam(nodeOldSeam);
		g->add_edge(seamNode[i].getStart(), nodeOldSeam, seamNode[i].getCapacity2(), seamNode[i].getCapacity2());
		g->add_edge(nodeOldSeam, seamNode[i].getEnd(), seamNode[i].getCapacity3(), seamNode[i].getCapacity3());
		g->add_tweights(nodeOldSeam, 0.0, seamNode[i].getCapacity1());
	}

#ifdef _DEBUG
	std::cout << "Graph created\n";
#endif

	// Assignments to source node
	for (int i = 0; i < croppedInput_w; i++)
	{
		int j = 0;
		// Node position in globalNode
		nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
		if (globalNode[nodeNbGlobal].notEmpty())
		{
			g->add_tweights(getNodeNbLocal(i, j), infiniteCap, 0.0);
		}

		j = croppedInput_h - 1;
		// Node position in globalNode
		nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
		if (globalNode[nodeNbGlobal].notEmpty())
		{
			g->add_tweights(getNodeNbLocal(i, j), infiniteCap, 0.0);
		}
	}

	for (int j = 0; j < croppedInput_h; j++)
	{
		int i = 0;
		// Node position in globalNode
		nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
		if (globalNode[nodeNbGlobal].notEmpty())
		{
			g->add_tweights(getNodeNbLocal(i, j), infiniteCap, 0.0);
		}

		i = croppedInput_w - 1;
		// Node position in globalNode
		nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
		if (globalNode[nodeNbGlobal].notEmpty())
		{
			g->add_tweights(getNodeNbLocal(i, j), infiniteCap, 0.0);
		}
	}

	int nodeNbGlobalLeft, nodeNbGlobalTop;
	int nbSink = 0;

	if (filling)
	{
		// Assignments to sink node
		for (int j = 0; j < croppedInput_h; j++)
		{
			for (int i = 0; i < croppedInput_w; i++)
			{
				// Node position in globalNode
				nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
				if (globalNode[nodeNbGlobal].notEmpty())
				{
					nodeNbGlobalLeft = getNodeNbGlobal(x, y, i - 1, j);
					nodeNbGlobalRight = getNodeNbGlobal(x, y, i + 1, j);
					nodeNbGlobalTop = getNodeNbGlobal(x, y, i, j - 1);
					nodeNbGlobalBottom = getNodeNbGlobal(x, y, i, j + 1);

					if ((nodeNbGlobalLeft != -1) && globalNode[nodeNbGlobalLeft].empty())
					{
						g->add_tweights(getNodeNbLocal(i, j), 0.0, infiniteCap);
						nbSink++;
					}
					else if ((nodeNbGlobalTop != -1) && globalNode[nodeNbGlobalTop].empty())
					{
						g->add_tweights(getNodeNbLocal(i, j), 0.0, infiniteCap);
						nbSink++;
					}
					else if ((nodeNbGlobalRight != -1) && globalNode[nodeNbGlobalRight].empty())
					{
						g->add_tweights(getNodeNbLocal(i, j), 0.0, infiniteCap);
						nbSink++;
					}
					else if ((nodeNbGlobalBottom != -1) && globalNode[nodeNbGlobalBottom].empty())
					{
						g->add_tweights(getNodeNbLocal(i, j), 0.0, infiniteCap);
						nbSink++;
					}
				}
			}
		}
	}
	else // Refinement case
	{
		// Add edge from the center of the new patch to the sink node
		// Note: improve this.
		//nodeNbLocal = getNodeNbLocal(croppedInput_w / 2, croppedInput_h / 2);
		//g->add_tweights(nodeNbLocal, 0.0, infiniteCap);

		// Crop subpatch
		int outputX = inputX + tX;
		int outputY = inputY + tY;
		v1X = std::max(0, outputX);
		v1Y = std::max(0, outputY);
		v2X = std::min(output_w - 1, outputX + radius - 1);
		v2Y = std::min(output_h - 1, outputY + radius - 1);
		// Subpatch origin coordinates in output space
		outputX = v1X;
		outputY = v1Y;
		// Subpatch origin coordinates in input space
		inputX = outputX - tX;
		inputY = outputY - tY;

		//std::cout << inputX << ", " << inputY << std::endl;
		//std::cout << outputX << ", " << outputY << " Size: " << (v2X - v1X + 1) << "*" << (v2Y - v1Y + 1) << std::endl;
		if ((inputX < 0) || (inputY < 0))
		{
			std::cout << inputX << ", " << inputY << std::endl;
			std::cout << outputX << ", " << outputY << " Size: " << (v2X - v1X + 1) << "*" << (v2Y - v1Y + 1) << std::endl;
			getchar();
		}

		for (int j = v1Y, cj = 0; j <= v2Y; j++, cj++)
		{
			for (int i = v1X, ci = 0; i <= v2X; i++, ci++)
			{
				if ((inputX + ci > 0) && (inputY + cj > 0))
				{
					nodeNbLocal = getNodeNbLocal(inputX + ci, inputY + cj);
					g->add_tweights(nodeNbLocal, 0.0, infiniteCap);
					nbSink++;
				}
			}
		}
	}

#ifdef _DEBUG
	std::cout << "Assignments done\n";
#endif

#ifdef _DEBUG
	std::cout << "Number of sink nodes : " << nbSink << std::endl;
#endif

	assert(nbSink);

	// Compute max flow
	double maxFlow = g->maxflow();

#ifdef _DEBUG
	std::cout << "Max flow computed (" << maxFlow << ")\n";
#endif

	// Store seams
	// Note: add the second parameter to what_segment so that the graph nodes outside the overlap region
	// have SINK assignment by default. This is required to not include the boundary of the overlap
	// region in the seam.
	// Note: edges to source and sink nodes are not taken into account into the node set.

	int currentSeamNode, currentSeamNodeEnd;
	unsigned int k = 0;
	for (int j = 0; j < croppedInput_h; j++)
	{
		for (int i = 0; i < croppedInput_w; i++)
		{
			nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
			nodeNbLocal = getNodeNbLocal(i, j);

			if (i < croppedInput_w - 1)
			{
				if (g->what_segment(nodeNbLocal, GraphType::SINK) != g->what_segment(getNodeNbLocal(i + 1, j), GraphType::SINK))
				{
					globalNode[nodeNbGlobal].setNewSeam();
				}
			}

			if (j < croppedInput_h - 1)
			{
				if (g->what_segment(nodeNbLocal, GraphType::SINK) != g->what_segment(getNodeNbLocal(i, j + 1), GraphType::SINK))
				{
					globalNode[nodeNbGlobal].setNewSeam();
				}
			}

			
			if (seamNode.size() && (k < seamNode.size()) && (nodeNbLocal == seamNode[k].getStart()))
			{
				// Process old seam

				currentSeamNode = seamNode[k].getSeam();
				currentSeamNodeEnd = seamNode[k].getEnd();

				if ((g->what_segment(nodeNbLocal, GraphType::SINK) == GraphType::SOURCE)
					&& (g->what_segment(currentSeamNodeEnd, GraphType::SINK) == GraphType::SINK))
				{
					
					if (g->what_segment(currentSeamNode, GraphType::SINK) == GraphType::SOURCE)
					{
						// Old seam remains with new seam cost
						if (seamNode[k].getOrientation() == 0)
						{
							// Right
							globalNode[nodeNbGlobal].setRightCost(seamNode[k].getCapacity3());
							globalNode[nodeNbGlobal].setSeamRight(maxFlow);
						}
						else
						{
							// Bottom
							globalNode[nodeNbGlobal].setBottomCost(seamNode[k].getCapacity3());
							globalNode[nodeNbGlobal].setSeamBottom(maxFlow);
						}
					}
					else
					{
						// Old seam remains with new seam cost
						if (seamNode[k].getOrientation() == 0)
						{
							// Right
							globalNode[nodeNbGlobal].setRightCost(seamNode[k].getCapacity2());
							globalNode[nodeNbGlobal].setSeamRight(maxFlow);
						}
						else
						{
							// Bottom
							globalNode[nodeNbGlobal].setBottomCost(seamNode[k].getCapacity2());
							globalNode[nodeNbGlobal].setSeamBottom(maxFlow);
						}
					}
					
				}
				else if ((g->what_segment(nodeNbLocal, GraphType::SINK) == GraphType::SINK)
					&& (g->what_segment(currentSeamNodeEnd, GraphType::SINK) == GraphType::SOURCE))
				{
					
					if (g->what_segment(currentSeamNode, GraphType::SINK) == GraphType::SOURCE)
					{
						// Old seam remains with new seam cost
						if (seamNode[k].getOrientation() == 0)
						{
							// Right
							globalNode[nodeNbGlobal].setRightCost(seamNode[k].getCapacity2());
							globalNode[nodeNbGlobal].setSeamRight(maxFlow);
						}
						else
						{
							// Bottom
							globalNode[nodeNbGlobal].setBottomCost(seamNode[k].getCapacity2());
							globalNode[nodeNbGlobal].setSeamBottom(maxFlow);
						}
					}
					else
					{
						// Old seam remains with new seam cost
						if (seamNode[k].getOrientation() == 0)
						{
							// Right
							globalNode[nodeNbGlobal].setRightCost(seamNode[k].getCapacity3());
							globalNode[nodeNbGlobal].setSeamRight(maxFlow);
						}
						else
						{
							// Bottom
							globalNode[nodeNbGlobal].setBottomCost(seamNode[k].getCapacity3());
							globalNode[nodeNbGlobal].setSeamBottom(maxFlow);
						}
					}
					
				}
				else if (g->what_segment(currentSeamNode, GraphType::SINK) == GraphType::SOURCE)
				{
					// Old seam with old seam cost
					if (seamNode[k].getOrientation() == 0)
					{
						// Right
						globalNode[nodeNbGlobal].setRightCost(seamNode[k].getCapacity1());
						globalNode[nodeNbGlobal].setSeamRight(maxFlow);
					}
					else
					{
						// Bottom
						globalNode[nodeNbGlobal].setBottomCost(seamNode[k].getCapacity1());
						globalNode[nodeNbGlobal].setSeamBottom(maxFlow);
					}
				}
				else
				{
					
					/*
					// No seam to set, but we want to keep the old seam, so we restore the
					// seam and its cost

					if (seamNode[k].getOrientation() == 0)
					{
						// Right
						globalNode[nodeNbGlobal].setRightCost(seamNode[k].getCapacity1());
						globalNode[nodeNbGlobal].setSeamRight(maxFlow);
					}
					else
					{
						// Bottom
						globalNode[nodeNbGlobal].setBottomCost(seamNode[k].getCapacity1());
						globalNode[nodeNbGlobal].setSeamBottom(maxFlow);
					}
					*/

				}

				k++;
			}
			else
			{
				// New seam
				// Note: costs are already assigned to edges in globalNode.
				if (i < croppedInput_w - 1)
				{
					if (g->what_segment(nodeNbLocal, GraphType::SINK) != g->what_segment(getNodeNbLocal(i + 1, j), GraphType::SINK))
					{
						globalNode[nodeNbGlobal].setSeamRight(maxFlow);
					}
				}

				if (j < croppedInput_h - 1)
				{
					if (g->what_segment(nodeNbLocal, GraphType::SINK) != g->what_segment(getNodeNbLocal(i, j + 1), GraphType::SINK))
					{
						globalNode[nodeNbGlobal].setSeamBottom(maxFlow);
					}
				}
			}
		}
	}

	if (blending)
	{
		std::vector<Point> seamPoints;

		int nodeNbGlobal;
		for (int j = 0; j < croppedInput_h; j++)
		{
			for (int i = 0; i < croppedInput_w; i++)
			{
				nodeNbGlobal = getNodeNbGlobal(x, y, i, j);

				if (globalNode[nodeNbGlobal].onNewSeam()
				    && (globalNode[nodeNbGlobal].onSeamRight()
					|| globalNode[nodeNbGlobal].onSeamBottom()))
				{
					seamPoints.push_back(Point(i, j));

					if (globalNode[nodeNbGlobal].onSeamRight())
						seamPoints.push_back(Point(i + 1, j));

					if (globalNode[nodeNbGlobal].onSeamBottom())
						seamPoints.push_back(Point(i, j + 1));
				}
			}
		}

		float dist, distMin, alpha;

		for (int j = 0; j < croppedInput_h; j++)
		{
			for (int i = 0; i < croppedInput_w; i++)
			{
				nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
				if (globalNode[nodeNbGlobal].empty())
				{
					// New pixel insertion
					globalNode[nodeNbGlobal].setColor(croppedInputImage.pixelAbsolute(i, j));
				}
				else
				{
					if (g->what_segment(getNodeNbLocal(i, j), GraphType::SINK) == GraphType::SOURCE)
					{
						// Source side: pixels from the old patch
						//globalNode[nodeNbGlobal].setColorOtherPatch(croppedInputImage.pixel(i, j));
						globalNode[nodeNbGlobal].setGradXOtherPatch(croppedInputImageGX.pixelAbsolute(i, j));
						globalNode[nodeNbGlobal].setGradYOtherPatch(croppedInputImageGY.pixelAbsolute(i, j));

						distMin = std::numeric_limits<float>::max();
						for (std::size_t k = 0; k < seamPoints.size(); k++)
						{
							dist = distPoint(Point(i, j), seamPoints[k]);
							if (dist < distMin)
								distMin = dist;
						}

						alpha = alphaFunction(distMin, blradius, bln);

						int rr = (int)(alpha * globalNode[nodeNbGlobal].getColor().GetRed() + (1 - alpha) * croppedInputImage.pixelAbsolute(i, j).GetRed());
						if (rr < 0)	rr = 0;
						if (rr > 255) rr = 255;
						int gg = (int)(alpha * globalNode[nodeNbGlobal].getColor().GetGreen() + (1 - alpha) * croppedInputImage.pixelAbsolute(i, j).GetGreen());
						if (gg < 0) gg = 0;
						if (gg > 255) gg = 255;
						int bb = (int)(alpha * globalNode[nodeNbGlobal].getColor().GetBlue() + (1 - alpha) * croppedInputImage.pixelAbsolute(i, j).GetBlue());
						if (bb < 0) bb = 0;
						if (bb > 255) bb = 255;

						globalNode[nodeNbGlobal].setColorOtherPatch(RGB<uint8_t >(rr, gg, bb));
					}
					else
					{
						// Sink side: pixels from the new patch
						globalNode[nodeNbGlobal].setColorOtherPatch(globalNode[nodeNbGlobal].getColor());
						globalNode[nodeNbGlobal].setGradXOtherPatch(imageGX.pixelAbsolute(i, j));
						globalNode[nodeNbGlobal].setGradYOtherPatch(imageGY.pixelAbsolute(i, j));

						// Blending on the new patch side
						distMin = std::numeric_limits<float>::max();
						for (std::size_t k = 0; k < seamPoints.size(); k++)
						{
							dist = distPoint(Point(i, j), seamPoints[k]);
							if (dist < distMin)
								distMin = dist;
						}

						alpha = alphaFunction(distMin, blradius, bln);

						int rr = (int)(alpha * globalNode[nodeNbGlobal].getColor().GetRed() + (1 - alpha) * croppedInputImage.pixelAbsolute(i, j).GetRed());
						if (rr < 0) rr = 0;
						if (rr > 255) rr = 255;
						int gg = (int)(alpha * globalNode[nodeNbGlobal].getColor().GetGreen() + (1 - alpha) * croppedInputImage.pixelAbsolute(i, j).GetGreen());
						if (gg < 0) gg = 0;
						if (gg > 255) gg = 255;
						int bb = (int)(alpha * globalNode[nodeNbGlobal].getColor().GetBlue() + (1 - alpha) * croppedInputImage.pixelAbsolute(i, j).GetBlue());
						if (bb < 0) bb = 0;
						if (bb > 255) bb = 255;

						globalNode[nodeNbGlobal].setColor(RGB<uint8_t >(rr, gg, bb));

						if (!globalNode[nodeNbGlobal].onNewSeam())
						{
							if (globalNode[nodeNbGlobal].onSeamRight())
								globalNode[nodeNbGlobal].unsetSeamRight();

							if (globalNode[nodeNbGlobal].onSeamBottom())
								globalNode[nodeNbGlobal].unsetSeamBottom();
						}
					}
				}

				if (globalNode[nodeNbGlobal].onNewSeam())
					globalNode[nodeNbGlobal].unsetNewSeam();
			}
		}

		seamPoints.clear();
	}
	else
	{
		for (int j = 0; j < croppedInput_h; j++)
		{
			for (int i = 0; i < croppedInput_w; i++)
			{
				// Remove seam in new patch region where there is no new seam
				// If new seam, unset new seam
				
				nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
				if (globalNode[nodeNbGlobal].empty())
				{
					// New pixel insertion
					globalNode[nodeNbGlobal].setColor(croppedInputImage.pixelAbsolute(i, j));
				}
				else
				{
					if (g->what_segment(getNodeNbLocal(i, j), GraphType::SINK) == GraphType::SOURCE)
					{
						// Source side: pixels from the old patch
						globalNode[nodeNbGlobal].setColorOtherPatch(croppedInputImage.pixelAbsolute(i, j));
						globalNode[nodeNbGlobal].setGradXOtherPatch(croppedInputImageGX.pixelAbsolute(i, j));
						globalNode[nodeNbGlobal].setGradYOtherPatch(croppedInputImageGY.pixelAbsolute(i, j));
					}
					else
					{
						// Sink side: pixels from the new patch
						globalNode[nodeNbGlobal].setColorOtherPatch(globalNode[nodeNbGlobal].getColor());
						globalNode[nodeNbGlobal].setGradXOtherPatch(imageGX.pixelAbsolute(i, j));
						globalNode[nodeNbGlobal].setGradYOtherPatch(imageGY.pixelAbsolute(i, j));
						globalNode[nodeNbGlobal].setColor(croppedInputImage.pixelAbsolute(i, j));

						if (!globalNode[nodeNbGlobal].onNewSeam())
						{
							if (globalNode[nodeNbGlobal].onSeamRight())
								globalNode[nodeNbGlobal].unsetSeamRight();

							if (globalNode[nodeNbGlobal].onSeamBottom())
								globalNode[nodeNbGlobal].unsetSeamBottom();
						}
					}
				}

				if (globalNode[nodeNbGlobal].onNewSeam())
					globalNode[nodeNbGlobal].unsetNewSeam();
			}
		}
	}

	g->reset();
	delete g;
#ifdef _DEBUG
	std::cout << "Graph deleted\n";
#endif

	seamNode.clear();
#ifdef _DEBUG
	std::cout << "Seam node cleared\n";
#endif

	return overlap;
}

// ---- Constructor and destructor ----

GCTexture::GCTexture()
{
	isInitialized = false;
	isFilled = false;
	globalNode = NULL;
}

GCTexture::~GCTexture()
{
	std::cout << "Cleaning\n";

	if (globalNode != NULL)
		delete[] globalNode;
}

// ---- Public member functions ----

// Initialization of the synthesizer
// - fileNameAbs: absolute path of the input image file
// - textureName: name of the texture (no extension)
// - n_w, n_h: multiplication factors that define the size of
// the output texture based on the size of the input image
void GCTexture::initialize(const ImageRGBu8& iImage, const std::string& tName, int o_w, int o_h)
{
	// Initialize texture name
	textureName = tName;

	// Assign input image
	inputImage = iImage;
	input_w = iImage.width();
	input_h = iImage.height();

	// Compute gradient images
	computeGradientImage(inputImage, inputImageGX, inputImageGY);

	// Create output image
	output_w = o_w;
	output_h = o_h;
	std::cout << "Output image size: " << output_w << "*" << output_h << std::endl;
	// Output image allocation with all pixels set to (0, 0, 0)
	outputImage.initItk(output_w, output_h, true);

	// Global graph nodes
	globalNode = new GlobalNode[output_w * output_h];

	// Infinite capacity
	// Note: not set to std::numeric_limits<double>::max()
	// because of possibility of augmentation when running
	// maxflow (to check).
	infiniteCap = std::numeric_limits<double>::max() / 2.0;

	// Minimum capacity
	// Note: used to avoid unreachable nodes in the graph
	// when not desired.
	minCap = std::numeric_limits<double>::min();

	// Initialize patch count
	patchCount = 0;

	// Initialize seam size (for the rendering of seams)
	seamSize = 2;

	isInitialized = true;
	isFilled = false;

	// Initialize max error node
	maxErrNodeNbGlobal = -1;

	// Blending radius
	blradius = 5;
	// Blending smoothness
	bln = 2;

	// SSD computation step size
	ssdStep = 3;

	borderSize = 32;
}

// Initial texture synthesis function (initialize has to be run first!)
// The overlap between the input image and the output is
// at least bandX along the horizontal axis and bandY along the
// vertical axis. bandX, bandY also controls the amount of overlap
// between patches. Typical values: input_w / 3 for bandX and
// input_h / 3 for bandY.
void GCTexture::topDownRandomFilling(int bandX, int bandY, bool blending)
{
	if (!isInitialized)
	{
		std::cout << "Synthesizer not initialized!\n";
		return;
	}

	std::cout << "---- Initial texture synthesis ----\n";
	
	int offsetY = 0;

	int x, y;

	// Random Y coordinate in output image
	y = offsetY - (bandY + (std::rand() % bandY));

	do
	{
		std::cout << "---- Next row ----\n";

		// Random X coordinate in output image
		x = -(bandX + (std::rand() % bandX));

		do
		{
			std::cout << "x=" << x << " y=" << y << std::endl;
			if (y < output_h)
			{
				if (insertPatch(x, y, true, blending) != -1)
				{
					updateSeamsMaxError();

					// Patch contributes to the output texture
					patchCount++;

//#ifdef _DEBUG
					std::cout << "Patch " << patchCount << " inserted\n";

					// Fill output image with current synthesis result
					fillOutputImage();

					// Write output image in a file
					writeOutputImagePatch();
//#endif
				}
			}

			x = x + (bandX + (std::rand() % bandX));
			y = offsetY - (bandY + (std::rand() % bandY));

		} while (x < output_w);

		offsetY += bandY;
		y = offsetY - (bandY + (std::rand() % bandY));

	} while (y < output_h);

	isFilled = true;
}

// step: step for SSD computations
// Test all possible positions, with constraint that the new patch has to overlap
// the previously inserted patch (guarantee that the output texture will be filled)
void GCTexture::topDownEntireMatchingMinErrorFilling(int bandX, int bandY, bool blending)
{
	if (!isInitialized)
	{
		std::cout << "Synthesizer not initialized!\n";
		return;
	}

	std::cout << "---- Initial texture synthesis ----\n";

	double err, minErr;

	int offsetY = 0;

	int x, y;
	int testX, testY;
	int minX, minY;

	// Random Y coordinate in output image
	y = offsetY - (bandY + (std::rand() % bandY));
	// Random X coordinate in output image
	x = -(bandX + (std::rand() % bandX));

	do
	{
		std::cout << "---- Next row ----\n";

		do
		{
			std::cout << "x=" << x << " y=" << y << std::endl;
			if (y < output_h)
			{
				if (insertPatch(x, y, true, blending) != -1)
				{
					updateSeamsMaxError();

					// Patch contributes to the output texture
					patchCount++;

#ifdef _DEBUG
					std::cout << "Patch " << patchCount << " inserted\n";

					// Fill output image with current synthesis result
					fillOutputImage();

					// Write output image in a file
					writeOutputImagePatch();
#endif
				}
			}

			minX = x + (bandX + (std::rand() % bandX));
			minY = offsetY - (bandY + (std::rand() % bandY));
			minErr = computeSSD(minX, minY);
			for (int j = 0; j < bandY; j += ssdStep)
			{
				for (int i = 0; i < bandX; i += ssdStep)
				{
					testX = x + bandX + i;
					testY = offsetY - bandY - j;
					err = computeSSD(testX, testY);
					if (err < minErr)
					{
						minErr = err;
						minX = testX;
						minY = testY;
					}
				}
			}
			x = minX;
			y = minY;

		} while (x < output_w);

		offsetY += bandY;

		minX = -(bandX + (std::rand() % bandX));
		minY = offsetY - (bandY + (std::rand() % bandY));
		minErr = computeSSD(minX, minY);
		for (int j = 0; j < bandY; j += ssdStep)
		{
			for (int i = 0; i < bandX; i += ssdStep)
			{
				testX = -bandX - i;
				testY = offsetY - bandY - j;
				err = computeSSD(testX, testY);
				if (err < minErr)
				{
					minErr = err;
					minX = testX;
					minY = testY;
				}
			}
		}
		x = minX;
		y = minY;

	} while (y < output_h);

	isFilled = true;
}

struct lessPos {
	bool operator() (const std::pair<int, int>& p1, const std::pair<int, int>& p2) const
	{
		if (p1.first < p2.first)
			return true;
		else
		{
			if (p1.first == p2.first)
			{
				if (p1.second < p2.second)
					return true;
			}
		}

		return false;
	}
};

// step: step for SSD computations
// Test all possible positions, with constraint that the new patch has to overlap
// the previously inserted patch (guarantee that the output texture will be filled)
void GCTexture::topDownEntireMatchingRandomnessFilling(int bandX, int bandY, double k, bool blending)
{
	if (!isInitialized)
	{
		std::cout << "Synthesizer not initialized!\n";
		return;
	}

	std::cout << "---- Initial texture synthesis ----\n";

	double err;
	int l, n;

	int offsetY = 0;

	int x, y;
	int testX, testY;

	// Random Y coordinate in output image
	y = offsetY - (bandY + (std::rand() % bandY));
	// Random X coordinate in output image
	x = -(bandX + (std::rand() % bandX));

	do
	{
		std::cout << "---- Next row ----\n";

		do
		{
			std::cout << "x=" << x << " y=" << y << std::endl;
			if (y < output_h)
			{
				if (insertPatch(x, y, true, blending) != -1)
				{
					updateSeamsMaxError();

					// Patch contributes to the output texture
					patchCount++;

#ifdef _DEBUG
					std::cout << "Patch " << patchCount << " inserted\n";

					// Fill output image with current synthesis result
					fillOutputImage();

					// Write output image in a file
					writeOutputImagePatch();
#endif
				}
			}
			
			std::multimap< double, std::pair<int, int> > errMap;
			std::multimap< double, std::pair<int, int> >::const_iterator emhit;
			for (int j = 0; j < bandY; j += ssdStep)
			{
				for (int i = 0; i < bandX; i += ssdStep)
				{
					testX = x + bandX + i;
					testY = offsetY - bandY - j;
					err = computeSSD(testX, testY);
					errMap.insert(std::pair< double, std::pair<int, int> >(err, std::pair<int, int>(testX, testY)));
				}
			}			
			n = rand() % ((int)(k * errMap.size()));
			emhit = errMap.begin();
			l = 0;
			while (l < n)
			{
				l++;
				++emhit;
			}
			x = (*emhit).second.first;
			y = (*emhit).second.second;

		} while (x < output_w);

		offsetY += bandY;

		std::multimap< double, std::pair<int, int> > errMap;
		std::multimap< double, std::pair<int, int> >::const_iterator emhit;
		for (int j = 0; j < bandY; j += ssdStep)
		{
			for (int i = 0; i < bandX; i += ssdStep)
			{
				testX = -bandX - i;
				testY = offsetY - bandY - j;
				err = computeSSD(testX, testY);
				errMap.insert(std::pair< double, std::pair<int, int> >(err, std::pair<int, int>(testX, testY)));
			}
		}
		n = rand() % ((int)(k * errMap.size()));
		emhit = errMap.begin();
		l = 0;
		while (l < n)
		{
			l++;
			++emhit;
		}
		x = (*emhit).second.first;
		y = (*emhit).second.second;

	} while (y < output_h);

	isFilled = true;
}

// Refinement with random placement (filling has to be run first!)
// bandX, bandY: define the part of the input patch that can lie outside
// the output texture.
// iter: number of iterations.
// blending: if true, blending enabled.
void GCTexture::randomRefinement(int bandX, int bandY, int iter, bool blending)
{
	if (!isFilled)
	{
		std::cout << "Initial texture not filled!\n";
		return;
	}

	std::cout << "---- Texture refinement ----\n";

	int x, y;

	// Insert new patches at random positions
	for (int k = 0; k < iter; k++)
	{
		x = -(2 * bandX) + (std::rand() % (output_w + bandX));
		y = -(2 * bandY) + (std::rand() % (output_h + bandY));
		
		insertPatch(x, y, false, blending);
		updateSeamsMaxError();

		patchCount++;

#ifdef _DEBUG
		std::cout << "Patch " << patchCount << " inserted\n";

		// Fill output image with current synthesis result
		fillOutputImage();

		// Write output image in a file
		writeOutputImagePatch();
#endif
	}
}

// Refinement based on seam errors. The best location for the
// input image is chosen based on SSD error. The center of the
// whole input image is translated in the region around the
// max error seam point. This region is forced to be inserted into
// the output texture. The size of the region around the max error
// seam point should not exceed the size of the input image.
// maxIter: maximum number of iterations.
// radius: radius around the max error seam point.
// blending: if true, blending enabled.
void GCTexture::seamsErrorCenterRefinement(int maxIter, int radius, bool blending)
{
	if (!isFilled)
	{
		std::cout << "Initial texture not filled!\n";
		return;
	}

//	std::cout << "---- Texture refinement ----\n";

	int maxErrX=0, maxErrY=0;
	int x, y;
	double err;
	double minErr;
	int minErrX=0, minErrY=0;

	updateSeamsMaxError();

	for (int k = 0; k < maxIter; k++)
	{
		if (maxErrNodeNbGlobal != -1)
		{
			minErr = std::numeric_limits<double>::max();
			minErrX = -1;
			minErrY = -1;

			getPixelGlobal(maxErrNodeNbGlobal, maxErrX, maxErrY);

			for (int ry = -(radius / 2); ry < (radius / 2); ry += ssdStep)
			{
				for (int rx = -(radius / 2); rx < (radius / 2); rx += ssdStep)
				{
					// Compute patch position
					y = (maxErrY + ry) - (input_h / 2);
					x = (maxErrX + rx) - (input_w / 2);

					// Compute the SSD for the current position
					err = computeSSD(x, y);

					if (err < minErr)
					{
						minErr = err;
						minErrY = y;
						minErrX = x;
					}
				}
			}
		}

#ifdef _DEBUG
		std::cout << "Patch location (" << minErrX << ", " << minErrY << ")\n";
#endif

		// Only a subpart of the input image is inserted.
		// Origin in input image sapce: minErrX, minErrY
		// Size: 3*radius * 3*radius

		// Fixed subpatch insertion

		// Cropped input image origin in input image space
		int croppedX = std::max(0, minErrX) - minErrX;
		int croppedY = std::max(0, minErrY) - minErrY;

		// Origin of fixed supatch in cropped input image space
		int inputX = maxErrX - (radius / 2) - minErrX - croppedX;
		int inputY = maxErrY - (radius / 2) - minErrY - croppedY;

		// Updated input image origin in output texture space
		minErrX = std::max(0, minErrX);
		minErrY = std::max(0, minErrY);

		// Insert the new patch at min SSD location
		//insertPatch(minErrX, minErrY, false, blending, radius, inputX, inputY, minErrX, minErrY);
		insertPatch(minErrX, minErrY, false, blending, radius, inputX, inputY, minErrX, minErrY);
		updateSeamsMaxError();

		patchCount++;

#ifdef _DEBUG
		std::cout << "Patch " << patchCount << " inserted\n";

		// Fill output image with current synthesis result
		fillOutputImage();

		// Write output image in a file
		writeOutputImagePatch();
#endif
	}
}

void GCTexture::seamsErrorFullPatchRefinement(int maxIter, int radius, bool blending)
{
	if (!isFilled)
	{
		std::cout << "Initial texture not filled!\n";
		return;
	}

	std::cout << "---- Texture refinement ----\n";

	int maxErrX, maxErrY;
	int x, y;
	double err;
	double minErr;
	int minErrI, minErrJ;

	updateSeamsMaxError();

	for (int k = 0; k < maxIter; k++)
	{
		if (maxErrNodeNbGlobal != -1)
		{
			minErr = std::numeric_limits<double>::max();
			minErrI = -1;
			minErrJ = -1;

			getPixelGlobal(maxErrNodeNbGlobal, maxErrX, maxErrY);

			// Fixed region position
			x = maxErrX - (radius / 2);
			y = maxErrY - (radius / 2);

			// Cropped fixed region
			int v1X = std::max(0, x);
			int v1Y = std::max(0, y);
			int v2X = std::min(output_w - 1, x + radius - 1);
			int v2Y = std::min(output_h - 1, y + radius - 1);

			int cropped_w = v2X - v1X + 1;
			int cropped_h = v2Y - v1Y + 1;

			// Sweeping of the input image with the fixed region in the output
			// texture
			for (int j = cropped_h; j < (input_h - 2 * cropped_h); j += ssdStep)
			{
				for (int i = cropped_w; i < (input_w - 2 * cropped_w); i += ssdStep)
				{
					// Compute the SSD for the current position
					err = computeLocalSSD(v1X, v1Y, i, j, cropped_w, cropped_h);

					if (err < minErr)
					{
						minErr = err;
						minErrI = i;
						minErrJ = j;
					}
				}
			}

#ifdef _DEBUG
			std::cout << "Patch location in input image (" << minErrI << ", " << minErrJ << ")\n";
#endif

			// Fixed subpatch insertion
			int minErrX = v1X - minErrI;
			int minErrY = v1Y - minErrJ;

			// Cropped input image origin in input image space
			int croppedX = std::max(0, minErrX);
			int croppedY = std::max(0, minErrY);

			// Origin in cropped input image space
			int inputX = maxErrX - (radius / 2) - croppedX;
			int inputY = maxErrY - (radius / 2) - croppedY;

			// Updated input image origin in output texture space
			minErrX = std::max(0, minErrX);
			minErrY = std::max(0, minErrY);

			// Insert the new patch at min SSD location
			insertPatch(minErrX, minErrY, false, blending, radius, inputX, inputY, minErrX, minErrY);

			updateSeamsMaxError();

			patchCount++;

#ifdef _DEBUG
			std::cout << "Patch " << patchCount << " inserted\n";

			// Fill output image with current synthesis result
			fillOutputImage();

			// Write output image in a file
			writeOutputImagePatch();
#endif
		}
	}
}

void GCTexture::seamsErrorSubPatchRefinementA(int maxIter, int radius, bool blending)
{
	if (!isFilled)
	{
		std::cout << "Initial texture not filled!\n";
		return;
	}

	std::cout << "---- Texture refinement ----\n";

	int maxErrX, maxErrY;
	int x, y;
	double err;
	double minErr;
	int minErrI, minErrJ;

	updateSeamsMaxError();

	for (int k = 0; k < maxIter; k++)
	{
		if (maxErrNodeNbGlobal != -1)
		{
			minErr = std::numeric_limits<double>::max();
			minErrI = -1;
			minErrJ = -1;

			getPixelGlobal(maxErrNodeNbGlobal, maxErrX, maxErrY);

			// Fixed region position
			x = maxErrX - (radius / 2);
			y = maxErrY - (radius / 2);

			// Cropped fixed region
			int v1X = std::max(0, x);
			int v1Y = std::max(0, y);
			int v2X = std::min(output_w - 1, x + radius - 1);
			int v2Y = std::min(output_h - 1, y + radius - 1);

			int cropped_w = v2X - v1X + 1;
			int cropped_h = v2Y - v1Y + 1;

			// Sweeping of the input image with the fixed region in the output
			// texture
			for (int j = cropped_h; j < (input_h - 2 * cropped_h); j += ssdStep)
			{
				for (int i = cropped_w; i < (input_w - 2 * cropped_w); i += ssdStep)
				{

					// Compute the SSD for the current position
					err = computeLocalSSD(v1X, v1Y, i, j, cropped_w, cropped_h);

					if (err < minErr)
					{
						minErr = err;
						minErrI = i;
						minErrJ = j;
					}
				}
			}

#ifdef _DEBUG
			std::cout << "Patch location in input image (" << minErrI << ", " << minErrJ << ")\n";
#endif

			int minErrX = v1X - radius;
			int minErrY = v1Y - radius;

			int minSX = 2 * radius + cropped_w;
			int minSY = 2 * radius + cropped_h;

			int croppedX = std::max(0, minErrX);
			int croppedY = std::max(0, minErrY);

			// Origin in cropped input image space
			int inputX = maxErrX - (radius / 2) - croppedX;
			int inputY = maxErrY - (radius / 2) - croppedY;

			// Updated input image origin in output texture space
			minErrX = std::max(0, minErrX);
			minErrY = std::max(0, minErrY);

			// Insert the new patch at min SSD location
			insertPatch(minErrX, minErrY, minErrI - radius, minErrJ - radius,
				        minSX, minSY, false, blending, radius,
						inputX, inputY, minErrX, minErrY);

			updateSeamsMaxError();

			patchCount++;

#ifdef _DEBUG
			std::cout << "Patch " << patchCount << " inserted\n";

			// Fill output image with current synthesis result
			fillOutputImage();

			// Write output image in a file
			writeOutputImagePatch();
#endif
		}
	}

	// Fill output image with current synthesis result
	fillOutputImage();

	// Write output image in a file
	writeOutputImagePatch();
}

void GCTexture::seamsErrorSubPatchRefinement(int maxIter, int radius, bool blending)
{
	if (!isFilled)
	{
		std::cout << "Initial texture not filled!\n";
		return;
	}

	std::cout << "---- Texture refinement ----\n";

	int inputRadius = radius;

	int maxErrX, maxErrY;
	int x, y;
	double err;
	double minErr;
	int minErrI, minErrJ;

	updateSeamsMaxError();

	for (int k = 0; k < maxIter; k++)
	{
		if (maxErrNodeNbGlobal != -1)
		{
			// Perturbation
			radius = (inputRadius - 2) + (rand() % 5);
			std::cout << "radius: " << radius << std::endl;
			
			minErr = std::numeric_limits<double>::max();
			minErrI = -1;
			minErrJ = -1;

			getPixelGlobal(maxErrNodeNbGlobal, maxErrX, maxErrY);

			// Fixed region position
			x = maxErrX - (radius / 2);
			y = maxErrY - (radius / 2);

			// Cropped fixed region
			int v1X = std::max(0, x);
			int v1Y = std::max(0, y);
			int v2X = std::min(output_w - 1, x + radius - 1);
			int v2Y = std::min(output_h - 1, y + radius - 1);

			int cropped_w = v2X - v1X + 1;
			int cropped_h = v2Y - v1Y + 1;

			std::cout << v1X - radius << ", " << v1Y - radius << std::endl;

			// Sweeping of the input image with the fixed region in the output
			// texture
			
			for (int j = 0; j < (input_h - 3 * radius); j += ssdStep)
			//for (int j = radius; j < (input_h - 2 * radius); j += ssdStep)
			{
				for (int i = 0; i < (input_w - 3 * radius); i += ssdStep)
				//for (int i = radius; i < (input_w - 2 * radius); i += ssdStep)
				{
					// Compute the SSD for the current position
					err = computeLocalSSD(v1X - radius, v1Y - radius, i, j, 2 * radius + cropped_w, 2 * radius + cropped_h);

					if (err < minErr)
					{
						minErr = err;
						minErrI = i;
						minErrJ = j;
					}
				}
			}

			minErrI += radius;
			minErrJ += radius;

#ifdef _DEBUG
			std::cout << "Patch location in input image (" << minErrI << ", " << minErrJ << ")\n";
#endif

			int minErrX = v1X - radius;
			int minErrY = v1Y - radius;

			int minSX = 2 * radius + cropped_w;
			int minSY = 2 * radius + cropped_h;

			int croppedX = std::max(0, minErrX);
			int croppedY = std::max(0, minErrY);

			// Origin in cropped input image space
			int inputX = maxErrX - (radius / 2) - croppedX;
			int inputY = maxErrY - (radius / 2) - croppedY;

			// Updated input image origin in output texture space
			minErrX = std::max(0, minErrX);
			minErrY = std::max(0, minErrY);

			// Insert the new patch at min SSD location
			insertPatch(minErrX, minErrY, minErrI - radius, minErrJ - radius,
				minSX, minSY, false, blending, radius,
				inputX, inputY, minErrX, minErrY);

			updateSeamsMaxError();

			patchCount++;

#ifdef _DEBUG
			std::cout << "Patch " << patchCount << " inserted\n";

			// Fill output image with current synthesis result
			fillOutputImage();

			// Write output image in a file
			writeOutputImagePatch();
#endif
		}
	}

	// Fill output image with current synthesis result
	fillOutputImage();

	// Write output image in a file
	writeOutputImagePatch();
}

// Draws the seams in red
void GCTexture::revealSeams()
{
	int nodeNbGlobal;

	for (int j = 0; j < output_h; j++)
	{
		for (int i = 0; i < output_w; i++)
		{
			nodeNbGlobal = getNodeNbGlobal(0, 0, i, j);
			if (globalNode[nodeNbGlobal].notEmpty())
			{
				if (globalNode[nodeNbGlobal].onSeamRight()
					|| globalNode[nodeNbGlobal].onSeamBottom())
				{
					for (int sj = -seamSize; sj < seamSize; sj++)
					{
						for (int si = -seamSize; si < seamSize; si++)
						{
							int tj = j + sj;
							int ti = i + si;
							if ((ti >= 0 && ti < output_w) && (tj >= 0 && tj < output_h))
							{
								outputImage.pixelAbsolute(ti, tj) = RGBu8(255, 0, 0);
							}
						}
					}
				}
			}
		}
	}
}

// Draws the seams with error depicted by colors ranging from
// green (low error) to red (high error)
void GCTexture::revealSeamsError()
{
	int nodeNbGlobal;

	for (int j = 0; j < output_h; j++)
	{
		for (int i = 0; i < output_w; i++)
		{
			nodeNbGlobal = getNodeNbGlobal(0, 0, i, j);
			if (globalNode[nodeNbGlobal].notEmpty())
			{
				if (globalNode[nodeNbGlobal].onSeamRight()
					|| globalNode[nodeNbGlobal].onSeamBottom())
				{
					for (int sj = -seamSize; sj < seamSize; sj++)
					{
						for (int si = -seamSize; si < seamSize; si++)
						{
							int tj = j + sj;
							int ti = i + si;
							if ((ti >= 0 && ti < output_w) && (tj >= 0 && tj < output_h))
							{
								if (globalNode[nodeNbGlobal].onSeamRight())
								{
									if (globalNode[nodeNbGlobal].onSeamBottom())
									{
										if (globalNode[nodeNbGlobal].getRightCost() < globalNode[nodeNbGlobal].getBottomCost())
										{
											int c = (int)((globalNode[nodeNbGlobal].getBottomCost() / two_sqrt3) * 255.0 * 6.0);
											if (c > 255) c = 255;
											outputImage.pixelAbsolute(ti, tj) = RGBu8(c, 255 - c, 0);
										}
										else
										{
											int c = (int)((globalNode[nodeNbGlobal].getRightCost() / two_sqrt3) * 255.0 * 6.0);
											if (c > 255) c = 255;
											outputImage.pixelAbsolute(ti, tj) = RGBu8(c, 255 - c, 0);
										}
									}
									else
									{
										int c = (int)((globalNode[nodeNbGlobal].getRightCost() / two_sqrt3) * 255.0 * 6.0);
										if (c > 255) c = 255;
										outputImage.pixelAbsolute(ti, tj) = RGBu8(c, 255 - c, 0);
									}
								}
								else
								{
									int c = (int)((globalNode[nodeNbGlobal].getBottomCost() / two_sqrt3) * 255.0 * 6.0);
									if (c > 255) c = 255;
									outputImage.pixelAbsolute(ti, tj) = RGBu8(c, 255 - c, 0);
								}
							}
						}
					}
				}
			}
		}
	}
}

void GCTexture::writeOutputImage(const std::string& fileName, const std::string& fileNameS)
{
	fillOutputImage();
	outputImage.save(fileName);
	revealSeamsError();
	outputImage.save(fileNameS);
}

// Normalized SSD in the overlapping region when inserting the
// input image at position (x, y)
double GCTexture::computeSSD(int x, int y)
{
	// Bounding box of the region of the new patch included
	// in the output texture (coordinates in input image
	// space, i.e. with origin at input image top left corner)
	int v1X = std::max(0, x) - x;
	int v1Y = std::max(0, y) - y;
	int v2X = std::min(output_w - 1, x + input_w - 1) - x;
	int v2Y = std::min(output_h - 1, y + input_h - 1) - y;

	croppedInput_w = v2X - v1X + 1;
	croppedInput_h = v2Y - v1Y + 1;

	// Update origin coordinates
	x = std::max(0, x);
	y = std::max(0, y);

	double err = 0.0;
	int nbpix = 0;
	int nodeNbGlobal;

	for (int j = v1Y, cj = 0; j <= v2Y; j++, cj++)
	{
		for (int i = v1X, ci = 0; i <= v2X; i++, ci++)
		{
			nodeNbGlobal = getNodeNbGlobal(x, y, ci, cj);
			if (!globalNode[nodeNbGlobal].empty())
			{
				err += sqDistRGB(globalNode[nodeNbGlobal].getColor(), inputImage.pixelAbsolute(i, j));
				nbpix++;
			}
		}
	}

	return err / nbpix;
}

// (x, y): position in output texture
// (inputX, inputY): position in input image
// size: size of the subpatch
double GCTexture::computeLocalSSD(int x, int y, int inputX, int inputY, int sizeX, int sizeY)
{
	// Region ton consider in the output texture:
	// from (x, y) to (x+sizeX, y+sizeY)
	// This region has been cropped to match the output
	// space

	double err = 0.0;
	int nbpix = 0;
	int nodeNbGlobal;

	for (int j = 0; j < sizeY; j++)
	{
		for (int i = 0; i < sizeX; i++)
		{
			nodeNbGlobal = getNodeNbGlobal(x, y, i, j);
			if (!globalNode[nodeNbGlobal].empty())
			{
				err += sqDistRGB(globalNode[nodeNbGlobal].getColor(), inputImage.pixelAbsolute(inputX + i, inputY + j));
				nbpix++;
			}
		}
	}

	if (nbpix == 0)
	{
		std::cout << "No pixel!\n";
		getchar();
	}

	return err / nbpix;
}
