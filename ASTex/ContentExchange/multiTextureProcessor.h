#ifndef __CTEXCH_MULTI_TEXTURE_PROCESSOR_H__
#define __CTEXCH_MULTI_TEXTURE_PROCESSOR_H__

#include "patchProcessor.h"
#include "ASTex/Stamping/sampler.h"

namespace ASTex
{

namespace ContentExchange
{

template<typename I>
/**
 * @brief The MultiTextureProcessor class is a mod made for
 * processing and storing multi-texture content exchange elements.
 * Feed it a tile, 0 to n textures (precise which are periodic textures)
 * and 0 to m transformations (Eigen matrices),
 * choose how you want patches and contents to be computed, then use
 * generate() to get an output tiling, or export the class for online rendering.
 */
class MultiTextureProcessor : public ASTex::ContentExchange::PatchProcessor<I>
{
public:

	MultiTextureProcessor();
	MultiTextureProcessor(const I& image);

	~MultiTextureProcessor();

	//get

	void addTextureToPool(const I& texture, bool periodicity);
	void addTransformationToPool(const Eigen::Matrix2d transformation);
	void contents_explorePools(bool gradientManipulation=false);
	void contents_randomPools(bool gradientManipulation=false);

private:

	std::vector<std::pair<I, bool>> m_texturePool;
	std::vector<Eigen::Matrix2d> m_transformationPool;
};

template<typename I>
MultiTextureProcessor<I>::MultiTextureProcessor():
	PatchProcessor<I>(),
	m_texturePool(),
	m_transformationPool()
{}

template<typename I>
MultiTextureProcessor<I>::MultiTextureProcessor(const I& tile):
	PatchProcessor<I>(tile),
	m_texturePool(),
	m_transformationPool()
{}

template<typename I>
MultiTextureProcessor<I>::~MultiTextureProcessor()
{}

template<typename I>
void MultiTextureProcessor<I>::addTextureToPool(const I& texture, bool periodicity)
{
	m_texturePool.push_back(std::make_pair(texture, periodicity));
}

template<typename I>
void MultiTextureProcessor<I>::addTransformationToPool(const Eigen::Matrix2d transformation)
{
	m_transformationPool.push_back(transformation);
}

template<typename I>
void MultiTextureProcessor<I>::contents_explorePools(bool gradientManipulation)
{
	assert(this->m_patches.size()>0 &&
		   "PatchProcessor::contents_explorePools: a patch initializer must be called before being able to choose contents");
	assert(this->m_mipmapMode != NO_FILTER && "PatchProcessor::contents_explorePools: Mipmap mode cannot be NO_FILTER even when ignoring mip-maps for this function");
	srand(this->m_seed);
	I selectedTexture, shiftedTile;
	bool textureIsPeriodic = true;
	Eigen::Matrix3d selectedTransformation;
	Mipmap<I> inputTileMipmap(this->m_tile);
	inputTileMipmap.generate();

	//some definitions...
	typedef struct
	{
		Mipmap<I> mipmap;
		MipmapBooleanImage oobMap; //< out of bounds map
		bool periodicity;
	} expTexture_t; //Like "exploratory texture" but shorter

	typedef struct
	{
		PixelPos origin; //< the origin in the transformed tile, for various resolutions...
		unsigned index; //< index of the transformed tile in the mipmap vector
		double score;
	} contentComparator_t;

	class CompareContents
	{
	public:
		CompareContents() {}
		bool operator()(const contentComparator_t &object, const contentComparator_t &other)
		{
			return object.score > other.score;
		}
	};

	std::vector<expTexture_t> expTexturesVector;
	using PriorityQueueComparator = std::priority_queue<contentComparator_t, std::vector<contentComparator_t>, CompareContents>;

	auto compareHigherResRegion = [&] (const Patch<I> &patch, expTexture_t expTexture, PixelPos origin, unsigned level) -> double
	{
		PixelPos patchOrigin = patch.originAt(level, level);
		const I& variableResContent = expTexture.mipmap.mipmap(level, level);
		unsigned width = patch.mipmap(level, level).width(), height = patch.mipmap(level, level).height();
		bool suited = true;
		double score = 0;
		for(unsigned x=0; x<width && suited; ++x)
			for(unsigned y=0; y<height && suited; ++y)
			{
				typename ImageAlphad::PixelType pix = patch.mipmap(level, level).pixelAbsolute(x, y);
				//alpha map space: x, y. Orinal tile space: patchOrigin+(x, y). expTile space: origin+(x, y).
				if(pix>0)
				{
					if(origin[0]+x >= variableResContent.width() || origin[1]+y >= variableResContent.height()
						|| (!expTexture.periodicity &&
							expTexture.oobMap.mipmap(level, level).pixelAbsolute(origin[0]+x, origin[1]+y))
					)
						suited = false;
					else
					{
						if(pix<1) //on the border of the patch.
						{
							const I &levelResInputTile = inputTileMipmap.mipmap(level, level);
							score += mse(	levelResInputTile, expTexture.mipmap.mipmap(level, level),
											(patchOrigin[0]+x)%levelResInputTile.width(), (patchOrigin[1]+y)%levelResInputTile.height(),
											origin[0]+x, origin[1]+y, 0, false);
						}
					}
				}
			}
		if(!suited)
			score=std::numeric_limits<double>::infinity();
		return score;
	};

	for(unsigned textureIndex=0; textureIndex<=m_texturePool.size(); ++textureIndex)
	{
		for(unsigned transformIndex=0; transformIndex<=m_transformationPool.size(); ++transformIndex)
		{
			if(textureIndex==m_texturePool.size())
			{
				selectedTexture.initItk(this->m_tile.width(), this->m_tile.height());
				selectedTexture.copy_pixels(this->m_tile);
				textureIsPeriodic=true;
			}
			else
			{
				selectedTexture.initItk(m_texturePool[textureIndex].first.width(),
										m_texturePool[textureIndex].first.height());
				selectedTexture.copy_pixels(m_texturePool[textureIndex].first);
				textureIsPeriodic = m_texturePool[textureIndex].second;
			}
			Eigen::Matrix2d rotationMatrix2d;
			Eigen::Matrix3d rotationMatrix3d;
			if(transformIndex==m_transformationPool.size())
			{
				selectedTransformation=Eigen::Matrix3d::Identity();
			}
			else
			{
				Eigen::Matrix2d transformation2d = m_transformationPool[transformIndex];
				selectedTransformation.row(0) << transformation2d(0, 0), transformation2d(0, 1), 0;
				selectedTransformation.row(1) << transformation2d(1, 0), transformation2d(1, 1), 0;
				selectedTransformation.row(2) << 0, 0, 1;
				selectedTransformation=Eigen::Affine2d(Eigen::Translation2d(selectedTexture.width()/2.0, selectedTexture.height()/2.0)).matrix()
						* selectedTransformation
						* Eigen::Affine2d(Eigen::Translation2d(-selectedTexture.width()/2.0, -selectedTexture.height()/2.0)).matrix();
				if(gradientManipulation)
				{
					Eigen::Affine2d transform;
					transform.matrix() = selectedTransformation;
					rotationMatrix2d = transform.rotation();
					rotationMatrix3d.row(0) << rotationMatrix2d(0, 0), rotationMatrix2d(0, 1), 0;
					rotationMatrix3d.row(1) << rotationMatrix2d(1, 0), rotationMatrix2d(1, 1), 0;
					rotationMatrix3d.row(2) << 0, 0, 1;
				}
			}
			Eigen::Matrix3d invTransformation = selectedTransformation.inverse();
			//we need a relatively large tile because rotations and scaling can cause the tile to go out of bounds.
			//I put a factor of 1.6 for the size of the tile because
			//it can encompasse a pi/2-rotated image with a small upscaling on top of it.
			shiftedTile.initItk(1.6*selectedTexture.width(), 1.6*selectedTexture.height(), true);
			nullify_image(shiftedTile); //because this thing ^ does not fucking work
			expTexture_t exploratoryTexture;
			if(textureIsPeriodic)
			{
				shiftedTile.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
				{
					Eigen::Vector3d v;
					v[0]=x;
					v[1]=y;
					v[2]=1;
					Eigen::Vector3d vt = invTransformation * v;
					//interpolation?
					/*pix = bilinear_interpolation(selectedTexture,
												 vt[0] + randomShiftX + selectedTexture.width()*4,
												 vt[1] + randomShiftY + selectedTexture.height()*4,
												 textureIsPeriodic);*/
					pix = selectedTexture.pixelAbsolute(int(vt[0] + selectedTexture.width()*12)%selectedTexture.width(),
														int(vt[1] + selectedTexture.height()*12)%selectedTexture.height());
					if(gradientManipulation)
					{
						Eigen::Vector3d rotatedGradientPixel = rotationMatrix3d*(Eigen::Vector3d(pix[0], pix[1], 1) - Eigen::Vector3d(0.5, 0.5, 0))
																+ Eigen::Vector3d(0.5, 0.5, 0);
						pix[0] = rotatedGradientPixel[0];
						pix[1] = rotatedGradientPixel[1];
					}
				});

				exploratoryTexture.periodicity = true;
			}
			else
			{
				PixelPos posMin, posMax;
				posMin[0] = shiftedTile.width()/2;
				posMax[0] = shiftedTile.width()/2;
				posMin[1] = shiftedTile.height()/2;
				posMax[1] = shiftedTile.height()/2;
				//first compute the transformed texture explicitely
				ImageGrayb outOfBoundsMap;
				outOfBoundsMap.initItk(shiftedTile.width(), shiftedTile.height());
				shiftedTile.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
				{
					Eigen::Vector3d v, vt;
					v[0]=x-(shiftedTile.width()-selectedTexture.width())/2.0;
					v[1]=y-(shiftedTile.height()-selectedTexture.height())/2.0;
					v[2]=1;
					vt = invTransformation * v;
					if(vt[0]<0 || vt[0]>=selectedTexture.width() || vt[1]<0 || vt[1]>=selectedTexture.height())
					{
						pix = this->ms_zero;
						outOfBoundsMap.pixelAbsolute(x, y)=true;
					}
					else
					{
						pix = selectedTexture.pixelAbsolute(vt[0], vt[1]);
						outOfBoundsMap.pixelAbsolute(x, y)=false;
						posMin[0] = std::min(PixelPos::IndexValueType(x), posMin[0]);
						posMin[1] = std::min(PixelPos::IndexValueType(y), posMin[1]);
						posMax[0] = std::max(PixelPos::IndexValueType(x), posMax[0]);
						posMax[1] = std::max(PixelPos::IndexValueType(y), posMax[1]);
						if(gradientManipulation)
						{
							Eigen::Vector3d rotatedGradientPixel = rotationMatrix3d*(Eigen::Vector3d(pix[0], pix[1], 1)
																	- Eigen::Vector3d(0.5, 0.5, 0))
																	+ Eigen::Vector3d(0.5, 0.5, 0);
							pix[0] = rotatedGradientPixel[0];
							pix[1] = rotatedGradientPixel[1];
						}
					}
				});
				I croppedTile;
				crop_image(shiftedTile, croppedTile, posMin[0], posMax[0], posMin[1], posMax[1]);
				shiftedTile = croppedTile;
				ImageGrayb croppedOutOfBoundsMap;
				crop_image(outOfBoundsMap, croppedOutOfBoundsMap, posMin[0], posMax[0], posMin[1], posMax[1]);
				outOfBoundsMap = croppedOutOfBoundsMap;
				exploratoryTexture.periodicity = false;
				exploratoryTexture.oobMap.setTexture(outOfBoundsMap);
				exploratoryTexture.oobMap.generate();
			}
			IO::save01_in_u8(shiftedTile, std::string("/home/nlutz/shiftedTile_tex") + std::to_string(textureIndex)
																			+ "_tr" + std::to_string(transformIndex) + ".png");
			Mipmap<I> tileMipmap(shiftedTile);
			tileMipmap.setMode(ISOTROPIC);
			tileMipmap.setMaxPowReductionLevel(0); //<optimize if the bounding box of the transformed image is computed
			tileMipmap.generate();
			exploratoryTexture.mipmap = tileMipmap;
			expTexturesVector.push_back(exploratoryTexture);
		}
	}
	for(unsigned p=0; p<this->m_patches.size(); ++p)
	{
		std::cout << p << std::endl;
		Patch<I> &patch=this->m_patches[p];
		PriorityQueueComparator contentComparatorPQ;

		unsigned nbKeptContents = std::max(nbContents(), expTexturesVector.size()); //TODO I guess, it feels wrong to hard code this
		std::cout << "nbKeptContents=" << nbKeptContents << std::endl;
		int maxK = std::max(0, int(std::max(expTexturesVector.back().mipmap.numberMipmapsWidth(), expTexturesVector.back().mipmap.numberMipmapsHeight()) - 5)); //TODO, again.
		for(int k=maxK; k>=0; --k) //^ back is where the tile with identity transform is (rember if modified)
		{
			PriorityQueueComparator ccPQCopy;
			if(k!=maxK)
			{
				while(!contentComparatorPQ.empty())
				{
					contentComparator_t contentComparator;
					const contentComparator_t &cc = contentComparatorPQ.top();
					PixelPos newOrigin = cc.origin;
					double scores[4], bestScore=std::numeric_limits<double>::infinity(), secondBestScore=std::numeric_limits<double>::infinity();
					unsigned bestIndex=0, secondBestIndex=0;
					PixelPos origins[4];
					//std::cout << "Score before: " << contentComparator.score << std::endl;
					newOrigin[0] *= 2;
					newOrigin[1] *= 2;
					origins[0] = newOrigin;
					scores[0] = compareHigherResRegion(patch, expTexturesVector[cc.index], newOrigin, k);
					++newOrigin[0];
					origins[1] = newOrigin;
					scores[1] = compareHigherResRegion(patch, expTexturesVector[cc.index], newOrigin, k);
					++newOrigin[1];
					origins[2] = newOrigin;
					scores[2] = compareHigherResRegion(patch, expTexturesVector[cc.index], newOrigin, k);
					--newOrigin[0];
					origins[3] = newOrigin;
					scores[3] = compareHigherResRegion(patch, expTexturesVector[cc.index], newOrigin, k);
					//std::cout << "Scores after: " << scores[0] << ' ' << scores[1] << ' ' << scores[2] << ' ' << scores[3] << std::endl;


					//find the two best indices. This is hard coded, might want to change later if absolutely necessary
					for(unsigned i=0; i<4; ++i)
					{
						if(scores[i] < bestScore)
						{
							bestIndex = i;
							bestScore = scores[i];
						}
					}
					for(unsigned i=0; i<4; ++i)
					{
						if(i!=bestIndex && scores[i] < secondBestScore)
						{
							secondBestIndex = i;
							secondBestScore = scores[i];
						}
					}
					contentComparatorPQ.pop();
					contentComparator.index = cc.index;
					contentComparator.score = bestScore;
					contentComparator.origin = origins[bestIndex];
					ccPQCopy.push(contentComparator);
					contentComparator.score = secondBestScore;
					contentComparator.origin = origins[secondBestIndex];
					//ccPQCopy.push(contentComparator);
				}
				contentComparatorPQ = ccPQCopy;
			}
			else
			{	//first step: a proportional number of contents per texture.
				unsigned nbContentsPerTexture = (nbContents()/expTexturesVector.size())*2 + 1;
				for(unsigned i=0; i<expTexturesVector.size(); ++i)
				{
					Stamping::SamplerPoissonGrid sampler;
					sampler.setNbPoints(1024); //hard-coded again... TODO
					std::vector<Eigen::Vector2f> samples = sampler.generate();
					for(unsigned j=0; j<samples.size(); ++j)
					{
						PixelPos origin;
						double score;
						contentComparator_t contentComparator;
						origin[0] = samples[j][0]*(expTexturesVector[i].mipmap.mipmap(k, k).width()-patch.mipmap(k, k).width()-1);
						origin[1] = samples[j][1]*(expTexturesVector[i].mipmap.mipmap(k, k).height()-patch.mipmap(k, k).height()-1);
						score = compareHigherResRegion(patch, expTexturesVector[i], origin, k);
						//std::cout << origin << std::endl;
						contentComparator.index = i;
						contentComparator.origin = origin;
						contentComparator.score = score;
						ccPQCopy.push(contentComparator);
					}
					while(ccPQCopy.size()>0)
					{
						if(ccPQCopy.size()<=nbContentsPerTexture)
						{
							contentComparatorPQ.push(ccPQCopy.top());
						}
						ccPQCopy.pop();
					}
				}
			}
			while(contentComparatorPQ.size() >= nbKeptContents)
				contentComparatorPQ.pop();
			nbKeptContents *= 1.5; //hard-coded. TODO. Also check if this doesn't generate garbage.
		}
		while(contentComparatorPQ.size() > nbContents())
			contentComparatorPQ.pop();
		//Now we can build the contents out of the remains!
		while(contentComparatorPQ.size() > 0)
		{
			contentComparator_t contentComparator = contentComparatorPQ.top();
			//IO::save01_in_u8(expTexturesVector[contentComparator.index].mipmap.texture(), std::string("/home/nlutz/shiftedTileEX_p") + std::to_string(p) + "_vec" + std::to_string(contentComparatorPQ.size()) + ".png");
			std::cout << contentComparator.score << ", " << contentComparator.origin << std::endl;
			I contentTile;
			contentTile.initItk(this->m_tile.width(), this->m_tile.height());
			for(int x=0; x<patch.alphaMap().width(); ++x)
				for(int y=0; y<patch.alphaMap().height(); ++y)
				{
					unsigned xPatch = (x+patch.originAt(0, 0)[0])%this->m_tile.width();
					unsigned yPatch = (y+patch.originAt(0, 0)[1])%this->m_tile.height();
					unsigned xContent = x+contentComparator.origin[0];
					unsigned yContent = y+contentComparator.origin[1];
					if(patch.alphaMap().pixelAbsolute(x, y)>0) //very important because some pixels may be out of bounds
						contentTile.pixelAbsolute(xPatch, yPatch) = expTexturesVector[contentComparator.index].mipmap.texture().pixelAbsolute(xContent, yContent);
				}
			contentComparatorPQ.pop();
			Content<I> content(contentTile, patch);
			IO::save01_in_u8(contentTile, std::string("/home/nlutz/content_p") + std::to_string(p) + "_c" + std::to_string(contentComparatorPQ.size()) + ".png");
			patch.addContent(content);
		}
	}
}

template<typename I>
void MultiTextureProcessor<I>::contents_randomPools(bool gradientManipulation)
{
	assert(this->m_patches.size()>0 &&
		   "PatchProcessor::contents_initRandom: a patch initializer must be called before being able to choose contents");
	srand(this->m_seed);
	unsigned randomTextureIndex=0;
	unsigned randomTransformationIndex=0;
	unsigned i,j;
	I selectedTexture;
	bool textureIsPeriodic = true;
	Eigen::Matrix3d selectedTransformation;
	I shiftedTile;
	shiftedTile.initItk(this->m_tile.width(), this->m_tile.height());
	for(j=0; j<this->m_patches.size(); ++j)
	{
		Patch<I> &patch=this->m_patches[j]; //we establish a bounding box of the patch to avoid out of bounds later on
		PixelPos patchOrigin = patch.originAt(0, 0);

		unsigned patchBBWidth=patch.alphaMap().width(), patchBBHeight=patch.alphaMap().height();

		for(i=0; i<this->m_nbContentsPerPatch-1; ++i)
		{
			randomTextureIndex = rand()%(m_texturePool.size()+1);
			if(randomTextureIndex==m_texturePool.size())
			{
				selectedTexture.initItk(this->m_tile.width(), this->m_tile.height());
				selectedTexture.copy_pixels(this->m_tile);
				textureIsPeriodic=true;
			}
			else
			{
				selectedTexture.initItk(m_texturePool[randomTextureIndex].first.width(),
										m_texturePool[randomTextureIndex].first.height());
				selectedTexture.copy_pixels(m_texturePool[randomTextureIndex].first);
				textureIsPeriodic = m_texturePool[randomTextureIndex].second;
			}
			randomTransformationIndex = rand()%(m_transformationPool.size()+1);
			Eigen::Matrix2d rotationMatrix2d;
			Eigen::Matrix3d rotationMatrix3d;
			if(randomTransformationIndex==m_transformationPool.size())
			{
				selectedTransformation=Eigen::Matrix3d::Identity();
			}
			else
			{
				Eigen::Matrix3d randomTransformation;
				Eigen::Matrix2d transformation2d = m_transformationPool[randomTransformationIndex];
				randomTransformation.row(0) << transformation2d(0, 0), transformation2d(0, 1), 0;
				randomTransformation.row(1) << transformation2d(1, 0), transformation2d(1, 1), 0;
				randomTransformation.row(2) << 0, 0, 1;
				selectedTransformation=Eigen::Affine2d(Eigen::Translation2d(selectedTexture.width()/2.0, selectedTexture.height()/2.0)).matrix()
						* randomTransformation
						* Eigen::Affine2d(Eigen::Translation2d(-selectedTexture.width()/2.0, -selectedTexture.height()/2.0)).matrix();
				if(gradientManipulation)
				{
					Eigen::Affine2d transform;
					transform.matrix() = selectedTransformation;
					rotationMatrix2d = transform.rotation();
					rotationMatrix3d.row(0) << rotationMatrix2d(0, 0), rotationMatrix2d(0, 1), 0;
					rotationMatrix3d.row(1) << rotationMatrix2d(1, 0), rotationMatrix2d(1, 1), 0;
					rotationMatrix3d.row(2) << 0, 0, 1;
				}
			}
			Eigen::Matrix3d invTransformation = selectedTransformation.inverse();
			if(textureIsPeriodic)
			{
				unsigned randomShiftX = rand()%selectedTexture.width(); //each tile is pre-translated for later.
				unsigned randomShiftY = rand()%selectedTexture.height();
				shiftedTile.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
				{
					Eigen::Vector3d v;
					v[0]=x;
					v[1]=y;
					v[2]=1;
					Eigen::Vector3d vt = invTransformation * v;
					//interpolation?
					/*pix = bilinear_interpolation(selectedTexture,
												 vt[0] + randomShiftX + selectedTexture.width()*4,
												 vt[1] + randomShiftY + selectedTexture.height()*4,
												 textureIsPeriodic);*/
					pix = selectedTexture.pixelAbsolute(int(vt[0] + randomShiftX + selectedTexture.width()*4)%selectedTexture.width(),
														int(vt[1] + randomShiftY + selectedTexture.height()*4)%selectedTexture.height());
					if(gradientManipulation)
					{
						Eigen::Vector3d rotatedGradientPixel = rotationMatrix3d*(Eigen::Vector3d(pix[0], pix[1], 1) - Eigen::Vector3d(0.5, 0.5, 0))
																+ Eigen::Vector3d(0.5, 0.5, 0);
						pix[0] = rotatedGradientPixel[0];
						pix[1] = rotatedGradientPixel[1];
					}
				});
				IO::save01_in_u8(shiftedTile, std::string("/home/nlutz/shiftedTile_") + std::to_string(randomTransformationIndex) + ".png");
				I correctionOfTheShiftedTile;//while the tile is periodic, the shifted tile is not always.
				//We solve this by shifting the origin of the previously translated tile to the origin of the patch.
				correctionOfTheShiftedTile.initItk(shiftedTile.width(), shiftedTile.height());
				correctionOfTheShiftedTile.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
				{
					pix = shiftedTile.pixelAbsolute((x-patchOrigin[0] + shiftedTile.width())%shiftedTile.width(),
													(y-patchOrigin[1] + shiftedTile.height())%shiftedTile.height());
				});
				Content<I> c(correctionOfTheShiftedTile, patch);
				c.setTranslationTag(Eigen::Vector2i(0, 0));
				patch.addContent(c);
			}
			else
			{
				//first compute the transformed texture explicitely
				I transformedTexture;
				using ImageGrayb = ImageGray8;
				ImageGrayb inBoundsMap;
				transformedTexture.initItk( sqrt(	selectedTexture.width()*selectedTexture.width()+
													selectedTexture.height()*selectedTexture.height()),
											sqrt(	selectedTexture.width()*selectedTexture.width()+
													selectedTexture.height()*selectedTexture.height()));
				inBoundsMap.initItk(transformedTexture.width(), transformedTexture.height());
				transformedTexture.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
				{
					Eigen::Vector3d v;
					v[0]=x-(transformedTexture.width()-selectedTexture.width())/2.0;
					v[1]=y-(transformedTexture.height()-selectedTexture.height())/2.0;
					v[2]=1;
					Eigen::Vector3d vt = invTransformation * v;
					if(vt[0]<0 || vt[0]>=selectedTexture.width() || vt[1]<0 || vt[1]>=selectedTexture.height())
					{
						pix = this->ms_zero;
						inBoundsMap.pixelAbsolute(x, y)=0;
					}
					else
					{
						pix = selectedTexture.pixelAbsolute(vt[0], vt[1]);
						inBoundsMap.pixelAbsolute(x, y)=1;
					}
					if(gradientManipulation)
					{
						Eigen::Vector3d rotatedGradientPixel = rotationMatrix3d*(Eigen::Vector3d(pix[0], pix[1], 1) - Eigen::Vector3d(0.5, 0.5, 0))
																+ Eigen::Vector3d(0.5, 0.5, 0);
						pix[0] = rotatedGradientPixel[0];
						pix[1] = rotatedGradientPixel[1];
					}
				});
				IO::save01_in_u8(transformedTexture, std::string("/home/nlutz/shiftedTexture_") + std::to_string(randomTransformationIndex) + ".png");
				//then choose a portion such that no pixel goes out of bound.
				//We assume it is possible because patches are supposed to be much smaller than the images in the pool
				Stamping::SamplerPoissonGrid sampler;
				sampler.setNbPoints(64);
				std::vector<Eigen::Vector2f> samples = sampler.generate();
				bool foundSuitableRegion=false;
				PixelPos contentOrigin;
				for(unsigned s=0; !foundSuitableRegion && s<samples.size(); ++s)
				{
					contentOrigin[0]=samples[s][0]*transformedTexture.width() - patchBBWidth/2;
					contentOrigin[1]=samples[s][1]*transformedTexture.height() - patchBBHeight/2;
					bool outOfBounds=false;
					for(int x=contentOrigin[0]; x<contentOrigin[0]+patchBBWidth && !outOfBounds; ++x)
					{
						if(x<0 || x>=transformedTexture.width())
							outOfBounds=true;
						else
							for(int y=contentOrigin[1]; y<contentOrigin[1]+patchBBHeight && !outOfBounds; ++y)
							{
								if(y<0 || y>=transformedTexture.height() ||
										(patch.alphaMap().pixelAbsolute(x-contentOrigin[0], y-contentOrigin[1])!=0
										 && inBoundsMap.pixelAbsolute(x, y)==0)) //meaning that pixel cannot be used
								outOfBounds=true;
							}

					}
					if(!outOfBounds)
						foundSuitableRegion=true;
				}
				if(!foundSuitableRegion)	//in case you get a crash here, the sampler may be broken (unlikely),
					exit(EXIT_FAILURE);		//the patch of this content too large (likely),
				else						//or the image you're trying to get a content from too small (very likely).
				{ //finally, build a tile texture such that the content you want is centered on the patch.
					I shiftedTile;
					shiftedTile.initItk(this->m_tile.width(), this->m_tile.height());
					unsigned countX=0;
					for(unsigned x=patchOrigin[0];
						x!=(patchOrigin[0]+patchBBWidth)%this->m_tile.width();
						x=(x+1)%this->m_tile.width(), ++countX)
					{
						unsigned countY=0;
						for(unsigned y=patchOrigin[1];
							y!=(patchOrigin[1]+patchBBHeight)%this->m_tile.height();
							y=(y+1)%this->m_tile.height(), ++countY)
						{
							PixelPos contentIndex = contentOrigin;
							contentIndex[0]=contentOrigin[0]+countX;
							contentIndex[1]=contentOrigin[1]+countY;
							shiftedTile.pixelAbsolute(x, y) = transformedTexture.pixelAbsolute(contentIndex);
						}
					}
					Content<I> c(shiftedTile, patch);
					c.setTranslationTag(Eigen::Vector2i(0, 0));
					patch.addContent(c);
				}
			}
		}
	}
}

}//namespace

}//namespace


#endif
