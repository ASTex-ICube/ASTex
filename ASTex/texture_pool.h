#include <Eigen/Eigen>
#include "mipmap.h"
#include "image_gray.h"
#include "image_common.h"
#include "easy_io.h"

namespace ASTex
{

template<typename I, typename MASK_TYPE>
class TexturePool
{
public:

	using PixelPosType = itk::Index<2>;

	using ImageType = I;
	using PixelType = typename I::PixelType;
	using DataType = typename I::DataType;
	using FuncType = std::function<void(ImageType &)>;

	using MatrixType = Eigen::Matrix3d;

	using ImageMaskType = ASTex::ImageCommon<ASTex::ImageGrayBase<MASK_TYPE>, false>;
	using PixelMaskType = typename ImageMaskType::PixelType;

	typedef struct
	{
		ImageType texture;
		ImageMaskType boundaries;
		bool periodicity;
	} TransformedTexture;

	TexturePool();

	void setTexturesAreGradient(bool b);
	bool texturesAreGradient() const;

	size_t size() const;
	const TransformedTexture &operator[](size_t i) const;

	void addTexture(const ImageType &texture, bool periodicity=false);
	void addTransformationMatrix(const MatrixType &transformation);
	void addFunction(const FuncType &function);

	void generate();
	bool isGenerated() const;

private:

	std::vector<std::pair<ImageType, bool>>		m_texturePool;
	std::vector<MatrixType>						m_transformationPool;
	std::vector<TransformedTexture>				m_pool;
	std::vector<FuncType>						m_functionPool;
	bool										m_generated;
	bool										m_texturesAreGradient;
};

template<typename I, typename MASK_TYPE>
TexturePool<I, MASK_TYPE>::TexturePool() :
	m_texturePool(),
	m_transformationPool(),
	m_pool(),
	m_functionPool(),
	m_generated(false),
	m_texturesAreGradient(false)
{}

template<typename I, typename MASK_TYPE>
void TexturePool<I, MASK_TYPE>::setTexturesAreGradient(bool b)
{
	m_texturesAreGradient = b;
}

template<typename I, typename MASK_TYPE>
bool TexturePool<I, MASK_TYPE>::texturesAreGradient() const
{
	return m_texturesAreGradient;
}

template<typename I, typename MASK_TYPE>
size_t TexturePool<I, MASK_TYPE>::size() const
{
	return m_pool.size();
}

template<typename I, typename MASK_TYPE>
const typename TexturePool<I, MASK_TYPE>::TransformedTexture &TexturePool<I, MASK_TYPE>::operator[](size_t i) const
{
	assert(i<m_pool.size());
	static bool warned = false;
	if(!warned && !m_generated)
	{
		std::cerr << "Warning: TexturePool::operator[] called at least once in a non-generated state." << std::endl;
		warned = true;
	}
	return m_pool[i];
}

template<typename I, typename MASK_TYPE>
void TexturePool<I, MASK_TYPE>::addTexture(const ImageType &texture, bool periodicity)
{
	m_generated = false;
	m_texturePool.push_back(std::make_pair(texture, periodicity));
}

template<typename I, typename MASK_TYPE>
void TexturePool<I, MASK_TYPE>::addTransformationMatrix(const MatrixType &transformation)
{
	m_generated = false;
	m_transformationPool.push_back(transformation);
}

template<typename I, typename MASK_TYPE>
void TexturePool<I, MASK_TYPE>::addFunction(const FuncType &function)
{
	m_generated = false;
	m_functionPool.push_back(function);
}

template<typename I, typename MASK_TYPE>
void TexturePool<I, MASK_TYPE>::generate()
{
	if(m_texturePool.size()==0)
	{
		std::cerr << "Warning: TexturePool::generate() called with 0 textures in the pool." << std::endl;
	}
	ImageType selectedTexture, shiftedTile;
	Eigen::Matrix3d selectedTransformation;
	FuncType selectedFunction;
	bool periodicity = true;
	TransformedTexture expTexture;

	m_pool.clear();
	for(unsigned textureIndex=0; textureIndex<m_texturePool.size(); ++textureIndex)
	{
		for(unsigned transformIndex=0; transformIndex<=m_transformationPool.size(); ++transformIndex)
		{
			for(unsigned functionIndex=0; functionIndex<=m_functionPool.size(); ++functionIndex)
			{

				selectedTexture.initItk(m_texturePool[textureIndex].first.width(),
										m_texturePool[textureIndex].first.height());
				selectedTexture.copy_pixels(m_texturePool[textureIndex].first);
				if(functionIndex!=m_functionPool.size())
				{
					m_functionPool[functionIndex](selectedTexture);
				}
				periodicity = m_texturePool[textureIndex].second;
				Eigen::Matrix2d rotationMatrix2d;
				Eigen::Matrix3d rotationMatrix3d;
				if(transformIndex==m_transformationPool.size())
				{
					selectedTransformation=Eigen::Matrix3d::Identity();
					rotationMatrix3d=Eigen::Matrix3d::Identity();
				}
				else
				{
					selectedTransformation=m_transformationPool[transformIndex];
					selectedTransformation=Eigen::Affine2d(Eigen::Translation2d(selectedTexture.width()/2.0, selectedTexture.height()/2.0)).matrix()
							* selectedTransformation
							* Eigen::Affine2d(Eigen::Translation2d(-selectedTexture.width()/2.0, -selectedTexture.height()/2.0)).matrix();
					if(m_texturesAreGradient)
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
				TransformedTexture transformedTexture;
				if(periodicity)
				{
					Eigen::Vector3d v;
					if(transformIndex==m_transformationPool.size())
					{
						shiftedTile.initItk(selectedTexture.width(), selectedTexture.height(), true);
						shiftedTile.copy_pixels(selectedTexture);
						transformedTexture.periodicity = true;
					}
					else
					{
						//we need a larger tile because rotations and scaling can cause the tile to go out of bounds
						//if the tile is non periodic, and lose information if the tile is periodic.
						//I put a factor of 1.6 for the size of the tile because
						//it can encompasse a pi/2-rotated image with a small upscaling on top of it.
						shiftedTile.initItk(1.6*selectedTexture.width(), 1.6*selectedTexture.height(), true);
						shiftedTile.for_all_pixels([&] (PixelType &pix, int x, int y)
						{
							v[0]=x;
							v[1]=y;
							v[2]=1;
							Eigen::Vector3d vt = invTransformation * v;
							pix = selectedTexture.pixelAbsolute(int(vt[0] + selectedTexture.width()*8)%selectedTexture.width(),
																int(vt[1] + selectedTexture.height()*8)%selectedTexture.height());
							if(m_texturesAreGradient)
							{
								/*
								Eigen::Vector3d rotatedGradientPixel = rotationMatrix3d*(Eigen::Vector3d(pix[0], pix[1], 1)
																						- Eigen::Vector3d(0.5, 0.5, 0))
																		+ Eigen::Vector3d(0.5, 0.5, 0); //in case the gradient is between 0 and 1
								//DataType *gradientData = reinterpret_cast<DataType*>(&pix);
								//gradientData[0] = std::min(std::max(rotatedGradientPixel[0], 0.0), 1.0);
								//gradientData[1] = std::min(std::max(rotatedGradientPixel[1], 0.0), 1.0);
								*/
								DataType *gradientData = reinterpret_cast<DataType*>(&pix);
								Eigen::Vector3d rotatedGradientPixel = rotationMatrix3d*(Eigen::Vector3d(gradientData[0], gradientData[1], 1));
								gradientData[0] = rotatedGradientPixel[0];
								gradientData[1] = rotatedGradientPixel[1];
							}
						});
						transformedTexture.periodicity = false; //periodicity is lost for some transformations.
						//process could be mildly optimized by discriminating transformations that keep periodicity
						//and transformations that do not.
					}
					//in any case, the entire domain is valid when the texture is periodic, no matter the transformation.
					transformedTexture.boundaries.initItk(shiftedTile.width(), shiftedTile.height());
					transformedTexture.boundaries.for_all_pixels([&] (PixelMaskType &pix)
					{
						pix = 1.0;
					});
				}
				else
				{
					shiftedTile.initItk(1.6*selectedTexture.width(), 1.6*selectedTexture.height(), true);
					PixelPosType posMin, posMax;
					ImageType croppedTile;
					ImageMaskType boundaries;
					ImageMaskType croppedBoundaries;
					posMin[0] = shiftedTile.width()/2;
					posMax[0] = shiftedTile.width()/2;
					posMin[1] = shiftedTile.height()/2;
					posMax[1] = shiftedTile.height()/2;
					//First we compute the transformed texture explicitely.
					boundaries.initItk(shiftedTile.width(), shiftedTile.height(), true);
					shiftedTile.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
					{
						Eigen::Vector3d v, vt;
						v[0]=x-(shiftedTile.width()-selectedTexture.width())/2.0;
						v[1]=y-(shiftedTile.height()-selectedTexture.height())/2.0;
						v[2]=1;
						vt = invTransformation * v;
						if(vt[0]<0 || vt[0]>=selectedTexture.width() || vt[1]<0 || vt[1]>=selectedTexture.height())
						{
							boundaries.pixelAbsolute(x, y)=0;
						}
						else
						{
							pix = selectedTexture.pixelAbsolute(vt[0], vt[1]);
							boundaries.pixelAbsolute(x, y)=1;
							posMin[0] = std::min(PixelPosType::IndexValueType(x), posMin[0]);
							posMin[1] = std::min(PixelPosType::IndexValueType(y), posMin[1]);
							posMax[0] = std::max(PixelPosType::IndexValueType(x), posMax[0]);
							posMax[1] = std::max(PixelPosType::IndexValueType(y), posMax[1]);
							if(m_texturesAreGradient)
							{
//								Eigen::Vector3d rotatedGradientPixel = rotationMatrix3d*(Eigen::Vector3d(pix[0], pix[1], 1)
//																		- Eigen::Vector3d(0.5, 0.5, 0))
//																		+ Eigen::Vector3d(0.5, 0.5, 0);
								DataType *gradientData = reinterpret_cast<DataType*>(&pix);
								Eigen::Vector3d rotatedGradientPixel = rotationMatrix3d*(Eigen::Vector3d(gradientData[0], gradientData[1], 1));
								gradientData[0] = rotatedGradientPixel[0];
								gradientData[1] = rotatedGradientPixel[1];
							}
						}
					});
					//Then, since transformations may leave holes, we attempt to cut most of those holes.
					//We also compute a boundary map that gives us the position of those holes (as 0).
					crop_image(shiftedTile, croppedTile, posMin[0], posMax[0], posMin[1], posMax[1]);
					shiftedTile = croppedTile;
					crop_image(boundaries, croppedBoundaries, posMin[0], posMax[0], posMin[1], posMax[1]);
					boundaries = croppedBoundaries;

					transformedTexture.periodicity = false;
					transformedTexture.boundaries = boundaries;
				}
				transformedTexture.texture = shiftedTile;
				m_pool.push_back(transformedTexture);
			}
		}
	}
	m_generated = true;
}

template<typename I, typename MASK_TYPE>
bool TexturePool<I, MASK_TYPE>::isGenerated() const
{
	return m_generated;
}

}
