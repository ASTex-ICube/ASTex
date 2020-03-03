#ifndef CSN_TEXTURE_H_
#define CSN_TEXTURE_H_

#include <Eigen/Eigen>

namespace CSN
{

template<typename I>
class CSN_Texture
{
public:

	using ImageType = I;
	using PixelType = typename ImageType::PixelType;
	using DataType = typename ImageType::DataType;

private:

	ImageType m_texture;
	Eigen::Vector2f m_cycles[2];
};

}

#endif
