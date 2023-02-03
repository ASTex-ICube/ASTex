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

namespace ASTex
{

/**
 Simple Tiling n blending algo 
 Need Gaussian imagea as input
 */
template<typename IMG>
class Tiling_n_Blending
{
    using PIXT = typename IMG::PixelType;
    using EPIXT = typename IMG::DoublePixelEigen;

    const IMG& img_input_;

    EPIXT img_average_;

protected:

    inline void set_zero(Eigen::Vector2d& v)
    {
        v.setZero();
    }

    inline void set_zero(Eigen::Vector3d& v)
    {
        v.setZero();
    }

    inline void set_zero(Eigen::Vector4d& v)
    {
        v.setZero();
    }

    inline void set_zero(double& v)
    {
        v = 0.0;
    }

	inline void clamp_channel(double& c)
	{
		if (std::is_floating_point<typename IMG::DataType>::value)
			c = std::max(0.0,std::min(1.0,c));
		else
			c = std::max(0.0,std::min(255.0,c));
	}

    inline void clamp(Eigen::Vector2d& v)
    {
        clamp_channel(v[0]);
		clamp_channel(v[1]);
    }

    inline void clamp(Eigen::Vector3d& v)
    {
        clamp_channel(v[0]);
		clamp_channel(v[1]);
        clamp_channel(v[2]);
    }

    inline void clamp(Eigen::Vector4d& v)
    {
        clamp_channel(v[0]);
		clamp_channel(v[1]);
        clamp_channel(v[2]);
		clamp_channel(v[3]);
    }

    inline void clamp(double& v)
    {
        clamp_channel(v);
    }


	void TriangleGrid(const Eigen::Vector2d& p_uv, Eigen::Vector3d& Bi, Eigen::Vector2i& vertex1, Eigen::Vector2i& vertex2, Eigen::Vector2i& vertex3)
	{
		const Eigen::Vector2d uv = p_uv * 2.0 * std::sqrt(3.0);

		Eigen::Matrix2d gridToSkewedGrid;
		gridToSkewedGrid << 1.0, -0.57735027,
							0.0, 01.15470054;

		Eigen::Vector2d skewedCoord = gridToSkewedGrid * uv;
		Eigen::Vector2d baseId{ std::floor(skewedCoord[0]), std::floor(skewedCoord[1]) };
		Eigen::Vector3d temp{ skewedCoord[0] - baseId[0], skewedCoord[1] - baseId[1], 0.0 };
		temp[2] = 1.0 - temp[0] - temp[1];

		if (temp[2] > 0.0)
		{
			Bi = Eigen::Vector3d(temp[2], temp[1], temp[0]);
			Eigen::Vector2i ibaseId = baseId.cast<int>();
			vertex1 = ibaseId;
			vertex2 = ibaseId + Eigen::Vector2i(0, 1);
			vertex3 = ibaseId + Eigen::Vector2i(1, 0);
		}
		else
		{
			Bi = Eigen::Vector3d(-temp[2], 1.0 - temp[1], 1.0 - temp[0]);
			Eigen::Vector2i ibaseId = baseId.cast<int>();
			vertex1 = ibaseId + Eigen::Vector2i(1, 1);
			vertex2 = ibaseId + Eigen::Vector2i(1, 0);
			vertex3 = ibaseId + Eigen::Vector2i(0, 1);
		}
	}


	//original hash version
	Eigen::Vector2d hash(const Eigen::Vector2i& p)
	{
		Eigen::Matrix2d hashMat;
		hashMat << 127.1, 269.5,
				311.7, 183.3;

		Eigen::Vector2d q = hashMat * p.cast<double>();
		q[0] = std::sin(q[0]);
		q[1] = std::sin(q[1]);
		q *= 43758.5453;
		return Eigen::Vector2d(q[0] - std::floor(q[0]), q[1] - std::floor(q[1]));

	}

public:
	Tiling_n_Blending(const IMG& in):
    img_input_(in)
	{
		// compute average img value
        EPIXT sum;
        set_zero(sum);
        img_input_.for_all_pixels([&] ( const PIXT& P)
        {
            EPIXT lv = eigenPixel<double>(P);
            sum += lv;
        });

        img_average_ = sum / double(img_input_.width()*img_input_.height());
        
    }

    EPIXT fetch(const Eigen::Vector2d& uv)
    {
		// take fract part of uv mult by image size
		Eigen::Vector2d uvd = Eigen::Vector2d{-0.5 + (uv[0] - std::floor(uv[0]))*img_input_.width(),
		 -0.5 + (uv[1] - std::floor(uv[1]))*img_input_.height()};
		Eigen::Vector2d uvfl {std::floor(uvd[0]), std::floor(uvd[1])};
		//  a = coef for linear interpolation
		Eigen::Vector2d a = uvd - uvfl;
		// c = integer texel coord
		Eigen::Vector2i c = uvfl.cast<int>();

		auto acces_repeat = [&] (int xp, int yp)
		{
            int xx = xp < 0 ? img_input_.width() - 1 : ( xp >= img_input_.width() ? 0 : xp );
			int yy = yp < 0 ? img_input_.height() - 1 : ( yp >= img_input_.height() ? 0 : yp );
			return eigenPixel<double>(img_input_.pixelAbsolute(xx, yy));
		};


        return acces_repeat(c[0],c[1]);

		EPIXT V1 = (1.0 - a[0]) * acces_repeat(c[0],c[1])	+ a[0] * acces_repeat(c[0]+1,c[1]);
		EPIXT V2 = (1.0 - a[0]) * acces_repeat(c[0],c[1]+1) + a[0] * acces_repeat(c[0]+1,c[1]+1);
		EPIXT V = (1.0 - a[1]) * V1 + a[1] * V2;

        return V;
    }

	    
    PIXT tile_pixel(const Eigen::Vector2d& uv )
	{
		Eigen::Vector3d B;
		Eigen::Vector2i  vertex1, vertex2, vertex3;
		TriangleGrid(uv, B,	vertex1, vertex2, vertex3);

		// Assign random offset to each triangle vertex
    	Eigen::Vector2d uv1 = uv + hash(vertex1);
		Eigen::Vector2d uv2 = uv + hash(vertex2);
		Eigen::Vector2d uv3 = uv + hash(vertex3);


		EPIXT t1 = fetch(uv1) - img_average_;
        EPIXT t2 = fetch(uv2) - img_average_;
        EPIXT t3 = fetch(uv3) - img_average_;

        auto W = B.normalized();

        EPIXT P = W[0] * t1 + W[1] * t2 + W[2] * t3 + img_average_ ;

		clamp(P);

        return IMG::itkPixel(P);
	}


	void tile_img(IMG& img_out)
	{
     	img_out.parallel_for_all_pixels([&] (typename IMG::PixelType& P, int x, int y) 
		{
			Eigen::Vector2d uv{ double(x) / (img_input_.width()), double(y) / (img_input_.height()) };
			P = tile_pixel(uv);
		});
	}

};

template<typename IMG>
Tiling_n_Blending<IMG> make_Tiling_n_Blending(const IMG& img)
{
    return Tiling_n_Blending<IMG>(img);
}

}
