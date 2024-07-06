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



#ifndef __ASTEX__NEIGHBORHOOD_H__
#define __ASTEX__NEIGHBORHOOD_H__


#include <ASTex/image_common.h>
#include <array>

namespace ASTex
{

/** Type representing the position of a pixel in an image.
 */
using PixelAbsolute = itk::Index<2>;


/** Type representing the position of a pixel in an image, in absolute and relative coordinates with respect to another pixel.
 *  In the case of images with periodic boundaries, a neighborhood query provides an absolute index which is always inside the
 *  image frame and a relative index which corresponds to contiguous coordinates wrt. the query pixel.
 */
struct PixelRelative
{
    itk::Index<2> absolute;
    itk::Index<2> relative;

    inline PixelRelative() {}
    inline PixelRelative( const itk::Index<2> &_absolute ) : absolute(_absolute), relative(_absolute) {}
    inline PixelRelative( const itk::Index<2> &_absolute, const itk::Index<2> &_relative ) : absolute(_absolute), relative(_relative) {}

    inline operator itk::Index<2>() const { return absolute; }
};


/** Neighborhood operator that considers images without periodic boudaries.
 *  When 4 or 8-neighborhoods are queried, pixels crossing boundaries are simply rejected.
 */
class NeighborhoodClamped
{
    int  w, h;

public:
    using PixelType      = PixelAbsolute;

    using NeighborsList  = std::vector<PixelType>;

    /** Creates a neighborhood operator on the basis of the provided image size.
     *  \param img          Image from which neighborhood queries have to be done.
     */
    template <class T>
    inline NeighborhoodClamped( const T& img ) :
        w( img.width () ),
        h( img.height() )
    {
    }

    /** Creates a neighborhood operator on the basis of the provided image size.
     *  \param size         Size of the image from which neighborhood queries have to be done.
     */
    inline NeighborhoodClamped( const ASTex::Size& size ) :
        w( size[0] ),
        h( size[1] )
    {
    }

    /** Creates a neighborhood operator on the basis of the provided image size.
     *  \param width        Width of the image from which neighborhood queries have to be done.
     *  \param height       Height of the image from which neighborhood queries have to be done.
     */
    inline NeighborhoodClamped( itk::SizeValueType width, itk::SizeValueType height ) :
        w( width  ),
        h( height )
    {
    }

    inline static bool isPeriodic() { return false; }

    /** Recovers the 4-neighbors of the given pixel.
     *  \param[in]  pixel       Pixel for which neighborhood is queried.
     *  \param[out] neighbors   List where the recovered neighbors are stored.
     */
    inline void get4( const PixelType &pixel, NeighborsList &neighbors ) const
    {
        neighbors.clear();

        if( pixel[0] > 0 )
            neighbors.push_back( PixelType({pixel[0]-1,pixel[1]}) );

        if( pixel[0] < w-1 )
            neighbors.push_back( PixelType({pixel[0]+1,pixel[1]}) );

        if( pixel[1] > 0 )
            neighbors.push_back( PixelType({pixel[0],pixel[1]-1}) );

        if( pixel[1] < h-1 )
            neighbors.push_back( PixelType({pixel[0],pixel[1]+1}) );
    }

    /** Recovers the 8-neighbors of the given pixel.
     *  \param[in]  pixel       Pixel for which neighborhood is queried.
     *  \param[out] neighbors   List where the recovered neighbors are stored.
     */
    inline void get8( const PixelType &pixel, NeighborsList &neighbors ) const
    {
        neighbors.clear();

        if( pixel[0] > 0 )
        {
            neighbors.push_back( PixelType({pixel[0]-1,pixel[1]}) );
            if( pixel[1] > 0 )
                neighbors.push_back( PixelType({pixel[0]-1,pixel[1]-1}) );
            if( pixel[1] < h-1 )
                neighbors.push_back( PixelType({pixel[0]-1,pixel[1]+1}) );
        }

        if( pixel[0] < w-1 )
        {
            neighbors.push_back( PixelType({pixel[0]+1,pixel[1]}) );
            if( pixel[1] > 0 )
                neighbors.push_back( PixelType({pixel[0]+1,pixel[1]-1}) );
            if( pixel[1] < h-1 )
                neighbors.push_back( PixelType({pixel[0]+1,pixel[1]+1}) );
        }

        if( pixel[1] > 0 )
            neighbors.push_back( PixelType({pixel[0],pixel[1]-1}) );

        if( pixel[1] < h-1 )
            neighbors.push_back( PixelType({pixel[0],pixel[1]+1}) );
    }

    /** Recovers the neighbor of the given pixel located at (pixel.x+dx, pixel.y+dy).
     *  \param[in]  pixel       Pixel for which neighborhood is queried.
     *  \param[in]  dx          Relative X offset of the pixel wrt. the query pixel.
     *  \param[in]  dy          Relative Y offset of the pixel wrt. the query pixel.
     *  \param[out] neighbor    Where the recovered neighbor is stored.
     */
    inline bool get( const PixelType &pixel, int dx, int dy, PixelType &neighbor ) const
    {
        neighbor = { pixel[0]+dx, pixel[1]+dy };
        return neighbor[0] >= 0  &&  neighbor[0] < w  &&
               neighbor[1] >= 0  &&  neighbor[1] < h;
    }
};


/** Neighborhood operator that considers images with periodic boudaries.
    *  When 4 or 8-neighborhoods are queried, coordinates of pixels that cross the image boundaries are cyclically wrapped to the opposite side.
    */
class NeighborhoodPeriodic
{
public:
    using PixelType      = PixelAbsolute;

    using NeighborsList  = std::vector<PixelType>;
    using NeighborsList4 = std::array<PixelType,4>;
    using NeighborsList8 = std::array<PixelType,8>;

private:
    inline itk::IndexValueType clamp( itk::IndexValueType v, itk::IndexValueType vmax ) const
    {
        return (v + vmax) % vmax;
    }

    int  w, h;

    inline void get4( const PixelType &pixel, PixelType *neighbors ) const
    {
        neighbors[0] = { clamp(pixel[0]-1,w), pixel[1] };
        neighbors[1] = { clamp(pixel[0]+1,w), pixel[1] };
        neighbors[2] = { pixel[0], clamp(pixel[1]-1,h) };
        neighbors[3] = { pixel[0], clamp(pixel[1]+1,h) };
    }

    inline void get8( const PixelType &pixel, PixelType *neighbors ) const
    {
        get4( pixel, neighbors );
        neighbors[4] = { neighbors[0][0], neighbors[2][1] };
        neighbors[5] = { neighbors[1][0], neighbors[2][1] };
        neighbors[6] = { neighbors[0][0], neighbors[3][1] };
        neighbors[7] = { neighbors[1][0], neighbors[3][1] };
    }

public:
    /** Creates a neighborhood operator on the basis of the provided image size.
     *  \param img          Image from which neighborhood queries have to be done.
     */
    template <class T>
    inline NeighborhoodPeriodic( const T& img ) :
        w( img.width () ),
        h( img.height() )
    {
    }

    /** Creates a neighborhood operator on the basis of the provided image size.
     *  \param size         Size of the image from which neighborhood queries have to be done.
     */
    inline NeighborhoodPeriodic( const ASTex::Size& size ) :
        w( size[0] ),
        h( size[1] )
    {
    }

    /** Creates a neighborhood operator on the basis of the provided image size.
     *  \param width        Width of the image from which neighborhood queries have to be done.
     *  \param height       Height of the image from which neighborhood queries have to be done.
     */
    inline NeighborhoodPeriodic( itk::SizeValueType width, itk::SizeValueType height ) :
        w( width  ),
        h( height )
    {
    }

    inline static bool isPeriodic() { return true; }

#if 0
    /** Recovers the 4-neighbors of the given pixel.
     *  \param[in]  pixel       Pixel for which neighborhood is queried.
     *  \param[out] neighbors   List where the recovered neighbors are stored.
     */
    inline void get4( const PixelType &pixel, NeighborsList4 &neighbors ) const
    {
        get4( pixel, &neighbors[0] );
    }
#endif

    /** Recovers the 4-neighbors of the given pixel.
     *  \param[in]  pixel       Pixel for which neighborhood is queried.
     *  \param[out] neighbors   List where the recovered neighbors are stored.
     */
    inline void get4( const PixelType &pixel, NeighborsList  &neighbors ) const
    {
        neighbors.resize( 4 );
        get4( pixel, &neighbors[0] );
    }

#if 0
    /** Recovers the 8-neighbors of the given pixel.
     *  \param[in]  pixel       Pixel for which neighborhood is queried.
     *  \param[out] neighbors   List where the recovered neighbors are stored.
     */
    inline void get8( const PixelType &pixel, NeighborsList8 &neighbors ) const
    {
        get8( pixel, &neighbors[0] );
    }
#endif

    /** Recovers the 8-neighbors of the given pixel.
     *  \param[in]  pixel       Pixel for which neighborhood is queried.
     *  \param[out] neighbors   List where the recovered neighbors are stored.
     */
    inline void get8( const PixelType &pixel, NeighborsList  &neighbors ) const
    {
        neighbors.resize( 8 );
        get8( pixel, &neighbors[0] );
    }

    /** Recovers the neighbor of the given pixel located at (pixel.x+dx, pixel.y+dy).
     *  \param[in]  pixel       Pixel for which neighborhood is queried.
     *  \param[in]  dx          Relative X offset of the pixel wrt. the query pixel.
     *  \param[in]  dy          Relative Y offset of the pixel wrt. the query pixel.
     *  \param[out] neighbor    Where the recovered neighbor is stored.
     */
    inline bool get( const PixelType &pixel, int dx, int dy, PixelType &neighbor ) const
    {
        neighbor = { clamp(pixel[0]+dx,w), clamp(pixel[1]+dy,h) };
        return true;
    }
};


/** Neighborhood operator that considers images with periodic boudaries.
    *  When 4 or 8-neighborhoods are queried, coordinates of pixels that cross the image boundaries are cyclically wrapped to the opposite side.
    */
class NeighborhoodPeriodicRelative
{
public:
    using PixelType      = PixelRelative;

    using NeighborsList  = std::vector<PixelType>;
    using NeighborsList4 = std::array<PixelType,4>;
    using NeighborsList8 = std::array<PixelType,8>;

private:
    inline itk::IndexValueType clamp( itk::IndexValueType v, itk::IndexValueType vmax ) const
    {
        return (v + vmax) % vmax;
    }

    int  w, h;

    inline void get4( const PixelType &pixel, PixelType *neighbors ) const
    {
        neighbors[0].relative = { pixel.relative[0]-1, pixel.relative[1]   };
        neighbors[1].relative = { pixel.relative[0]+1, pixel.relative[1]   };
        neighbors[2].relative = { pixel.relative[0]  , pixel.relative[1]-1 };
        neighbors[3].relative = { pixel.relative[0]  , pixel.relative[1]+1 };

        for( int i=0; i<4; ++i )
        {
            neighbors[i].absolute[0] = clamp( neighbors[i].relative[0], w );
            neighbors[i].absolute[1] = clamp( neighbors[i].relative[1], h );
        }
    }
    inline void get8( const PixelType &pixel, PixelType *neighbors ) const
    {
        get4( pixel, neighbors );

        neighbors[4].absolute = { neighbors[0].absolute[0], neighbors[2].absolute[1] };
        neighbors[5].absolute = { neighbors[1].absolute[0], neighbors[2].absolute[1] };
        neighbors[6].absolute = { neighbors[0].absolute[0], neighbors[3].absolute[1] };
        neighbors[7].absolute = { neighbors[1].absolute[0], neighbors[3].absolute[1] };

        neighbors[4].relative = { neighbors[0].relative[0], neighbors[2].relative[1] };
        neighbors[5].relative = { neighbors[1].relative[0], neighbors[2].relative[1] };
        neighbors[6].relative = { neighbors[0].relative[0], neighbors[3].relative[1] };
        neighbors[7].relative = { neighbors[1].relative[0], neighbors[3].relative[1] };
    }

public:
    /** Creates a neighborhood operator on the basis of the provided image size.
     *  \param img          Image from which neighborhood queries have to be done.
     */
    template <class T>
    inline NeighborhoodPeriodicRelative( const T& img ) :
        w( img.width () ),
        h( img.height() )
    {
    }

    /** Creates a neighborhood operator on the basis of the provided image size.
        *  \param size         Size of the image from which neighborhood queries have to be done.
        */
    inline NeighborhoodPeriodicRelative( const ASTex::Size& size ) :
        w( size[0] ),
        h( size[1] )
    {
    }

    /** Creates a neighborhood operator on the basis of the provided image size.
        *  \param width        Width of the image from which neighborhood queries have to be done.
        *  \param height       Height of the image from which neighborhood queries have to be done.
        */
    inline NeighborhoodPeriodicRelative( itk::SizeValueType width, itk::SizeValueType height ) :
        w( width  ),
        h( height )
    {
    }

    inline static bool isPeriodic() { return true; }

    /** Recovers the 4-neighbors of the given pixel.
        *  \param[in]  pixel       Pixel for which neighborhood is queried.
        *  \param[out] neighbors   List where the recovered neighbors are stored.
        */
    inline void get4( const PixelType &pixel, NeighborsList4 &neighbors ) const
    {
        get4( pixel, &neighbors[0] );
    }

    /** Recovers the 4-neighbors of the given pixel.
        *  \param[in]  pixel       Pixel for which neighborhood is queried.
        *  \param[out] neighbors   List where the recovered neighbors are stored.
        */
    inline void get4( const PixelType &pixel, NeighborsList &neighbors ) const
    {
        neighbors.resize( 4 );
        get4( pixel, &neighbors[0] );
    }

    /** Recovers the 8-neighbors of the given pixel.
        *  \param[in]  pixel       Pixel for which neighborhood is queried.
        *  \param[out] neighbors   List where the recovered neighbors are stored.
        */
    inline void get8( const PixelType &pixel, NeighborsList8 &neighbors ) const
    {
        get8( pixel, &neighbors[0] );
    }

    /** Recovers the 8-neighbors of the given pixel.
        *  \param[in]  pixel       Pixel for which neighborhood is queried.
        *  \param[out] neighbors   List where the recovered neighbors are stored.
        */
    inline void get8( const PixelType &pixel, NeighborsList &neighbors ) const
    {
        neighbors.resize( 8 );
        get8( pixel, &neighbors[0] );
    }

    /** Recovers the neighbor of the given pixel located at (pixel.x+dx, pixel.y+dy).
        *  \param[in]  pixel       Pixel for which neighborhood is queried.
        *  \param[in]  dx          Relative X offset of the pixel wrt. the query pixel.
        *  \param[in]  dy          Relative Y offset of the pixel wrt. the query pixel.
        *  \param[out] neighbor    Where the recovered neighbor is stored.
        */
    inline bool get( const PixelType &pixel, int dx, int dy, PixelType &neighbor ) const
    {
        neighbor.relative = { pixel.relative[0]+dx, pixel.relative[1]+dy };
        neighbor.absolute[0] = clamp( neighbor.relative[0], w );
        neighbor.absolute[1] = clamp( neighbor.relative[1], h );
        return true;
    }
};


} // namespace ASTex




#endif // __ASTEX__NEIGHBORHOOD_H__
