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



#ifndef __PERMUTOHEDRAL_LATTICE_H__
#define __PERMUTOHEDRAL_LATTICE_H__




#include <unordered_map>
#include "Eigen/Eigen"

template <typename ARRAY>
void assign_array( ARRAY& ar, typename ARRAY::value_type x)
{
	for ( auto& v : ar )
		v = x;
}

template < typename TScalar, int D, int DValue >
class PermutohedralLattice
{
    using Point         = Eigen::Matrix< TScalar, D  , 1 >;
    using Embedding     = Eigen::Matrix< TScalar, D+1, 1 >;
    using Vertex        = Eigen::Matrix< int32_t, D+1, 1 >;
    using Permutation   = std::array< int, D+1 >;
    using BaryCoords    = std::array< TScalar, D+/*1*/2 >;
    using TValue        = Eigen::Matrix< TScalar, DValue, 1 >;


    class Hasher
    {
    public:
        inline size_t operator()( const Vertex& key ) const
        {
            size_t hash = key[0] * 2531011;
            for( int i=1; i<=D; ++i )
                hash = (hash + key[i]) * 2531011;
            return hash;
        }
    };

    using HashMap       = std::unordered_map< Vertex, uint32_t, Hasher >;

    struct Simplex
    {
        Vertex      vertex0;
        Permutation permutation;
    };

    struct Recording
    {
        BaryCoords  barycentric;
        uint32_t    vertexIds[D+1];
    };

    Vertex                  canonical_[D+1];
    Vertex                  latticeNeighborhood_[D+1];
    TScalar                 embeddingScaleFactor_[D];
    TValue                  zeroValue_;
    HashMap                 latticeHash_;
    std::vector<TValue>     latticeValues_;
    std::vector<Recording>  vertexRecording_;

    inline void buildCanonicalSimplex()
    {
        for( int k=0; k<=D; ++k )
        {
            for( int i=0; i<(D+1)-k; ++i )
                canonical_[k][i] = k;
            for( int i=(D+1)-k; i<=D; ++i )
                canonical_[k][i] = k - (D+1);
        }
    }

    inline void getHyperplaneCoordinates( const Point &position, Embedding &embedding ) const
    {
		TScalar sm( 0 );

        for( int i=D; i>0; i-- )
        {
			TScalar cf = position[i-1] * embeddingScaleFactor_[i-1];
			embedding[i] = sm - i*cf;
			sm += cf;
		}

        embedding[0] = sm;
    }

    inline void getEnclosingSimplex( const Embedding &x, Simplex &simplex ) const
    {
        Embedding differential;
        int32_t h = 0;

        for( int i=0; i<=D; ++i )
        {
            TScalar sx = x[i] / (D+1);

            TScalar lower = std::floor(sx) * (D+1);
            TScalar upper = std::ceil (sx) * (D+1);

            if( upper - x[i] < x[i] - lower )
                simplex.vertex0[i] = (int32_t) upper;
            else
                simplex.vertex0[i] = (int32_t) lower;

            differential[i] = x[i] - simplex.vertex0[i];

            h += simplex.vertex0[i];
        }

        h /= D+1;


//        simplex.permutation.assign( 0 );
		assign_array(simplex.permutation,0);

        for( int i=0; i<D; ++i )
            for( int j=i+1; j<=D; ++j )
                if( differential[i] < differential[j] )
                    simplex.permutation[i] ++;
                else
                    simplex.permutation[j] ++;


        // See lemma 2.9 page 4.
        if( h < 0 )
        {
            for( int i=0; i<=D; ++i )
                if( simplex.permutation[i] < -h )
                {
                    simplex.vertex0[i] += (D+1);
                    simplex.permutation[i] += (D+1) + h;
                }
                else
                    simplex.permutation[i] += h;
        }
        else if( h > 0 )
        {
            for( int i=0; i<=D; ++i )
                if( simplex.permutation[i] >= (D+1) - h )
                {
                    simplex.vertex0[i] -= (D+1);
                    simplex.permutation[i] += h - (D+1);
                }
                else
                    simplex.permutation[i] += h;
        }
    }

    void getSimplexVertex( const Simplex &simplex, int n, Vertex &v ) const
    {
        for( int i=0; i<=D; ++i )
            v[i] = simplex.vertex0[i] + canonical_[n][ simplex.permutation[i] ];
    }

    void getBarycentricCoordinates( const Embedding &x, const Simplex &simplex, BaryCoords &barycentric ) const
    {
        const TScalar scale = 1 / TScalar(D+1);

//        barycentric.assign( TScalar(0) );
		assign_array(barycentric, TScalar(0));

        for( int i=0; i<=D; ++i )
        {
            int pi = simplex.permutation[i];
            TScalar yi = (x[i] - simplex.vertex0[i]) * scale;
            barycentric[D  -pi] += yi;
            barycentric[D+1-pi] -= yi;
        }

        barycentric[0] += 1 + barycentric[D+1];
    }

    template < typename TCreateRecordFunc,
               typename TRecordingFunc   ,
               typename TSplattingFunc   >
    inline void splatCore( const Point &point,
                           const TValue &value,
                           const TCreateRecordFunc &createRecord,
                           const TRecordingFunc &doRecording,
                           const TSplattingFunc &doSplatting )
    {
        Embedding y;
        getHyperplaneCoordinates( point, y );

        Simplex simplex;
        getEnclosingSimplex( y, simplex );

        Recording recording;
        getBarycentricCoordinates( y, simplex, recording.barycentric );

        for( int i=0; i<=D; ++i )
        {
            Vertex vertex;
            getSimplexVertex( simplex, i, vertex );

            auto entry = latticeHash_.find( vertex );
            if( entry == latticeHash_.end() )
            {
                doRecording( recording.vertexIds[i], (uint32_t) latticeValues_.size() );
                latticeHash_[vertex] = (uint32_t) latticeValues_.size();
                latticeValues_.push_back( zeroValue_ );
                doSplatting( latticeValues_.back(), value * recording.barycentric[i] );
            }
            else
            {
                doRecording( recording.vertexIds[i], entry->second );
                doSplatting( latticeValues_[entry->second], value * recording.barycentric[i] );
            }
        }

        createRecord( recording );
    }

public:
    inline PermutohedralLattice()
    {
        TScalar invStdDev = std::sqrt(TScalar(2.0/3.0)) * (D+1);
        for( int i=0; i<D; ++i )
            embeddingScaleFactor_[i] = invStdDev / std::sqrt( (TScalar)(i+1)*(i+2) );

        for( int n=0; n<=D; ++n )
        {
            for( int i=0; i<=D; ++i )
                latticeNeighborhood_[n][i] = -1;
            latticeNeighborhood_[n][n] = D;
        }

        for( int i=0; i<DValue; ++i )
            zeroValue_[i] = TScalar(0);

        buildCanonicalSimplex();
    }

    inline void clear()
    {
        latticeHash_.clear();
        latticeValues_.clear();
        vertexRecording_.clear();
    }

    inline void reserveRecordings( unsigned int nPoints )
    {
        vertexRecording_.reserve( nPoints );
    }

    inline void splatAndDeclare( const Point &point, const TValue &value )
    {
        splatCore( point, value,
            [&]( const Recording &rec ){ vertexRecording_.push_back( rec ); },
            [&]( uint32_t &dst, uint32_t src ) { dst = src; },
            [&]( TValue &dst, const TValue &src ) { dst += src; }
        );
    }

    inline void splatOnly( const Point &point, const TValue &value )
    {
        splatCore( point, value,
			[&]( const Recording &/*rec*/ ){},
			[&]( uint32_t &/*dst*/, uint32_t /*src*/ ) {},
            [&]( TValue &dst, const TValue &src ) { dst += src; }
        );
    }

    inline void declareOnly( const Point &point )
    {
        splatCore( point, zeroValue_,
            [&]( const Recording &rec ){ vertexRecording_.push_back( rec ); },
            [&]( uint32_t &dst, uint32_t src ) { dst = src; },
			[&]( TValue &/*dst*/, const TValue &/*src*/ ) {}
        );
    }

    //inline void buildLatticeOnly( const Point &point );
    //inline TValue slice( const Point &point );

    inline void blur()
    {
#if 0
size_t nonZero = 0;
size_t maxBucket = 0;
for( size_t i=0; i<latticeHash_.bucket_count(); ++i )
    if( latticeHash_.bucket_size(i) )
    {
        maxBucket = std::max( maxBucket, latticeHash_.bucket_size(i) );
        ++ nonZero;
    }
std::cout << latticeValues_.size() << std::endl;
std::cout << latticeHash_.bucket_count() << "  " << maxBucket << "  " << nonZero << std::endl;
#endif
        std::vector<TValue> valuesTmp( latticeHash_.size() );
        TValue *valuesNew = &valuesTmp[0];
        TValue *valuesOld = &latticeValues_[0];

        // Apply at each lattice vertex a 1D blurring kernel of radius 1 along each axis.

        for( int n=0; n<=D; ++n )
        {
            for( auto entry = latticeHash_.begin(); entry != latticeHash_.end(); ++entry )
            {
                const Vertex &vertex = entry->first;
                uint32_t vertexId = entry->second;

                TValue &blurredValue = valuesNew[vertexId];
                blurredValue = /*TScalar(0.5) * */valuesOld[vertexId] + valuesOld[vertexId];

                auto neighborVertex = latticeHash_.find( vertex - latticeNeighborhood_[n] );
                if( neighborVertex != latticeHash_.end() )
                    blurredValue += /*TScalar(0.25) * */valuesOld[ neighborVertex->second ];

                neighborVertex = latticeHash_.find( vertex + latticeNeighborhood_[n] );
                if( neighborVertex != latticeHash_.end() )
                    blurredValue += /*TScalar(0.25) * */valuesOld[ neighborVertex->second ];
            }

            std::swap( valuesOld, valuesNew );
        }

        if( valuesNew == &latticeValues_[0] )
            latticeValues_ = valuesTmp;
    }

    inline TValue slice( unsigned int pointId )
    {
        auto &recording = vertexRecording_[pointId];

        TValue slicedValue = zeroValue_;
        for( int i=0; i<=D; ++i )
            slicedValue += latticeValues_[ recording.vertexIds[i] ] * recording.barycentric[i];

        return slicedValue;
    }
};




#endif // __PERMUTOHEDRAL_LATTICE_H__
