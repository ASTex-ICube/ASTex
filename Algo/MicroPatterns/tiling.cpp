
#include <Algo/MicroPatterns/tiling.h>

namespace ASTex
{

void Tiling_n_Blendinng::TriangleGrid(const Eigen::Vector2d& p_uv, Eigen::Vector3d& Bi, Eigen::Vector2i& vertex1, Eigen::Vector2i& vertex2, Eigen::Vector2i& vertex3)
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


Eigen::Vector2d Tiling_n_Blendinng::hash(const Eigen::Vector2i& p)
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

Tiling_n_Blendinng::Tiling_n_Blendinng()
{}

void Tiling_n_Blendinng::operator()(const MicroPatternsInput& noises, double level, const Eigen::Vector2d& uv, Eigen::Vector2d& mean, Eigen::Vector2d& variance)
{
    Eigen::Vector3d B;
    Eigen::Vector2i  vertex1, vertex2, vertex3;
    TriangleGrid(uv, B,	vertex1, vertex2, vertex3);

    // Assign random offset to each triangle vertex
    Eigen::Vector2d uv1 = uv + hash(vertex1);
    Eigen::Vector2d uv2 = uv + hash(vertex2);
    Eigen::Vector2d uv3 = uv + hash(vertex3);

    Eigen::Vector4d n1 = noises.fetch(uv1, level);
    Eigen::Vector4d n2 = noises.fetch(uv2, level);
    Eigen::Vector4d n3 = noises.fetch(uv3, level);

    Eigen::Vector2d nu = noises.fetch_average();

    B.normalize();
    Eigen::Matrix<double, 2, 3> M;
    M << n1[0], n2[0], n3[0],
            n1[1], n2[1], n3[1];
    mean = M * B +nu;
    mean[0] = std::max(0.0,std::min(1.0,mean[0]));
    mean[1] = std::max(0.0,std::min(1.0,mean[1]));

    B = B.cwiseProduct(B);
    Eigen::Matrix<double, 2, 3> S; 
    S << n1[2], n2[2], n3[2],
            n1[3], n2[3], n3[3];
    variance = S * B;
}

} // namespace
