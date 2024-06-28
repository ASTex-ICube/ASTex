#ifndef _MP_TILING_
#define _MP_TILING_

#include <Eigen/Dense>

#include <Algo/MicroPatterns/mp_inputs.h>
#include <Algo/MicroPatterns/vec_ops.h>

namespace ASTex
{

class Tiling_n_Blendinng
{
	void TriangleGrid(const Eigen::Vector2d& p_uv, Eigen::Vector3d& Bi, Eigen::Vector2i& vertex1, Eigen::Vector2i& vertex2, Eigen::Vector2i& vertex3);
	//original hash version
	Eigen::Vector2d hash(const Eigen::Vector2i& p);

public:
	Tiling_n_Blendinng();
	void operator()(const MicroPatternsInput& noises, double level, const Eigen::Vector2d& uv, Eigen::Vector2d& mean, Eigen::Vector2d& variance);
};

} // namespace

#endif
