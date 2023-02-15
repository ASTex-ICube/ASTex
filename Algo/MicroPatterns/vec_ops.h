#ifndef _MP_VEC_OPS_
#define _MP_VEC_OPS_

#include <Eigen/Dense>

namespace ASTex
{

template <typename VEC, typename FUNC>
inline VEC applyFuncCompo(const VEC& u , const FUNC& f)
{ 
	VEC v;
	for (int i = 0; i < v.size(); ++i)
		v[i] = f(u[i]);
	return v;
}


template <typename VEC>
inline VEC vec_floor(const VEC& u)
{
	return applyFuncCompo(u, [](double d) { return std::floor(d); });
}

template <typename VEC>
inline VEC vec_fract(const VEC& u)
{
	return applyFuncCompo(u, [](double d) { return d - std::floor(d); });
}


template <typename VEC>
inline VEC vec_clamp01(const VEC& u)
{
	return applyFuncCompo(u, [](double d) { return std::max(0.0,std::min(1.0,d)); });
}

} // namespace

#endif
