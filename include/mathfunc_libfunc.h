#pragma once
#include "macro.h"
#include <cmath>


#define INT_ABS abs


#if TINKER_DOUBLE_PRECISION
#   define REAL_SQRT sqrt
#   define REAL_EXP exp
#   define REAL_FLOOR floor
#   define REAL_ABS fabs
#   define REAL_POW pow
#   define REAL_RECIP(x) (1 / static_cast<double>(x))
#   define REAL_RSQRT(x) (1 / sqrt(x))
#   define REAL_COS cos
#   define REAL_SIN sin
#   define REAL_ACOS acos
#   define REAL_ASIN asin
#   define REAL_ERF erf
// #   define REAL_ERFC erfc
#   define REAL_ERFC(x) (1 - erf(x))
#   define REAL_MIN fmin
#   define REAL_MAX fmax
#   define REAL_SIGN copysign
#endif


#if TINKER_SINGLE_PRECISION
#   define REAL_SQRT sqrtf
#   define REAL_EXP expf
#   define REAL_FLOOR floorf
#   define REAL_ABS fabsf
#   define REAL_POW powf
#   define REAL_RECIP(x) (1 / static_cast<float>(x))
#   define REAL_RSQRT(x) (1 / sqrtf(x))
#   define REAL_COS cosf
#   define REAL_SIN sinf
#   define REAL_ACOS acosf
#   define REAL_ASIN asinf
#   define REAL_ERF erff
// #   define REAL_ERFC erfcf
#   define REAL_ERFC(x) (1 - erff(x))
#   define REAL_MIN fminf
#   define REAL_MAX fmaxf
#   define REAL_SIGN copysignf
#endif


#define REAL_SQ(x) ((x) * (x))
#define REAL_CUBE(x) ((x) * (x) * (x))


#ifdef _OPENACC
#   include <accelmath.h>


#pragma acc routine(abs) seq


#pragma acc routine(sqrt) seq
#pragma acc routine(exp) seq
#pragma acc routine(floor) seq
#pragma acc routine(fabs) seq
#pragma acc routine(pow) seq
#pragma acc routine(cos) seq
#pragma acc routine(sin) seq
#pragma acc routine(acos) seq
#pragma acc routine(asin) seq
#pragma acc routine(erf) seq
#pragma acc routine(erfc) seq
#pragma acc routine(fmin) seq
#pragma acc routine(fmax) seq
#pragma acc routine(copysign) seq


#pragma acc routine(sqrtf) seq
#pragma acc routine(expf) seq
#pragma acc routine(floorf) seq
#pragma acc routine(fabsf) seq
#pragma acc routine(powf) seq
#pragma acc routine(cosf) seq
#pragma acc routine(sinf) seq
#pragma acc routine(acosf) seq
#pragma acc routine(asinf) seq
#pragma acc routine(erff) seq
#pragma acc routine(erfcf) seq
#pragma acc routine(fminf) seq
#pragma acc routine(fmaxf) seq
#pragma acc routine(copysignf) seq
#endif


TINKER_NAMESPACE_BEGIN
#ifndef __CUDACC__
inline float2 make_float2(float x, float y)
{
   float2 f{.x = x, .y = y};
   return f;
}

inline float3 make_float3(float x, float y, float z)
{
   float3 f{.x = x, .y = y, .z = z};
   return f;
}

inline float4 make_float4(float x, float y, float z, float w)
{
   float4 f{.x = x, .y = y, .z = z, .w = w};
   return f;
}

inline double2 make_double2(double x, double y)
{
   double2 f{.x = x, .y = y};
   return f;
}

inline double3 make_double3(double x, double y, double z)
{
   double3 f{.x = x, .y = y, .z = z};
   return f;
}

inline double4 make_double4(double x, double y, double z, double w)
{
   double4 f{.x = x, .y = y, .z = z, .w = w};
   return f;
}
#endif


#if TINKER_SINGLE_PRECISION
#   define make_real2(x, y) make_float2((x), (y))
#   define make_real3(x, y, z) make_float3((x), (y), (z))
#   define make_real4(x, y, z, w) make_float4((x), (y), (z), (w))
#endif


#if TINKER_DOUBLE_PRECISION
#   define make_real2(x, y) make_double2((x), (y))
#   define make_real3(x, y, z) make_double3((x), (y), (z))
#   define make_real4(x, y, z, w) make_double4((x), (y), (z), (w))
#endif
TINKER_NAMESPACE_END
