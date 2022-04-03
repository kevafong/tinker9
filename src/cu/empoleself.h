#pragma once
#include "add.h"
#include "ff/amoeba/mpole.h"
#include "tool/energybuffer.h"

namespace tinker {
template <bool do_a>
__global__
void empoleSelf_cu(CountBuffer restrict nem, EnergyBuffer restrict em,
   const real (*restrict rpole)[10], int n, real f, real aewald)
{
   real aewald_sq_2 = 2 * aewald * aewald;
   real fterm = -f * aewald * 0.5f * (real)(M_2_SQRTPI);

   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n; i += blockDim.x * gridDim.x) {
      real ci = rpole[i][MPL_PME_0];
      real dix = rpole[i][MPL_PME_X];
      real diy = rpole[i][MPL_PME_Y];
      real diz = rpole[i][MPL_PME_Z];
      real qixx = rpole[i][MPL_PME_XX];
      real qixy = rpole[i][MPL_PME_XY];
      real qixz = rpole[i][MPL_PME_XZ];
      real qiyy = rpole[i][MPL_PME_YY];
      real qiyz = rpole[i][MPL_PME_YZ];
      real qizz = rpole[i][MPL_PME_ZZ];

      real cii = ci * ci;
      real dii = dix * dix + diy * diy + diz * diz;
      real qii =
         2 * (qixy * qixy + qixz * qixz + qiyz * qiyz) + qixx * qixx + qiyy * qiyy + qizz * qizz;

      int offset = threadIdx.x + blockIdx.x * blockDim.x;
      real e = fterm * (cii + aewald_sq_2 * (dii / 3 + 2 * aewald_sq_2 * qii * (real)0.2));
      atomic_add(e, em, offset);
      if CONSTEXPR (do_a)
         atomic_add(1, nem, offset);
   }
}
}
