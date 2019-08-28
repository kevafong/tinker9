#include "acc_add.h"
#include "e_polar.h"
#include "md.h"

TINKER_NAMESPACE_BEGIN
void epolar0_dotprod(const real (*gpu_uind)[3], const real (*gpu_udirp)[3]) {
  const real f = -0.5 * electric / dielec;

  auto* ep = ep_handle.e()->buffer();
  auto bufsize = ep_handle.buffer_size();

  #pragma acc parallel num_gangs(bufsize)\
              deviceptr(ep,gpu_uind,gpu_udirp,polarity_inv)
  #pragma acc loop gang independent
  for (int i = 0; i < n; ++i) {
    int offset = i & (bufsize - 1);
    real e = polarity_inv[i] *
        (gpu_uind[i][0] * gpu_udirp[i][0] + gpu_uind[i][1] * gpu_udirp[i][1] +
         gpu_uind[i][2] * gpu_udirp[i][2]);
    atomic_add_value(f * e, ep, offset);
  }
}
TINKER_NAMESPACE_END
