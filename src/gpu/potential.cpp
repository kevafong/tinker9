#include "gpu/potential.h"
#include "gpu/acc.h"
#include "gpu/nblist.h"
#include "util/format.print.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void potential_data(int op) {
  e_bond_data(op);
  e_vdw_data(op);
}
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_gradient1() {
  m_tinker_using_namespace;
  const char* title = " Energy Component Breakdown :{:>20s}{:>20s}\n\n";
  const char* fmt = " {:28s}{:>20.4f}{:>17d}        {}\n";

  gpu::async_launches_begin(&gpu::queue_b);
  if (gpu::use_ebond()) {
    tinker_gpu_ebond_harmonic0();
    tinker_gpu_ebond_harmonic4();
    tinker_gpu_ebond_harmonic5();
    tinker_gpu_ebond_harmonic6();
    tinker_gpu_ebond_harmonic1();
  }

  gpu::async_launches_begin(&gpu::queue_nb);

  if (gpu::use_evdw()) {
    tinker_gpu_evdw0();
    tinker_gpu_evdw3();
    tinker_gpu_evdw4();
    tinker_gpu_evdw5();
    tinker_gpu_evdw6();
    tinker_gpu_evdw1();
  }

  tinker_gpu_vlist_update();

  gpu::async_launches_end();

  print(stdout, title, "Kcal/mole", "Interactions");

  print(stdout, fmt, "Bond Stretching", gpu::get_ebond(), gpu::count_ebond(),
        gpu::bndtyp_str);

  print(stdout, fmt, "Van der Waals", gpu::get_evdw(), gpu::count_evdw(),
        gpu::vdwtyp_str);
}
}
