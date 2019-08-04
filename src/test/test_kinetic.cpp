#include "util_files.h"
#include "util_md.h"
#include "util_test.h"
#include "util_test_rt.h"

using namespace TINKER_NAMESPACE;

static const char* verlet_intg = "integrator  verlet\n";
static int usage_ = calc::xyz | calc::vel | calc::mass | calc::vmask | calc::md;

TEST_CASE("Kinetic-ArBox", "[ff][kinetic][arbox]") {
  const char* k = "test_arbox.key";
  const char* d = "test_arbox.dyn";
  const char* x = "test_arbox.xyz";
  const char* p = "amoeba09.prm";

  std::string k0 = arbox_key;
  k0 += verlet_intg;
  TestFile fke(k, k0);

  TestFile fd(d, arbox_dyn2);
  TestFile fx(x, arbox_xyz);
  TestFile fp(p, amoeba09_prm);

  const char* argv[] = {"dummy", x};
  int argc = 2;
  test_begin_with_args(argc, argv);
  test_mdinit(0, 0);

  use_data = usage_;
  initialize();

  real temp;
  kinetic(temp);

  const double ref_eksum = 100446.40376;
  const double ref_temp = 156008.001336;
  const double eps_e = 0.0001;

  REQUIRE(eksum == Approx(ref_eksum).margin(eps_e));
  REQUIRE(temp == Approx(ref_temp).margin(eps_e));

  finish();
  test_end();
}
