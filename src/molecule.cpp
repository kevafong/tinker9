#include "molecule.h"
#include "dev_array.h"
#include "md.h"
#include <tinker/detail/molcul.hh>


TINKER_NAMESPACE_BEGIN
Molecule::~Molecule()
{
   device_array::deallocate(imol, kmol, molecule, molmass);
   imol = nullptr;
   kmol = nullptr;
   molecule = nullptr;
   molmass = nullptr;
}


void molecule_data(rc_op op)
{
   if (op & rc_dealloc) {
      molecule.~Molecule();
   }


   if (op & rc_alloc) {
      auto& st = molecule;
      device_array::allocate(n, &st.imol, &st.kmol, &st.molecule, &st.molmass);
   }


   if (op & rc_init) {
      auto& st = molecule;

      std::vector<int> buf(2 * n);
      st.nmol = molcul::nmol;
      for (int i = 0; i < st.nmol; ++i) {
         int j = 2 * i;
         buf[j] = molcul::imol[j] - 1;
         buf[j + 1] = molcul::imol[j + 1];
      }
      device_array::copyin(st.nmol, st.imol, buf.data());
      for (int i = 0; i < n; ++i) {
         buf[i] = molcul::kmol[i] - 1;
      }
      device_array::copyin(n, st.kmol, buf.data());
      for (int i = 0; i < n; ++i) {
         buf[i] = molcul::molcule[i] - 1;
      }
      device_array::copyin(n, st.molecule, buf.data());
      st.totmass = molcul::totmass;
      device_array::copyin(st.nmol, st.molmass, molcul::molmass);
   }
}
TINKER_NAMESPACE_END