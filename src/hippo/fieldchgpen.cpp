#include "ff/amoeba/induce.h"
#include "ff/elec.h"
#include "ff/energy.h"
#include "ff/hippo/inducechgpen.h"
#include "ff/nblist.h"

namespace tinker {
void dfieldEwaldRecipSelf_acc(real (*field)[3]);

void dfield_chgpen_nonewald(real (*field)[3]);
void dfield_chgpen_ewald(real (*field)[3]);
void dfield_chgpen_ewald_recip_self(real (*field)[3]);
void dfield_chgpen_ewald_real(real (*field)[3]);

void dfield_chgpen_nonewald_acc(real (*field)[3]);
void dfield_chgpen_ewald_real_acc(real (*field)[3]);
void dfield_chgpen_nonewald_cu(real (*field)[3]);
void dfield_chgpen_ewald_real_cu(real (*field)[3]);

void ufield_chgpen_nonewald(const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald(const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald_recip_self(const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald_real(const real (*uind)[3], real (*field)[3]);

void ufield_chgpen_nonewald_acc(const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald_recip_self_acc(const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald_real_acc(const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_nonewald_cu(const real (*uind)[3], real (*field)[3]);
void ufield_chgpen_ewald_real_cu(const real (*uind)[3], real (*field)[3]);
}

namespace tinker {
void dfieldChgpen(real (*field)[3])
{
   if (useEwald())
      dfield_chgpen_ewald(field);
   else
      dfield_chgpen_nonewald(field);
}

void dfield_chgpen_nonewald(real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      dfield_chgpen_nonewald_cu(field);
   else
#endif
      dfield_chgpen_nonewald_acc(field);
}

void dfield_chgpen_ewald(real (*field)[3])
{
   dfield_chgpen_ewald_recip_self(field);
   dfield_chgpen_ewald_real(field);
}

void dfield_chgpen_ewald_recip_self(real (*field)[3])
{
   dfieldEwaldRecipSelf_acc(field);
}

void dfield_chgpen_ewald_real(real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      dfield_chgpen_ewald_real_cu(field);
   else
#endif
      dfield_chgpen_ewald_real_acc(field);
}

void ufieldChgpen(const real (*uind)[3], real (*field)[3])
{
   if (useEwald())
      ufield_chgpen_ewald(uind, field);
   else
      ufield_chgpen_nonewald(uind, field);
}

void ufield_chgpen_nonewald(const real (*uind)[3], real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      ufield_chgpen_nonewald_cu(uind, field);
   else
#endif
      ufield_chgpen_nonewald_acc(uind, field);
}

void ufield_chgpen_ewald(const real (*uind)[3], real (*field)[3])
{
   ufield_chgpen_ewald_recip_self(uind, field);
   ufield_chgpen_ewald_real(uind, field);
}

void ufield_chgpen_ewald_recip_self(const real (*uind)[3], real (*field)[3])
{
   ufield_chgpen_ewald_recip_self_acc(uind, field);
}

void ufield_chgpen_ewald_real(const real (*uind)[3], real (*field)[3])
{
#if TINKER_CUDART
   if (mlistVersion() & Nbl::SPATIAL)
      ufield_chgpen_ewald_real_cu(uind, field);
   else
#endif
      ufield_chgpen_ewald_real_acc(uind, field);
}
}