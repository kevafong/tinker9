#pragma once
#include "ff/elec.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup mpole
/// \{
void empoleData(RcOp);
void empole(int vers);
void torque(int vers, grad_prec* dx, grad_prec* dy, grad_prec* dz);
void mpoleInit(int vers);
/// \}
}
