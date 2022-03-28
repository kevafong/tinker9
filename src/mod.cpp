#define TINKER_EXTERN_DEFINITION_FILE 1

#include "mod/accasync.h"
#include "mod/disp.h"
#include "mod/elecamoeba.h"
#include "mod/elechippo.h"
#include "mod/elecpchg.h"
#include "mod/energy.h"
#include "mod/evalence.h"
#include "mod/md.h"
#include "mod/mutant.h"
#include "mod/repel.h"
#include "mod/vdw.h"

#include "ff/atom.h"
#include "ff/box.h"
#include "ff/elec.h"
#include "ff/molecule.h"
#include "ff/nblist.h"
#include "ff/pme.h"
#include "ff/rattle.h"
#include "ff/spatial.h"

#include "tool/cudalib.h"
#include "tool/gpucard.h"
#include "tool/rcman.h"

#include "platform.h"
