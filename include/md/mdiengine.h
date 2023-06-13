#pragma once
#include "tool/error.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>
#include "mdi.h"

namespace tinker {

class MDIEngine
{
protected:
   static std::FILE* outputFile;
   static int target_node;
   static MDI_Comm mdi_comm;
   static bool exit_flag;
   
public:
   static void initialize(std::FILE* o);
   static void mdi_commands();
   static void run_mdi(const char* node);

   // temporary function for MDI testing purposes
   template <class F, class... Ts>
   static void mdiprint(const F& f, Ts&&... args){
       print(outputFile, f, args...);
   };
};

}
