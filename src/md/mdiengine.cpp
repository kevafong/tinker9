#include "md/mdiengine.h"
#include "md/integrator.h"
#include "ff/energy.h"
#include "md/lflpiston.h"
#include "md/misc.h"
#include "md/pq.h"
#include "md/rattle.h"
#include "tool/error.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>
#include <string.h>
#include "mdi.h"

namespace tinker {

std::FILE* MDIEngine::outputFile;
MDI_Comm MDIEngine::mdi_comm;
bool MDIEngine::exit_flag = false;
int MDIEngine::target_node = 0;

void MDIEngine::initialize(std::FILE* o)
{
   outputFile = o;

   int ret;
  int is_initialized = 0;
  ret = MDI_Initialized(&is_initialized);
  if ( is_initialized ) {
    mdi_commands();
    ret = MDI_Accept_communicator(&mdi_comm);
    if ( ret ) {
      TINKER_THROW(format("MDI  --  Error in MDI_Accept_communicator\n"));
    }
    
  }
   // the call to MDI_Init and all other initialization should go here
}

void MDIEngine::mdi_commands()
{
  MDI_Register_node("@DEFAULT");
  MDI_Register_command("@DEFAULT", "<@");
  MDI_Register_command("@DEFAULT", "<NATOMS");
  MDI_Register_command("@DEFAULT", "<MASSES");
  MDI_Register_command("@DEFAULT", "@INIT_MD");

  MDI_Register_node("@INIT_MD");
  MDI_Register_command("@INIT_MD", "<@");
  MDI_Register_command("@INIT_MD", "@");
  MDI_Register_command("@INIT_MD", "@COORDS");
  MDI_Register_command("@INIT_MD", "@FORCES");
  MDI_Register_command("@INIT_MD", "<NATOMS");
  MDI_Register_command("@INIT_MD", "<MASSES");
  MDI_Register_command("@INIT_MD", "<ENERGY");
  MDI_Register_command("@INIT_MD", "<PE");
  MDI_Register_command("@INIT_MD", "<KE");

  MDI_Register_node("@COORDS");
  MDI_Register_command("@COORDS", "<@");
  MDI_Register_command("@COORDS", "@");
  MDI_Register_command("@COORDS", "@COORDS");
  MDI_Register_command("@COORDS", "@FORCES");
  MDI_Register_command("@COORDS", "<NATOMS");
  MDI_Register_command("@COORDS", "<MASSES");
  MDI_Register_command("@COORDS", "<COORDS");
  MDI_Register_command("@COORDS", ">COORDS");
  MDI_Register_command("@COORDS", "<ENERGY");
  MDI_Register_command("@COORDS", "<PE");
  MDI_Register_command("@COORDS", "<KE");
  

  MDI_Register_node("@FORCES");
  MDI_Register_command("@FORCES", "<@");
  MDI_Register_command("@FORCES", "@");
  MDI_Register_command("@FORCES", "@COORDS");
  MDI_Register_command("@FORCES", "@FORCES");
  MDI_Register_command("@FORCES", "<NATOMS");
  MDI_Register_command("@FORCES", "<MASSES");
  MDI_Register_command("@FORCES", "<FORCES");
  MDI_Register_command("@FORCES", ">FORCES");
  MDI_Register_command("@FORCES", "<ENERGY");
  MDI_Register_command("@FORCES", "<PE");
  MDI_Register_command("@FORCES", "<KE");
}

void MDIEngine::run_mdi(const char* node)
{
  int ret;
  int exists;
  int is_initialized = 0;
  double bohrA_conv = 0.529177249; // 1 bohr = 0.529177249 angstroms
  double hartreekcal_conv = 627.503; // 1 hartree = 627.503 kcal/mol
  
  // Check if MDI is initialized
  ret = MDI_Initialized(&is_initialized);
  if ( ret ) {
    TINKER_THROW(format("MDI  --  Error in MDI_Initialized\n"));
  }
  if (not is_initialized) { return; }

  /*Check if a target node has been set
  If current node matches target node, exit*/
  if (target_node)  {
    if ((target_node == 1) && ( strcmp(node, "@COORDS") != 0)) return;
    if ((target_node == 2) && ( strcmp(node, "@FORCES") != 0)) return;
  }

  /* MDI command from the driver */
  char* command = new char[MDI_COMMAND_LENGTH];
 
  /* Main MDI loop */
  while (not exit_flag) {
    /* Receive a command from the driver */
    ret = MDI_Recv_command(command, mdi_comm);
    if ( ret ) {
       TINKER_THROW(format("MDI  --  Error in MDI_Recv_command\n"));
    }



    /*check if command is support at node*/
    // ret = MDI_Check_command_exists(node, command, MDI_COMM_NULL, &exists);
    // if ( !exists ) {
    //    TINKER_THROW(format("MDI  --  Unsupported Command at Node\n"));
    // }

    /* Respond to the received command */
    if ( strcmp(command, "EXIT") == 0 ) {
      exit_flag = true;
    }
    else if ( strcmp(command, "@") == 0 ) {
      //if (strcmp(node, "@DEFAULT") == 0)  TINKER_THROW(format("MDI  --  Unsupported command @DEFAULT\n"));
      break;
    }
    else if ( strcmp(command, "<@") == 0 ) {
      ret = MDI_Send(node, MDI_NAME_LENGTH, MDI_CHAR, mdi_comm);
    }
    else if ( strcmp(command, "<NATOMS") == 0 ) {
      ret = MDI_Send(&n, 1, MDI_INT, mdi_comm);
    }
    else if ( strcmp(command, "<COORDS") == 0 ) {
      double coords[n*3];
      for (int i=0; i<n ; i++)  {
        coords[3*i]= x[i] / bohrA_conv;
        coords[3*i + 1]= y[i] / bohrA_conv;
        coords[3*i + 2]= z[i] / bohrA_conv;
      }


      // mdiprint("coords  in tinker: \n");
      // for (int i=0; i< 3*n; i++)  mdiprint("%f\n", coords[i]);


      // units should be in Bohr not Angstroms, 1.22something
      ret = MDI_Send(coords, n * 3, MDI_DOUBLE, mdi_comm);
    }
    else if ( strcmp(command, ">COORDS") == 0 ) {
      double recv_coords[n*3];
      ret = MDI_Recv(recv_coords, n * 3, MDI_DOUBLE, mdi_comm);
      for (int i=0; i<n ; i++)  {
        mdiprint("%f\n",x[i]);
        x[i]=recv_coords[3*i] * bohrA_conv;
        y[i]=recv_coords[3*i + 1] * bohrA_conv;
        z[i]=recv_coords[3*i + 2] * bohrA_conv;
      }
      // units converted back into angstroms frpm Bohr
    }
    else if ( strcmp(command, "<FORCES") == 0 ) {
      double forces[n*3];
      
      for (int i=0; i<n ; i++)  {
        forces[3*i]= - gx[i] * bohrA_conv / hartreekcal_conv;
        forces[3*i + 1]= - gy[i] * bohrA_conv / hartreekcal_conv;
        forces[3*i + 2]= - gz[i] * bohrA_conv / hartreekcal_conv;
      }
      // units should be in atomic units. (something/angstrom)
      ret = MDI_Send(forces, n * 3, MDI_DOUBLE, mdi_comm);
    }
    else if ( strcmp(command, ">FORCES") == 0 ) {
      double recv_forces[n*3];
      // units should be in atomic units. (something/angstrom)
      // presumable (energy/unit distance)
      //QE is in hartree for energy
      ret = MDI_Recv(recv_forces, n * 3, MDI_DOUBLE, mdi_comm);

      for (int i=0; i<n ; i++)  {
        gx[i] = - recv_forces[3*i] * hartreekcal_conv / bohrA_conv;
        gy[i] = - recv_forces[3*i + 1] * hartreekcal_conv / bohrA_conv;
        gz[i] = - recv_forces[3*i + 2] * hartreekcal_conv / bohrA_conv;
      }
    }
    else if ( strcmp(command, "<MASSES") == 0 ) {
      ret = MDI_Send(mass, n, MDI_DOUBLE, mdi_comm);
    }
    else if ( strcmp(command, ">MASSES") == 0 ) {
      ret = MDI_Recv(mass, n, MDI_DOUBLE, mdi_comm);
    }
    else if ( strcmp(command, "<PE") == 0 ) {
      //double potential = epot / hartreekcal_conv;
      ret = MDI_Send(&esum, 1, MDI_DOUBLE, mdi_comm);
    }
    else if ( strcmp(command, "<KE") == 0 ) {
      double kinetic = eksum / hartreekcal_conv;
      ret = MDI_Send(&eksum, 1, MDI_DOUBLE, mdi_comm);
    }
    else if ( strcmp(command, "<ENERGY") == 0 ) {
      double total_energy= esum + eksum;
      total_energy = total_energy / hartreekcal_conv;
      ret = MDI_Send(&total_energy, 1, MDI_DOUBLE, mdi_comm);
    }
    else if ( strcmp(command, "@INIT_MD") == 0 ) {
      break;
    }
    else if ( strcmp(command, "@COORDS") == 0 ) {
      target_node = 1;
      break;
    }
    else if ( strcmp(command, "@FORCES") == 0 ) {
      target_node = 2;
      break;
    }
    
    else {
      TINKER_THROW(format("MDI  --  Received unsupported command: %s\n",
         command));
    }
  }

  // Free any memory allocations
  delete [] command;

}

}
