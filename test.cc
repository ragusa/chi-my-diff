#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiConsole/chi_console.h"
#include "CFEMDiffSolver/lua/lua_utils.h"

#include "ChiMacros/lua_register_macro.h"

// #include "ChiMesh/MeshHandler/chi_meshhandler.h"
// #include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

// #include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
// #include "ChiMath/PETScUtils/petsc_utils.h"
 
// #include "ChiPhysics/FieldFunction2/fieldfunction2.h"

int main(int argc, char* argv[])
{
  chi::Initialize(argc,argv);

  chi::log.Log() << "\n------------------------\nCFEM Diffusion main code!\n------------------------\n";

  // Adding lua functiont hat requires console
  auto& console = chi_objects::ChiConsole::GetInstance();
  auto L = console.consoleState;
  RegisterFunction(chiPrintStatus);
  RegisterFunction(chiCFEMDiffusionSolverCreate);

  chi::RunBatch(argc, argv);

  // Adding rest of the code here
  chi::log.Log() << "Coding Tutorial CFEM MMS";
  chi::Finalize();
  return 0;

}