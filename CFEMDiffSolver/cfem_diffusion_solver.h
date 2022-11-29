#ifndef CHI_CFEM_DIFFUSION_SOLVER_H
#define CHI_CFEM_DIFFUSION_SOLVER_H

#include "ChiPhysics/SolverBase/chi_solver.h"
#include "ChiMath/PETScUtils/petsc_utils.h"

#include "cfem_diffusion_bndry.h"
#include "ChiTimer/chi_timer.h"

// forword declaration 
namespace chi_mesh
{
class MeshContinuum; 
typedef std::shared_ptr<MeshContinuum> MeshContinuumPtr;
};
namespace chi_math
{
class SpatialDiscretization; 
typedef std::shared_ptr<SpatialDiscretization> SDMPtr ;
};

namespace cfem_diffusion
{
/** My CFEM diffusion solver 
 * 
*/
class Solver : public chi_physics::Solver
{
private:
  chi_objects::ChiTimer t_assembly;
  chi_objects::ChiTimer t_solve;

  double time_assembly=0.0, time_solve=0.0;
  bool verbose_info=true;

public:
  chi_mesh::MeshContinuumPtr grid_ptr=nullptr;
  chi_math::SDMPtr sdm_ptr =nullptr;
  size_t num_local_dofs = 0;
  size_t num_globl_dofs = 0;

  Vec            x = nullptr;            // approx solution
  Vec            b = nullptr;            // RHS
  Mat            A = nullptr;            // linear system matrix

  typedef std::pair<BoundaryType,std::vector<double>> BoundaryInfo;
  typedef std::map<uint, BoundaryInfo> BoundaryPreferences;
  BoundaryPreferences      boundary_preferences;
  std::vector<Boundary*>   boundaries;

  chi_math::UnknownManager                 unknown_manager;

  explicit Solver(const std::string& in_solver_name);
  // void Initialize() override;
  void Initialize() override {Initialize(true);}
  void Initialize(bool verbose);
  void Execute() override;
};

}; // namespace cfem_diffusion


#endif // CHI_CFEM_DIFFUSION_SOLVER_H

