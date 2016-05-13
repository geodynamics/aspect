#include <aspect/postprocess/velocity_binary.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/lac/block_vector.h>

namespace aspect {
    namespace Postprocess {
        template <int dim>
        VelocityBinary<dim>::VelocityBinary()
        {}

        template <int dim>
       VelocityBinary<dim>::~VelocityBinary()
        {}

       template <int dim>
       void VelocityBinary<dim>::initialize() {
          // Interface::initialize();
       }

       template <int dim>
       std::pair<std::string, std::string> VelocityBinary<dim>::execute(TableHandler &statistics) {
          std::string fileName = "solution-" + Utilities::int_to_string(this->get_timestep_number(), 5) + ".mesh";
          parallel::distributed::SolutionTransfer<dim, LinearAlgebra::BlockVector> sol_trans(this->get_dof_handler());
          sol_trans.prepare_serialization (this->get_solution());
          this->get_triangulation().save(fileName.c_str());
           return std::make_pair("Writing binary output to: ", fileName);
       }
   }
}


namespace aspect
{
    namespace Postprocess
    {
        ASPECT_REGISTER_POSTPROCESSOR(VelocityBinary,
                                      "velocity binary",
                                      "A Postprocessor that output the velocity solution data per timestep.")
    }
}



