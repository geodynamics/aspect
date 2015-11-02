#include <aspect/material_model/simple.h>
#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>


namespace aspect
{
  using namespace dealii;

  template <int dim>
  class EVPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      /**
       * Generate graphical output from the current solution.
       */
      virtual
      std::pair<std::string,std::string>
      execute (TableHandler &statistics);
  };

  template <int dim>
  std::pair<std::string,std::string>
  EVPostprocessor<dim>::execute (TableHandler &statistics)
  {
    std::ostringstream os;
    Vector<float> ev(this->get_triangulation().n_active_cells());
    this->get_artificial_viscosity(ev);
    std::cout << "EV temperature: " << std::endl;
    ev.print(std::cout);
    os << ev.l2_norm() << " ";

    for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
      {
        ev = 0.0;
        this->get_artificial_viscosity_composition(ev, c);
        std::cout << "EV composition " << c << ": " << std::endl;
        ev.print(std::cout, 10);
        os << ev.l2_norm() << " ";
      }
    return std::make_pair("EV norms:", os.str());
  }

}



// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_POSTPROCESSOR(EVPostprocessor,
                                "EVPostprocessor",
                                "")
}
