/*
  Copyright (C) 2022 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/material_model/simple.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/postprocess/interface.h>
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
    os << ev.l2_norm() << ' ';

    for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
      {
        ev = 0.0;
        this->get_artificial_viscosity_composition(ev, c);
        std::cout << "EV composition " << c << ": " << std::endl;
        ev.print(std::cout, 10);
        os << ev.l2_norm() << ' ';
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
