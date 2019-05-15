/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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


#include <aspect/mesh_deformation/function.h>

#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    Function<dim>::Function()
      :
      function(dim)
    {}



    template <int dim>
    void
    Function<dim>::update ()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        function.set_time (this->get_time() / year_in_seconds);
      else
        function.set_time (this->get_time());
    }



    template <int dim>
    void
    Function<dim>::deformation_constraints(const DoFHandler<dim> &free_surface_dof_handler,
                                           ConstraintMatrix &mesh_constraints) const
    {
      for (std::set<types::boundary_id>::const_iterator p = this->get_parameters().free_surface_boundary_indicators.begin();
           p != this->get_parameters().free_surface_boundary_indicators.end(); ++p)
        {
          VectorTools::interpolate_boundary_values (free_surface_dof_handler,
                                                    *p,
                                                    function,
                                                    mesh_constraints);
        }
    }



    template <int dim>
    void Function<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void Function<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("Function");
        try
          {
            function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Mesh deformation.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }
        prm.leave_subsection();

      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(Function,
                                           "function",
                                           "A plugin, which prescribes the surface mesh to "
                                           "deform according to an analytically prescribed "
                                           "function. Note that the function prescribes a "
                                           "deformation velocity, i.e. the return value of "
                                           "this plugin is later multiplied by the time step length "
                                           "to compute the displacement increment in this time step. "
                                           "The format of the "
                                           "functions follows the syntax understood by the "
                                           "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
