/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/boundary_composition/initial_composition.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace BoundaryComposition
  {
// ------------------------------ InitialComposition -------------------

    template <int dim>
    double
    InitialComposition<dim>::
    composition (const GeometryModel::Interface<dim> &geometry_model,
                 const unsigned int                   boundary_indicator,
                 const Point<dim>                    &location,
                 const unsigned int                   compositional_field) const
    {
      return this->get_compositional_initial_conditions().initial_composition(location, compositional_field);
    }


    template <int dim>
    double
    InitialComposition<dim>::
    minimal_composition (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return min_composition;
    }



    template <int dim>
    double
    InitialComposition<dim>::
    maximal_composition (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return max_composition;
    }



    template <int dim>
    void
    InitialComposition<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Initial composition");
        {
          prm.declare_entry ("Minimal composition", "0",
                             Patterns::Double (),
                             "Minimal composition. Units: none.");
          prm.declare_entry ("Maximal composition", "1",
                             Patterns::Double (),
                             "Maximal composition. Units: none.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    InitialComposition<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Initial composition");
        {
          min_composition = prm.get_double ("Minimal composition");
          max_composition = prm.get_double ("Maximal composition");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryComposition
  {
    ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(InitialComposition,
                                               "initial composition",
                                               "A model in which the composition at the boundary "
                                               "is chosen to be the same as given in the initial "
                                               "conditions."
                                               "\n\n"
                                               "Because this class simply takes what the initial "
                                               "composition had described, this class can not "
                                               "know certain pieces of information such as the "
                                               "minimal and maximal composition on the boundary. "
                                               "For operations that require this, for example in "
                                               "postprocessing, this boundary composition model "
                                               "must therefore be told what the minimal and "
                                               "maximal values on the boundary are. This is done "
                                               "using parameters set in section ``Boundary composition model/Initial composition''.")
  }
}
