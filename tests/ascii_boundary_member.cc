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

#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>



namespace aspect
{
  /**
   * A boundary velocity plugin that uses an AsciiDataBoundary object as member
   */
  template <int dim>
  class AsciiBoundaryMember : public BoundaryVelocity::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      /**
       * Empty Constructor.
       */
      AsciiBoundaryMember ();

      /**
       * Initialization function. This function is called once at the
       * beginning of the program. Checks preconditions.
       */
      void
      initialize ();

      /**
       * A function that is called at the beginning of each time step. For
       * the current plugin, this function loads the next data files if
       * necessary and outputs a warning if the end of the set of data files
       * is reached.
       */
      void
      update ();

      /**
       * Return the boundary velocity as a function of position. For the
       * current class, this function returns value from the text files.
       */
      Tensor<1,dim>
      boundary_velocity (const types::boundary_id boundary_indicator,
                         const Point<dim> &position) const;

      /**
       * Declare the parameters this class takes through input files.
       */
      static
      void
      declare_parameters (ParameterHandler &prm);

      /**
       * Read the parameters this class declares from the parameter file.
       */
      void
      parse_parameters (ParameterHandler &prm);

      std::unique_ptr<Utilities::AsciiDataBoundary<dim>> member;

      std::set<types::boundary_id> boundary_ids;
  };

  template <int dim>
  AsciiBoundaryMember<dim>::AsciiBoundaryMember ()
  {}


  template <int dim>
  void
  AsciiBoundaryMember<dim>::initialize ()
  {
    unsigned int i=0;
    for (const auto &plugin : this->get_boundary_velocity_manager().get_active_plugins())
      {
        if (plugin.get() == this)
          boundary_ids.insert(this->get_boundary_velocity_manager().get_active_plugin_boundary_indicators()[i]);

        ++i;
      }

    AssertThrow(boundary_ids.empty() == false,
                ExcMessage("Did not find the boundary indicator for the velocity ascii data plugin."));

    member->initialize(boundary_ids,
                       dim);
  }

  template <int dim>
  void
  AsciiBoundaryMember<dim>::update ()
  {
    member->update();
  }


  template <int dim>
  Tensor<1,dim>
  AsciiBoundaryMember<dim>::
  boundary_velocity (const types::boundary_id ,
                     const Point<dim> &position) const
  {
    Tensor<1,dim> velocity;
    for (unsigned int i = 0; i < dim; i++)
      velocity[i] = member->get_data_component(*(boundary_ids.begin()),
                                               position,
                                               i);
    return velocity;
  }


  template <int dim>
  void
  AsciiBoundaryMember<dim>::declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection("Boundary velocity model");
    {
      Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                            "$ASPECT_SOURCE_DIR/data/boundary-velocity/ascii-data/test/",
                                                            "box_2d_%s.%d.txt");
    }
    prm.leave_subsection();
  }


  template <int dim>
  void
  AsciiBoundaryMember<dim>::parse_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection("Boundary velocity model");
    {
      member = std::make_unique<Utilities::AsciiDataBoundary<dim>>();
      member->initialize_simulator(this->get_simulator());

      member->parse_parameters(prm);

      if (this->convert_output_to_years() == true)
        {
          member->scale_factor               /= year_in_seconds;
        }

    }
    prm.leave_subsection();
  }

}



// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(AsciiBoundaryMember,
                                          "ascii boundary member",
                                          "Implementation of a model in which the boundary "
                                          "velocity is derived from files containing data "
                                          "in ascii format. Note the required format of the "
                                          "input data: The first lines may contain any number of comments "
                                          "if they begin with '#', but one of these lines needs to "
                                          "contain the number of grid points in each dimension as "
                                          "for example '# POINTS: 3 3'. "
                                          "The order of the data columns "
                                          "has to be 'x', 'velocity_x', 'velocity_y' in a 2d model "
                                          "or 'x', 'y', 'velocity_x', 'velocity_y', "
                                          "'velocity_z' in a 3d model. "
                                          "Note that the data in the input "
                                          "files need to be sorted in a specific order: "
                                          "the first coordinate needs to ascend first, "
                                          "followed by the second in order to "
                                          "assign the correct data to the prescribed coordinates."
                                          "If you use a spherical model, "
                                          "then the velocities will still be handled as Cartesian, "
                                          "however the assumed grid changes. 'x' will be replaced by "
                                          "the radial distance of the point to the bottom of the model, "
                                          "'y' by the azimuth angle and 'z' by the polar angle measured "
                                          "positive from the north pole. The grid will be assumed to be "
                                          "a latitude-longitude grid. Note that the order "
                                          "of spherical coordinates is 'r', 'phi', 'theta' "
                                          "and not 'r', 'theta', 'phi', since this allows "
                                          "for dimension independent expressions. "
                                          "No matter which geometry model is chosen, "
                                          "the unit of the velocities is assumed to be "
                                          "m/s or m/yr depending on the 'Use years in output "
                                          "instead of seconds' flag.")
}
