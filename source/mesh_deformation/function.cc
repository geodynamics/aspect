// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2018-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code
// and contributions can be found in the folder LICENSES
// and in CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/mesh_deformation/function.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    BoundaryFunction<dim>::BoundaryFunction()
      :
      function(dim)
    {}



    template <int dim>
    void
    BoundaryFunction<dim>::update ()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        function.set_time (this->get_time() / year_in_seconds);
      else
        function.set_time (this->get_time());
    }



    /**
     * A function that creates constraints for the velocity of certain mesh
     * vertices (e.g. the surface vertices) for a specific boundary.
     * The calling class will respect
     * these constraints when computing the new vertex positions.
     */
    template <int dim>
    void
    BoundaryFunction<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                                    AffineConstraints<double> &mesh_velocity_constraints,
                                                                    const std::set<types::boundary_id> &boundary_ids) const
    {
      // Loop over all boundary indicators to set the velocity constraints
      for (const auto boundary_id : boundary_ids)
        {
          Utilities::VectorFunctionFromVelocityFunctionObject<dim> vel
          (dim,
           [&] (const dealii::Point<dim> &x) -> Tensor<1,dim>
          {
            Tensor<1,dim> velocity;

            // convert the position into the selected coordinate system
            const Utilities::NaturalCoordinate<dim> point = this->get_geometry_model().cartesian_to_other_coordinates(x, coordinate_system);

            for (unsigned int d=0; d<dim; ++d)
              velocity[d] = function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()), d);

            if (this->convert_output_to_years())
              velocity /= year_in_seconds;

            return velocity;
          });

          VectorTools::interpolate_boundary_values (this->get_mapping(),
                                                    mesh_deformation_dof_handler,
                                                    boundary_id,
                                                    vel,
                                                    mesh_velocity_constraints);
        }
    }



    template <int dim>
    bool
    BoundaryFunction<dim>::
    needs_surface_stabilization () const
    {
      return false;
    }



    template <int dim>
    void BoundaryFunction<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Boundary function");
        {
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2d/3d "
                             "respectively with theta being the polar angle. `depth' "
                             "will create a function, in which only the first "
                             "parameter is non-zero, which is interpreted to "
                             "be the depth of the point.");
          Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void BoundaryFunction<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("Boundary function");
        {
          try
            {
              function.parse_parameters (prm);
            }
          catch (...)
            {
              std::cerr << "ERROR: FunctionParser failed to parse\n"
                        << "\t'Mesh deformation.BoundaryFunction'\n"
                        << "with expression\n"
                        << "\t'" << prm.get("Function expression") << "'"
                        << "More information about the cause of the parse error \n"
                        << "is shown below.\n";
              throw;
            }

          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
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
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(BoundaryFunction,
                                           "boundary function",
                                           "A plugin, which prescribes the surface mesh to "
                                           "deform according to an analytically prescribed "
                                           "function. Note that the function prescribes a "
                                           "deformation velocity, i.e. the return value of "
                                           "this plugin is later multiplied by the time step length "
                                           "to compute the displacement increment in this time step. "
                                           "The function's time variable is interpreted as "
                                           "years when Use years instead of seconds is set to true, "
                                           "as is the boundary deformation velocity. Otherwise it should be given "
                                           "in m/s. The format of the "
                                           "functions follows the syntax understood by the "
                                           "muparser library, see {ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`.")
  }
}
