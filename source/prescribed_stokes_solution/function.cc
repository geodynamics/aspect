// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/global.h>
#include <aspect/prescribed_stokes_solution/function.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace PrescribedStokesSolution
  {
    template <int dim>
    Function<dim>::Function ()
      :
      prescribed_velocity_function (dim),
      prescribed_pressure_function (1),
      prescribed_fluid_pressure_function (1),
      prescribed_compaction_pressure_function (1),
      prescribed_fluid_velocity_function (dim)
    {}

    template <int dim>
    void
    Function<dim>::stokes_solution (const Point<dim> &position, Vector<double> &value) const
    {
      // convert the position into the selected coordinate system
      const Utilities::NaturalCoordinate<dim> point = this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system);
      const Point<dim> converted_point = Utilities::convert_array_to_point<dim>(point.get_coordinates());

      // velocity
      for (unsigned int d=0; d<dim; ++d)
        {
          value[d] = prescribed_velocity_function.value(converted_point,d);

          if (this->convert_output_to_years())
            value[d] /= year_in_seconds;
        }

      if (this->get_parameters().include_melt_transport)
        {
          // fluid pressure
          value(dim) = prescribed_fluid_pressure_function.value(converted_point);

          // compaction pressure
          value(dim+1) = prescribed_compaction_pressure_function.value(converted_point);

          // fluid velocity
          for (unsigned int d=0; d<dim; ++d)
            {
              value[dim+2+d] = prescribed_fluid_velocity_function.value(converted_point,d);

              if (this->convert_output_to_years())
                value[dim+2+d] /= year_in_seconds;
            }

          // pressure
          value(2*dim+2) = prescribed_pressure_function.value(converted_point);
        }
      else
        {
          // pressure
          value(dim) = prescribed_pressure_function.value(converted_point);
        }
    }

    template <int dim>
    void
    Function<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        {
          prescribed_velocity_function.set_time (this->get_time() / year_in_seconds);
          prescribed_pressure_function.set_time (this->get_time() / year_in_seconds);
          prescribed_fluid_pressure_function.set_time (this->get_time() / year_in_seconds);
          prescribed_compaction_pressure_function.set_time (this->get_time() / year_in_seconds);
          prescribed_fluid_velocity_function.set_time (this->get_time() / year_in_seconds);
        }
      else
        {
          prescribed_velocity_function.set_time (this->get_time());
          prescribed_pressure_function.set_time (this->get_time());
          prescribed_fluid_pressure_function.set_time (this->get_time());
          prescribed_compaction_pressure_function.set_time (this->get_time());
          prescribed_fluid_velocity_function.set_time (this->get_time());
        }
    }

    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed Stokes solution");
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
        prm.enter_subsection("Velocity function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
        }
        prm.leave_subsection();
        prm.enter_subsection("Pressure function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
        prm.enter_subsection("Fluid pressure function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
        prm.enter_subsection("Compaction pressure function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
        prm.enter_subsection("Fluid velocity function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed Stokes solution");
      {
        coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        prm.enter_subsection("Velocity function");
        try
          {
            prescribed_velocity_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Prescribed Stokes solution.Function.Velocity function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }
        prm.leave_subsection();
        prm.enter_subsection("Pressure function");
        try
          {
            prescribed_pressure_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Prescribed Stokes solution.Function.Pressure function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }
        prm.leave_subsection();
        prm.enter_subsection("Fluid pressure function");
        try
          {
            prescribed_fluid_pressure_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Prescribed Stokes solution.Function.Fluid pressure function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }
        prm.leave_subsection();
        prm.enter_subsection("Compaction pressure function");
        try
          {
            prescribed_compaction_pressure_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Prescribed Stokes solution.Function.Compaction pressure function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }
        prm.leave_subsection();
        prm.enter_subsection("Fluid velocity function");
        try
          {
            prescribed_fluid_velocity_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Prescribed Stokes solution.Function.Fluid velocity function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace PrescribedStokesSolution
  {
    ASPECT_REGISTER_PRESCRIBED_STOKES_SOLUTION(Function,
                                               "function",
                                               "This plugin allows to prescribe the Stokes solution "
                                               "for the velocity and pressure field in terms of an "
                                               "explicit formula. The format of these "
                                               "functions follows the syntax understood by the "
                                               "muparser library, see "
                                               "{ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`.")
  }
}
