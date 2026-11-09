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


#include <aspect/prescribed_dilation/function.h>

#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace PrescribedDilation
  {
    template <int dim>
    Function<dim>::Function ()
      :
      function (1)
    {}





    template <int dim>
    void
    Function<dim>::
    evaluate ( const aspect::MaterialModel::MaterialModelInputs<dim> &in,
               PrescribedDilationOutputs<dim> &out) const
    {
      Assert (  out.n_evaluation_points() == in.n_evaluation_points(),
                ExcMessage("The number of points in the MaterialModelInputs and "
                           "PrescribedDilationOutputs do not match."));
      for ( unsigned int q = 0; q < out.n_evaluation_points(); ++q )
        {
          const Utilities::NaturalCoordinate<dim> point =
            this->get_geometry_model().cartesian_to_other_coordinates(in.position[q], coordinate_system);
          const Point<dim> point_in_coordinate_system = Utilities::convert_array_to_point<dim>(point.get_coordinates());
          out.dilation[q] = function.value(point_in_coordinate_system);
        }
    }



    template <int dim>
    void
    Function<dim>::update ()
    {
      const double time = this->get_time();
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        function.set_time (time / year_in_seconds);
      else
        function.set_time (time);
    }



    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed dilation");
      {
        prm.enter_subsection("Function");
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

          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed dilation");
      {
        prm.enter_subsection("Function");
        {
          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        }
        try
          {
            function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Prescribed dilation.Function'\n"
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
  namespace PrescribedDilation
  {
    ASPECT_REGISTER_DILATION_MODEL(Function,
                                   "function",
                                   "Prescribed dilation is given in terms of an explicit formula "
                                   "that is elaborated in the parameters in section "
                                   "``Prescribed dilation|Function''. The format of these "
                                   "functions follows the syntax understood by the "
                                   "muparser library, see "
                                   "{ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`."
                                   "\n\n"
                                   "Physical units: $\\si{\\per\\second}$.")
  }
}
