/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/postprocess/rotation_statistics.h>
#include <aspect/material_model/simple.h>
#include <aspect/global.h>
#include <aspect/simulator.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace
    {
      void add_scientific_column(const std::string &name,
                                 const double value,
                                 TableHandler &statistics)
      {
        statistics.add_value(name,value);
        statistics.set_precision (name, 8);
        statistics.set_scientific (name, true);
      }
    }



    template <int dim>
    std::pair<std::string,std::string>
    RotationStatistics<dim>::execute (TableHandler &statistics)
    {
      RotationProperties<dim> rotation = this->compute_net_angular_momentum(use_constant_density,
                                                                            this->get_solution(),
                                                                            false);

      RotationProperties<dim> surface_rotation = this->compute_net_angular_momentum(true,
                                                                                    this->get_solution(),
                                                                                    true);

      const std::vector<std::string> names = {"Angular momentum", "Moment of inertia", "Angular velocity", "Surface angular velocity"};
      std::vector<std::string> units;

      if (this->convert_output_to_years() == true)
        {
          units = {"kg*m^2/year", "kg*m^2", "1/year", "1/year"};

          if (dim == 2)
            {
              rotation.scalar_angular_momentum *= year_in_seconds;
              rotation.scalar_rotation *= year_in_seconds;
              surface_rotation.scalar_rotation *= year_in_seconds;
            }
          else if (dim == 3)
            {
              rotation.tensor_angular_momentum *= year_in_seconds;
              rotation.tensor_rotation *= year_in_seconds;
              surface_rotation.tensor_rotation *= year_in_seconds;
            }
        }
      else
        {
          units = {"kg*m^2/s", "kg*m^2", "1/s", "1/s"};
        }

      std::ostringstream output;
      output.precision(3);

      if (dim == 2)
        {
          add_scientific_column(names[0] + " (" + units[0] +")", rotation.scalar_angular_momentum, statistics);
          add_scientific_column(names[1] + " (" + units[1] +")", rotation.scalar_moment_of_inertia, statistics);
          add_scientific_column(names[2] + " (" + units[2] +")", rotation.scalar_rotation, statistics);
          add_scientific_column(names[3] + " (" + units[3] +")", surface_rotation.scalar_rotation, statistics);

          output << rotation.scalar_angular_momentum << ' ' << units[0] << ", "
                 << rotation.scalar_moment_of_inertia << ' ' << units[1] << ", "
                 << rotation.scalar_rotation << ' ' << units[2] << ", "
                 << surface_rotation.scalar_rotation << ' ' << units[3];
        }
      else if (dim == 3)
        {
          add_scientific_column(names[0] + " (" + units[0] +")", rotation.tensor_angular_momentum.norm(), statistics);

          const double scalar_moment_of_inertia = (rotation.tensor_moment_of_inertia * rotation.tensor_rotation / rotation.tensor_rotation.norm()).norm();

          if (output_full_tensor == false)
            {
              add_scientific_column(names[1] + " (" + units[1] +")", scalar_moment_of_inertia, statistics);
            }
          else
            {
              add_scientific_column(names[1] + "_xx (" + units[1] +")", rotation.tensor_moment_of_inertia[0][0], statistics);
              add_scientific_column(names[1] + "_yy (" + units[1] +")", rotation.tensor_moment_of_inertia[1][1], statistics);
              add_scientific_column(names[1] + "_zz (" + units[1] +")", rotation.tensor_moment_of_inertia[2][2], statistics);
              add_scientific_column(names[1] + "_xy (" + units[1] +")", rotation.tensor_moment_of_inertia[0][1], statistics);
              add_scientific_column(names[1] + "_xz (" + units[1] +")", rotation.tensor_moment_of_inertia[0][2], statistics);
              add_scientific_column(names[1] + "_yz (" + units[1] +")", rotation.tensor_moment_of_inertia[1][2], statistics);
            }

          add_scientific_column(names[2] + " (" + units[2] +")", rotation.tensor_rotation.norm(), statistics);
          add_scientific_column(names[3] + " (" + units[3] +")", surface_rotation.tensor_rotation.norm(), statistics);

          output << rotation.tensor_angular_momentum.norm() << ' ' << units[0] << ", "
                 << scalar_moment_of_inertia << ' ' << units[1] << ", "
                 << rotation.tensor_rotation.norm() << ' ' << units[2] << ", "
                 << surface_rotation.tensor_rotation.norm() << ' ' << units[3];
        }

      return std::pair<std::string, std::string> (names[0]+ ", " + names[1] + ", " + names[2] + ", " + names[3] + ":",
                                                  output.str());
    }



    template <int dim>
    void
    RotationStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Rotation statistics");
        {
          prm.declare_entry("Use constant density of one","false",
                            Patterns::Bool(),
                            "Whether to use a constant density of one for the computation of the "
                            "angular momentum and moment of inertia. This is an approximation "
                            "that assumes that the 'volumetric' rotation is equal to the 'mass' "
                            "rotation. If this parameter is true this postprocessor computes "
                            "'net rotation' instead of 'angular momentum'.");
          prm.declare_entry("Output full moment of inertia tensor","false",
                            Patterns::Bool(),
                            "Whether to write the full moment of inertia tensor into the "
                            "statistics output instead of its norm for the current rotation "
                            "axis. This is a second-order symmetric tensor with "
                            "6 components in 3d. In 2d this option has no effect, because "
                            "the rotation axis is fixed and thus the moment of inertia "
                            "is always a scalar.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    RotationStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Rotation statistics");
        {
          use_constant_density = prm.get_bool("Use constant density of one");
          output_full_tensor = prm.get_bool("Output full moment of inertia tensor");
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(RotationStatistics,
                                  "rotation statistics",
                                  "A postprocessor that computes some statistics about the "
                                  "rotational velocity of the model (i.e. integrated "
                                  "net rotation and angular momentum). In 2d we assume the "
                                  "model to be a cross-section through an infinite domain in "
                                  "z direction, with a zero z-velocity. Thus, the z-axis is "
                                  "the only possible rotation axis and both moment of inertia "
                                  "and angular momentum are scalar instead of tensor "
                                  "quantities.")
  }
}
