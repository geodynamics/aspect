/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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


#include "trench_location.h"
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    TrenchLocation<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the compositional element alone.
      AssertThrow(this->introspection().n_compositional_fields > 0,
                  ExcMessage("This postprocessor cannot be used without compositional fields."));
      const Quadrature<dim - 1> &quadrature_formula = this->introspection().face_quadratures.compositional_fields;

      FEFaceValues<dim>
      fe_face_values(this->get_mapping(),
                     this->get_fe(),
                     quadrature_formula,
                     update_values |
                     update_quadrature_points |
                     update_JxW_values);

      const unsigned int n_q_points = fe_face_values.n_quadrature_points;

      // Vector to store the values of the field that defines the trench.
      std::vector<double> compositional_values_trench(n_q_points);
      // Vector to store the positions of the quadrature points.
      std::vector<Point<dim>> position_values(n_q_points);
      // For each process, find the right most occurrence of the trench field.
      double local_trench_location = 0.;

      const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      // For cells with a face along the domain surface, compute the trench location
      // by looping over these cells and if the trench field is present at a location
      // further right than the current local_trench_location, set the local_trench_location
      // to this new x-coordinate.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          for (const unsigned int f : cell->face_indices())
            if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
              {
                fe_face_values.reinit(cell,f);

                position_values = fe_face_values.get_quadrature_points();

                fe_face_values[this->introspection().extractors.compositional_fields[selected_field]].get_function_values(this->get_solution(), compositional_values_trench);

                // Calculate the rightmost occurrence of the selected field per processor.
                // A field is counted as present when its value is equal or higher than 0.5.
                for (unsigned int q = 0; q < n_q_points; ++q)
                  {
                    if (compositional_values_trench[q] >= 0.5 && position_values[q][0] > local_trench_location)
                      {
                        local_trench_location = position_values[q][0];
                      }
                  }
              }

      // compute the rightmost point over all processors
      const  double global_trench_location =
        Utilities::MPI::max (local_trench_location, this->get_mpi_communicator());

      statistics.add_value ("Trench location [m]", global_trench_location);

      // also make sure that the other columns filled by the this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Trench location [m]"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }

      std::ostringstream output;
      output.precision(4);
      output << global_trench_location
             << " m";
      return std::pair<std::string, std::string> ("Trench location [m]:",
                                                  output.str());
    }

    template <int dim>
    void TrenchLocation<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Trench location");
        {
          prm.declare_entry("Name of trench compositional field", "",
                            Patterns::Anything(),
                            "The name of the trench compositional field that "
                            "you want to compute the rightmost position for.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void TrenchLocation<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Trench location");
        {
          selected_field = this->introspection().compositional_index_for_name(prm.get("Name of trench compositional field"));
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
    ASPECT_REGISTER_POSTPROCESSOR(TrenchLocation,
                                  "trench location",
                                  "A postprocessor that computes the trench location at the domain surface, "
                                  "i.e., the rightmost point of a given compositional field at the surface.")
  }
}
