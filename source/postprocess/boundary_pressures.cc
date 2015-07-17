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


#include <aspect/postprocess/boundary_pressures.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    BoundaryPressures<dim>::execute (TableHandler &statistics)
    {
      const QGauss<dim-1> quadrature_formula_face (this->get_fe()
                                                   .base_element(this->introspection().base_elements.pressure)
                                                   .degree+1);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula_face,
                                        update_values |
                                        update_gradients |
                                        update_q_points |
                                        update_JxW_values);

      double local_top_pressure = 0.;
      double local_bottom_pressure = 0.;
      double local_top_area = 0.;
      double local_bottom_area = 0.;

      std::vector<double> pressure_vals( fe_face_values.n_quadrature_points );

      // loop over all of the surface cells and if one less than h/3 away from
      // the top or bottom surface, evaluate the pressure on that face
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          if (cell->at_boundary())
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              {
                bool cell_at_top = false;
                bool cell_at_bottom = false;

                //Test for top or bottom surface cell faces
                if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center())
                    < cell->face(f)->minimum_vertex_distance()/3.)
                  cell_at_top = true;
                if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center())
                    > (this->get_geometry_model().maximal_depth() - cell->face(f)->minimum_vertex_distance()/3.))
                  cell_at_bottom = true;

                if ( cell_at_top || cell_at_bottom )
                  {
                    //evaluate the pressure on the face
                    fe_face_values.reinit (cell, f);
                    fe_face_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(), pressure_vals);

                    //calculate the top properties
                    if (cell_at_top)
                      for ( unsigned int q = 0; q < fe_face_values.n_quadrature_points; ++q)
                        {
                          local_top_pressure += pressure_vals[q] * fe_face_values.JxW(q);
                          local_top_area += fe_face_values.JxW(q);
                        }
                    if (cell_at_bottom)
                      for ( unsigned int q = 0; q < fe_face_values.n_quadrature_points; ++q)
                        {
                          local_bottom_pressure += pressure_vals[q] * fe_face_values.JxW(q);
                          local_bottom_area += fe_face_values.JxW(q);
                        }
                  }
              }

      //vector for packing local values before mpi summing them
      double values[4] = {local_bottom_area, local_top_area, local_bottom_pressure, local_top_pressure};

      Utilities::MPI::sum<double, 4>( values, this->get_mpi_communicator(), values );

      top_pressure = values[3] / values[1]; //density over area
      bottom_pressure = values[2] / values[0]; //density over area

      statistics.add_value ("Pressure at top (Pa)",
                            top_pressure);
      statistics.add_value ("Pressure at bottom (Pa)",
                            bottom_pressure);

      // also make sure that the other columns filled by the this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Pressure at top (Pa)",
                                  "Pressure at bottom (Pa)"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }

      std::ostringstream output;
      output.precision(4);
      output << top_pressure << " Pa, "
             << bottom_pressure << " Pa, ";

      return std::pair<std::string, std::string> ("Pressure at top/bottom of domain:",
                                                  output.str());
    }

    template <int dim>
    double
    BoundaryPressures<dim>::pressure_at_top() const
    {
      return top_pressure;
    }

    template <int dim>
    double
    BoundaryPressures<dim>::pressure_at_bottom() const
    {
      return bottom_pressure;
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(BoundaryPressures,
                                  "boundary pressures",
                                  "A postprocessor that computes the laterally averaged "
                                  "pressure at the top and bottom of the domain.")
  }
}
