/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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

#include "drip_depth.h"
#include <aspect/geometry_model/box.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    DripDepth<dim>::execute (TableHandler &statistics)
    {
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());
      AssertThrow (geometry != 0,
                   ExcMessage ("This initial condition can only be used if the geometry "
                               "is a box."));

      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(0).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points);

      std::vector<double> compositional_values_drip(n_q_points);
      std::vector<Point<dim> > position_values(n_q_points);

      // The minimal depth of the drip
      double local_drip_depth      = 0.0;

      // The maximal depth of the drip
      const double max_depth = geometry->maximal_depth();

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            // Retrieve position of quadrature points
            position_values = fe_values.get_quadrature_points();
            // Only look at the drip on the left side of the domain.
            if (position_values[0][0] > 250e3)
              continue;

            // Retrieve compositional values for the drip and trench
            fe_values[this->introspection().extractors.compositional_fields[0]].get_function_values
            (this->get_solution(), compositional_values_drip);

            // Obtain current drip tip depth and most left point of trench in a loop over the quadrature points
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                // Depth of quadrature point
                const double depth = geometry->depth(position_values[q]);

                // The compositional fields go from 0 to 1, so 0.5 is taken to represent their boundary
                if (compositional_values_drip[q] >= 0.5 && depth > local_drip_depth)
                  local_drip_depth = depth;
              }
          }


      // Compute the maximum drip tip depth over all processors
      const  double global_drip_depth =
        Utilities::MPI::max (local_drip_depth, this->get_mpi_communicator());

      AssertThrow(global_drip_depth<=max_depth+1000.,ExcMessage("Computed drip depth bigger than domain depth."));

      // Add output to the statistics file
      statistics.add_value ("Drip depth [km]", global_drip_depth/1000.0);

      // Also make sure that the other columns filled by the this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Drip depth [km]"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }

      // Add output to the log file
      std::ostringstream output;
      output.precision(3);
      output << global_drip_depth/1000.0
             << " km";
      return std::pair<std::string, std::string> ("Drip depth[km]:",
                                                  output.str());
    }

    template <int dim>
    void
    DripDepth<dim>::declare_parameters (ParameterHandler &)
    {
    }



    template <int dim>
    void
    DripDepth<dim>::parse_parameters (ParameterHandler &)
    {
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(DripDepth,
                                  "drip depth",
                                  "A postprocessor that computes the depth of the left drip "
                                  "in the Rayleigh Taylor instability benchmark described in Kaus et al. (2010) "
                                  "and Rose et al. (2017) ")
  }
}
