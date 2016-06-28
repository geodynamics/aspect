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

#include <aspect/postprocess/topography.h>

#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/fe/fe_values.h>

#include <cmath>
#include <limits>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    Topography<dim>::execute (TableHandler &statistics)
    {

      double reference_height = 0.0;
      bool vertical_gravity = false;
      types::boundary_id relevant_boundary = 0;

      if (GeometryModel::Box<dim> *gm = dynamic_cast<GeometryModel::Box<dim> *>
                                        (const_cast<GeometryModel::Interface<dim> *>(&this->get_geometry_model())))
        {
          Point<dim> extents = gm->get_extents();
          reference_height = extents[dim-1];
          vertical_gravity = true;
          relevant_boundary = (dim == 2 ? 3 : 5); //select top boundary
        }
      else if (GeometryModel::Sphere<dim> *gm = dynamic_cast<GeometryModel::Sphere<dim> *>
                                                (const_cast<GeometryModel::Interface<dim> *>(&this->get_geometry_model())))
        {
          reference_height = gm->radius();
          vertical_gravity = false;
          relevant_boundary = 0;  //select top boundary
        }
      else if (GeometryModel::SphericalShell<dim> *gm = dynamic_cast<GeometryModel::SphericalShell<dim> *>
                                                        (const_cast<GeometryModel::Interface<dim> *>(&this->get_geometry_model())))
        {
          reference_height = gm->outer_radius();
          vertical_gravity = false;
          relevant_boundary = 1;  //select top boundary
        }
      else
        {
          AssertThrow(false, ExcMessage("The topography postprocessor does not recognize the geometry model."
                                        "Consider using a box, a spherical shell, or a sphere.") );
        }

      //Get a quadrature rule that exists only on the corners
      QTrapez<dim-1> face_corners;
      FEFaceValues<dim> face_vals (this->get_mapping(), this->get_fe(), face_corners, update_quadrature_points);

      //get maximum surface topography
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = this->get_triangulation().begin_active(),
                                                                               endc = this->get_triangulation().end();

      //Choose stupidly large values for initialization
      double local_max_height = -std::numeric_limits<double>::max();
      double local_min_height = std::numeric_limits<double>::max();

      for (; cell != endc; ++cell)
        if (cell->is_locally_owned() && cell->at_boundary())
          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
            if (cell->face(face_no)->at_boundary())
              {
                if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                  continue;

                face_vals.reinit( cell, face_no);

                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                    Point<dim> vertex = face_vals.quadrature_point(corner);

                    double topography = (vertical_gravity ? vertex[dim-1] : vertex.norm()) - reference_height;
                    if ( topography > local_max_height)
                      local_max_height = topography;
                    if ( topography < local_min_height)
                      local_min_height = topography;

                  }
              }

      double max_topography = Utilities::MPI::max(local_max_height, this->get_mpi_communicator());
      double min_topography = Utilities::MPI::min(local_min_height, this->get_mpi_communicator());

      statistics.add_value ("Minimum topography (m)",
                            min_topography);
      statistics.add_value ("Maximum topography (m)",
                            max_topography);
      {
        const char *columns[] = { "Minimum topography (m)",
                                  "Maximum topography (m)"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }

      std::ostringstream output;
      output.precision(4);
      output << min_topography << " m, "
             << max_topography << " m, ";

      return std::pair<std::string, std::string> ("Topography min/max:",
                                                  output.str());
    }
  }
}


namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Topography,
                                  "topography",
                                  "A postprocessor intended for use with a free surface.  After every step, "
                                  "it loops over all the vertices on the top surface and determines the "
                                  "maximum and minimum topography relative to a reference datum (initial "
                                  "box height for a box geometry model or initial radius for a sphere/"
                                  "spherical shell geometry model).  Outputs topography in meters")
  }
}
