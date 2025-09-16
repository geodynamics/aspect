/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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


#include <aspect/postprocess/current_surface.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>

#include <cmath>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    CurrentSurface<dim>::execute (TableHandler &)
    {
      AssertThrow(Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("The current surface postprocessor is only implemented with a 2D box model."));

      AssertThrow(dim==2,
                  ExcMessage("Depth with mesh deformation currently only works with a 2D model."));

      this->get_computing_timer().enter_subsection("Geometry model surface update");

      // loop over all of the surface cells and save the elevation to a stored value.
      // This needs to be sent to 1 processor, sorted, and broadcast so that every processor knows the entire surface.
      std::vector<std::vector<double>> local_surface_height;
      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
      const QTrapezoid<dim-1> face_corners;
      FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                       this->get_fe(),
                                       face_corners,
                                       update_quadrature_points);

      // Loop over all corners at the surface and save their position.
      // TODO: Update this to work in 3D and spherical geometries
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int face_no : cell->face_indices())
            if (cell->face(face_no)->at_boundary())
              if ( cell->face(face_no)->boundary_id() == relevant_boundary)
                {
                  fe_face_values.reinit(cell, face_no);

                  for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                    {
                      const Point<dim> vertex = fe_face_values.quadrature_point(corner);

                      // We can't push back a point so we convert it into a vector.
                      // This is needed later to keep the vertices together when sorting.
                      std::vector<double> vertex_row;
                      for (unsigned int i=0; i<dim; ++i)
                        vertex_row.push_back(vertex[i]);

                      local_surface_height.push_back(vertex_row);
                    }
                }

      // Combine all local_surfaces and broadcast back.
      std::vector<std::vector<double>> temp_surface =
        Utilities::MPI::compute_set_union(local_surface_height, this->get_mpi_communicator());

      // Sort the vector so that it ascends in X.
      std::sort(temp_surface.begin(), temp_surface.end(), [](const std::vector<double> &a, const std::vector<double> &b)
      {
        return a[0] < b[0];
      });

      // Define a comparison to remove duplicate surface points.
      const auto compareRows = [](const std::vector<double> &row1, const std::vector<double> &row2)
      {
        return row1 == row2;
      };

      // Remove non-unique rows from the sorted 2D vector
      const auto last = std::unique(temp_surface.begin(), temp_surface.end(), compareRows);
      temp_surface.erase(last, temp_surface.end());

      // Resize data table.
      TableIndices<dim-1> size_idx;
      for (unsigned int d=0; d<dim-1; ++d)
        size_idx[d] = temp_surface.size();

      data_table.TableBase<dim-1,double>::reinit(size_idx);
      TableIndices<dim-1> idx;

      // Fill the data table with the y-values that correspond to the surface.
      // This only works in 2D and will need to be updated for 3D models.
      for (unsigned int x=0; x<data_table.size()[0]; ++x)
        {
          idx[0] = x;
          data_table(idx) = temp_surface[x][1];
        }

      // Fill the coordinates with the x-values used for the data table.
      // This only works in 2D, and will need to be updated for 3D.
      coordinates[0].clear();
      for (unsigned int i=0; i<temp_surface.size(); ++i)
        coordinates[0].push_back(temp_surface[i][0]);

      // Create a surface function for the elevations.
      surface_function = std::make_unique<Functions::InterpolatedTensorProductGridData<dim-1>>(coordinates, data_table);

      this->get_computing_timer().leave_subsection("Geometry model surface update");

      return std::make_pair ("Storing deformed surface topography: ",
                             "Done");
    }



    template <int dim>
    double
    CurrentSurface<dim>::depth_including_mesh_deformation(const Point<dim> &position) const
    {
      AssertThrow(Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("Depth with mesh deformation currently only works with a 2D box geometry model."));

      AssertThrow(dim==2,
                  ExcMessage("Depth with mesh deformation currently only works with a 2D geometry model."));

      AssertThrow (this->get_parameters().mesh_deformation_enabled,
                   ExcMessage("Depth with mesh deformation must be used with mesh deformation activated."));

      // Convert the point to dim-1, as we aren't interested in the vertical component
      // for the function.
      Point<dim-1> p;
      for (unsigned int d=0; d<dim-1; ++d)
        p[d] = position[d];

      double depth_from_surface = surface_function->value(p) - position[dim-1];
      return std::max (depth_from_surface, 0.);
    }



    template <int dim>
    template <class Archive>
    void CurrentSurface<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &coordinates
      & data_table;
    }



    template <int dim>
    void
    CurrentSurface<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      // Serialize into a stringstream. Put the following into a code
      // block of its own to ensure the destruction of the 'oa'
      // archive triggers a flush() on the stringstream so we can
      // query the completed string below.
      std::ostringstream os;
      {
        aspect::oarchive oa (os);
        oa << (*this);
      }

      status_strings["CurrentSurface"] = os.str();
    }


    template <int dim>
    void
    CurrentSurface<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("CurrentSurface") != status_strings.end())
        {
          std::istringstream is (status_strings.find("CurrentSurface")->second);
          aspect::iarchive ia (is);
          ia >> (*this);
        }

      // Recreate function after loading data.
      surface_function = std::make_unique<Functions::InterpolatedTensorProductGridData<dim-1>>(coordinates, data_table);
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(CurrentSurface,
                                  "current surface",
                                  "A postprocessor that computes a function "
                                  "of the surface that includes the mesh deformation. "
                                  "This postprocessor has a function that can be called from other "
                                  "plugins to get the depth.")
  }
}
