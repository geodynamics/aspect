// TODO:
// - make R_values a parameter so that one specifies
//    . R0
//    . R1
//    . any intermediate values (if any)
// - make the lateral refinement a parameter
// - write a 2d test
// - write a 3d test


/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#include <aspect/geometry_model/custom_spherical_shell.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>

#include <deal.II/grid/grid_generator.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    CustomSphericalShell<dim>::CustomSphericalShell()
      :
      spherical_manifold()
    {}

    template <int dim>
    void
    CustomSphericalShell<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      Triangulation<dim-1,dim> sphere_mesh;
      GridGenerator::hyper_sphere (sphere_mesh);
      sphere_mesh.refine_global (2);  // define as a run-time parameter

      std::vector<double> R_values (n_slices+1);  // also define as a run-time parameter
      for (unsigned int s=0; s<n_slices+1; ++s)
        {
        R_values[s] = R0 + (R1-R0)/n_slices * s;
        std::cout << R_values[s] << std::endl;
        }
      std::vector<Point<dim>>    points(R_values.size() * sphere_mesh.n_vertices());

      // copy the array of points as many times as there will be slices,
      // one slice at a time. The z-axis value are defined in slices_coordinates
      for (unsigned int point_layer = 0; point_layer < R_values.size(); ++point_layer)
        {
          for (unsigned int vertex_n = 0; vertex_n < sphere_mesh.n_vertices();
               ++vertex_n)
            {
              const Point<dim> vertex = sphere_mesh.get_vertices()[vertex_n];
              points[point_layer * sphere_mesh.n_vertices() + vertex_n] =
                vertex * R_values[point_layer];
            }
        }

      // then create the cells of each of the slices, one stack at a
      // time
      std::vector<CellData<dim>> cells;
      cells.reserve((R_values.size() - 1) * sphere_mesh.n_active_cells());

      SubCellData               subcell_data;

      for (const auto &cell : sphere_mesh.active_cell_iterators())
        {
          for (unsigned int cell_layer = 0; cell_layer < R_values.size() - 1; ++cell_layer)
            {
              CellData<dim> this_cell;
              for (unsigned int vertex_n = 0;
                   vertex_n < GeometryInfo<dim-1>::vertices_per_cell;
                   ++vertex_n)
                {
                  this_cell.vertices[vertex_n] =
                    cell->vertex_index(vertex_n) + cell_layer * sphere_mesh.n_vertices();
                  this_cell.vertices[vertex_n + GeometryInfo<dim-1>::vertices_per_cell] =
                    cell->vertex_index(vertex_n) +
                    (cell_layer + 1) * sphere_mesh.n_vertices();
                }

              this_cell.material_id = cell->material_id();

              cells.push_back(this_cell);


              // Mark the bottom face of the cell as boundary 0 if we are in
              // the bottom layer of cells
              if (cell_layer == 0)
                {
                  CellData<dim-1> face;
                  for (unsigned int vertex_n = 0;
                       vertex_n < GeometryInfo<dim-1>::vertices_per_cell;
                       ++vertex_n)
                    face.vertices[vertex_n] =
                      cell->vertex_index(vertex_n) + cell_layer * sphere_mesh.n_vertices();
                  face.boundary_id = 0;

                  if (dim == 2)
                    subcell_data.boundary_lines.push_back(reinterpret_cast<CellData<1>&>(face));
                  else
                    subcell_data.boundary_quads.push_back(reinterpret_cast<CellData<2>&>(face));
                }

              // Mark the top face of the cell as boundary 1 if we are in
              // the top layer of cells
              if (cell_layer == R_values.size()-2)
                {
                  CellData<dim-1> face;
                  for (unsigned int vertex_n = 0;
                       vertex_n < GeometryInfo<dim-1>::vertices_per_cell;
                       ++vertex_n)
                    face.vertices[vertex_n] =
                      cell->vertex_index(vertex_n) +
                      (cell_layer + 1) * sphere_mesh.n_vertices();
                  face.boundary_id = 1;

                  if (dim == 2)
                    subcell_data.boundary_lines.push_back(reinterpret_cast<CellData<1>&>(face));
                  else
                    subcell_data.boundary_quads.push_back(reinterpret_cast<CellData<2>&>(face));
                }

            }
        }

      // Then create the actual mesh:
      coarse_grid.create_triangulation(points, cells, subcell_data);

      // Use a manifold description for all cells.
      coarse_grid.set_manifold (99, spherical_manifold);
      set_manifold_ids(coarse_grid);
    }



    template <int dim>
    void
    CustomSphericalShell<dim>::set_manifold_ids (parallel::distributed::Triangulation<dim> &triangulation) const
    {
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active(); cell != triangulation.end(); ++cell)
        cell->set_all_manifold_ids (99);
    }


    template <int dim>
    std::set<types::boundary_id>
    CustomSphericalShell<dim>::
    get_used_boundary_indicators () const
    {
      // follow what is described in the documentation of this class.
      // see the documentation of the various GridGenerator::*hyper_shell
      // functions for a description of which boundary indicators are
      // set and how they correlate to what's used below
      const types::boundary_id s[] = { 0, 1 };
      return std::set<types::boundary_id>(&s[0],
                                          &s[sizeof(s)/sizeof(s[0])]);
    }


    template <int dim>
    std::map<std::string,types::boundary_id>
    CustomSphericalShell<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id> ("bottom", 0),
                  std::pair<std::string,types::boundary_id> ("top", 1)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[2]);
          }

          case 3:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("bottom", 0),
                  std::pair<std::string,types::boundary_id>("top",    1)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[2]);
          }
        }

      Assert (false, ExcNotImplemented());
      return std::map<std::string,types::boundary_id>();
    }


    template <int dim>
    double
    CustomSphericalShell<dim>::
    length_scale () const
    {
      // As described in the first ASPECT paper, a length scale of
      // 10km = 1e4m works well for the pressure scaling for earth
      // sized spherical shells. So use a formulation that yields this
      // value for the R0,R1 corresponding to earth, and that more
      // generally scales with (R1-R0), which we can get by calling
      // maximal_depth(). In essence, the factor in front of the call
      // to maximal_depth() is just a magic number that has turned out
      // to work well.
      return (1e4 / (6336000.-3481000.)) * maximal_depth();
    }



    template <int dim>
    double
    CustomSphericalShell<dim>::depth(const Point<dim> &position) const
    {
      return std::min (std::max (R1-position.norm(), 0.), maximal_depth());
    }

    template <int dim>
    double
    CustomSphericalShell<dim>::height_above_reference_surface(const Point<dim> &position) const
    {
      return position.norm()-outer_radius();
    }


    template <int dim>
    Point<dim>
    CustomSphericalShell<dim>::representative_point(const double depth) const
    {
      Point<dim> p;
      p(dim-1) = std::min (std::max(R1 - depth, R0), R1);
      return p;
    }



    template <int dim>
    double
    CustomSphericalShell<dim>::maximal_depth() const
    {
      return R1-R0;
    }

    template <int dim>
    double CustomSphericalShell<dim>::inner_radius () const
    {
      return R0;
    }



    template <int dim>
    double CustomSphericalShell<dim>::outer_radius () const
    {
      return R1;
    }



    template <int dim>
    bool
    CustomSphericalShell<dim>::has_curved_elements () const
    {
      return true;
    }


    template <int dim>
    bool
    CustomSphericalShell<dim>::point_is_in_domain(const Point<dim> &point) const
    {
      AssertThrow(this->get_free_surface_boundary_indicators().size() == 0 ||
                  this->get_timestep_number() == 0,
                  ExcMessage("After displacement of the free surface, this function can no longer be used to determine whether a point lies in the domain or not."));

      AssertThrow(dynamic_cast<const InitialTopographyModel::ZeroTopography<dim>*>(&this->get_initial_topography_model()) != nullptr,
                  ExcMessage("After adding topography, this function can no longer be used to determine whether a point lies in the domain or not."));

      const std::array<double, dim> spherical_point = Utilities::Coordinates::cartesian_to_spherical_coordinates(point);

      std::array<double, dim> point1, point2;
      point1[0] = R0;
      point2[0] = R1;
      point1[1] = 0.0;
      point2[1] = 2 * numbers::PI;
      for (unsigned int d = 0; d < dim; d++)
        if (spherical_point[d] > point2[d]+std::numeric_limits<double>::epsilon()*std::abs(point2[d]) ||
            spherical_point[d] < point1[d]-std::numeric_limits<double>::epsilon()*std::abs(point2[d]))
          return false;

      return true;
    }


    template <int dim>
    aspect::Utilities::Coordinates::CoordinateSystem
    CustomSphericalShell<dim>::natural_coordinate_system() const
    {
      return aspect::Utilities::Coordinates::CoordinateSystem::spherical;
    }


    template <int dim>
    void
    CustomSphericalShell<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Custom spherical shell");
        {
        	// instead just provide a list of R values
          prm.declare_entry ("Inner radius", "3481000",  // 6371-2890 in km
                             Patterns::Double (0),
                             "Inner radius of the spherical shell. Units: m. "
                             "\n\n"
                             "\\note{The default value of 3,481,000 m equals the "
                             "radius of a sphere with equal volume as Earth (i.e., "
                             "6371 km) minus the average depth of the core-mantle "
                             "boundary (i.e., 2890 km).}");
          prm.declare_entry ("Outer radius", "6336000",  // 6371-35 in km
                             Patterns::Double (0),
                             "Outer radius of the spherical shell. Units: m. "
                             "\n\n"
                             "\\note{The default value of 6,336,000 m equals the "
                             "radius of a sphere with equal volume as Earth (i.e., "
                             "6371 km) minus the average depth of the mantle-crust "
                             "interface (i.e., 35 km).}");
          prm.declare_entry ("Cells along circumference", "0",
                             Patterns::Integer (0),
                             "...");

          // remove: duplicative of R_values
          prm.declare_entry ("Number of slices", "2",
                             Patterns::Integer (0),
                             "Number of slices for the tailored spherical shell case."
                             "\n\n"
                             "The number of slices is used for extrusion and must be "
                             "at least 2");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    CustomSphericalShell<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Custom spherical shell");
        {
          R0  = prm.get_double ("Inner radius");
          R1  = prm.get_double ("Outer radius");
          n_cells_along_circumference = prm.get_integer ("Cells along circumference");
          n_slices = prm.get_integer ("Number of slices");

          AssertThrow (R0 < R1,
                       ExcMessage ("Inner radius must be less than outer radius."));
          AssertThrow (n_slices > 0, ExcMessage("You must set a positive number of slices for extrusion. "));
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
  namespace GeometryModel
  {
    ASPECT_REGISTER_GEOMETRY_MODEL(CustomSphericalShell,
                                   "custom spherical shell",
                                   "...")
  }
}
