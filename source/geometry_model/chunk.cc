/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>
#include <aspect/geometry_model/initial_topography_model/ascii_data.h>

#include <aspect/simulator_signals.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>


namespace aspect
{
  namespace GeometryModel
  {
    namespace internal
    {
      template <int dim>
      ChunkGeometry<dim>::ChunkGeometry(const InitialTopographyModel::Interface<dim> &topo,
                                        const double min_longitude,
                                        const double min_radius,
                                        const double max_depth)
        :
        topo (&topo),
        point1_lon(min_longitude),
        inner_radius(min_radius),
        max_depth(max_depth)
      {}



      template <int dim>
      DerivativeForm<1,dim,dim>
      ChunkGeometry<dim>::
      push_forward_gradient(const Point<dim> &chart_point) const
      {
        const double R = chart_point[0]; // Radius

        Assert (R > 0.0, ExcMessage("Negative radius for given point."));

        Tensor<2,dim> DX;

        // We only have access to the topography gradients
        // (dtopo/dphi and dtopo/dtheta) for the AsciiData
        // initial topography model. We assume 0.0 otherwise;
        // while this is valid for the ZeroTopography model,
        // this will result in incorrect gradients for
        // other initial topography models. Hence only AsciiData
        // and ZeroTopography are allowed for now in Chunk<dim>::initialize().
        Tensor<1,dim-1> topo_derivatives;
        if (const InitialTopographyModel::AsciiData<dim> *itm = dynamic_cast<const InitialTopographyModel::AsciiData<dim> *> (topo))
          topo_derivatives = itm->vector_gradient(push_forward_sphere(chart_point));
        else if (dynamic_cast<const InitialTopographyModel::ZeroTopography<dim> *> (topo))
          {
            // Gradient is zero (which it is already initialized to from before)
          }
        else
          Assert (false, ExcNotImplemented());

        // Construct surface point in lon(,lat) coordinates
        Point<dim-1> surface_point;
        for (unsigned int d=0; d<dim-1; ++d)
          surface_point[d] = chart_point[d+1];

        // Convert latitude to colatitude
        if (dim == 3)
          surface_point[1] = 0.5*numbers::PI - surface_point[1];

        // get the maximum topography at the surface point
        const double d_topo = topo->value(surface_point);

        // get the spherical point including topography
        const Point<dim> topo_point = push_forward_topo(chart_point);
        const double R_topo = topo_point[0];
        const double phi_topo = topo_point[1];

        // The derivatives of topo_point to chart_point
        Tensor<2, dim> Dtopo;
        // The derivatives of the cartesian point to chart_point
        DerivativeForm<1, dim, dim> Dtotal;


        switch (dim)
          {
            case 2:
            {
              // R_topo = R + topo(phi) * ((R-R_0)/(R_1-R_0)) = R + topo(phi) * ((R-R_0)/max_depth)
              // phi_topo = phi
              //dR_topo/dR
              Dtopo[0][0] = (d_topo / max_depth) + 1.;
              //dR_topo/dphi = dR_topo/dtopo * dtopo/dphi
              Dtopo[0][1] = (R-inner_radius)/max_depth * topo_derivatives[0];
              //dphi_topo/dR
              Dtopo[1][0] = 0.;
              //dphi_topo/dphi
              Dtopo[1][1] = 1.;

              //dx/dR_topo
              DX[0][0] =           std::cos(phi_topo);
              //dx/dphi_topo
              DX[0][1] = -R_topo * std::sin(phi_topo);
              //dy/dR_topo
              DX[1][0] =           std::sin(phi_topo);
              //dy/dphi_topo
              DX[1][1] =  R_topo * std::cos(phi_topo);

              break;
            }
            case 3:
            {
              // R_topo = R + topo(phi,theta) * ((R-R_0)/(R_1-R_0))
              // phi_topo = phi
              // theta_topo = theta
              //dR_topo/dR
              Dtopo[0][0] = (d_topo / max_depth) + 1.;
              //dR_topo/dphi
              Dtopo[0][1] = (R - inner_radius) / max_depth * topo_derivatives[0];
              //dR_topo/dtheta
              Dtopo[0][2] = (R - inner_radius) / max_depth * topo_derivatives[1];
              //dphi_topo/dR
              Dtopo[1][0] = 0.;
              //dphi_topo/dphi
              Dtopo[1][1] = 1.;
              //dphi_topo/dtheta
              Dtopo[1][2] = 0.;
              //dtheta_topo/dR
              Dtopo[2][0] = 0.;
              //dtheta_topo/dphi
              Dtopo[2][1] = 0.;
              //dtheta_topo/dtheta
              Dtopo[2][2] = 1.;

              const double theta_topo = topo_point[2];

              // The derivatives of the cartesian points to topo_point
              DX[0][0] =      std::cos(theta_topo) * std::cos(phi_topo);
              DX[0][1] = -R_topo * std::cos(theta_topo) * std::sin(phi_topo);
              DX[0][2] = -R_topo * std::sin(theta_topo) * std::cos(phi_topo); //reorder
              DX[1][0] =      std::cos(theta_topo) * std::sin(phi_topo);
              DX[1][1] =  R_topo * std::cos(theta_topo) * std::cos(phi_topo);
              DX[1][2] = -R_topo * std::sin(theta_topo) * std::sin(phi_topo);
              DX[2][0] =      std::sin(theta_topo);
              DX[2][1] = 0;
              DX[2][2] =  R_topo * std::cos(theta_topo);

              break;
            }
            default:
              Assert (false, ExcNotImplemented ());
          }

        Dtotal = DerivativeForm<1,dim,dim>(DX * Dtopo);

        return Dtotal;
      }



      template <int dim>
      Point<dim>
      ChunkGeometry<dim>::
      push_forward(const Point<dim> &r_phi_theta) const
      {
        // Only take into account topography when we're not using the ZeroTopography plugin
        if (dynamic_cast<const InitialTopographyModel::ZeroTopography<dim>*>(topo) != nullptr)
          return push_forward_sphere(r_phi_theta);
        else
          return push_forward_sphere(push_forward_topo(r_phi_theta));

      }



      template <int dim>
      Point<dim>
      ChunkGeometry<dim>::
      pull_back(const Point<dim> &x_y_z) const
      {
        // Only take into account topography when we're not using the ZeroTopography plugin
        if (dynamic_cast<const InitialTopographyModel::ZeroTopography<dim>*>(topo) != nullptr)
          return pull_back_sphere(x_y_z);
        else
          return pull_back_topo(pull_back_sphere(x_y_z));
      }



      template <int dim>
      Point<dim>
      ChunkGeometry<dim>::
      push_forward_sphere(const Point<dim> &input_vertex) const
      {
        Point<dim> output_vertex;
        switch (dim)
          {
            case 2:
            {
              output_vertex[0] = input_vertex[0]*std::cos(input_vertex[1]); // x=rcosphi
              output_vertex[1] = input_vertex[0]*std::sin(input_vertex[1]); // z=rsinphi
              break;
            }
            case 3:
            {
              output_vertex[0] = input_vertex[0]*std::cos(input_vertex[2])*std::cos(input_vertex[1]); // x=rsinthetacosphi
              output_vertex[1] = input_vertex[0]*std::cos(input_vertex[2])*std::sin(input_vertex[1]); // y=rsinthetasinphi
              output_vertex[2] = input_vertex[0]*std::sin(input_vertex[2]); // z=rcostheta
              break;
            }
            default:
              Assert (false, ExcNotImplemented ());
          }
        return output_vertex;
      }



      template <int dim>
      Tensor<1, dim>
      ChunkGeometry<dim>::
      normal_vector(const typename Triangulation<dim>::face_iterator &face,
                    const Point<dim> &p) const
      {
        // Let us first test whether we are on a "horizontal" face
        // (tangential to the sphere).  In this case, the normal vector is
        // easy to compute since it is proportional to the vector from the
        // center to the point 'p'.
        //
        // We test whether a face is horizontal by checking that the vertices
        // all have roughly the same distance from the center: If the
        // maximum deviation for the distances from the vertices to the
        // center is less than 1.e-5 of the distance between vertices (as
        // measured by the minimum distance from any of the other vertices
        // to the first vertex), then we call this a horizontal face.
        constexpr unsigned int max_n_vertices_per_face = (dim==2 ? 2 : 4);
        std::array<double, max_n_vertices_per_face>     distances_to_center {};
        std::array<double, max_n_vertices_per_face - 1> distances_to_first_vertex {};
        distances_to_center[0] = face->vertex(0).norm_square();
        for (unsigned int i = 1; i < face->n_vertices(); ++i)
          {
            AssertIndexRange (i, distances_to_center.size());
            AssertIndexRange (i-1, distances_to_first_vertex.size());

            distances_to_center[i] = face->vertex(i).norm_square();
            distances_to_first_vertex[i - 1] =
              (face->vertex(i) - face->vertex(0)).norm_square();
          }
        const auto minmax_distance =
          std::minmax_element(distances_to_center.begin(),
                              distances_to_center.begin()+face->n_vertices());
        const auto min_distance_to_first_vertex =
          std::min_element(distances_to_first_vertex.begin(),
                           distances_to_first_vertex.begin()+face->n_vertices()-1);

        // So, if this is a "horizontal" face, then just compute the normal
        // vector as the one from the center to the point 'p', adequately
        // scaled.
        if (*minmax_distance.second - *minmax_distance.first <
            1.e-10 * *min_distance_to_first_vertex)
          {
            const Tensor<1, dim> unnormalized_spherical_normal = p;
            const Tensor<1, dim> normalized_spherical_normal =
              unnormalized_spherical_normal / unnormalized_spherical_normal.norm();
            return normalized_spherical_normal;
          }

        // If it is not a horizontal face, just use the machinery of the
        // base class.
        return Manifold<dim>::normal_vector(face, p);
      }



      template <int dim>
      Point<dim>
      ChunkGeometry<dim>::
      pull_back_sphere(const Point<dim> &v) const
      {
        Point<dim> output_vertex;
        switch (dim)
          {
            case 2:
            {
              output_vertex[1] = std::atan2(v[1], v[0]);
              output_vertex[0] = v.norm();

              // We must guarantee that all returned points have a longitude coordinate
              // value that is larger than (or equal to) the longitude of point1.
              // For example:
              // If the domain runs from longitude -10 to 200 degrees,
              // atan2 will also return a negative value (-180 to -160) for the points
              // with longitude 180 to 200. These values must be corrected
              // so that they are larger than the minimum longitude value of -10,
              // by adding 360 degrees.
              // A 100*epsilon ensures we catch all cases.
              if (output_vertex[1] < 0.0)
                if (output_vertex[1] < point1_lon - 100 * std::abs(point1_lon)*std::numeric_limits<double>::epsilon())
                  output_vertex[1] += 2.0 * numbers::PI;
              break;
            }
            case 3:
            {
              const double radius=v.norm();
              output_vertex[0] = radius;
              output_vertex[1] = std::atan2(v[1], v[0]);
              // See 2d case
              if (output_vertex[1] < 0.0)
                if (output_vertex[1] < point1_lon - 100 * std::abs(point1_lon)*std::numeric_limits<double>::epsilon())
                  output_vertex[1] += 2.0 * numbers::PI;
              output_vertex[2] = std::asin(v[2]/radius);
              break;
            }
            default:
              Assert (false, ExcNotImplemented ());
          }
        return output_vertex;
      }



      template <int dim>
      std::unique_ptr<Manifold<dim,dim>>
      ChunkGeometry<dim>::
      clone() const
      {
        return std::make_unique<ChunkGeometry>(*this);
      }



      template <int dim>
      Point<dim>
      ChunkGeometry<dim>::
      push_forward_topo(const Point<dim> &r_phi_theta) const
      {
        // the radius of the current point without topography
        const double radius = r_phi_theta[0];

        // Grab lon,lat coordinates
        Point<dim-1> surface_point;
        for (unsigned int d=0; d<dim-1; ++d)
          surface_point[d] = r_phi_theta[d+1];

        // While throughout ASPECT we use co-latitude as a convention when
        // representing points in spherical coordinates (see for example
        // the Utilities::Coordinates::cartesian_to_spherical_coordinates()
        // function), the current class uses latitude in its pull back and
        // push forward functions. As a consequence, the argument provided
        // to this function has latitude as one coordinate (at least in 3d).
        // On the other hand, the InitialTopography::Interface::value()
        // function expects colatitude (=standard form spherical coordinates)
        // for geometry models that are (parts of) spherical objects.
        //
        // So do the conversion:
        if (dim == 3)
          surface_point[1] = 0.5*numbers::PI - surface_point[1];

        // Then query the topography at this point:
        const double topography = topo->value(surface_point);

        // adjust the radius based on the radius of the point
        // through a linear interpolation between 0 at max depth and
        // "topography" at the surface
        const double topo_radius = std::max(inner_radius,radius + (radius-inner_radius)/max_depth*topography);

        // return the point with adjusted radius
        Point<dim> topo_r_phi_theta = r_phi_theta;
        topo_r_phi_theta[0] = topo_radius;

        return topo_r_phi_theta;
      }



      template <int dim>
      Point<dim>
      ChunkGeometry<dim>::
      pull_back_topo(const Point<dim> &topo_r_phi_theta) const
      {
        // the radius of the point with topography
        const double topo_radius = topo_r_phi_theta[0];

        // Grab lon,lat coordinates
        Point<dim-1> surface_point;
        for (unsigned int d=0; d<dim-1; ++d)
          surface_point[d] = topo_r_phi_theta[d+1];

        // While throughout ASPECT we use co-latitude as a convention when
        // representing points in spherical coordinates (see for example
        // the Utilities::Coordinates::cartesian_to_spherical_coordinates()
        // function), the current class uses latitude in its pull back and
        // push forward functions. As a consequence, the argument provided
        // to this function has latitude as one coordinate (at least in 3d).
        // On the other hand, the InitialTopography::Interface::value()
        // function expects colatitude (=standard form spherical coordinates)
        // for geometry models that are (parts of) spherical objects.
        //
        // So do the conversion:
        if (dim == 3)
          surface_point[1] = 0.5*numbers::PI - surface_point[1];

        // Then query the topography at this point:
        const double topography = topo->value(surface_point);

        // remove the topography (which scales with radius)
        const double radius = std::max(inner_radius,
                                       (topo_radius*max_depth+inner_radius*topography)/(max_depth+topography));

        // return the point without topography
        Point<dim> r_phi_theta = topo_r_phi_theta;
        r_phi_theta[0] = radius;

        return r_phi_theta;
      }



      template <int dim>
      double
      ChunkGeometry<dim>::
      get_radius(const Point<dim> &x_y_z) const
      {
        const Point<dim> r_phi_theta = pull_back(x_y_z);
        Point<dim-1> surface_point;
        for (unsigned int d=0; d<dim-1; ++d)
          surface_point[d] = r_phi_theta[d+1];
        // Convert latitude to colatitude
        if (dim == 3)
          surface_point[1] = 0.5*numbers::PI - surface_point[1];
        const double topography = topo->value(surface_point);

        // return the outer radius at this phi, theta point including topography
        return topography + inner_radius + max_depth;
      }
    }



    template <int dim>
    void
    Chunk<dim>::initialize ()
    {
      AssertThrow(Plugins::plugin_type_matches<const InitialTopographyModel::ZeroTopography<dim>>(this->get_initial_topography_model()) ||
                  Plugins::plugin_type_matches<const InitialTopographyModel::AsciiData<dim>>(this->get_initial_topography_model()),
                  ExcMessage("At the moment, only the Zero or AsciiData initial topography model can be used with the Chunk geometry model."));

      manifold = std::make_unique<internal::ChunkGeometry<dim>>(this->get_initial_topography_model(),
                                                                 point1[1],
                                                                 point1[0],
                                                                 point2[0]-point1[0]);
    }


    template <int dim>
    void
    Chunk<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      // First create a box in Cartesian coordinates:
      const std::vector<unsigned int> rep_vec(repetitions.begin(), repetitions.end());
      GridGenerator::subdivided_hyper_rectangle (coarse_grid,
                                                 rep_vec,
                                                 point1,
                                                 point2,
                                                 true);

      // Then transform this box into a spherical chunk (possibly with
      // topography -- the 'manifold' has that built in):
      GridTools::transform (
        [&](const Point<dim> &p) -> Point<dim>
      {
        return manifold->push_forward(p);
      },
      coarse_grid);

      // Deal with a curved mesh by assigning a manifold. We arbitrarily
      // choose manifold_id 15 for this.
      coarse_grid.set_manifold (my_manifold_id, *manifold);
      for (const auto &cell : coarse_grid.active_cell_iterators())
        cell->set_all_manifold_ids (my_manifold_id);
    }



    template <int dim>
    std::set<types::boundary_id>
    Chunk<dim>::
    get_used_boundary_indicators () const
    {
      // boundary indicators are zero through 2*dim-1
      std::set<types::boundary_id> s;
      for (unsigned int i=0; i<2*dim; ++i)
        s.insert (i);
      return s;
    }



    template <int dim>
    std::map<std::string,types::boundary_id>
    Chunk<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            return
            {
              {"bottom", 0},
              {"top",    1},
              {"west",   2},
              {"east",   3}
            };
          }

          case 3:
          {
            return
            {
              {"bottom", 0},
              {"top",    1},
              {"west",   2},
              {"east",   3},
              {"south",  4},
              {"north",  5}
            };
          }
        }

      Assert (false, ExcNotImplemented());
      return {};
    }



    template <int dim>
    double
    Chunk<dim>::
    length_scale () const
    {
      // As described in the first ASPECT paper, a length scale of
      // 10km = 1e4m works well for the pressure scaling for earth
      // sized spherical shells. use a length scale that
      // yields this value for the R0,R1 corresponding to earth
      // but otherwise scales like (R1-R0)
      return 1e4 * maximal_depth() / (6336000.-3481000.);
    }



    template <int dim>
    double
    Chunk<dim>::depth(const Point<dim> &position) const
    {
      // depth is defined wrt the reference surface point2[0]
      // negative depth is not allowed
      return std::max (0., std::min (point2[0]-position.norm(), maximal_depth()));
    }



    template <int dim>
    double
    Chunk<dim>::height_above_reference_surface(const Point<dim> &position) const
    {
      return position.norm()-point2[0];
    }



    template <int dim>
    Point<dim>
    Chunk<dim>::representative_point(const double depth) const
    {
      Assert (depth >= 0,
              ExcMessage ("Given depth must be positive or zero."));
      Assert (depth <= maximal_depth(),
              ExcMessage ("Given depth must be less than or equal to the maximal depth of this geometry."));

      // Choose a point at the mean longitude (and latitude)
      Point<dim> p = 0.5*(point2+point1);
      // at a depth beneath the top surface
      p[0] = point2[0]-depth;

      // Now convert to Cartesian coordinates. This ignores the surface topography,
      // but that is as documented.
      return manifold->push_forward_sphere(p);
    }



    template <int dim>
    double
    Chunk<dim>::west_longitude () const
    {
      return point1[1];
    }



    template <int dim>
    double
    Chunk<dim>::east_longitude () const
    {
      return point2[1];
    }



    template <int dim>
    double
    Chunk<dim>::longitude_range () const
    {
      return point2[1] - point1[1];
    }



    template <int dim>
    double
    Chunk<dim>::south_latitude () const
    {
      if (dim == 3)
        return point1[2];
      else
        return 0;
    }



    template <int dim>
    double
    Chunk<dim>::north_latitude () const
    {
      if (dim==3)
        return point2[2];
      else
        return 0;
    }



    template <int dim>
    double
    Chunk<dim>::latitude_range () const
    {
      if (dim==3)
        return point2[2] - point1[2];
      else
        return 0;
    }



    template <int dim>
    double
    Chunk<dim>::maximal_depth() const
    {
      // The depth is defined as relative to a reference surface (without
      // topography) and since we don't apply topography on the CMB,
      // the maximal depth really is the formula below, unless one applies a
      // topography that is always strictly below zero (i.e., where the
      // actual surface lies strictly below the reference surface).
      return point2[0]-point1[0];
    }



    template <int dim>
    double
    Chunk<dim>::inner_radius () const
    {
      return point1[0];
    }



    template <int dim>
    double
    Chunk<dim>::outer_radius () const
    {
      return point2[0];
    }



    template <int dim>
    bool
    Chunk<dim>::has_curved_elements() const
    {
      return true;
    }



    template <int dim>
    bool
    Chunk<dim>::point_is_in_domain(const Point<dim> &point) const
    {
      // If mesh deformation is enabled, we have to loop over all the current
      // grid cells to see if the given point lies in the domain.
      // If mesh deformation is not enabled, or has not happened yet,
      // we can use the global extents of the model domain with or without
      // initial topography instead.
      // This function only checks if the given point lies in the domain
      // in its current shape at the current time. It can be called before
      // mesh deformation is applied in the first timestep (e.g., by the boundary
      // traction plugins), and therefore there is no guarantee
      // that the point will still lie in the domain after initial mesh deformation.
      if (this->get_parameters().mesh_deformation_enabled &&
          this->simulator_is_past_initialization())
        {
          return Utilities::point_is_in_triangulation(this->get_mapping(), this->get_triangulation(), point, this->get_mpi_communicator());
        }
      // Without mesh deformation enabled, it is much cheaper to check whether the point lies in the domain.
      else
        {
          const Point<dim> spherical_point = manifold->pull_back(point);

          for (unsigned int d = 0; d < dim; ++d)
            if (spherical_point[d] > point2[d]+std::numeric_limits<double>::epsilon()*std::abs(point2[d]) ||
                spherical_point[d] < point1[d]-std::numeric_limits<double>::epsilon()*std::abs(point2[d]))
              return false;

          return true;
        }
    }



    template <int dim>
    std::array<double,dim>
    Chunk<dim>::cartesian_to_natural_coordinates(const Point<dim> &position_point) const
    {
      // The chunk manifold uses (radius, longitude, latitude).
      // This is exactly what we need.

      // Ignore the topography to avoid a loop when calling the
      // AsciiDataBoundary for topography which uses this function....
      const Point<dim> transformed_point = manifold->pull_back_sphere(position_point);
      std::array<double,dim> position_array;
      for (unsigned int i = 0; i < dim; ++i)
        position_array[i] = transformed_point(i);

      // Internally, the Chunk geometry uses longitudes in the range -pi...pi,
      // and that is what pull_back_sphere() produces. But externally, we need
      // to use 0...2*pi to match what Utilities::Coordinates::cartesian_to_spherical_coordinates()
      // returns, for example, and what we document the AsciiBoundaryData class
      // wants to see from input files.
      if (position_array[1] < 0)
        position_array[1] += 2*numbers::PI;

      return position_array;
    }



    template <int dim>
    aspect::Utilities::Coordinates::CoordinateSystem
    Chunk<dim>::natural_coordinate_system() const
    {
      // TODO This will give problems somewhere down the line
      // if geometry models are asked for their coordinate system,
      // chunk returns spherical and then Utilitiess::Coordinates::cartesian_to_spherical
      // is used
      return aspect::Utilities::Coordinates::CoordinateSystem::spherical;
    }



    template <int dim>
    Point<dim>
    Chunk<dim>::natural_to_cartesian_coordinates(const std::array<double,dim> &position_tensor) const
    {
      // Ignore the topography to avoid a loop when calling the
      // AsciiDataBoundary for topography which uses this function....
      Point<dim> position_point;
      for (unsigned int i = 0; i < dim; ++i)
        position_point[i] = position_tensor[i];
      const Point<dim> transformed_point = manifold->push_forward_sphere(position_point);

      return transformed_point;
    }



    template <int dim>
    void
    Chunk<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Chunk");
        {
          prm.declare_entry ("Chunk inner radius", "0.",
                             Patterns::Double (0.),
                             "Radius at the bottom surface of the chunk. Units: \\si{\\meter}.");
          prm.declare_entry ("Chunk outer radius", "1.",
                             Patterns::Double (0.),
                             "Radius at the top surface of the chunk. Units: \\si{\\meter}.");

          prm.declare_entry ("Chunk minimum longitude", "0.",
                             Patterns::Double (-180., 360.), // enables crossing of either hemisphere
                             "Minimum longitude of the chunk. Units: degrees.");
          prm.declare_entry ("Chunk maximum longitude", "1.",
                             Patterns::Double (-180., 360.), // enables crossing of either hemisphere
                             "Maximum longitude of the chunk. Units: degrees.");

          prm.declare_entry ("Chunk minimum latitude", "0.",
                             Patterns::Double (-90., 90.),
                             "Minimum latitude of the chunk. This value is ignored "
                             "if the simulation is in 2d. Units: degrees.");
          prm.declare_entry ("Chunk maximum latitude", "1.",
                             Patterns::Double (-90., 90.),
                             "Maximum latitude of the chunk. This value is ignored "
                             "if the simulation is in 2d. Units: degrees.");

          prm.declare_entry ("Radius repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in radius.");
          prm.declare_entry ("Longitude repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in longitude.");
          prm.declare_entry ("Latitude repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in latitude. This value is ignored "
                             "if the simulation is in 2d");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Chunk<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Chunk");
        {
          point1[0] = prm.get_double ("Chunk inner radius");
          point2[0] = prm.get_double ("Chunk outer radius");
          repetitions[0] = prm.get_integer ("Radius repetitions");
          point1[1] = prm.get_double ("Chunk minimum longitude") * constants::degree_to_radians;
          point2[1] = prm.get_double ("Chunk maximum longitude") * constants::degree_to_radians;
          repetitions[1] = prm.get_integer ("Longitude repetitions");

          AssertThrow (point1[0] < point2[0],
                       ExcMessage ("Inner radius must be less than outer radius."));
          AssertThrow (point1[1] < point2[1],
                       ExcMessage ("Minimum longitude must be less than maximum longitude."));
          AssertThrow (point2[1] - point1[1] < 2.*numbers::PI,
                       ExcMessage ("Maximum - minimum longitude should be less than 360 degrees."));

          if (dim == 3)
            {
              point1[2] = prm.get_double ("Chunk minimum latitude") * constants::degree_to_radians;
              point2[2] = prm.get_double ("Chunk maximum latitude") * constants::degree_to_radians;
              repetitions[2] = prm.get_integer ("Latitude repetitions");

              AssertThrow (point1[2] < point2[2],
                           ExcMessage ("Minimum latitude must be less than maximum latitude."));
              AssertThrow (point1[2] > -0.5*numbers::PI,
                           ExcMessage ("Minimum latitude needs to be larger than -90 degrees."));
              AssertThrow (point2[2] < 0.5*numbers::PI,
                           ExcMessage ("Maximum latitude needs to be less than 90 degrees."));
            }

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
    namespace internal
    {
#define INSTANTIATE(dim) \
  template class ChunkGeometry<dim>;
      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }

    ASPECT_REGISTER_GEOMETRY_MODEL(Chunk,
                                   "chunk",
                                   "A geometry which can be described as a chunk of a spherical shell, "
                                   "bounded by lines of longitude, latitude and radius. "
                                   "The minimum and maximum longitude, latitude (if in 3d) and depth of the chunk "
                                   "is set in the parameter file. The chunk geometry labels its "
                                   "2*dim sides as follows: ``west'' and ``east'': minimum and maximum "
                                   "longitude, ``south'' and ``north'': minimum and maximum latitude, "
                                   "``inner'' and ``outer'': minimum and maximum radii. "
                                   "\n\n"
                                   "The dimensions of the model are specified by parameters "
                                   "of the following form: "
                                   "Chunk (minimum || maximum) (longitude || latitude): "
                                   "edges of geographical quadrangle (in degrees)"
                                   "Chunk (inner || outer) radius: Radii at bottom and top of chunk"
                                   "(Longitude || Latitude || Radius) repetitions: "
                                   "number of cells in each coordinate direction."
                                   "\n\n"
                                   "When used in 2d, this geometry does not imply the use of "
                                   "a spherical coordinate system. Indeed, "
                                   "in 2d the geometry is simply a sector of an annulus in a Cartesian "
                                   "coordinate system and consequently would correspond to "
                                   "a sector of a cross section of the fluid filled space between two "
                                   "infinite cylinders where one has made the assumption that "
                                   "the velocity in direction of the cylinder axes is zero. "
                                   "This is consistent with the definition of what we consider "
                                   "the two-dimension case given in "
                                   "Section~\\ref{sec:methods:2d-models}. "
                                   "It is also possible to add initial topography to the chunk geometry, "
                                   "based on an ascii data file. ")
  }
}
