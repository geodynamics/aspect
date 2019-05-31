/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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
    template <int dim>
    Chunk<dim>::ChunkGeometry::ChunkGeometry()
      :
      point1_lon(0.0),
      inner_radius(3471e3),
      max_depth(2900e3)
    {}



    template <int dim>
    Chunk<dim>::ChunkGeometry::ChunkGeometry(const ChunkGeometry &other)
      :
      ChartManifold<dim,dim>(other),
      point1_lon(other.point1_lon),
      inner_radius(other.inner_radius),
      max_depth(other.max_depth)
    {
      this->initialize(other.topo);
    }



    template <int dim>
    void
    Chunk<dim>::ChunkGeometry::initialize(const InitialTopographyModel::Interface<dim> *topo_pointer)
    {
      topo = topo_pointer;
    }



    template <int dim>
    DerivativeForm<1,dim,dim>
    Chunk<dim>::ChunkGeometry::
    push_forward_gradient(const Point<dim> &chart_point) const
    {
      const double R = chart_point[0]; // Radius

      Assert (R > 0.0, ExcMessage("Negative radius for given point."));

      Tensor<2,dim> DX;

      // The initial topography derivatives (dtopo/dphi and dtopo/dtheta)
      // They are zero for the ZeroTopography model
      Tensor<1,dim-1> topo_derivatives;
      if (const InitialTopographyModel::AsciiData<dim> *itm = dynamic_cast<const InitialTopographyModel::AsciiData<dim> *> (topo))
        topo_derivatives = itm->vector_gradient(push_forward_sphere(chart_point));

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
    Chunk<dim>::ChunkGeometry::
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
    Chunk<dim>::ChunkGeometry::
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
    Chunk<dim>::ChunkGeometry::
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
    Point<dim>
    Chunk<dim>::ChunkGeometry::
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
            // See 2D case
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
    void
    Chunk<dim>::ChunkGeometry::
    set_min_longitude(const double p1_lon)
    {
      point1_lon = p1_lon;
    }



    template <int dim>
    std::unique_ptr<Manifold<dim,dim> >
    Chunk<dim>::ChunkGeometry::
    clone() const
    {
      return std_cxx14::make_unique<ChunkGeometry>(*this);
    }



    template <int dim>
    Point<dim>
    Chunk<dim>::ChunkGeometry::
    push_forward_topo(const Point<dim> &r_phi_theta) const
    {
      // the radius of the current point without topography
      const double radius = r_phi_theta[0];

      // Grab lon,lat coordinates
      Point<dim-1> surface_point;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_point[d] = r_phi_theta[d+1];
      // Convert latitude to colatitude
      if (dim == 3)
        surface_point[1] = 0.5*numbers::PI - surface_point[1];
      const double topography = topo->value(surface_point);

      // adjust the radius based on the radius of the point
      // through a linear interpolation between 0 at max depth and
      // "topography" at the surface
      const double topo_radius = std::max(inner_radius,radius + (radius-inner_radius)/max_depth*topography);

      // return the point with adjusted radius
      Point<dim> topor_phi_theta = r_phi_theta;
      topor_phi_theta[0] = topo_radius;

      return topor_phi_theta;
    }



    template <int dim>
    Point<dim>
    Chunk<dim>::ChunkGeometry::
    pull_back_topo(const Point<dim> &topor_phi_theta) const
    {
      // the radius of the point with topography
      const double topo_radius = topor_phi_theta[0];

      // Grab lon,lat coordinates
      Point<dim-1> surface_point;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_point[d] = topor_phi_theta[d+1];
      // Convert latitude to colatitude
      if (dim == 3)
        surface_point[1] = 0.5*numbers::PI - surface_point[1];
      const double topography = topo->value(surface_point);

      // remove the topography (which scales with radius)
      const double radius = std::max(inner_radius,(topo_radius*max_depth+inner_radius*topography)/(max_depth+topography));

      // return the point without topography
      Point<dim> r_phi_theta = topor_phi_theta;
      r_phi_theta[0] = radius;
      return r_phi_theta;
    }



    template <int dim>
    double
    Chunk<dim>::ChunkGeometry::
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



    template <int dim>
    void
    Chunk<dim>::ChunkGeometry::
    set_min_radius(const double p1_rad)
    {
      inner_radius = p1_rad;
    }



    template <int dim>
    void
    Chunk<dim>::ChunkGeometry::
    set_max_depth(const double p2_p1_rad)
    {
      max_depth = p2_p1_rad;
    }



    template <int dim>
    void
    Chunk<dim>::initialize ()
    {
      AssertThrow(dynamic_cast<const InitialTopographyModel::ZeroTopography<dim>*>(&this->get_initial_topography_model()) != nullptr ||
                  dynamic_cast<const InitialTopographyModel::AsciiData<dim>*>(&this->get_initial_topography_model()) != nullptr,
                  ExcMessage("At the moment, only the Zero or AsciiData initial topography model can be used."));

      manifold.initialize(&(this->get_initial_topography_model()));
    }


    template <int dim>
    void
    Chunk<dim>::initialize_for_test (const InitialTopographyModel::Interface<dim> *topo_pointer)
    {
      manifold.initialize(topo_pointer);
    }



    template <int dim>
    void
    Chunk<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      std::vector<unsigned int> rep_vec(repetitions, repetitions+dim);
      GridGenerator::subdivided_hyper_rectangle (coarse_grid,
                                                 rep_vec,
                                                 point1,
                                                 point2,
                                                 true);

      // Transform box into spherical chunk
      GridTools::transform (
        [&](const Point<dim> &p) -> Point<dim>
      {
        return manifold.push_forward(p);
      },
      coarse_grid);

      // Deal with a curved mesh
      // Attach the real manifold to slot 15.
      coarse_grid.set_manifold (15, manifold);
      for (const auto &cell : coarse_grid.active_cell_iterators())
        cell->set_all_manifold_ids (15);
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
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("bottom", 0),
                  std::pair<std::string,types::boundary_id>("top",    1),
                  std::pair<std::string,types::boundary_id>("west",   2),
                  std::pair<std::string,types::boundary_id>("east",   3)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }

          case 3:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("bottom", 0),
                  std::pair<std::string,types::boundary_id>("top",    1),
                  std::pair<std::string,types::boundary_id>("west",   2),
                  std::pair<std::string,types::boundary_id>("east",   3),
                  std::pair<std::string,types::boundary_id>("south",  4),
                  std::pair<std::string,types::boundary_id>("north",  5)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }
        }

      Assert (false, ExcNotImplemented());
      return std::map<std::string,types::boundary_id>();
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
    Chunk<dim>::depth_wrt_topo(const Point<dim> &position) const
    {
      // depth is defined wrt the reference surface point2[0] + the topography
      // depth is therefore always positive
      const double outer_radius = manifold.get_radius(position);
      const Point<dim> rtopo_phi_theta = manifold.pull_back_sphere(position);
      Assert (rtopo_phi_theta[0] <= outer_radius, ExcMessage("The radius is bigger than the maximum radius."));
      return std::max(0.0, outer_radius - rtopo_phi_theta[0]);
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

      // Now convert to Cartesian coordinates
      return manifold.push_forward_sphere(p);
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
      AssertThrow(!this->get_parameters().mesh_deformation_enabled ||
                  // we are still before the first time step has started
                  this->get_timestep_number() == 0 ||
                  this->get_timestep_number() == numbers::invalid_unsigned_int,
                  ExcMessage("After displacement of the mesh, this function can no longer be used to determine whether a point lies in the domain or not."));

      AssertThrow(dynamic_cast<const InitialTopographyModel::ZeroTopography<dim>*>(&this->get_initial_topography_model()) != nullptr,
                  ExcMessage("After adding topography, this function can no longer be used to determine whether a point lies in the domain or not."));

      const Point<dim> spherical_point = manifold.pull_back(point);

      for (unsigned int d = 0; d < dim; ++d)
        if (spherical_point[d] > point2[d]+std::numeric_limits<double>::epsilon()*std::abs(point2[d]) ||
            spherical_point[d] < point1[d]-std::numeric_limits<double>::epsilon()*std::abs(point2[d]))
          return false;

      return true;
    }



    template <int dim>
    std::array<double,dim>
    Chunk<dim>::cartesian_to_natural_coordinates(const Point<dim> &position_point) const
    {
      // the chunk manifold has a order of radius, longitude, latitude.
      // This is exactly what we need.
      // Ignore the topography to avoid a loop when calling the
      // AsciiDataBoundary for topography which uses this function....
      const Point<dim> transformed_point = manifold.pull_back_sphere(position_point);
      std::array<double,dim> position_array;
      for (unsigned int i = 0; i < dim; i++)
        position_array[i] = transformed_point(i);

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
      for (unsigned int i = 0; i < dim; i++)
        position_point[i] = position_tensor[i];
      const Point<dim> transformed_point = manifold.push_forward_sphere(position_point);

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
          prm.declare_entry ("Chunk inner radius", "0",
                             Patterns::Double (0),
                             "Radius at the bottom surface of the chunk. Units: $\\si{m}$.");
          prm.declare_entry ("Chunk outer radius", "1",
                             Patterns::Double (0),
                             "Radius at the top surface of the chunk. Units: $\\si{m}$.");

          prm.declare_entry ("Chunk minimum longitude", "0",
                             Patterns::Double (-180, 360), // enables crossing of either hemisphere
                             "Minimum longitude of the chunk. Units: degrees.");
          prm.declare_entry ("Chunk maximum longitude", "1",
                             Patterns::Double (-180, 360), // enables crossing of either hemisphere
                             "Maximum longitude of the chunk. Units: degrees.");

          prm.declare_entry ("Chunk minimum latitude", "0",
                             Patterns::Double (-90, 90),
                             "Minimum latitude of the chunk. This value is ignored "
                             "if the simulation is in 2d. Units: degrees.");
          prm.declare_entry ("Chunk maximum latitude", "1",
                             Patterns::Double (-90, 90),
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

          const double degtorad = dealii::numbers::PI/180;

          Assert (dim >= 2, ExcInternalError());
          Assert (dim <= 3, ExcInternalError());

          if (dim >= 2)
            {
              point1[0] = prm.get_double ("Chunk inner radius");
              point2[0] = prm.get_double ("Chunk outer radius");
              repetitions[0] = prm.get_integer ("Radius repetitions");
              point1[1] = prm.get_double ("Chunk minimum longitude") * degtorad;
              point2[1] = prm.get_double ("Chunk maximum longitude") * degtorad;
              repetitions[1] = prm.get_integer ("Longitude repetitions");

              AssertThrow (point1[0] < point2[0],
                           ExcMessage ("Inner radius must be less than outer radius."));
              AssertThrow (point1[1] < point2[1],
                           ExcMessage ("Minimum longitude must be less than maximum longitude."));
              AssertThrow (point2[1] - point1[1] < 2.*numbers::PI,
                           ExcMessage ("Maximum - minimum longitude should be less than 360 degrees."));
            }

          // Inform the manifold about the minimum longitude
          manifold.set_min_longitude(point1[1]);
          // Inform the manifold about the minimum radius
          manifold.set_min_radius(point1[0]);
          // Inform the manifold about the maximum depth (without topo)
          manifold.set_max_depth(point2[0]-point1[0]);

          if (dim == 3)
            {
              point1[2] = prm.get_double ("Chunk minimum latitude") * degtorad;
              point2[2] = prm.get_double ("Chunk maximum latitude") * degtorad;
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
                                   "Section~\\ref{sec:meaning-of-2d}. "
                                   "It is also possible to add initial topography to the chunk geometry, "
                                   "based on an ascii data file. ")
  }
}

