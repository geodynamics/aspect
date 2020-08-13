/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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


#include <aspect/boundary_traction/initial_lithostatic_pressure.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <array>

#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/two_merged_boxes.h>

namespace aspect
{
  namespace BoundaryTraction
  {

    template <int dim>
    void
    InitialLithostaticPressure<dim>::initialize()
    {
      // Ensure the initial lithostatic pressure traction boundary conditions are used,
      // and register for which boundary indicators these conditions are set.
      std::set<types::boundary_id> traction_bi;
      for (const auto &p : this->get_boundary_traction())
        {
          if (p.second.get() == this)
            traction_bi.insert(p.first);
        }
      AssertThrow(*(traction_bi.begin()) != numbers::invalid_boundary_id,
                  ExcMessage("Did not find any boundary indicators for the initial lithostatic pressure plugin."));

      // The below is adapted from adiabatic_conditions/initial_profile.cc,
      // but we use the initial temperature and composition and only calculate
      // a pressure profile with depth.

      // The spacing of the depth profile
      delta_z = this->get_geometry_model().maximal_depth() / (n_points-1);

      // The number of compositional fields
      const unsigned int n_compositional_fields = this->n_compositional_fields();

      // The pressure at the surface
      pressure[0]    = this->get_surface_pressure();

      // For spherical(-like) domains, modify the representative point:
      // go from degrees to radians...
      const double degrees_to_radians = dealii::numbers::PI/180.0;
      std::array<double, dim> spherical_representative_point;
      for (unsigned int d=0; d<dim; d++)
        spherical_representative_point[d] = representative_point[d];
      spherical_representative_point[1] *= degrees_to_radians;
      // and go from latitude to colatitude.
      if (dim == 3)
        {
          spherical_representative_point[2] = 90.0 - spherical_representative_point[2];
          spherical_representative_point[2] *= degrees_to_radians;
        }

      // Check that the representative point lies in the domain.
      AssertThrow(((this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::CoordinateSystem::cartesian) ?
                   (this->get_geometry_model().point_is_in_domain(representative_point)) :
                   (this->get_geometry_model().point_is_in_domain(Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(spherical_representative_point)))),
                  ExcMessage("The reference point does not lie with the domain."));

      // Set the radius of the representative point to the surface radius for spherical domains
      // or set the vertical coordinate to the surface value for box domains.
      if (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model()))
        spherical_representative_point[0] = Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model()).outer_radius();
      else if (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>> (this->get_geometry_model()))
        spherical_representative_point[0] =  Plugins::get_plugin_as_type<const GeometryModel::Chunk<dim>>(this->get_geometry_model()).outer_radius();
      else if (Plugins::plugin_type_matches<const GeometryModel::EllipsoidalChunk<dim>> (this->get_geometry_model()))
        {
          const GeometryModel::EllipsoidalChunk<dim> &gm = Plugins::get_plugin_as_type<const GeometryModel::EllipsoidalChunk<dim>> (this->get_geometry_model());
          // TODO
          // If the eccentricity of the EllipsoidalChunk is non-zero, the radius can vary along a boundary,
          // but the maximal depth is the same everywhere and we could calculate a representative pressure
          // profile. However, it requires some extra logic with ellipsoidal
          // coordinates, so for now we only allow eccentricity zero.
          // Using the EllipsoidalChunk with eccentricity zero can still be useful,
          // because the domain can be non-coordinate parallel.
          AssertThrow(gm.get_eccentricity() == 0.0, ExcMessage("This initial lithospheric pressure plugin cannot be used with a non-zero eccentricity. "));

          spherical_representative_point[0] = gm.get_semi_major_axis_a();
        }
      else if (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>> (this->get_geometry_model()))
        spherical_representative_point[0] =  Plugins::get_plugin_as_type<const GeometryModel::Sphere<dim>>(this->get_geometry_model()).radius();
      else if (Plugins::plugin_type_matches<const GeometryModel::Box<dim>> (this->get_geometry_model()))
        representative_point[dim-1]=  Plugins::get_plugin_as_type<const GeometryModel::Box<dim>>(this->get_geometry_model()).get_extents()[dim-1];
      else if (Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model()))
        representative_point[dim-1]=  Plugins::get_plugin_as_type<const GeometryModel::TwoMergedBoxes<dim>>(this->get_geometry_model()).get_extents()[dim-1];
      else
        AssertThrow(false, ExcNotImplemented());

      // Set up the input for the density function of the material model.
      typename MaterialModel::Interface<dim>::MaterialModelInputs in0(1, n_compositional_fields);
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out0(1, n_compositional_fields);

      // Where to calculate the density
      // for cartesian domains
      if (Plugins::plugin_type_matches<const GeometryModel::Box<dim>> (this->get_geometry_model()) ||
          Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model()))
        in0.position[0] = representative_point;
      // and for spherical domains
      else
        in0.position[0] = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(spherical_representative_point);

      // We need the initial temperature at this point
      in0.temperature[0] = this->get_initial_temperature_manager().initial_temperature(in0.position[0]);

      // and the surface pressure.
      in0.pressure[0] = pressure[0];

      // Then the compositions at this point.
      for (unsigned int c=0; c<n_compositional_fields; ++c)
        in0.composition[0][c] = this->get_initial_composition_manager().initial_composition(in0.position[0], c);

      // We do not need the viscosity.
      in0.strain_rate.resize(0);

      // Evaluate the material model to get the density.
      this->get_material_model().evaluate(in0, out0);
      const double density0 = out0.densities[0];

      // Get the magnitude of gravity. We assume
      // that gravity always points along the depth direction. This
      // may not strictly be true always but is likely a good enough
      // approximation here.
      const double gravity0 = this->get_gravity_model().gravity_vector(in0.position[0]).norm();

      // Now integrate pressure downward using trapezoidal integration
      // p'(z) = rho(p,c,T) * |g| * delta_z
      double sum = delta_z * 0.5 * density0 * gravity0;

      for (unsigned int i=1; i<n_points; ++i)
        {
          AssertThrow (i < pressure.size(), ExcMessage(std::string("The current index ")
                                                       + dealii::Utilities::int_to_string(i)
                                                       + std::string(" is bigger than the size of the pressure vector ")
                                                       + dealii::Utilities::int_to_string(pressure.size())));

          // Set up the input for the density function of the material model
          typename MaterialModel::Interface<dim>::MaterialModelInputs in(1, n_compositional_fields);
          typename MaterialModel::Interface<dim>::MaterialModelOutputs out(1, n_compositional_fields);

          // Where to calculate the density:
          // for cartesian domains
          if (Plugins::plugin_type_matches<const GeometryModel::Box<dim>> (this->get_geometry_model()) ||
              Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model()))
            {
              // decrease z coordinate with depth increment
              representative_point[dim-1] -= delta_z;
              in.position[0] = representative_point;
            }
          // and for spherical domains
          else
            {
              // decrease radius with depth increment
              spherical_representative_point[0] -= delta_z;
              in.position[0] = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(spherical_representative_point);
            }

          // Retrieve the initial temperature at this point.
          in.temperature[0] = this->get_initial_temperature_manager().initial_temperature(in.position[0]);
          // and use the previous pressure
          in.pressure[0] = pressure[i-1];

          // Retrieve the compositions at this point.
          for (unsigned int c=0; c<n_compositional_fields; ++c)
            in.composition[0][c] = this->get_initial_composition_manager().initial_composition(in.position[0], c);

          // We do not need the viscosity.
          in.strain_rate.resize(0);

          // Evaluate the material model to get the density at the current point.
          this->get_material_model().evaluate(in, out);
          const double density = out.densities[0];

          // Get the magnitude of gravity.
          const double gravity = this->get_gravity_model().gravity_vector(in.position[0]).norm();

          // Trapezoid integration
          pressure[i] = sum + delta_z * 0.5 * density * gravity;
          sum += delta_z * density * gravity;
        }

      Assert (*std::min_element (pressure.begin(), pressure.end()) >=
              -std::numeric_limits<double>::epsilon() * pressure.size(),
              ExcInternalError());

    }

    template <int dim>
    Tensor<1,dim>
    InitialLithostaticPressure<dim>::
    traction (const Point<dim> &p,
              const Tensor<1,dim> &normal) const
    {
      // We want to set the normal component to the vertical boundary
      // to the lithostatic pressure, the rest of the traction
      // components are left set to zero. We get the lithostatic pressure
      // from a linear interpolation of the calculated profile.
      Tensor<1,dim> traction;
      traction = -interpolate_pressure(p) * normal;

      return traction;
    }

    template <int dim>
    double
    InitialLithostaticPressure<dim>::
    interpolate_pressure (const Point<dim> &p) const
    {
      // The depth at which we need the pressure.
      const double z = this->get_geometry_model().depth(p);

      // Check that depth does not exceed the maximal depth.
      if (z >= this->get_geometry_model().maximal_depth())
        {
          Assert (z <= this->get_geometry_model().maximal_depth() + delta_z,
                  ExcInternalError());
          // Return deepest (last) pressure
          return pressure.back();
        }

      const unsigned int i = static_cast<unsigned int>(z/delta_z);
      Assert ((z/delta_z) >= 0, ExcInternalError());
      Assert (i+1 < pressure.size(), ExcInternalError());

      // Now do the linear interpolation.
      const double d=1.0+i-z/delta_z;
      Assert ((d>=0) && (d<=1), ExcInternalError());

      return d*pressure[i]+(1-d)*pressure[i+1];
    }


    template <int dim>
    void
    InitialLithostaticPressure<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary traction model");
      {
        prm.enter_subsection("Initial lithostatic pressure");
        {
          prm.declare_entry ("Representative point", "",
                             Patterns::List(Patterns::Double()),
                             "The point where the pressure profile will be calculated. "
                             "Cartesian coordinates $(x,y,z)$ when geometry is a box, otherwise enter radius, "
                             "longitude, and in 3D latitude. Note that the coordinate related to the depth "
                             "($y$ in 2D cartesian, $z$ in 3D cartesian and radius in spherical coordinates) is "
                             "not used. "
                             "Units: \\si{\\meter} or degrees.");
          prm.declare_entry("Number of integration points", "1000",
                            Patterns::Integer(0),
                            "The number of integration points over which we integrate the lithostatic pressure "
                            "downwards.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    InitialLithostaticPressure<dim>::parse_parameters (ParameterHandler &prm)
    {
      unsigned int refinement;
      prm.enter_subsection("Mesh refinement");
      {
        refinement = prm.get_integer("Initial adaptive refinement") + prm.get_integer("Initial global refinement");
      }
      prm.leave_subsection();

      prm.enter_subsection("Boundary traction model");
      {
        prm.enter_subsection("Initial lithostatic pressure");
        {
          // The representative point where to calculate the depth profile.
          const std::vector<double> rep_point =
            dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Representative point")));
          AssertThrow(rep_point.size() == dim, ExcMessage("Representative point does not have the right dimensions."));
          for (unsigned int d = 0; d<dim; d++)
            representative_point[d] = rep_point[d];
          n_points = prm.get_integer("Number of integration points");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Check that we have enough integration points for this mesh.
      AssertThrow(std::pow(2.0,refinement) <= n_points, ExcMessage("Not enough integration points for this resolution."));

      pressure.resize(n_points,-1);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTraction
  {
    ASPECT_REGISTER_BOUNDARY_TRACTION_MODEL(InitialLithostaticPressure,
                                            "initial lithostatic pressure",
                                            "Implementation of a model in which the boundary "
                                            "traction is given in terms of a normal traction component "
                                            "set to the lithostatic pressure "
                                            "calculated according to the parameters in section "
                                            "``Boundary traction model|Lithostatic pressure''. "
                                            "\n\n"
                                            "The lithostatic pressure is calculated by integrating "
                                            "the pressure downward based on the initial composition "
                                            "and temperature along the user-specified depth profile. "
                                            "The user-specified profile is given in terms of a point "
                                            "in cartesian coordinates for box geometries and "
                                            "in spherical coordinates for all "
                                            "other geometries (radius, longitude, latitude), and "
                                            "the number of integration points. "
                                            "The lateral coordinates of the point are used to calculate "
                                            "the lithostatic pressure profile with depth. This means that "
                                            "the depth coordinate is not used."
                                            "\n\n"
                                            "Gravity is expected to point along the depth direction. ")
  }
}
