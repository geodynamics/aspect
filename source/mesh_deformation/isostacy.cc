/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE. If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/global.h>

#include <aspect/mesh_deformation/isostacy.h>

#include <aspect/geometry_model/box.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/material_model/interface.h>

#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_composition/interface.h>

#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    Isostacy<dim>::Isostacy ()
    {}



    template <int dim>
    void
    Isostacy<dim>::initialize ()
    {
      AssertThrow(Plugins::plugin_type_matches<GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("The isostatic topography model currently only supports the box geometry model."));

      const auto &geometry_model = dynamic_cast<const GeometryModel::Box<dim> &>(this->get_geometry_model());

      AssertThrow(balanced_depth >= 0.0 &&
                  balanced_depth <= geometry_model.maximal_depth(),
                  ExcMessage("The balanced depth must be between zero and the maximal model depth."));

      const double x_min = 0.0;
      const double x_max = geometry_model.get_extents()[0];

      const double y_top = geometry_model.get_extents()[1];


      const double grid_spacing = (x_max - x_min) / (n_lateral_points - 1.0);
      const double delta_z = this->get_geometry_model().maximal_depth() / (n_vertical_points-1);

      // Determine the reference density at the balanced depth beneath
      // the center of the model domain.
      double ref_density = std::numeric_limits<double>::lowest();
      const unsigned int ref_column_index = static_cast<unsigned int>(std::round((x_max - x_min) / 2.0 / grid_spacing));
      const unsigned int ref_depth_index = static_cast<unsigned int>(std::round(balanced_depth/delta_z));

      topography.resize(n_lateral_points);
      // Compute the mass of each vertical column and accumulate the
      // total column mass.
      std::vector<double> column_masses(n_lateral_points);
      double total_column_mass = 0.0;

      for (unsigned int i = 0; i < n_lateral_points; ++i)
        {
          const double x = x_min + i * grid_spacing;
          double column_mass = 0.0;

          // Compute pressure by vertically integrating the density,
          // following the same approach used by the adiabatic conditions.
          double pressure;
          double density;

          MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());
          MaterialModel::MaterialModelOutputs<dim> out(1, this->n_compositional_fields());

          in.requested_properties = MaterialModel::MaterialProperties::equation_of_state_properties;
          in.velocity[0] = Tensor <1,dim> ();

          const Tensor <1,dim> g = this->get_gravity_model().gravity_vector(this->get_geometry_model().representative_point(0));
          const Point<dim> point_surf = this->get_geometry_model().representative_point(0);
          const Point<dim> point_bot = this->get_geometry_model().representative_point(this->get_geometry_model().maximal_depth());
          const int gravity_direction =  (g * (point_bot - point_surf) >= 0) ?
                                         1 :
                                         -1;

          for (unsigned int j=0; j<n_vertical_points; ++j)
            {
              Point<dim> position;
              position[0] = x;
              position[1] = y_top - j * delta_z;

              if (j==0)
                {
                  // Use the surface pressure of the adiabatic condition
                  pressure = this->get_adiabatic_conditions().pressure(position);
                }
              else
                {
                  // Update the pressure using the density evaluated at the
                  // previous vertical sampling point.
                  const double gravity = gravity_direction * this->get_gravity_model().gravity_vector(position).norm();

                  pressure = pressure + density * gravity * delta_z;
                }

              // Use the initial temperature and composition at the local position
              std::vector<double> composition(this->n_compositional_fields());

              for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
                composition[c] =
                  this->get_initial_composition_manager().initial_composition(position, c);

              in.position[0] = position;
              in.temperature[0] =  this->get_initial_temperature_manager().initial_temperature(position);
              in.pressure[0] = pressure;
              in.composition[0] = composition;

              this->get_material_model().evaluate(in, out);
              density = out.densities[0];

              // Store the reference density at the balanced depth beneath
              // the center of the model domain.
              if (i==ref_column_index and j==ref_depth_index)
                ref_density = density;

              const double integration_weight = (j == 0 || j == n_vertical_points-1) ? 0.5 : 1.0;
              column_mass += integration_weight * density * delta_z;
            }

          column_masses[i] = column_mass;
          total_column_mass += column_mass;

        }

      AssertThrow(ref_density > 0.0,
                  ExcMessage("The reference density is not set properly by the isostatic topography model."));

      const double average_column_mass = total_column_mass / n_lateral_points;

      for (unsigned int i = 0; i < n_lateral_points; ++i)
        {
          topography[i] = (average_column_mass - column_masses[i]) / ref_density;
          topography[i] = std::clamp(topography[i], -max_isostatic_topography, max_isostatic_topography);
        }
    }



    template <int dim>
    Tensor<1,dim>
    Isostacy<dim>::
    compute_initial_deformation_on_boundary (
      const types::boundary_id boundary_indicator,
      const Point<dim> &position) const
    {
      AssertThrow(boundary_indicator == this->get_geometry_model().translate_symbolic_boundary_name_to_id("top"),
                  ExcMessage("The isostaty mesh deformation option should only be applied to the top boundary."));

      AssertThrow(!topography.empty(),
                  ExcMessage("The isostatic topography profile has not been properly initialized."));

      const auto &geometry_model = dynamic_cast<const GeometryModel::Box<dim> &>(this->get_geometry_model());

      const double x_min = 0.0;
      const double x_max = geometry_model.get_extents()[0];
      const double x = position[0];

      const double grid_spacing = (x_max - x_min) / (n_lateral_points - 1.0);
      const double index = (x - x_min) / grid_spacing;
      const unsigned int i = static_cast<unsigned int>(std::floor(index));
      const double alpha = index - i;

      const double topo =  (1.0 - alpha) * topography[i] + alpha * topography[i+1];

      const Tensor<1,dim> gravity =
        this->get_gravity_model().gravity_vector(position);

      Tensor<1,dim> deformation_direction;

      if (gravity.norm() > 0.0)
        deformation_direction = -gravity / gravity.norm();

      return topo * deformation_direction;
    }



    template <int dim>
    bool
    Isostacy<dim>::needs_surface_stabilization () const
    {
      return false;
    }



    template <int dim>
    void
    Isostacy<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh deformation");
      {
        prm.enter_subsection("Isostacy");
        {
          prm.declare_entry("Number of lateral points",
                            "2000",
                            Patterns::Integer(2),
                            "Number of equally spaced lateral sampling points.");

          prm.declare_entry("Number of vertical points",
                            "2000",
                            Patterns::Integer(2),
                            "Number of equally spaced vertical sampling points.");

          prm.declare_entry("Balanced depth",
                            "500000",
                            Patterns::Double(0),
                            "Depth at which the reference density is evaluated.");

          prm.declare_entry("Maximum isostatic topography",
                            "10000",
                            Patterns::Double(0),
                            "Maximum magnitude of the computed initial topography.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Isostacy<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh deformation");
      {
        prm.enter_subsection("Isostacy");
        {
          n_lateral_points =
            prm.get_integer("Number of lateral points");

          n_vertical_points =
            prm.get_integer("Number of vertical points");

          balanced_depth =
            prm.get_double("Balanced depth");

          max_isostatic_topography =
            prm.get_double("Maximum isostatic topography");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(
      Isostacy,
      "isostacy",
      "A mesh deformation model that computes an initial isostatic surface "
      "deformation by balancing the mass of vertical columns.")
  }
}
