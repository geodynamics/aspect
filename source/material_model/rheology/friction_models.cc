/*
  Copyright (C) 2019 - 2024 by the authors of the ASPECT code.

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

#include <aspect/material_model/rheology/friction_models.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_signals.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <boost/geometry/algorithms/comparable_distance.hpp>
#include <boost/geometry/index/predicates.hpp>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      double
      FrictionModels<dim>::
      compute_friction_angle(const double current_edot_ii,
                             const unsigned int volume_fraction_index,
                             const double static_friction_angle,
                             const Point<dim> &position) const
      {

        switch (friction_mechanism)
          {
            case static_friction:
            {
              return static_friction_angle;
            }
            case dynamic_friction:
            {
              // Calculate effective steady-state friction coefficient.
              // This is based on the former material model dynamic friction.
              // The formula below is equivalent to the equation 13 in \cite{van_dinther_seismic_2013}.
              // Although here the dynamic friction coefficient is directly specified. In addition,
              // we also use a reference strain rate in place of a characteristic
              // velocity divided by local element size. This reference strain rate is called
              // the dynamic characteristic strain rate and is used to compute what value between
              // dynamic and static angle of internal friction should be used.
              // Furthermore a smoothness exponent X is added, which determines whether the
              // friction vs strain rate curve is rather step-like or more gradual.
              // mu  = mu_d + (mu_s - mu_d) / ( (1 + strain_rate_dev_inv2/dynamic_characteristic_strain_rate)^X );
              // Angles of friction are used in radians within ASPECT. The coefficient
              // of friction (mu) is the tangent of the internal angle of friction, hence convergence is needed.
              // The incoming variable static_friction_angle holds the static friction angle. It's value is
              // then updated with the strain rate dependent friction angle which is returned by the function.
              const double mu = (std::tan(dynamic_angles_of_internal_friction[volume_fraction_index])
                                 + (std::tan(static_friction_angle) - std::tan(dynamic_angles_of_internal_friction[volume_fraction_index]))
                                 / (1. + std::pow((current_edot_ii / dynamic_characteristic_strain_rate),
                                                  dynamic_friction_smoothness_exponent)));
              const double dynamic_friction_angle = std::atan (mu);
              Assert((mu < 1) && (0 < dynamic_friction_angle) && (dynamic_friction_angle <= 1.6), ExcMessage(
                       "The friction coefficient should be larger than zero and smaller than 1. "
                       "The friction angle should be smaller than 1.6 rad."));
              return dynamic_friction_angle;
            }
            case differential_dynamic_friction:
            {
              // Select a convergent, divergent, or default target angle from
              // the projected tangential surface-velocity divergence. The full
              // velocity divergence can vanish in an incompressible model while
              // this surface indicator is nonzero.
              // The default dynamic angle represents transform/neutral deformation.
              double target_angle = dynamic_angles_of_internal_friction[volume_fraction_index];
              const TectonicRegime tectonic_regime = compute_tectonic_regime(position);
              if (tectonic_regime == convergent_regime)
                target_angle = dynamic_angles_of_internal_friction_for_convergence[volume_fraction_index];
              else if (tectonic_regime == divergent_regime)
                target_angle = dynamic_angles_of_internal_friction_for_divergence[volume_fraction_index];

              // Compute the strain-rate dependent friction coefficient:
              // μ = μ_d + (μ_s - μ_d) / [1 + (ε̇ / ε̇_c)^x]
              // where μ_d is derived from target_angle and μ_s from static friction angle.
              const double mu = std::tan(target_angle)
                                + (std::tan(static_friction_angle) - std::tan(target_angle))
                                  / (1. + std::pow((current_edot_ii / dynamic_characteristic_strain_rate),
                                                  dynamic_friction_smoothness_exponent));

              const double dynamic_friction_angle = std::atan(mu);
              Assert((mu < 1) && (0 < dynamic_friction_angle) && (dynamic_friction_angle <= 1.6),
                     ExcMessage("The friction coefficient should be larger than zero and smaller than 1. "
                                "The friction angle should be smaller than 1.6 rad."));
              return dynamic_friction_angle;
            }
            case function:
            {
              // Use a given function input per composition to get the friction angle
              Utilities::NaturalCoordinate<dim> point =
                this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system_friction_function);

              // we get time passed as seconds (always) but may want
              // to reinterpret it in years
              if (this->convert_output_to_years())
                friction_function->set_time (this->get_time() / year_in_seconds);
              else
                friction_function->set_time (this->get_time());

              // determine the friction angle based on position and composition
              // This plugin reads a function for background material and every chemical compositional field.
              // We assume the order of the functions is the same as the order of the volume fractions.
              double friction_from_function =
                friction_function->value(Utilities::convert_array_to_point<dim>(point.get_coordinates()), volume_fraction_index);

              // Convert angles from degrees to radians
              friction_from_function *= constants::degree_to_radians;

              return friction_from_function;
            }
          }
        // we should never get here, return something anyway, so the compiler does not complain...
        AssertThrow (false, ExcMessage("Unknown friction model."));
        return static_friction_angle;
      }

      template <int dim>
      FrictionMechanism
      FrictionModels<dim>::
      get_friction_mechanism() const
      {
        return friction_mechanism;
      }



      template <int dim>
      double
      FrictionModels<dim>::
      compute_tectonic_divergence_indicator(const Point<dim> &position) const
      {
        if (this->get_geometry_model().depth(position) > surface_regime_projection_depth)
          return 0.0;

        // No velocity field exists before the first completed Stokes solve, so
        // the initial solve is neutral. Never fall back to a volume strain-rate
        // diagnostic: surface velocity divergence is the only regime source.
        if (surface_divergence_index.empty())
          return 0.0;

        return interpolate_surface_velocity_divergence(position);
      }



      template <int dim>
      TectonicRegime
      FrictionModels<dim>::
      compute_tectonic_regime(const Point<dim> &position) const
      {
        const double divergence_indicator = compute_tectonic_divergence_indicator(position);
        if (divergence_indicator < -convergence_threshold)
          return convergent_regime;
        if (divergence_indicator > divergence_threshold)
          return divergent_regime;
        return transform_or_neutral_regime;
      }



      template <int dim>
      typename FrictionModels<dim>::SurfaceIndexPoint
      FrictionModels<dim>::
      make_surface_index_point(const Point<dim> &position) const
      {
        SurfaceIndexPoint index_point;
        boost::geometry::set<0>(index_point, 0.0);
        boost::geometry::set<1>(index_point, 0.0);
        boost::geometry::set<2>(index_point, 0.0);

        if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical)
          {
            const double radius = position.norm();
            AssertThrow(radius > 0.0, ExcMessage("Cannot project the model origin onto the surface."));
            boost::geometry::set<0>(index_point, position[0] / radius);
            boost::geometry::set<1>(index_point, position[1] / radius);
            if (dim == 3)
              boost::geometry::set<2>(index_point, position[2] / radius);
          }
        else
          {
            const double length_scale = this->get_geometry_model().length_scale();
            boost::geometry::set<0>(index_point, position[0] / length_scale);
            if (dim == 3)
              boost::geometry::set<1>(index_point, position[1] / length_scale);
          }
        return index_point;
      }



      template <int dim>
      double
      FrictionModels<dim>::
      interpolate_surface_velocity_divergence(const Point<dim> &position) const
      {
        // Surface divergence is stored at surface-face quadrature points, but
        // projected interior material points generally do not coincide with
        // those samples (especially after refinement or MPI partitioning).
        // Interpolate the signed divergence before classifying the regime to
        // avoid discontinuous, mesh-dependent nearest-neighbor patches.
        namespace bgi = boost::geometry::index;
        const SurfaceIndexPoint query_point = make_surface_index_point(position);
        std::vector<SurfaceIndexValue> nearest_samples;
        surface_divergence_index.query(bgi::nearest(query_point, 4),
                                       std::back_inserter(nearest_samples));
        AssertThrow(!nearest_samples.empty(), ExcInternalError());

        double weighted_value = 0.0;
        double weight_sum = 0.0;
        for (const SurfaceIndexValue &sample : nearest_samples)
          {
            const double distance_squared = boost::geometry::comparable_distance(query_point, sample.first);
            if (distance_squared < 1e-30)
              return sample.second;
            const double weight = 1.0 / distance_squared;
            weighted_value += weight * sample.second;
            weight_sum += weight;
          }
        return weighted_value / weight_sum;
      }



      template <int dim>
      void
      FrictionModels<dim>::
      update_surface_velocity_divergence()
      {
        const Quadrature<dim-1> &quadrature = this->introspection().face_quadratures.velocities;
        FEFaceValues<dim> face_values(this->get_mapping(),
                                      this->get_fe(),
                                      quadrature,
                                      update_values | update_gradients |
                                      update_quadrature_points | update_normal_vectors);
        std::vector<Tensor<1,dim>> velocities(quadrature.size());
        std::vector<Tensor<2,dim>> velocity_gradients(quadrature.size());
        std::vector<std::vector<double>> local_samples;
        const types::boundary_id surface_boundary
          = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned())
            for (const unsigned int face_number : cell->face_indices())
              if (cell->face(face_number)->at_boundary()
                  && cell->face(face_number)->boundary_id() == surface_boundary)
                {
                  face_values.reinit(cell, face_number);
                  face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), velocities);
                  face_values[this->introspection().extractors.velocities].get_function_gradients(this->get_solution(), velocity_gradients);

                  for (unsigned int q = 0; q < quadrature.size(); ++q)
                    {
                      const Point<dim> &surface_position = face_values.quadrature_point(q);
                      const Tensor<1,dim> normal = face_values.normal_vector(q);
                      const SymmetricTensor<2,dim> strain_rate = symmetrize(velocity_gradients[q]);
                      double surface_divergence = trace(strain_rate)
                                                  - normal * (strain_rate * normal);

                      // Remove the apparent surface-area change caused by radial
                      // motion so that the spherical case uses tangential velocity.
                      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical)
                        surface_divergence -= velocities[q] * normal * (dim-1) / surface_position.norm();

                      const SurfaceIndexPoint index_point = make_surface_index_point(surface_position);
                      local_samples.push_back({boost::geometry::get<0>(index_point),
                                               boost::geometry::get<1>(index_point),
                                               boost::geometry::get<2>(index_point),
                                               surface_divergence});
                    }
                }

        const std::vector<std::vector<std::vector<double>>> gathered_samples
          = Utilities::MPI::all_gather(this->get_mpi_communicator(), local_samples);
        std::vector<SurfaceIndexValue> index_values;
        for (const std::vector<std::vector<double>> &process_samples : gathered_samples)
          for (const std::vector<double> &sample : process_samples)
            {
              SurfaceIndexPoint point;
              boost::geometry::set<0>(point, sample[0]);
              boost::geometry::set<1>(point, sample[1]);
              boost::geometry::set<2>(point, sample[2]);
              index_values.emplace_back(point, sample[3]);
            }
        // Differentiating the discrete surface velocity amplifies cell-scale
        // oscillations. Average each sample with its nearest surface neighbors
        // before projection so regimes represent coherent plate-scale motion,
        // rather than alternating finite-element patches.
        namespace bgi = boost::geometry::index;
        const SurfaceIndex raw_index(index_values.begin(), index_values.end());
        std::vector<SurfaceIndexValue> smoothed_values;
        smoothed_values.reserve(index_values.size());
        for (const SurfaceIndexValue &value : index_values)
          {
            std::vector<SurfaceIndexValue> neighbors;
            raw_index.query(bgi::nearest(value.first, std::min<std::size_t>(8, index_values.size())),
                            std::back_inserter(neighbors));
            double sum = 0.0;
            for (const SurfaceIndexValue &neighbor : neighbors)
              sum += neighbor.second;
            smoothed_values.emplace_back(value.first, sum / neighbors.size());
          }
        surface_divergence_index = SurfaceIndex(smoothed_values.begin(), smoothed_values.end());
      }



      template <int dim>
      void
      FrictionModels<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Friction mechanism", "none",
                           Patterns::Selection("none|dynamic friction|differential dynamic friction|function"),
                           "Whether to make the friction angle dependent on strain rate or not. This rheology "
                           "is intended to be used together with the visco-plastic rheology model."
                           "\n\n"
                           "\\item ``none'': No dependence of the friction angle is applied. "
                           "\n\n"
                           "\\item ``dynamic friction'': The friction angle is rate dependent."
                           "When 'dynamic angles of internal friction' are specified, "
                           "the friction angle will be weakened for high strain rates with: "
                           "$\\mu = \\mu_d + \\frac{\\mu_s-\\mu_d}{1+\\frac{\\dot{\\epsilon}_{ii}}{\\dot{\\epsilon}_C}}^x$  "
                           "where $\\mu_s$ and $\\mu_d$ are the friction angles at low and high strain rates, "
                           "respectively. $\\dot{\\epsilon}_{ii}$ is the second invariant of the strain rate and "
                           "$\\dot{\\epsilon}_C$ is the 'dynamic characteristic strain rate' where $\\mu = (\\mu_s+\\mu_d)/2$. "
                           "The 'dynamic friction smoothness exponent' x controls how "
                           "smooth or step-like the change from $\\mu_s$ to $\\mu_d$ is. "
                           "The equation is modified after Equation (13) in \\cite{van_dinther_seismic_2013}. "
                           "$\\mu_s$ and $\\mu_d$ can be specified by setting 'Angles of internal friction' and "
                           "'Dynamic angles of internal friction', respectively. "
                           "This relationship is similar to rate-and-state friction constitutive relationships, which "
                           "are applicable to the strength of rocks during earthquakes."
                           "\n\n"
                           "\\item ``differential dynamic friction'': As for ``dynamic friction'', but the dynamic "
                           "angle is selected from convergent, divergent, and transform/neutral values using the "
                           "divergence of tangential velocity on the top boundary. The surface field is smoothed "
                           "and projected down to the specified depth. This indicator remains meaningful in an "
                           "incompressible model because it is not the divergence of the full velocity field."
                           "\n\n"
                           "\\item ``function'': Specify the friction angle as a function of space and time "
                           "for background material and compositional fields, for a total of N$+$1 values, "
                           "where N is the number of all compositional fields corresponding to chemical compositions.");

        // Dynamic friction parameters
        prm.declare_entry ("Dynamic characteristic strain rate", "1e-12",
                           Patterns::Double (0),
                           "The characteristic strain rate value at which the angle of friction is "
                           "equal to $\\mu = (\\mu_s+\\mu_d)/2$. When the effective strain rate "
                           "is very high, the dynamic angle of friction is taken, when it is very low, "
                           "the static angle of internal friction is used. Around the dynamic characteristic "
                           "strain rate, there is a smooth gradient from the static to the dynamic angle "
                           "of internal friction. "
                           "Units: \\si{\\per\\second}.");

        prm.declare_entry ("Dynamic angles of internal friction", "2",
                           Patterns::List(Patterns::Double(0)),
                           "List of dynamic angles of internal friction, $\\phi$, for background material and compositional "
                           "fields, for a total of N$+$1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "Dynamic angles of friction are used as the current friction angle when the effective "
                           "strain rate is well above the 'dynamic characteristic strain rate'. For the "
                           "'differential dynamic friction' mechanism these values apply to transform or neutral flow. "
                           "Units: \\si{\\degree}.");

        prm.declare_entry ("Dynamic friction smoothness exponent", "1",
                           Patterns::Double(0),
                           "An exponential factor in the equation for the calculation of the friction angle when a "
                           "static and a dynamic angle of internal friction are specified. A factor of 1 returns the equation "
                           "to Equation (13) in \\cite{van_dinther_seismic_2013}. A factor between 0 and 1 makes the "
                           "curve of the friction angle vs. the strain rate smoother, while a factor $>$ 1 makes "
                           "the change between static and dynamic friction angle more steplike. "
                           "Units: none.");

        prm.declare_entry ("Convergence threshold", "1e-15",
                          Patterns::Double(0),
                          "Magnitude of negative tangential surface-velocity divergence required "
                          "to classify flow as convergent. Units: \\si{\\per\\second}.");

        prm.declare_entry ("Divergence threshold", "1e-15",
                          Patterns::Double(0),
                          "Positive tangential surface-velocity divergence required to classify "
                          "flow as divergent. Units: \\si{\\per\\second}.");

        prm.declare_entry ("Surface regime projection depth", "200e3",
                          Patterns::Double(0),
                          "Maximum depth to which the surface-velocity-divergence classification "
                          "is projected. Below this depth the default dynamic friction angle is "
                          "selected. Units: \\si{\\meter}.");

        prm.declare_entry ("Dynamic angles of internal friction for convergent flow", "2",
          Patterns::List(Patterns::Double(0)),
          "Dynamic friction angles for convergent flow. Units: \\si{\\degree}.");

        prm.declare_entry ("Dynamic angles of internal friction for divergent flow", "4",
          Patterns::List(Patterns::Double(0)),
          "Dynamic friction angles for divergent flow. Units: \\si{\\degree}.");

        /**
         * If friction is specified as a function input.
         */
        prm.enter_subsection("Friction function");
        {
          /**
           * The function to specify the friction angle per composition can be declared in dependence
           * of depth, cartesian coordinates or spherical coordinates. Note that the order
           * of spherical coordinates is r,phi,theta and not r,theta,phi, since
           * this allows for dimension independent expressions.
           */
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2d/3d "
                             "respectively with theta being the polar angle. `depth' "
                             "will create a function, in which only the first "
                             "parameter is non-zero, which is interpreted to "
                             "be the depth of the point.");

          Functions::ParsedFunction<dim>::declare_parameters(prm,1);
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      FrictionModels<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Friction dependence parameters
        if (prm.get ("Friction mechanism") == "none")
          friction_mechanism = static_friction;
        else if (prm.get ("Friction mechanism") == "dynamic friction")
          friction_mechanism = dynamic_friction;
        else if (prm.get ("Friction mechanism") == "differential dynamic friction")
          friction_mechanism = differential_dynamic_friction;
        else if (prm.get ("Friction mechanism") == "function")
          friction_mechanism = function;
        else
          AssertThrow(false, ExcMessage("Not a valid friction mechanism option!"));

        // Dynamic friction parameters
        dynamic_characteristic_strain_rate = prm.get_double("Dynamic characteristic strain rate");

        // Differential dynamic friction parameters
        convergence_threshold = prm.get_double("Convergence threshold");
        divergence_threshold = prm.get_double("Divergence threshold");
        surface_regime_projection_depth = prm.get_double("Surface regime projection depth");
        if (friction_mechanism == differential_dynamic_friction)
          this->get_signals().post_nonlinear_solver.connect(
            [this](const SolverControl &)
            {
              this->update_surface_velocity_divergence();
            });


        // Retrieve the list of composition names
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
        // plastic strain
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(), "background");
        chemical_field_names.insert(chemical_field_names.begin(), "background");

        Utilities::MapParsing::Options options_default(chemical_field_names, "Dynamic angles of internal friction");
        options_default.list_of_allowed_keys = compositional_field_names;

        dynamic_angles_of_internal_friction = Utilities::MapParsing::parse_map_to_double_array (prm.get("Dynamic angles of internal friction"),
                                              options_default);

        // --- Convergent dynamic angles ---
        Utilities::MapParsing::Options options_convergence(chemical_field_names, "Dynamic angles of internal friction for convergent flow");
        options_convergence.list_of_allowed_keys = compositional_field_names;

        dynamic_angles_of_internal_friction_for_convergence = Utilities::MapParsing::parse_map_to_double_array (prm.get("Dynamic angles of internal friction for convergent flow"),
                                              options_convergence);

        // --- Divergent dynamic angles ---
        Utilities::MapParsing::Options options_divergence(chemical_field_names, "Dynamic angles of internal friction for divergent flow");
        options_divergence.list_of_allowed_keys = compositional_field_names;

        dynamic_angles_of_internal_friction_for_divergence = Utilities::MapParsing::parse_map_to_double_array (prm.get("Dynamic angles of internal friction for divergent flow"),
                                              options_divergence);


        // Convert angles from degrees to radians
        auto convert_angles_to_radians = [](std::vector<double> &angles, const std::string &description)
        {
          for (double &angle : angles)
            {
              AssertThrow(angle <= 90,
                          ExcMessage(description + " must be <= 90 degrees"));
              angle *= constants::degree_to_radians;
            }
        };

        convert_angles_to_radians(dynamic_angles_of_internal_friction, "Dynamic angles of internal friction");
        convert_angles_to_radians(dynamic_angles_of_internal_friction_for_convergence, "Dynamic angles of internal friction for convergent flow");
        convert_angles_to_radians(dynamic_angles_of_internal_friction_for_divergence, "Dynamic angles of internal friction for divergent flow");


        dynamic_friction_smoothness_exponent = prm.get_double("Dynamic friction smoothness exponent");



        // Get the number of fields for composition-dependent material properties
        // including the background field.
        // Parse the functions for the background material plus all chemical composition fields,
        // if friction is specified as a function.
        if (friction_mechanism == function)
          {
            prm.enter_subsection("Friction function");
            {
              coordinate_system_friction_function = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
              try
                {
                  friction_function
                    = std::make_unique<Functions::ParsedFunction<dim>>(chemical_field_names.size());
                  friction_function->parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "FunctionParser failed to parse\n"
                            << "\t friction function\n"
                            << "with expression \n"
                            << "\t' " << prm.get("Function expression") << "'";
                  throw;
                }
            }
            prm.leave_subsection();
          }

      }
    }
  }
}



// Explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
#define INSTANTIATE(dim) \
  template class FrictionModels<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
