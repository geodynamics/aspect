/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/initial_temperature/adiabatic.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/two_merged_boxes.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>

#include <cmath>

namespace aspect
{
  namespace InitialTemperature
  {

    template <int dim>
    Adiabatic<dim>::Adiabatic ()
      :
      surface_boundary_id(numbers::invalid_unsigned_int)
    {}

    template <int dim>
    void
    Adiabatic<dim>::initialize ()
    {
      // Find the boundary indicator that represents the surface
      surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      std::set<types::boundary_id> surface_boundary_set;
      surface_boundary_set.insert(surface_boundary_id);

      // The input ascii table contains one data column (LAB depths(m)) in addition to the coordinate columns.
      Utilities::AsciiDataBoundary<dim>::initialize(surface_boundary_set,
                                                    1);
    }

    template <int dim>
    double
    Adiabatic<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      double age_top = 0;
      const double age_bottom = (this->convert_output_to_years() ? age_bottom_boundary_layer * year_in_seconds
                                 : age_bottom_boundary_layer);
      if (read_from_ascii_file)
        {
          // The input ascii contains the age of the seafloor. User must provide ages in seconds
          age_top = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id,
                                                                          position,
                                                                          0);
        }
      else
        {
          // convert input ages to seconds
          age_top =    (this->convert_output_to_years() ? age_top_boundary_layer * year_in_seconds
                        : age_top_boundary_layer);
        }
      // First, get the temperature of the adiabatic profile at a representative
      // point at the top and bottom boundary of the model
      // if adiabatic heating is switched off, assume a constant profile
      const Point<dim> surface_point = this->get_geometry_model().representative_point(0.0);
      const Point<dim> bottom_point = this->get_geometry_model().representative_point(this->get_geometry_model().maximal_depth());
      const double adiabatic_surface_temperature = this->get_adiabatic_conditions().temperature(surface_point);
      const double adiabatic_bottom_temperature = (this->include_adiabatic_heating())
                                                  ?
                                                  this->get_adiabatic_conditions().temperature(bottom_point)
                                                  :
                                                  adiabatic_surface_temperature;

      // then, get the temperature at the top and bottom boundary of the model
      // if no boundary temperature is prescribed simply use the adiabatic.
      // This implementation assumes that the top and bottom boundaries have
      // prescribed temperatures and minimal_temperature() returns the value
      // at the surface and maximal_temperature() the value at the bottom.
      const double T_surface = (this->has_boundary_temperature()
                                ?
                                this->get_boundary_temperature_manager().minimal_temperature(
                                  this->get_fixed_temperature_boundary_indicators())
                                :
                                adiabatic_surface_temperature);
      const double T_bottom = (this->has_boundary_temperature()
                               ?
                               this->get_boundary_temperature_manager().maximal_temperature(
                                 this->get_fixed_temperature_boundary_indicators())
                               :
                               adiabatic_bottom_temperature);
      const double depth = this->get_geometry_model().depth(position);

      // look up material properties
      MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(1, this->n_compositional_fields());

      in.position[0]=position;
      in.temperature[0]=this->get_adiabatic_conditions().temperature(position);
      in.pressure[0]=this->get_adiabatic_conditions().pressure(position);
      in.velocity[0]= Tensor<1,dim> ();
      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        in.composition[0][c] = function->value(Point<1>(depth),c);
      in.strain_rate.resize(0); // adiabat has strain=0.
      this->get_material_model().evaluate(in, out);

      const double kappa = ( (this->get_parameters().formulation_temperature_equation ==
                              Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
                             ?
                             out.thermal_conductivities[0] /
                             (this->get_adiabatic_conditions().density(in.position[0]) * out.specific_heat[0])
                             :
                             out.thermal_conductivities[0] / (out.densities[0] * out.specific_heat[0])
                           );

      double surface_cooling_temperature = 0;
      double bottom_heating_temperature = 0;
      if (cooling_model == "half-space cooling")
        {
          // analytical solution for the thermal boundary layer from half-space cooling model
          surface_cooling_temperature = age_top > 0.0 ?
                                        (T_surface - adiabatic_surface_temperature) *
                                        erfc(this->get_geometry_model().depth(position) /
                                             (2 * sqrt(kappa * age_top)))
                                        : 0.0;
          bottom_heating_temperature = (age_bottom > 0.0 && this->get_adiabatic_conditions().is_initialized()) ?
                                       (T_bottom - adiabatic_bottom_temperature + subadiabaticity)
                                       * erfc((this->get_geometry_model().maximal_depth()
                                               - this->get_geometry_model().depth(position)) /
                                              (2 * sqrt(kappa * age_bottom)))
                                       : 0.0;
        }

      if (cooling_model == "plate cooling")
        {
          if (depth > lithosphere_thickness)
            {
              surface_cooling_temperature = 0;
            }
          else
            {
              const double exponential = -kappa * std::pow(numbers::PI, 2) * age_top / std::pow(lithosphere_thickness, 2);
              double sum_terms = 0;
              for (unsigned int n=1; n<11; ++n)
                {
                  sum_terms += 1/(double)n * std::exp(std::pow((double)n, 2) * exponential) * std::sin((double)n * depth * numbers::PI / lithosphere_thickness);
                  surface_cooling_temperature = T_surface - adiabatic_surface_temperature + (adiabatic_surface_temperature - T_surface) * (depth / lithosphere_thickness + 2 / numbers::PI * sum_terms);
                }
            }
        }

      // set the initial temperature perturbation
      // first: get the center of the perturbation, then check the distance to the
      // evaluation point. the center is supposed to lie at the center of the bottom
      // surface.
      Point<dim> mid_point;
      if (perturbation_position == "center")
        {
          if (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model()))
            {
              const GeometryModel::SphericalShell<dim> &shell_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model());

              const double inner_radius = shell_geometry_model.inner_radius();
              const double half_opening_angle = numbers::PI/180.0 * 0.5 * shell_geometry_model.opening_angle();
              if (dim==2)
                {
                  // choose the center of the perturbation at half angle along the inner radius
                  mid_point(0) = inner_radius * std::sin(half_opening_angle),
                  mid_point(1) = inner_radius * std::cos(half_opening_angle);
                }
              else if (dim==3)
                {
                  // if the opening angle is 90 degrees (an eighth of a full spherical
                  // shell, then choose the point on the inner surface along the first
                  // diagonal
                  if (shell_geometry_model.opening_angle() == 90)
                    {
                      mid_point(0) = inner_radius*std::sqrt(1./3),
                      mid_point(1) = inner_radius*std::sqrt(1./3),
                      mid_point(2) = inner_radius*std::sqrt(1./3);
                    }
                  else
                    {
                      // otherwise do the same as in 2d
                      mid_point(0) = inner_radius * std::sin(half_opening_angle) * std::cos(half_opening_angle),
                      mid_point(1) = inner_radius * std::sin(half_opening_angle) * std::sin(half_opening_angle),
                      mid_point(2) = inner_radius * std::cos(half_opening_angle);
                    }
                }
            }
          else if (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>>(this->get_geometry_model()))
            {
              const GeometryModel::Chunk<dim> &chunk_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::Chunk<dim>> (this->get_geometry_model());

              const double inner_radius = chunk_geometry_model.inner_radius();

              const double west_longitude = chunk_geometry_model.west_longitude(); // in radians
              const double longitude_range = chunk_geometry_model.longitude_range(); // in radians
              const double longitude_midpoint = west_longitude + 0.5 * longitude_range;

              if (dim==2)
                {
                  // choose the center of the perturbation at half angle along the inner radius
                  mid_point(0) = inner_radius * std::cos(longitude_midpoint),
                  mid_point(1) = inner_radius * std::sin(longitude_midpoint);
                }
              else if (dim==3)
                {

                  const double south_latitude = chunk_geometry_model.south_latitude(); // in radians
                  const double latitude_range = chunk_geometry_model.latitude_range(); // in radians
                  const double latitude_midpoint = south_latitude + 0.5 * latitude_range;
                  mid_point(0) = inner_radius * std::cos(latitude_midpoint) * std::cos(longitude_midpoint);
                  mid_point(1) = inner_radius * std::cos(latitude_midpoint) * std::sin(longitude_midpoint);
                  mid_point(2) = inner_radius * std::sin(latitude_midpoint);
                }
            }
          else if (Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()))
            {
              const GeometryModel::Box<dim> &box_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::Box<dim>> (this->get_geometry_model());

              // for the box geometry, choose a point at the center of the bottom face.
              // (note that the loop only runs over the first dim-1 coordinates, leaving
              // the depth variable at zero)
              mid_point = box_geometry_model.get_origin();
              for (unsigned int i=0; i<dim-1; ++i)
                mid_point(i) += 0.5 * box_geometry_model.get_extents()[i];
            }
          else if (Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model()))
            {
              const GeometryModel::TwoMergedBoxes<dim> &two_merged_boxes_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model());

              // for the box geometry, choose a point at the center of the bottom face.
              // (note that the loop only runs over the first dim-1 coordinates, leaving
              // the depth variable at zero)
              mid_point = two_merged_boxes_geometry_model.get_origin();
              for (unsigned int i=0; i<dim-1; ++i)
                mid_point(i) += 0.5 * two_merged_boxes_geometry_model.get_extents()[i];
            }
          else
            AssertThrow (false,
                         ExcMessage ("Not a valid geometry model for the initial temperature model"
                                     "adiabatic."));
        }

      const double perturbation = (mid_point.distance(position) < radius) ? amplitude
                                  : 0.0;


      // add the subadiabaticity
      const double zero_depth = 0.174;
      const double nondimensional_depth = (this->get_geometry_model().depth(position) / this->get_geometry_model().maximal_depth() - zero_depth)
                                          / (1.0 - zero_depth);
      double subadiabatic_T = 0.0;
      if (nondimensional_depth > 0)
        subadiabatic_T = -subadiabaticity * nondimensional_depth * nondimensional_depth;

      // If adiabatic heating is disabled, apply all perturbations to
      // constant adiabatic surface temperature instead of adiabatic profile.
      const double temperature_profile = (this->include_adiabatic_heating())
                                         ?
                                         this->get_adiabatic_conditions().temperature(position)
                                         :
                                         adiabatic_surface_temperature;

      // return sum of the adiabatic profile, the boundary layer temperatures and the initial
      // temperature perturbation.
      return temperature_profile + surface_cooling_temperature
             + (perturbation > 0.0 ? std::max(bottom_heating_temperature + subadiabatic_T,perturbation)
                : bottom_heating_temperature + subadiabatic_T);
    }

    template <int dim>
    void
    Adiabatic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-temperature/adiabatic/",
                                                          "adiabatic.txt",
                                                          "Adiabatic");
        prm.enter_subsection("Adiabatic");
        {
          prm.declare_entry ("Age top boundary layer", "0.",
                             Patterns::Double (0.),
                             "The age of the upper thermal boundary layer, used for the calculation "
                             "of the half-space cooling model temperature. Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Age bottom boundary layer", "0.",
                             Patterns::Double (0.),
                             "The age of the lower thermal boundary layer, used for the calculation "
                             "of the half-space cooling model temperature. Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Radius", "0.",
                             Patterns::Double (0.),
                             "The Radius (in m) of the initial spherical temperature perturbation "
                             "at the bottom of the model domain.");
          prm.declare_entry ("Amplitude", "0.",
                             Patterns::Double (0.),
                             "The amplitude (in K) of the initial spherical temperature perturbation "
                             "at the bottom of the model domain. This perturbation will be added to "
                             "the adiabatic temperature profile, but not to the bottom thermal "
                             "boundary layer. Instead, the maximum of the perturbation and the bottom "
                             "boundary layer temperature will be used.");
          prm.declare_entry ("Position", "center",
                             Patterns::Selection ("center"),
                             "Where the initial temperature perturbation should be placed. If `center' is "
                             "given, then the perturbation will be centered along a `midpoint' of some "
                             "sort of the bottom boundary. For example, in the case of a box geometry, "
                             "this is the center of the bottom face; in the case of a spherical shell "
                             "geometry, it is along the inner surface halfway between the bounding "
                             "radial lines.");
          prm.declare_entry ("Subadiabaticity", "0.",
                             Patterns::Double (0.),
                             "If this value is larger than 0, the initial temperature profile will "
                             "not be adiabatic, but subadiabatic. This value gives the maximal "
                             "deviation from adiabaticity. Set to 0 for an adiabatic temperature "
                             "profile. Units: \\si{\\kelvin}.\n\n"
                             "The function object in the Function subsection "
                             "represents the compositional fields that will be used as a reference "
                             "profile for calculating the thermal diffusivity. "
                             "This function is one-dimensional and depends only on depth. The format of this "
                             "functions follows the syntax understood by the "
                             "muparser library, see Section~\\ref{sec:muparser-format}.");
          prm.declare_entry ("Use ASCII file for seafloor age", "false",
                             Patterns::Bool (),
                             "Whether to define seafloor ages with an ASCII data file.");
          prm.declare_entry ("Cooling model", "half-space cooling",
                             Patterns::Selection ("half-space cooling|plate cooling"),
                             "Whether to use the half space cooling model or the plate cooling model");
          prm.declare_entry ("Lithosphere thickness", "125e3",
                             Patterns::Double (0.),
                             "Thickness of the lithosphere for plate cooling model.");
          prm.enter_subsection("Function");
          {
            Functions::ParsedFunction<1>::declare_parameters (prm, 1);
          }
          prm.leave_subsection();
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Adiabatic<dim>::parse_parameters (ParameterHandler &prm)
    {
      // we need to get the number of compositional fields here to
      // initialize the function parser. unfortunately, we can't get it
      // via SimulatorAccess from the simulator itself because at the
      // current point the SimulatorAccess hasn't been initialized
      // yet. so get it from the parameter file directly.
      prm.enter_subsection ("Compositional fields");
      const unsigned int n_compositional_fields = prm.get_integer ("Number of fields");
      prm.leave_subsection ();

      prm.enter_subsection ("Initial temperature model");
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm, "Adiabatic");
        prm.enter_subsection("Adiabatic");
        {
          age_top_boundary_layer = prm.get_double ("Age top boundary layer");
          age_bottom_boundary_layer = prm.get_double ("Age bottom boundary layer");
          radius = prm.get_double ("Radius");
          amplitude = prm.get_double ("Amplitude");
          perturbation_position = prm.get("Position");
          subadiabaticity = prm.get_double ("Subadiabaticity");
          read_from_ascii_file = prm.get_bool ("Use ASCII file for seafloor age");
          cooling_model = prm.get ("Cooling model");
          lithosphere_thickness = prm.get_double ("Lithosphere thickness");
          if (n_compositional_fields > 0)
            {
              prm.enter_subsection("Function");
              try
                {
                  function
                    = std_cxx14::make_unique<Functions::ParsedFunction<1>>(n_compositional_fields);
                  function->parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "ERROR: FunctionParser failed to parse\n"
                            << "\t'Initial temperature model.Adiabatic.Function'\n"
                            << "with expression\n"
                            << "\t'" << prm.get("Function expression") << "'"
                            << "More information about the cause of the parse error \n"
                            << "is shown below.\n";
                  throw;
                }

              prm.leave_subsection();
            }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(Adiabatic,
                                              "adiabatic",
                                              "Temperature is prescribed as an adiabatic "
                                              "profile with upper and lower thermal boundary layers, "
                                              "whose ages are given as input parameters.")
  }
}
