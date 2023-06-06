/*
  Copyright (C) 2018 - 2021 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.


   Note that the empirical model used to define how Tmin increases with depth
   and how the position of Tmin shift with depth is expected to change somewhat
   after better calibrating with further tests.
*/

#include "world_builder/features/subducting_plate_models/temperature/mass_conserving.h"

#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/bool.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/point.h"
#include "world_builder/utilities.h"
#include "world_builder/world.h"

namespace WorldBuilder
{
  using namespace Utilities;

  namespace Features
  {
    namespace SubductingPlateModels
    {
      namespace Temperature
      {
        MassConserving::MassConserving(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          density(NaN::DSNAN),
          plate_velocity(NaN::DSNAN),
          mantle_coupling_depth(NaN::DSNAN),
          forearc_cooling_factor(NaN::DSNAN),
          thermal_conductivity(NaN::DSNAN),
          thermal_expansion_coefficient(NaN::DSNAN),
          specific_heat(NaN::DSNAN),
          thermal_diffusivity(NaN::DSNAN),
          potential_mantle_temperature(NaN::DSNAN),
          surface_temperature(NaN::DSNAN),
          taper_distance(NaN::DSNAN),
          adiabatic_heating(true),
          operation(Operations::REPLACE)
        {
          this->world = world_;
          this->name = "mass conserving";
        }

        MassConserving::~MassConserving() = default;

        void
        MassConserving::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add `plate velocity` to the required parameters.
          prm.declare_entry("", Types::Object({"plate velocity"}),
                            "Mass conserving temperature model. The temperature "
                            "model uses the heat content (proportional to to thermal mass anomaly) to "
                            "define a smooth temperature profile that conserves mass along the slab length. "
                            "An empirical model, using error functions for smooth transitions, is used to "
                            " define how the minimum temperature increases with depth and how the location of "
                            "the minimum temperature shifts into the slab interior. The slab is divided into top "
                            "and bottom parts, which meet at the location where the minimum temperature occurs in the slab. "
                            "For the bottom slab, the temperature is defined by a half-space cooling model. "
                            "For the top of the slab the temperature is defined by one side of a 1D infinite "
                            "space cooling model: this function was chosen to have a smoother temperature function across "
                            "the minimum temperature position. The age of the overriding plate is used so the slab temperature "
                            "at shallow depth smoothly transitions to the temperature of the overriding plate: "
                            "this is not perfect, and is affected by the value of \"top truncation\" parameter "
                            "subducting plate. Notes:"
                            "1) the parameter \"thickness\" for the subducting plate segments needs to be defined but is not used. "
                            "2) because we use a negative truncation for distance above the slab, it is recommended to use"
                            "depth method:begin at end segment, in the main part of the world-builder file."
                            "Other methods may lead to gpas in temperatures at the segment boundaries."
                            "3)the empirical model used to define how Tmin increases with depth "
                            "and how the position of Tmin shift with depth is expected to change somewhat "
                            "after better calibrating with further tests.");

          // Declare entries of this plugin
          prm.declare_entry("min distance slab top", Types::Double(0),
                            "The distance in meters from the top surface of the slab over which the temperature is "
                            "determined by this feature. This parameter should be negative and should be 1.5-2 times "
                            "larger than the nominal slab thickness to allow the diffusion of cold "
                            "temperatures from in the slab into the mantle above the slab surface. "
                            "Also note that the top truncation value for the slab segment needs to have a value "
                            "of -1, otherwise the temperature above the slab will be cut off at a distance less than "
                            "the value set here.");

          prm.declare_entry("max distance slab top", Types::Double(std::numeric_limits<double>::max()),
                            "The distance in meters from the top surface of the slab over which the temperature is "
                            "determined by this feature. This parameter should be positive and approximately 2.5-3.0 times "
                            "larger than the nominal slab thickness to allow the diffusion of cold"
                            "temperatures from in the slab into the mantle below the slab surface."
                            "For example if the slab starts with cold temperatures over a 100 km wide region, this"
                            "parameters should be about 250 km.");

          prm.declare_entry("density", Types::Double(3300),
                            "The reference density of the subducting plate in $kg/m^3$");

          prm.declare_entry("plate velocity", Types::Double(0.05),
                            "The velocity with which the plate subducts in meters per year. Default is 5 cm/yr");

          prm.declare_entry("coupling depth", Types::Double(100e3),
                            "The depth at which the slab surface first comes in contact with the hot mantle wedge "
                            "in meters. Default is 100 km.");

          prm.declare_entry("forearc cooling factor", Types::Double(1.0),
                            "Increase the value to create thin (~2 km) cold thermal boundary layer above the slab."
                            "Any value greater than 1 does NOT meet the instantaneous conservation of mass, but does allow "
                            "one to account for the history of insulating the forearc from heating up to this point in time. "
                            "Note younger subducting lithosphere provides less insulation, while thicker, older slabs "
                            "provide more insulation. Values up to 10 to 30 have been tested and don't cause any other "
                            "extraneous effects. The larger th value the more you are not meeting the mass conserving criteria, "
                            "so you don't want to see this affecting the temperature beyond the coupling depth as it will "
                            "increase the mass of the slab and affect how it sinks.  If you use higher values, you will start to "
                            "see that this creates a very thick cool layer above the entire slab - if you see this extending beyond "
                            "the coupling zone reduce the value. You should use a value of 1 first and then "
                            "only increase as little as possible to cool just the forearc region. "
                            "Please examine the output temperature carefully. ");

          prm.declare_entry("thermal conductivity", Types::Double(3.3),
                            "The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.");

          prm.declare_entry("thermal expansion coefficient", Types::Double(-1),
                            "The thermal expansivity of the subducting plate material in $K^{-1}$. If smaller than zero, the global value is used.");

          prm.declare_entry("specific heat", Types::Double(-1),
                            "The specific heat of the subducting plate material in $J kg^{-1} K^{-1}$. If smaller than zero, the global value is used.");

          prm.declare_entry("thermal diffusivity", Types::Double(-1),
                            "The thermal conductivity of the subducting plate material in $W m^{-1} K^{-1}$.");

          prm.declare_entry("adiabatic heating", Types::Bool(true),
                            "Whether adiabatic heating should be used for the slab.");

          prm.declare_entry("taper distance", Types::Double(100e3),
                            "Distance over which to taper the slab tip."
                            "tapers the initial heat content to zero and the minimum temperature to the background temperature.");

          prm.declare_entry("potential mantle temperature", Types::Double(-1),
                            "The potential temperature of the mantle at the surface in Kelvin. If smaller than zero, the global value is used.");

          prm.declare_entry("ridge coordinates", Types::Array(Types::Array(Types::Point<2>(), 2),1),
                            "An list of ridges. Each ridge is a lists of at least 2 2d points which "
                            "define the location of the ridge. You need to define at least one ridge."
                            "So the an example with two ridges is "
                            "[[[10,20],[20,30],[10,40]],[[50,10],[60,10]]].");
        }

        void
        MassConserving::parse_entries(Parameters &prm)
        {

          min_depth = prm.get<double>("min distance slab top");
          max_depth = prm.get<double>("max distance slab top");
          operation = string_operations_to_enum(prm.get<std::string>("operation"));

          density = prm.get<double>("density");
          thermal_conductivity = prm.get<double>("thermal conductivity");
          plate_velocity = prm.get<double>("plate velocity");

          mantle_coupling_depth = prm.get<double>("coupling depth");
          forearc_cooling_factor = prm.get<double>("forearc cooling factor");

          taper_distance = prm.get<double>("taper distance");

          thermal_expansion_coefficient = prm.get<double>("thermal expansion coefficient");
          if (thermal_expansion_coefficient < 0)
            thermal_expansion_coefficient = this->world->thermal_expansion_coefficient;

          specific_heat = prm.get<double>("specific heat");
          if (specific_heat < 0)
            specific_heat = this->world->specific_heat;

          thermal_diffusivity = prm.get<double>("thermal diffusivity");
          if (thermal_diffusivity < 0)
            thermal_diffusivity = this->world->thermal_diffusivity;

          adiabatic_heating = prm.get<bool>("adiabatic heating");

          potential_mantle_temperature = this->world->potential_mantle_temperature >= 0
                                         ?
                                         this->world->potential_mantle_temperature
                                         :
                                         prm.get<double>("potential mantle temperature");

          surface_temperature = this->world->surface_temperature;

          mid_oceanic_ridges = prm.get_vector<std::vector<Point<2>>>("ridge coordinates");
          const double dtr = prm.coordinate_system->natural_coordinate_system() == spherical ? Consts::PI / 180.0 : 1.0;
          for (auto &ridge_coordinates : mid_oceanic_ridges)
            for (auto &ridge_coordinate : ridge_coordinates)
              {
                ridge_coordinate *= dtr;
              }
        }

        double
        MassConserving::get_temperature(const Point<3> & /*position_in_cartesian_coordinates*/,
                                        const double depth,
                                        const double gravity_norm,
                                        double temperature_,
                                        const double /*feature_min_depth*/,
                                        const double /*feature_max_depth*/,
                                        const WorldBuilder::Utilities::PointDistanceFromCurvedPlanes &distance_from_planes,
                                        const AdditionalParameters &additional_parameters) const
        {

          const double distance_from_plane = distance_from_planes.distance_from_plane;

          if (distance_from_plane <= max_depth && distance_from_plane >= min_depth)
            {

              const CoordinateSystem coordinate_system = world->parameters.coordinate_system->natural_coordinate_system();
              double distance_ridge = std::numeric_limits<double>::max();
              const Point<3> trench_point = distance_from_planes.closest_trench_point;
              const Objects::NaturalCoordinate trench_point_natural = Objects::NaturalCoordinate(trench_point,
                                                                      *(world->parameters.coordinate_system));
              const Point<2> trench_point_2d(trench_point_natural.get_surface_coordinates(),trench_point_natural.get_coordinate_system());
              // find the distance between the trench and ridge


              // first find if the coordinate is on this side of a ridge
              unsigned int relevant_ridge = 0;


              // if there is only one ridge, there is no transform
              if (mid_oceanic_ridges.size() > 1)
                {
                  // There are more than one ridge, so there are transform faults
                  // Find the first which is on the same side
                  for (relevant_ridge = 0; relevant_ridge < mid_oceanic_ridges.size()-1; relevant_ridge++)
                    {
                      const Point<2> transform_point_0 = mid_oceanic_ridges[relevant_ridge+1][0];
                      const Point<2> transform_point_1 = mid_oceanic_ridges[relevant_ridge][mid_oceanic_ridges[relevant_ridge].size()-1];
                      const Point<2> reference_point   = mid_oceanic_ridges[relevant_ridge][0];

                      const bool reference_on_side_of_line = (transform_point_1[0] - transform_point_0[0])
                                                             * (reference_point[1] - transform_point_0[1])
                                                             - (transform_point_1[1] - transform_point_0[1])
                                                             * (reference_point[0] - transform_point_0[0])
                                                             < 0;
                      const bool checkpoint_on_side_of_line = (transform_point_1[0] - transform_point_0[0])
                                                              * (trench_point_2d[1] - transform_point_0[1])
                                                              - (transform_point_1[1] - transform_point_0[1])
                                                              * (trench_point_2d[0] - transform_point_0[0])
                                                              < 0;

                      if (reference_on_side_of_line == checkpoint_on_side_of_line)
                        {
                          break;
                        }

                    }
                }

              for (unsigned int i_coordinate = 0; i_coordinate < mid_oceanic_ridges[relevant_ridge].size() - 1; i_coordinate++)
                {
                  const Point<2> segment_point0 = mid_oceanic_ridges[relevant_ridge][i_coordinate];
                  const Point<2> segment_point1 = mid_oceanic_ridges[relevant_ridge][i_coordinate + 1];

                  // based on http://geomalgorithms.com/a02-_lines.html
                  const Point<2> v = segment_point1 - segment_point0;
                  const Point<2> w = trench_point_2d - segment_point0;

                  const double c1 = (w[0] * v[0] + w[1] * v[1]);
                  const double c2 = (v[0] * v[0] + v[1] * v[1]);

                  Point<2> Pb(coordinate_system);
                  // This part is needed when we want to consider segments instead of lines
                  // If you want to have infinite lines, use only the else statement.

                  if (c1 <= 0)
                    Pb = segment_point0;
                  else if (c2 <= c1)
                    Pb = segment_point1;
                  else
                    Pb = segment_point0 + (c1 / c2) * v;

                  Point<3> compare_point(coordinate_system);

                  compare_point[0] = coordinate_system == cartesian ? Pb[0] : trench_point_natural.get_depth_coordinate();
                  compare_point[1] = coordinate_system == cartesian ? Pb[1] : Pb[0];
                  compare_point[2] = coordinate_system == cartesian ? trench_point_natural.get_depth_coordinate() : Pb[1];

                  distance_ridge = std::min(distance_ridge, this->world->parameters.coordinate_system->distance_between_points_at_same_depth(Point<3>(trench_point_natural.get_coordinates(),trench_point_natural.get_coordinate_system()), compare_point));
                }

              const double km2m = 1.0e3; // 1000 m/km
              const double cm2m = 100; // 100 cm/m
              const double my = 1.0e6;  // 1e6 y/my
              const double seconds_in_year = 60.0 * 60.0 * 24.0 * 365.25;  // sec/y

              const double age_at_trench = distance_ridge / plate_velocity; // m/(m/y) = yr
              const double plate_age_sec = age_at_trench * seconds_in_year; // y --> seconds

              /* information about nearest point on the slab segment */
              const double distance_along_plane = distance_from_planes.distance_along_plane;
              const double depth_to_reference_surface = distance_from_planes.depth_reference_surface;
              const double total_segment_length = additional_parameters.total_local_segment_length;
              const double average_angle = distance_from_planes.average_angle;

              // 1. Determine initial heat content of the slab based on age of plate at trench
              //    This uses the integral of the half-space temperature profile
              double initial_heat_content = 2 * thermal_conductivity * (surface_temperature - potential_mantle_temperature) *
                                            std::sqrt(plate_age_sec / (thermal_diffusivity * Consts::PI));

              // Plate age increases with distance along the slab in the mantle
              double effective_plate_age = plate_age_sec + (distance_along_plane / plate_velocity) * seconds_in_year; // m/(m/y) = y(seconds_in_year)


              // Need adiabatic temperature at position of grid point
              const double background_temperature = adiabatic_heating ? potential_mantle_temperature *
                                                    std::exp(thermal_expansion_coefficient * gravity_norm * depth / specific_heat)
                                                    : potential_mantle_temperature;

              WBAssert(!std::isnan(background_temperature), "Internal error: temp is not a number: " << background_temperature << ". In exponent: "
                       << std::exp(((thermal_expansion_coefficient * gravity_norm) / specific_heat) * depth)
                       << ", thermal_expansion_coefficient = " << thermal_expansion_coefficient << ", gravity_norm = " << gravity_norm
                       << ", specific_heat = " << specific_heat << ", depth = " << depth);

              const double adiabatic_gradient = adiabatic_heating ? background_temperature - potential_mantle_temperature : 0;

              //  2. Get Tmin and offset as a function of depth: these depend on plate velocity and plate age_at_trench.
              //     shallow-dipping slabs will take longer to reach the same depth - this leads to larger effective age at a given depth
              //     causing these slabs to be broader than steeper dipping slabs
              //     These equations are empirical based on fitting the temperature profiles from dynamic subduction models.
              //     and published kinematic models for specific subduction zones.

              // increases Tmin slope for slower relative to slope for maximum plate velocity
              // will be between 0.1 (fast) and 0.35 (slow)
              const double max_plate_vel = 20/cm2m;  // e.g., 20 cm/yr -> 0.2 m/yr
              const double vsubfact = std::min( std::max( 0.35 + ((0.1-0.35) / max_plate_vel) * plate_velocity, 0.1), 0.35);

              // increases Tmin slope for younger plate relative to slope for old place
              // will be between 0.1 (old) and 0.35 *(young)
              const double max_plate_age = 100*my; // years
              const double max_age_fact = 1.0;
              const double agefact = std::min( std::max( max_age_fact + ((0.1-max_age_fact) / max_plate_age) * age_at_trench, 0.1), max_age_fact);

              // increases offset  relative to slab surface for older plates
              // will be between 0.1 (young) and 0.35 (old)
              const double agefact2 = std::max( std::min( 0.1 + ((0.35-0.1) / max_plate_age) * age_at_trench, 0.35), 0.1);

              // ranges between 0.5 (old and fast) and 1.0 (young and slow)
              const double subfact = 0.3 + vsubfact + agefact;  // this has a minimum value of 0.5 for the erfc
              // ranges between 0.5 (young and fast) and 1 (old and slow)
              const double subfact2 = 0.3 + vsubfact + agefact2;

              // Minimum Temperature
              const double min_Tcoup = 10; // 0.3*mantle_coupling_depth/km2m; // deg C: from 0.3deg/km * 100 km = 30 deg
              const double max_Tcoup = 350; // deg C: from Cascadia subduction models by Gao and Wang, Gcubd 2017
              const double Tcoup = min_Tcoup + (subfact - 0.5)*(max_Tcoup - min_Tcoup);

              const double min_Tmin660 = 300; //110; // based on syracruse models (with 0.5 deg/km adiabat; (265-30)-0.5*240)
              const double max_Tmin660 = 900;
              const double Tmin660 =  min_Tmin660 + (subfact - 0.5)*(max_Tmin660 - min_Tmin660);

              // Temperature offset
              const double min_offset_coup = 2*km2m;  // values from syracuse and time-dep Citcom run
              const double max_offset_coup = 10*km2m;
              const double offset_coup = min_offset_coup + (subfact2)*(max_offset_coup - min_offset_coup);

              const double min_offset660 = 15*km2m;  // values from time-dep Citcom run
              const double max_offset660 = 25*km2m;
              const double offset660 =  min_offset660 + (subfact2)*(max_offset660 - min_offset660);

              // For tapering the slab tip temperature to the mantle temperature
              const double start_taper_distance = total_segment_length - taper_distance;

              const double upper_mantle_lengthscale = 660e3 - mantle_coupling_depth; // m


              const double taper_con = 0.8;  // controls how close taper gets to end of segment
              double theta = 0;
              double min_temperature = 0;
              double offset = 0;

              if (depth_to_reference_surface < mantle_coupling_depth)
                {
                  // above coupling depth
                  theta = (mantle_coupling_depth - depth_to_reference_surface) / (subfact * mantle_coupling_depth); // must be scaled depth_coupling to match Tcoup ad Tsurface.
                  min_temperature = Tcoup * std::erfc(theta);
                  offset = offset_coup * std::erfc(theta);
                }
              else if (distance_along_plane >= start_taper_distance)
                {
                  // beyond start taper distance to taper the slab tip

                  const double depth_start_taper = depth_to_reference_surface - (distance_along_plane - start_taper_distance)
                                                   * std::sin(average_angle * Consts::PI / 180.0);
                  const double theta_start = (mantle_coupling_depth - depth_start_taper) / (subfact * upper_mantle_lengthscale);
                  const double Tmin_start_taper = Tcoup + Tmin660 * (std::erfc(theta_start)) - Tmin660;
                  //keep the offset location constant in the taper
                  const double offset_start_taper = offset_coup + (offset660) * std::erfc(theta_start) - offset660;

                  // taper Tmin to the mantle temperature
                  theta = (distance_along_plane - start_taper_distance) / (taper_distance);
                  min_temperature =  Tmin_start_taper + (potential_mantle_temperature - Tmin_start_taper) * (1-std::erfc(taper_con*theta));
                  offset =  offset_start_taper + (2*max_offset660 - offset_start_taper) * (1-std::erfc(taper_con*theta));

                  // Also taper the initial heat content and effective plate age
                  initial_heat_content = initial_heat_content * std::erfc(1.5*taper_con*theta);
                  effective_plate_age = effective_plate_age * std::erfc(1.5*taper_con*theta);
                }

              else
                {
                  theta = (mantle_coupling_depth - depth_to_reference_surface) / (subfact * upper_mantle_lengthscale);
                  min_temperature = Tcoup + Tmin660 * std::erfc(theta) - Tmin660;
                  offset = offset_coup + offset660 * std::erfc(theta) - offset660;
                }

              min_temperature = min_temperature + adiabatic_gradient + surface_temperature;

              // Adjust distance for the offset of the minimum temperature from the top of the slab
              const double adjusted_distance = distance_from_plane - offset;

              // Base value chosen to insure that even young slabs have some cooling of the top of the slab
              // Use forearc_cooling_factor input variable to increase the amount of long-term cooling of the forearc.
              const double max_top_heat_content = -1.0e9 * forearc_cooling_factor * background_temperature; // need to multiply by background temperature  since it gets divided by something of this scale.

              double temperature = 0.0;  // temperature (to be determined) at this location in the mesh

              if (min_temperature < background_temperature)
                {

                  // 3. Determine the heat content for side 1 (bottom) of the slab
                  // Comes from integrating the half-space cooling model temperature

                  const double bottom_heat_content = 2 * thermal_conductivity * (min_temperature - potential_mantle_temperature) *
                                                     std::sqrt(effective_plate_age /(thermal_diffusivity * Consts::PI));

                  // 4. The difference in heat content goes into the temperature above where Tmin occurs.
                  // Should not be a positive value
                  //double top_heat_content = initial_heat_content - bottom_heat_content;

                  double top_heat_content =  std::min(max_top_heat_content, (initial_heat_content - bottom_heat_content) );

                  // Also need to taper the top_heat_content otherwise slab top will continue to thicken to the tip.
                  // Can't do this above because need min_temperature to get bottom_heat_content first
                  if (distance_along_plane > start_taper_distance)
                    {
                      top_heat_content = top_heat_content * std::erfc(taper_con*theta);
                    }

                  // Assign the temperature depending on whether distance is negative (above) or positive (below) the slab
                  if (adjusted_distance < 0)
                    {
                      // use 1D infinite space solution for top (side 2) of slab the slab
                      // 2 times the "top_heat_content" because all this heat needs to be on one side of the Gaussian
                      const double time_top_slab = (1/(Consts::PI*thermal_diffusivity)) *
                                                   pow(((2 * top_heat_content) / (2 * density * specific_heat * (min_temperature - temperature_ + 1e-16))),2) + 1e-16;

                      // for overriding plate region where plate temperature is less the minimum slab temperature
                      // need to set temperature = temperature_ otherwise end up with temperature less than surface temperature ;
                      if (temperature_ < min_temperature)
                        {
                          temperature = temperature_;
                        }
                      else
                        {
                          temperature  = temperature_ + (2 * top_heat_content /
                                                         (2*density*specific_heat * std::sqrt(Consts::PI * thermal_diffusivity * time_top_slab)))*
                                         std::exp(-(adjusted_distance*adjusted_distance)/(4*thermal_diffusivity*time_top_slab));
                        }
                    }
                  else
                    {
                      // use half-space cooling model for the bottom (side 1) of the slab
                      temperature = background_temperature + (min_temperature - background_temperature) *
                                    std::erfc(adjusted_distance / (2 * std::sqrt(thermal_diffusivity * effective_plate_age)));
                    }
                }
              else
                {
                  // slab temperature anomaly is gone.
                  temperature = temperature_;
                }

              WBAssert(!std::isnan(temperature), "Internal error: temperature is not a number: " << temperature << '.');
              WBAssert(std::isfinite(temperature), "Internal error: temperature is not finite: " << temperature << '.');

              return apply_operation(operation, temperature_, temperature);
            }

          return temperature_;
        }

        WB_REGISTER_FEATURE_SUBDUCTING_PLATE_TEMPERATURE_MODEL(MassConserving, mass conserving)
      } // namespace Temperature
    } // namespace SubductingPlateModels
  } // namespace Features
} // namespace WorldBuilder
