/*
  Copyright (C) 2018 - 2024 by the authors of the World Builder code.

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
*/

#include "world_builder/features/subducting_plate_models/composition/tian2019_water_content.h"

#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/unsigned_int.h"
#include "world_builder/utilities.h"
#include "world_builder/world.h"


namespace WorldBuilder
{

  using namespace Utilities;

  namespace Features
  {
    namespace SubductingPlateModels
    {
      namespace Composition
      {
        TianWaterContent::TianWaterContent(WorldBuilder::World *world_)
          :
          min_depth(NaN::DSNAN),
          max_depth(NaN::DSNAN),
          density(NaN::DSNAN)
        {
          this->world = world_;
          this->name = "water content";
        }

        TianWaterContent::~TianWaterContent()
          = default;

        void
        TianWaterContent::declare_entries(Parameters &prm, const std::string & /*unused*/)
        {
          // Document plugin and require entries if needed.
          // Add compositions to the required parameters.
          prm.declare_entry("", Types::Object({"compositions"}),
                            "TianWaterContent compositional model. Sets bound water content as a compositional field. The returned "
                            "water content is based on the the temperature and pressure at a point within the world. Currently, "
                            "the bound water content can be determined for four different lithologies: 'sediment', mid-ocean "
                            "ridge basalt ('MORB'), 'gabbro', and 'peridotite', using parameterized phase diagrams from Tian et al., 2019 "
                            "(https://doi.org/10.1029/2019GC008488). The pressure is lithostatic, calculated with a constant user defined "
                            "density, and is limited by a user defined cutoff pressure (in GPa) for each lithology. This is required because the "
                            "parameterization breaks down at large pressures. Recommended cutoff pressures are 10 GPa is used for 'peridotite', "
                            "26 GPa is used for 'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.");

          // Declare entries of this plugin
          prm.declare_entry("min distance slab top", Types::Double(0),
                            "todo The depth in meters from which the composition of this feature is present.");
          prm.declare_entry("max distance slab top", Types::Double(std::numeric_limits<double>::max()),
                            "todo The depth in meters to which the composition of this feature is present.");
          prm.declare_entry("density", Types::Double(3000.0),
                            "The reference density used for determining the lithostatic pressure for calculating "
                            "the bound water content.");
          prm.declare_entry("compositions", Types::Array(Types::UnsignedInt(),0),
                            "A list with the labels of the composition which are present there.");
          prm.declare_entry("lithology",  Types::String("peridotite"),
                            "The lithology used to determine which polynomials to use for calculating the water content. Valid options are: "
                            "'sediment', 'MORB', 'gabbro', and 'peridotite'.");
          prm.declare_entry("initial water content", Types::Double(5),
                            "The value of the initial water content (in wt%) for the lithology at the trench. This represents the "
                            "max value applied to this lithology.");
          prm.declare_entry("cutoff pressure", Types::Double(10),
                            "The upper bound for the pressure, in GPa, for the specified lithology in the Tian parameterization. This is necessary because "
                            "the parameterization breaks down for high pressures. It is recommended that 10 GPa is used for 'peridotite', 26 GPa is used for "
                            "'gabbro', 16 GPa is used for 'MORB', and 1 GPa is used for 'sediment'.");
          prm.declare_entry("operation", Types::String("replace", std::vector<std::string> {"replace", "replace defined only", "add", "subtract"}),
                            "Whether the value should replace any value previously defined at this location (replace) or "
                            "add the value to the previously define value. Replacing implies that all compositions not "
                            "explicitly defined are set to zero. To only replace the defined compositions use the replace only defined option.");

        }

        void
        TianWaterContent::parse_entries(Parameters &prm)
        {
          min_depth = prm.get<double>("min distance slab top");
          max_depth = prm.get<double>("max distance slab top");
          density = prm.get<double>("density");
          compositions = prm.get_vector<unsigned int>("compositions");
          max_water_content = prm.get<double>("initial water content");
          cutoff_pressure = prm.get<double>("cutoff pressure");
          operation = string_operations_to_enum(prm.get<std::string>("operation"));
          std::string lithology_string = prm.get<std::string>("lithology");

          if (lithology_string=="peridotite")
            lithology_type = peridotite;
          else if (lithology_string=="gabbro")
            lithology_type = gabbro;
          else if (lithology_string=="MORB")
            lithology_type = MORB;
          else if (lithology_string=="sediment")
            lithology_type = sediment;
        }


        double
        TianWaterContent::calculate_water_content(double pressure,
                                                  double temperature) const
        {
          double ln_LR_value = 0;
          double ln_c_sat_value = 0;
          double Td_value = 0;
          std::vector<double> LR_polynomial_coeffs;
          std::vector<double> c_sat_polynomial_coeffs;
          std::vector<double> Td_polynomial_coeffs;

          // Calculate the c_sat value from Tian et al., 2019
          if (lithology_type == sediment)
            {
              for (unsigned int c_sat_index = 0; c_sat_index < c_sat_poly[lithology_type].size(); ++c_sat_index)
                ln_c_sat_value += c_sat_poly[lithology_type][c_sat_index] * (std::pow(std::log10(pressure), c_sat_poly[lithology_type].size() - 1 - c_sat_index));
            }
          else
            {
              for (unsigned int c_sat_index = 0; c_sat_index < c_sat_poly[lithology_type].size(); ++c_sat_index)
                ln_c_sat_value += c_sat_poly[lithology_type][c_sat_index] * (std::pow(pressure, c_sat_poly[lithology_type].size() - 1 - c_sat_index));
            }

          // Calculate the LR value from Tian et al., 2019
          for (unsigned int LR_coeff_index = 0; LR_coeff_index < LR_poly[lithology_type].size(); ++LR_coeff_index)
            ln_LR_value += LR_poly[lithology_type][LR_coeff_index] * (std::pow(1/pressure, LR_poly[lithology_type].size() - 1 - LR_coeff_index));

          // Calculate the Td value from Tian et al., 2019
          for (unsigned int Td_coeff_index = 0; Td_coeff_index < Td_poly[lithology_type].size(); ++Td_coeff_index)
            Td_value += Td_poly[lithology_type][Td_coeff_index] * (std::pow(pressure, Td_poly[lithology_type].size() - 1 - Td_coeff_index));

          double partition_coeff = std::exp(ln_c_sat_value) * std::exp(std::exp(ln_LR_value) * (1/temperature - 1/Td_value));
          return partition_coeff;
        }


        double
        TianWaterContent::get_composition(const Point<3> &position_in_cartesian_coordinates,
                                          const double depth,
                                          const unsigned int composition_number,
                                          double composition,
                                          const double  /*feature_min_depth*/,
                                          const double  /*feature_max_depth*/,
                                          const WorldBuilder::Utilities::PointDistanceFromCurvedPlanes &distance_from_plane,
                                          const AdditionalParameters & /*additional_parameters*/) const
        {
          if (distance_from_plane.distance_from_plane <= max_depth && distance_from_plane.distance_from_plane >= min_depth)
            {
              // The polynomials break down for pressures less than 0.5 GPa, and for pressures above a user defined cutoff pressure
              // ensure that the pressure is never below 0.5 GPa
              const double lithostatic_pressure = std::max(0.5, std::min(density * 9.81 * depth / 1e9, cutoff_pressure)); // GPa
              const double slab_temperature = world->properties(position_in_cartesian_coordinates.get_array(), depth, {{{1,0,0}}})[0];
              double partition_coefficient = calculate_water_content(lithostatic_pressure,
                                                                     slab_temperature);

              partition_coefficient = std::min(max_water_content, partition_coefficient);

              for (unsigned int i = 0; i < compositions.size(); ++i)
                {
                  if (compositions[i] == composition_number)
                    {
                      return apply_operation(operation,composition, partition_coefficient);
                    }
                }

              if (operation == Operations::REPLACE)
                {
                  return 0.0;
                }
            }
          return composition;
        }
        WB_REGISTER_FEATURE_SUBDUCTING_PLATE_COMPOSITION_MODEL(TianWaterContent, tian water content)
      } // namespace Composition
    } // namespace SubductingPlateModels
  } // namespace Features
} // namespace WorldBuilder


