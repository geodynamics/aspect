/*
  Copyright (C) 2018-2024 by the authors of the World Builder code.

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

#include "world_builder/world.h"


#include "world_builder/config.h"
#include "world_builder/features/subducting_plate.h"
#include "world_builder/gravity_model/interface.h"
#include "world_builder/nan.h"
#include "world_builder/types/array.h"
#include "world_builder/types/bool.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/plugin_system.h"
#include "world_builder/types/point.h"
#include "world_builder/types/int.h"

#include <iostream>
#include <world_builder/coordinate_system.h>
#include <world_builder/objects/distance_from_surface.h>

#ifdef WB_WITH_MPI
// we don't need the c++ MPI wrappers
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif

#ifndef NDEBUG
#ifdef WB_USE_FP_EXCEPTIONS
#include <cfenv>
#endif
#endif

namespace WorldBuilder
{
  using namespace Utilities;

  World::World(std::string filename, bool has_output_dir, const std::string &output_dir, unsigned long random_number_seed, const bool limit_debug_consistency_checks_)
    :
    parameters(*this),
    surface_coord_conversions(invalid),
    dim(NaN::ISNAN),
    random_number_engine(random_number_seed),
    limit_debug_consistency_checks(limit_debug_consistency_checks_)
  {

#ifndef NDEBUG
#ifdef WB_USE_FP_EXCEPTIONS
    // Some implementations seem to not initialize the floating point exception
    // bits to zero. Make sure we start from a clean state.
    feclearexcept(FE_DIVBYZERO|FE_INVALID);

    // enable floating point exceptions
    feenableexcept(FE_DIVBYZERO|FE_INVALID);
#endif
#endif

#ifdef WB_WITH_MPI
    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized == 0)
      {
        MPI_RANK = 0;
        MPI_SIZE = 1;
      }
    else
      {
        MPI_Comm_rank(MPI_COMM_WORLD, &MPI_RANK);
        MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE);
      }
#else
    MPI_RANK = 0;
    MPI_SIZE = 1;
#endif

    WorldBuilder::World::declare_entries(parameters);

    parameters.initialize(filename, has_output_dir, output_dir);

    this->parse_entries(parameters);
  }

  World::~World()
    = default;

  void World::declare_entries(Parameters &prm)
  {
    prm.enter_subsection("properties");
    {
      prm.declare_entry("", Types::Object({"version", "features"}), "Root object");

      prm.declare_entry("version", Types::String(""),"The major and minor version number for which the input file was written.");

      prm.declare_entry("$schema", Types::String(""),"The optional filename or https address to a JSON schema file");

      prm.declare_entry("cross section", Types::Array(Types::Point<2>(),2,2),"This is an array of two points along where the cross section is taken");

      prm.declare_entry("potential mantle temperature", Types::Double(1600),
                        "The potential temperature of the mantle at the surface in Kelvin.");
      prm.declare_entry("surface temperature", Types::Double(293.15),
                        "The temperature at the surface in Kelvin.");
      prm.declare_entry("force surface temperature", Types::Bool(false),
                        "Force the provided surface temperature to be set at the surface");
      prm.declare_entry("thermal expansion coefficient", Types::Double(3.5e-5),
                        "The thermal expansion coefficient in $K^{-1}$.");
      prm.declare_entry("specific heat", Types::Double(1250),
                        "The specific heat in $J kg^{-1} K^{-1}.$");
      prm.declare_entry("thermal diffusivity", Types::Double(0.804e-6),
                        "The thermal diffusivity in $m^{2} s^{-1}$.");

      prm.declare_entry("maximum distance between coordinates",Types::Double(0),
                        "This enforces a maximum distance (in degree for spherical coordinates "
                        "or meter in cartesian coordinates) between coordinates in the model. "
                        "If the distance is larger, extra points are added by interpolation. "
                        "Requires interpolation to be not 'none'.");

      prm.declare_entry("interpolation",Types::String("continuous monotone spline"),
                        "What type of interpolation should be used to enforce the minimum points per "
                        "distance parameter. Options are none, linear, monotone spline and "
                        "continuous monotone spline interpolation.");


      prm.declare_entry("coordinate system", Types::PluginSystem("cartesian", CoordinateSystems::Interface::declare_entries, {"model"}, false),"A coordinate system. Cartesian or spherical.");

      prm.declare_entry("gravity model", Types::PluginSystem("uniform", GravityModel::Interface::declare_entries, {"model"}, false),"A gravity model for the world.");

      prm.declare_entry("features", Types::PluginSystem("",Features::Interface::declare_entries, {"model"}),"A list of features.");

      prm.declare_entry("random number seed", Types::Int(-1),
                        "This allows the input of a preferred random number seed to generate random numbers."
                        " If no input is given, this value is -1 and triggers the use of default seed = 1.");

    }
    prm.leave_subsection();

  }


  void World::parse_entries(Parameters &prm)
  {
    using namespace rapidjson;

    /**
     * First load the major version number in the file and check the major
     * version number of the program.
     */

    WBAssertThrow((prm.get<std::string>("version") == Version::MAJOR + "." + Version::MINOR),
                  "The major and minor version combination for which is input file was written "
                  "is not the same as the version of the World Builder you are running. This means "
                  "That there may have been incompatible changes made between the versions. \n\n"
                  "Verify those changes and whether they affect your model. If this is not "
                  "the case, adjust the version number in the input file. \n\nThe provided version "
                  "number is \"" << prm.get<std::string>("version") << "\", while the used world builder "
                  "has  (major.minor) version \"" << Version::MAJOR << "." << Version::MINOR << "\". "
                  "If you created this file from scratch, fill set the version number to \"" <<
                  Version::MAJOR << "." << Version::MINOR << "\" to continue. If you got the world builder "
                  "file from somewhere, make sure that the output is what you expect it to be, because "
                  "backwards incompatible changes may have been made to the code.");

    /**
     * Secondly load the coordinate system parameters.
     */
    prm.coordinate_system = prm.get_unique_pointer<CoordinateSystems::Interface>("coordinate system");
    prm.coordinate_system->parse_entries(prm);

    /**
     * Thirdly load the gravity model parameters.
     */
    prm.gravity_model = prm.get_unique_pointer<GravityModel::Interface>("gravity model");

    prm.enter_subsection("gravity model");
    {
      prm.gravity_model->parse_entries(prm);
    }
    prm.leave_subsection();

    prm.get_unique_pointers<Features::Interface>("features",prm.features);

    const bool set_cross_section = prm.check_entry("cross section");

    const CoordinateSystem coordinate_system = prm.coordinate_system->natural_coordinate_system();

    if (set_cross_section)
      {
        dim = 2;
        const std::vector<Point<2> > cross_section_natural = prm.get_vector<Point<2> >("cross section");

        WBAssertThrow(cross_section_natural.size() == 2, "The cross section should contain two points, but it contains "
                      << cross_section.size() << " points.");

        for (const auto &it : cross_section_natural)
          cross_section.push_back(it *  (coordinate_system == spherical ? Consts::PI / 180.0 : 1.0));


        /**
         * pre-compute stuff for the cross section
         */
        surface_coord_conversions = cross_section[0]-cross_section[1];
        surface_coord_conversions *= -1/(surface_coord_conversions.norm());


      }
    else
      {
        dim = 3;
      }

    /**
     * Temperature parameters.
     */
    potential_mantle_temperature = prm.get<double>("potential mantle temperature");
    surface_temperature = prm.get<double>("surface temperature");
    surface_temperature = prm.get<double>("surface temperature");
    force_surface_temperature = prm.get<bool>("force surface temperature");
    thermal_expansion_coefficient = prm.get<double>("thermal expansion coefficient");
    specific_heat = prm.get<double>("specific heat");
    thermal_diffusivity = prm.get<double>("thermal diffusivity");

    /**
     * Model discretization parameters
     */
    maximum_distance_between_coordinates = prm.get<double>("maximum distance between coordinates");
    interpolation = prm.get<std::string>("interpolation");

    /**
     * Local random number seed parameter
    */
    const int local_seed = prm.get<int>("random number seed");

    if (local_seed>=0)
      random_number_engine.seed(static_cast<unsigned int>(local_seed+MPI_RANK));

    /**
     * Now load the features. Some features use for example temperature values,
     * so it is important that this is parsed the last.
     */
    prm.enter_subsection("features");
    {
      for (unsigned int i = 0; i < prm.features.size(); ++i)
        {
          prm.enter_subsection(std::to_string(i));
          {
            prm.features[i]->parse_entries(prm);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }



  std::vector<double>
  World::properties(const std::array<double, 2> &point,
                    const double depth,
                    const std::vector<std::array<unsigned int,3>> &properties) const
  {
    // turn it into a 3d coordinate and call the 3d temperature function
    WBAssertThrow(dim == 2, "This function can only be called when the cross section "
                  "variable in the world builder file has been set. Dim is "
                  << dim << '.');

    const CoordinateSystem coordinate_system = this->parameters.coordinate_system->natural_coordinate_system();

    Point<2> point_natural(point[0], point[1],coordinate_system);
    if (coordinate_system == spherical)
      {
        point_natural[1] = std::sqrt(point[0]*point[0]+point[1]*point[1]);
        point_natural[0] = std::atan2(point[1],point[0]);
      }

    Point<3> coord_3d(coordinate_system);
    if (coordinate_system == spherical)
      {
        coord_3d[0] = point_natural[1];
        coord_3d[1] = cross_section[0][0] + point_natural[0] * surface_coord_conversions[0];
        coord_3d[2] = cross_section[0][1] + point_natural[0] * surface_coord_conversions[1];
      }
    else
      {
        coord_3d[0] = cross_section[0][0] + point_natural[0] * surface_coord_conversions[0];
        coord_3d[1] = cross_section[0][1] + point_natural[0] * surface_coord_conversions[1];
        coord_3d[2] = point_natural[1];
      }


    const std::array<double, 3> point_3d_cartesian = this->parameters.coordinate_system->natural_to_cartesian_coordinates(coord_3d.get_array());

    return this->properties(point_3d_cartesian, depth, properties);
  }


  std::vector<double>
  World::properties(const std::array<double, 3> &point_,
                    const double depth,
                    const std::vector<std::array<unsigned int,3>> &properties) const
  {
    // We receive the cartesian points from the user.
    const Point<3> point(point_,cartesian);
    (void) this->limit_debug_consistency_checks;
    WBAssert(this->limit_debug_consistency_checks || this->parameters.coordinate_system->natural_coordinate_system() == cartesian
             || approx(depth, this->parameters.coordinate_system->max_model_depth()-sqrt(point_[0]*point_[0]+point_[1]*point_[1]+point_[2]*point_[2])),
             "Inconsistent input. Please check whether the radius in the spherical coordinates is consistent with the radius of the planet as defined "
             << "in the program that uses the Geodynamic World Builder. This is a debug check in GWB and can be disabled by setting "
             << "limit_debug_consistency_checks to true. "
             << "Depth = " << depth << ", radius = " << this->parameters.coordinate_system->max_model_depth()
             << ", point = " << point_[0] << " " << point_[1] << " " << point_[2]
             << ", radius-point.norm() = " << this->parameters.coordinate_system->max_model_depth()-sqrt(point_[0]*point_[0]+point_[1]*point_[1]+point_[2]*point_[2]));

    const Objects::NaturalCoordinate natural_coordinate = Objects::NaturalCoordinate(point,*(this->parameters.coordinate_system));

    // create output vector
    std::vector<double> output;
    std::vector<size_t> entry_in_output;
    std::vector<std::array<unsigned int,3>> properties_local;
    const double gravity_norm = this->parameters.gravity_model->gravity_norm(point);
    for (unsigned int i_property = 0; i_property < properties.size(); ++i_property)
      {
        switch (properties[i_property][0])
          {
            case 1: // Temperature
              if (std::fabs(depth) < 2.0 * std::numeric_limits<double>::epsilon() && force_surface_temperature)
                {
                  entry_in_output.emplace_back(output.size());
                  output.emplace_back(this->surface_temperature);
                  if (properties.size() == 1)
                    return output;
                }
              else
                {
                  entry_in_output.emplace_back(output.size());
                  output.emplace_back(potential_mantle_temperature * std::exp(((thermal_expansion_coefficient * gravity_norm) / specific_heat) * depth));
                }
              properties_local.emplace_back(properties[i_property]);
              break;
            case 2: // composition
              entry_in_output.emplace_back(output.size());
              output.emplace_back(0.);
              properties_local.emplace_back(properties[i_property]);
              break;
            case 3: // grains (10 entries per grain)
            {
              entry_in_output.emplace_back(output.size());
              const std::vector<double> tmp_vector(properties[i_property][2]*10,0.);
              output.insert(output.end(), tmp_vector.begin(), tmp_vector.end());
              properties_local.emplace_back(properties[i_property]);
              break;
            }
            case 4: // tag
            {
              entry_in_output.emplace_back(output.size());
              output.emplace_back(-1);
              properties_local.emplace_back(properties[i_property]);
              break;
            }
            default:
              WBAssertThrow(false,
                            "Internal error: Unimplemented property provided. " <<
                            "Only temperature (1), composition (2), grains (3) or tag (4) are allowed. "
                            "Provided property number was: " << properties[i_property][0]);
          }
      }
    for (auto &&it : parameters.features)
      {
        it->properties(point, natural_coordinate, depth, properties_local, gravity_norm, entry_in_output, output);
      }

    return output;
  }

  double
  World::temperature(const std::array<double,2> &point,
                     const double depth) const
  {
    return properties(point, depth, {{{1,0,0}}})[0];
  }

  double
  World::temperature(const std::array<double,2> &point,
                     const double depth,
                     const double /*gravity_norm*/) const
  {
    return properties(point, depth, {{{1,0,0}}})[0];
  }

  double
  World::temperature(const std::array<double,3> &point,
                     const double depth) const
  {
    return properties(point, depth, {{{1,0,0}}})[0];
  }

  double
  World::temperature(const std::array<double,3> &point,
                     const double depth,
                     const double /*gravity_norm*/) const
  {
    return properties(point, depth, {{{1,0,0}}})[0];
  }

  double
  World::composition(const std::array<double,2> &point,
                     const double depth,
                     const unsigned int composition_number) const
  {
    return properties(point, depth, {{{2,composition_number,0}}})[0];
  }

  double
  World::composition(const std::array<double,3> &point,
                     const double depth,
                     const unsigned int composition_number) const
  {
    return properties(point, depth, {{{2,composition_number,0}}})[0];
  }



  WorldBuilder::grains
  World::grains(const std::array<double,2> &point,
                const double depth,
                const unsigned int composition_number,
                size_t number_of_grains) const
  {
    return WorldBuilder::grains(properties(point, depth, {{{3,composition_number,static_cast<unsigned int>(number_of_grains)}}}),static_cast<unsigned int>(number_of_grains),0);
  }

  WorldBuilder::grains
  World::grains(const std::array<double,3> &point,
                const double depth,
                const unsigned int composition_number,
                size_t number_of_grains) const
  {
    return WorldBuilder::grains(properties(point, depth, {{{3,composition_number,static_cast<unsigned int>(number_of_grains)}}}),static_cast<unsigned int>(number_of_grains),0);
  }

  std::mt19937 &
  World::get_random_number_engine()
  {
    return random_number_engine;
  }

  Objects::PlaneDistances
  World::distance_to_plane(const std::array<double, 3> &point_,
                           const double depth,
                           const std::string &name) const
  {
    // We receive the cartesian points from the user.
    const Point<3> point(point_,cartesian);

    WBAssert(this->limit_debug_consistency_checks || this->parameters.coordinate_system->natural_coordinate_system() == cartesian
             || approx(depth, this->parameters.coordinate_system->max_model_depth()-sqrt(point_[0]*point_[0]+point_[1]*point_[1]+point_[2]*point_[2])),
             "Inconsistent input. Please check whether the radius in the spherical coordinates is consistent with the radius of the planet as defined "
             << "in the program that uses the Geodynamic World Builder. This is a debug check in GWB and can be disabled by setting "
             << "limit_debug_consistency_checks to true. "
             << "Depth = " << depth << ", radius = " << this->parameters.coordinate_system->max_model_depth()
             << ", point = " << point_[0] << " " << point_[1] << " " << point_[2]
             << ", radius-point.norm() = " << this->parameters.coordinate_system->max_model_depth()-sqrt(point_[0]*point_[0]+point_[1]*point_[1]+point_[2]*point_[2]));

    const Objects::NaturalCoordinate natural_coordinate = Objects::NaturalCoordinate(point,*(this->parameters.coordinate_system));

    Objects::PlaneDistances plane_distances(0.0, 0.0);
    for (auto &&it : this->parameters.features)
      {
        if (it->get_name() == name)
          {
            plane_distances = it->distance_to_feature_plane(point, natural_coordinate, depth);
            break;
          }
      }
    return plane_distances;
  }

} // namespace WorldBuilder

