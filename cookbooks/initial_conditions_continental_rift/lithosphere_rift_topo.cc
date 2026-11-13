/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#include "lithosphere_rift_topo.h"
#include "lithosphere_rift_IC.h"
#include <aspect/geometry_model/box.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/utilities.h>
#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace InitialTopographyModel
  {
    template <int dim>
    void
    LithosphereRift<dim>::
    initialize ()
    {
      // Compute the maximum topography amplitude based on isostasy.
      // Assume the reference density is representative for each layer (despite temperature dependence)

      // For now, we assume a 3-layer system with an upper crust, lower crust and lithospheric mantle
      const unsigned int id_upper = this->introspection().compositional_index_for_name("upper");
      const unsigned int id_lower = this->introspection().compositional_index_for_name("lower");
      const unsigned int id_mantle_L = this->introspection().compositional_index_for_name("mantle_L");

      // Assemble a list of phase densities for each composition.
      // Add 1 for background material.
      std::vector<std::vector<double>> densities_per_composition(this->n_compositional_fields()+1);
      unsigned int counter = 0;
      for (unsigned int i = 0; i < (*n_phases_for_each_composition).size(); ++i)
        {
          for (unsigned int j = 0; j < (*n_phases_for_each_composition)[i]; ++j)
            {
              densities_per_composition[i].push_back(temp_densities[counter]);
              ++counter;
            }
        }

      // Get the relevant densities for the lithosphere.
      // We take the reference density of the first phase.
      densities.push_back(densities_per_composition[0][0]);
      densities.push_back(densities_per_composition[id_upper+1][0]);
      densities.push_back(densities_per_composition[id_lower+1][0]);
      densities.push_back(densities_per_composition[id_mantle_L+1][0]);

      // The reference column
      ref_rgh = 0;
      // Assume constant gravity magnitude, so ignore
      for (unsigned int l=0; l<3; ++l)
        ref_rgh += densities[l+1] * reference_thicknesses[l];

      // The total lithosphere thickness
      const double sum_thicknesses = std::accumulate(reference_thicknesses.begin(), reference_thicknesses.end(), 0);

      // The column at the rift center
      double rift_rgh = 0;
      rift_thicknesses = reference_thicknesses;
      for (unsigned int l=0; l<rift_thicknesses.size(); ++l)
        rift_thicknesses[l] *= (1.-A_rift[l]);

      for (unsigned int l=0; l<3; ++l)
        rift_rgh += densities[l+1] * rift_thicknesses[l];

      // The total lithosphere thickness at the rift
      const double sum_rift_thicknesses = std::accumulate(rift_thicknesses.begin(), rift_thicknesses.end(),0);

      // The column at the polygon center
      const unsigned int n_polygons = polygon_thicknesses.size();
      std::vector<double> polygon_rgh(n_polygons);
      std::vector<double> sum_polygon_thicknesses(n_polygons);
      double max_sum_polygon_thicknesses = 0.;
      for (unsigned int i_polygons=0; i_polygons<n_polygons; ++i_polygons)
        {
          for (unsigned int l=0; l<3; ++l)
            {
              polygon_rgh[i_polygons] += densities[l+1] * polygon_thicknesses[i_polygons][l];
            }
          // The total lithosphere thickness
          sum_polygon_thicknesses[i_polygons] = std::accumulate(polygon_thicknesses[i_polygons].begin(), polygon_thicknesses[i_polygons].end(),0);
          max_sum_polygon_thicknesses = std::max(max_sum_polygon_thicknesses, sum_polygon_thicknesses[i_polygons]);

        }

      // Make sure the compensation depth is in the sublithospheric mantle
      compensation_depth = std::max(sum_thicknesses, std::max(sum_rift_thicknesses, max_sum_polygon_thicknesses)) + 5e3;

      // Add sublithospheric mantle part to the columns
      ref_rgh += (compensation_depth - sum_thicknesses) * densities[0];
      rift_rgh += (compensation_depth - sum_rift_thicknesses) * densities[0];
      for (unsigned int i_polygons = 0; i_polygons < n_polygons; ++i_polygons)
        {
          polygon_rgh[i_polygons] += (compensation_depth - sum_polygon_thicknesses[i_polygons]) * densities[0];
        }

      // Compute the maximum topography based on mass surplus/deficit
      topo_rift_amplitude = (ref_rgh - rift_rgh) / densities[0];
      for (unsigned int i_polygons = 0; i_polygons < n_polygons; ++i_polygons)
        topo_polygon_amplitude = std::max((ref_rgh - polygon_rgh[i_polygons]) / densities[0], topo_polygon_amplitude);

      // TODO: probably there are combinations of rift and polygon topography
      // that result in a higher topography
      maximum_topography = std::max(topo_rift_amplitude, topo_polygon_amplitude);

      this->get_pcout() << "   Maximum initial topography of rift: " << topo_rift_amplitude << " m" << std::endl;
      this->get_pcout() << "   Maximum initial topography of polygon: " << topo_polygon_amplitude << " m" << std::endl;
    }


    template <int dim>
    double
    LithosphereRift<dim>::
    value (const Point<dim-1> &position) const
    {
      // The simulator only keeps the initial conditions around for
      // the first time step. As a consequence, we have to save a
      // shared pointer to that object ourselves the first time we get
      // here.
      if (initial_composition_manager == nullptr)
        const_cast<std::shared_ptr<const aspect::InitialComposition::Manager<dim>>&>(initial_composition_manager) = this->get_initial_composition_manager_pointer();

      // Check that the required initial composition model is used
      // We have to do it here instead of in initialize() because
      // the names are not available upon initialization of the
      // initial topography model yet.
      const std::vector<std::string> active_initial_composition_models = initial_composition_manager->get_active_initial_composition_names();
      AssertThrow(initial_composition_manager->template has_matching_initial_composition_model<const InitialComposition::LithosphereRift<dim>>(),
                  ExcMessage("The 'lithosphere with rift' initial topography plugin requires the 'lithosphere with rift' initial composition plugin."));

      // When cartesian, position contains x(,y); when spherical, position contains lon(,lat) (in degrees).
      // Turn into a Point<dim-1>
      Point<dim-1> surface_position;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_position[d] = position[d];

      // Get the initial composition plugin
      const InitialComposition::LithosphereRift<dim> &ic = initial_composition_manager->template get_matching_initial_composition_model<const InitialComposition::LithosphereRift<dim>>();

      // Compute the topography based on distance to the rift and distance to the polygon
      std::vector<double> local_thicknesses = ic.compute_local_thicknesses(surface_position);

      // The local lithospheric column
      double local_rgh = 0;
      for (unsigned int l=0; l<3; ++l)
        local_rgh += densities[l+1] * local_thicknesses[l];
      // The total local lithosphere thickness
      const double sum_local_thicknesses = std::accumulate(local_thicknesses.begin(), local_thicknesses.end(),0);
      local_rgh += (compensation_depth - sum_local_thicknesses) * densities[0];

      return (ref_rgh - local_rgh) / densities[0];
    }


    template <int dim>
    double
    LithosphereRift<dim>::
    max_topography () const
    {
      return maximum_topography;
    }


    template <int dim>
    void
    LithosphereRift<dim>::
    declare_parameters (ParameterHandler &)
    {
    }



    template <int dim>
    void
    LithosphereRift<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial composition model");
      {
        prm.enter_subsection("Lithosphere with rift");
        {
          A_rift = Utilities::possibly_extend_from_1_to_N(Utilities::string_to_double(Utilities::split_string_list(prm.get("Amplitude of Gaussian rift geometry"))),
                                                          3,
                                                          "Amplitude of Gaussian rift geometry");
          reference_thicknesses = Utilities::possibly_extend_from_1_to_N(Utilities::string_to_double(Utilities::split_string_list(prm.get("Layer thicknesses"))),
                                                                         3,
                                                                         "Layer thicknesses");
          // Split the string into the separate polygons
          const std::vector<std::string> temp_thicknesses = Utilities::split_string_list(prm.get("Lithospheric polygon layer thicknesses"),';');
          const unsigned int n_polygons = temp_thicknesses.size();
          polygon_thicknesses.resize(n_polygons);
          for (unsigned int i_polygons = 0; i_polygons < n_polygons; ++i_polygons)
            {
              polygon_thicknesses[i_polygons] = Utilities::string_to_double(Utilities::split_string_list(temp_thicknesses[i_polygons],','));
              AssertThrow(polygon_thicknesses[i_polygons].size()==3, ExcMessage ("The number of layer thicknesses should be equal to 3."));
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection ("Compositional fields");
      {
        list_of_composition_names = Utilities::split_string_list (prm.get("Names of fields"));
      }
      prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Visco Plastic");
        {
          n_phases_for_each_composition = std::make_unique<std::vector<unsigned int>>();
          temp_densities = Utilities::parse_map_to_double_array (prm.get("Densities"),
                                                                 list_of_composition_names,
                                                                 /*has_background_field=*/true,
                                                                 "Densities",
                                                                 true,
                                                                 n_phases_for_each_composition);
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
  namespace InitialTopographyModel
  {
    ASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(LithosphereRift,
                                             "lithosphere with rift",
                                             "An initial topography model that defines the initial topography "
                                             "based on isostasy. It takes into account lithospheric thickness "
                                             "variations as specified in the InitialComposition model lithosphere "
                                             "with rift' as should only be used in conjunction with this model. ")
  }
}
