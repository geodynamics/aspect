/*
  Copyright (C) 2023 - 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/property/cpo_elastic_tensor.h>
#include <aspect/particle/manager.h>

#include <aspect/utilities.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {


      template <int dim>
      CpoElasticTensor<dim>::CpoElasticTensor ()
      {
        // The following values are directly from D-Rex.
        // Todo: make them a input parameter
        // Stiffness matrix for Olivine (GigaPascals)
        stiffness_matrix_olivine[0][0] = 320.71;
        stiffness_matrix_olivine[0][1] = 69.84;
        stiffness_matrix_olivine[0][2] = 71.22;
        stiffness_matrix_olivine[1][1] = 197.25;
        stiffness_matrix_olivine[1][2] = 74.8;
        stiffness_matrix_olivine[2][2] = 234.32;
        stiffness_matrix_olivine[3][3] = 63.77;
        stiffness_matrix_olivine[4][4] = 77.67;
        stiffness_matrix_olivine[5][5] = 78.36;


        // Stiffness matrix for Enstatite (GPa)
        stiffness_matrix_enstatite[0][0] = 236.9;
        stiffness_matrix_enstatite[0][1] = 79.6;
        stiffness_matrix_enstatite[0][2] = 63.2;
        stiffness_matrix_enstatite[1][1] = 180.5;
        stiffness_matrix_enstatite[1][2] = 56.8;
        stiffness_matrix_enstatite[2][2] = 230.4;
        stiffness_matrix_enstatite[3][3] = 84.3;
        stiffness_matrix_enstatite[4][4] = 79.4;
        stiffness_matrix_enstatite[5][5] = 80.1;
      }



      template <int dim>
      void
      CpoElasticTensor<dim>::initialize ()
      {
        const auto &manager = this->get_particle_manager(this->get_particle_manager_index()).get_property_manager();
        AssertThrow(manager.plugin_name_exists("crystal preferred orientation"),
                    ExcMessage("No crystal preferred orientation property plugin found."));

        AssertThrow(manager.check_plugin_order("crystal preferred orientation","cpo elastic tensor"),
                    ExcMessage("To use the cpo elastic tensor plugin, the cpo plugin needs to be defined before this plugin."));

        cpo_data_position = manager.get_data_info().get_position_by_plugin_index(manager.get_plugin_index_by_name("crystal preferred orientation"));
      }



      template <int dim>
      SymmetricTensor<2,6>
      CpoElasticTensor<dim>::voigt_average_elastic_tensor (const Particle::Property::CrystalPreferredOrientation<dim> &cpo_particle_property,
                                                           const unsigned int cpo_data_position,
                                                           const ArrayView<double> &data) const
      {
        SymmetricTensor<2,6> C_average;
        const SymmetricTensor<2,6> *stiffness_matrix = &stiffness_matrix_olivine;
        for (size_t mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
          {
            if (cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == DeformationType::olivine_a_fabric
                || cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == DeformationType::olivine_b_fabric
                || cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == DeformationType::olivine_d_fabric
                || cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == DeformationType::olivine_c_fabric
                || cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == DeformationType::olivine_e_fabric
               )
              {
                stiffness_matrix = &stiffness_matrix_olivine;
              }
            else if (cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == DeformationType::enstatite)
              {
                stiffness_matrix = &stiffness_matrix_enstatite;
              }
            else
              {
                AssertThrow(false, ExcMessage("Stiffness matrix not implemented for deformation type "
                                              + std::to_string(static_cast<unsigned int>(cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i)))));
              }

            for (size_t grain_i = 0; grain_i < n_grains; grain_i++)
              {
                const auto rotated_matrix = Utilities::Tensors::rotate_voigt_stiffness_matrix(transpose(cpo_particle_property.get_rotation_matrix_grains(cpo_data_position,data,mineral_i,grain_i)),*stiffness_matrix);
                C_average += cpo_particle_property.get_volume_fractions_grains(cpo_data_position,data,mineral_i,grain_i) * cpo_particle_property.get_volume_fraction_mineral(cpo_data_position,data,mineral_i) * rotated_matrix;
              }
          }

        return C_average;
      }



      template <int dim>
      void
      CpoElasticTensor<dim>::initialize_one_particle_property(const Point<dim> &,
                                                              std::vector<double> &data) const
      {

        // At initialization, the deformation type for cpo is initialized to -1.
        // Initialize with the stiffness matrix of olivine to avoid errors in the computation.

        for (unsigned int i = 0; i < SymmetricTensor<2,6>::n_independent_components ; ++i)
          {
            data.push_back(stiffness_matrix_olivine[SymmetricTensor<2,6>::unrolled_to_component_indices(i)]);
          }


      }



      template <int dim>
      void
      CpoElasticTensor<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &/*inputs*/,
                                                        typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        // Get a reference to the CPO particle property.
        const Particle::Property::CrystalPreferredOrientation<dim> &cpo_particle_property =
          this->get_particle_manager(this->get_particle_manager_index()).get_property_manager().template get_matching_active_plugin<Particle::Property::CrystalPreferredOrientation<dim>>();

        for (auto &particle: particles)
          {
            const SymmetricTensor<2,6> C_average = voigt_average_elastic_tensor(cpo_particle_property,
                                                                                cpo_data_position,
                                                                                particle.get_properties());

            Particle::Property::CpoElasticTensor<dim>::set_elastic_tensor(this->data_position,
                                                                          particle.get_properties(),
                                                                          C_average);
          }
      }



      template <int dim>
      SymmetricTensor<2,6>
      CpoElasticTensor<dim>::get_elastic_tensor(unsigned int cpo_data_position,
                                                const ArrayView<double> &data)
      {
        return Utilities::Tensors::to_symmetric_tensor<6>(&data[cpo_data_position],
                                                          &data[cpo_data_position]+SymmetricTensor<2,6>::n_independent_components);
      }



      template <int dim>
      void
      CpoElasticTensor<dim>::set_elastic_tensor(unsigned int cpo_data_position,
                                                const ArrayView<double> &data,
                                                const SymmetricTensor<2,6> &elastic_tensor)
      {
        Utilities::Tensors::unroll_symmetric_tensor_into_array(elastic_tensor,
                                                               &data[cpo_data_position],
                                                               &data[cpo_data_position]+SymmetricTensor<2,6>::n_independent_components);
      }



      template <int dim>
      UpdateTimeFlags
      CpoElasticTensor<dim>::need_update() const
      {
        return update_output_step;
      }



      template <int dim>
      UpdateFlags
      CpoElasticTensor<dim>::get_update_flags (const unsigned int /*component*/) const
      {
        return update_default;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      CpoElasticTensor<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int>> property_information;

        property_information.emplace_back("cpo_elastic_tensor", SymmetricTensor<2,6>::n_independent_components);
        return property_information;
      }



      template <int dim>
      void
      CpoElasticTensor<dim>::declare_parameters (ParameterHandler &)
      {}



      template <int dim>
      void
      CpoElasticTensor<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Crystal Preferred Orientation");
        {
          n_grains = prm.get_integer("Number of grains per particle");
          prm.enter_subsection("Initial grains");
          {
            n_minerals = dealii::Utilities::split_string_list(prm.get("Minerals")).size();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(CpoElasticTensor,
                                        "cpo elastic tensor",
                                        "A plugin in which the particle property tensor is "
                                        "defined as the Voigt average of the "
                                        "elastic tensors of the minerals in the textured rock."
                                        "Currently only Olivine and Enstatite are supported.")
    }
  }
}
