/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

//#include <cstdlib>
#include <aspect/particle/property/cpo_elastic_tensor.h>
#include <aspect/particle/world.h>

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
        permutation_operator_3d[0][1][2]  = 1;
        permutation_operator_3d[1][2][0]  = 1;
        permutation_operator_3d[2][0][1]  = 1;
        permutation_operator_3d[0][2][1]  = -1;
        permutation_operator_3d[1][0][2]  = -1;
        permutation_operator_3d[2][1][0]  = -1;

        // The following values are directly form D-Rex.
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
        const auto &manager = this->get_particle_world().get_property_manager();
        AssertThrow(manager.plugin_name_exists("crystal preferred orientation"),
                    ExcMessage("No crystal preferred orientation property plugin found."));
        Assert(manager.plugin_name_exists("cpo elastic tensor"),
               ExcMessage("No s wave anisotropy property plugin found."));

        AssertThrow(manager.check_plugin_order("crystal preferred orientation","cpo elastic tensor"),
                    ExcMessage("To use the cpo elastic tensor plugin, the cpo plugin need to be defined before this plugin."));

        cpo_data_position = manager.get_data_info().get_position_by_plugin_index(manager.get_plugin_index_by_name("cpo"));
      }





      template <int dim>
      SymmetricTensor<2,6>
      CpoElasticTensor<dim>::voigt_average_elastic_tensor (const Particle::Property::CrystalPreferredOrientation<dim> &cpo_particle_property,
                                                           const unsigned int cpo_data_position,
                                                           const ArrayView<double> &data) const
      {
        SymmetricTensor<2,6,double> Sav;
        const SymmetricTensor<2,6> *stiffness_matrix = &stiffness_matrix_olivine;
        for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
          {
            if (cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == (unsigned int)DeformationTypeSelector::olivine_a_fabric
                || cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == (unsigned int)DeformationTypeSelector::olivine_b_fabric
                || cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == (unsigned int)DeformationTypeSelector::olivine_c_fabric
                || cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == (unsigned int)DeformationTypeSelector::olivine_d_fabric
                || cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == (unsigned int)DeformationTypeSelector::olivine_e_fabric
                || cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == (unsigned int)DeformationTypeSelector::olivine_karato_2008
               )
              {
                stiffness_matrix = &stiffness_matrix_olivine;
              }
            else if (cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i) == (unsigned int)DeformationTypeSelector::enstatite)
              {
                stiffness_matrix = &stiffness_matrix_enstatite;
              }
            else
              {
                AssertThrow(false, ExcMessage("Stiffness matrix not implemented for deformation type " + std::to_string(cpo_particle_property.get_deformation_type(cpo_data_position,data,mineral_i))));
              }

            for (size_t grain_i = 0; grain_i < n_grains; grain_i++)
              {
                auto rotated_matrix = rotate_6x6_matrix(*stiffness_matrix, transpose(cpo_particle_property.get_rotation_matrix_grains(cpo_data_position,data,mineral_i,grain_i)));
                Sav += cpo_particle_property.get_volume_fractions_grains(cpo_data_position,data,mineral_i,grain_i) * cpo_particle_property.get_volume_fraction_mineral(cpo_data_position,data,mineral_i) * rotated_matrix;
              }
          }

        return Sav;
      }



      template <int dim>
      void
      CpoElasticTensor<dim>::initialize_one_particle_property(const Point<dim> &,
                                                              std::vector<double> &data) const
      {

        // Get a reference to the CPO particle property.
        const Particle::Property::CrystalPreferredOrientation<dim> &cpo_particle_property =
          this->get_particle_world().get_property_manager().template get_matching_property<Particle::Property::CrystalPreferredOrientation<dim>>();

        SymmetricTensor<2,6> S_average = voigt_average_elastic_tensor(cpo_particle_property,
                                                                      cpo_data_position,
                                                                      data);

        for (unsigned int i = 0; i < SymmetricTensor<2,6>::n_independent_components ; ++i)
          {
            data.push_back(S_average[SymmetricTensor<2,6>::unrolled_to_component_indices(i)]);
          }


      }

      template <int dim>
      void
      CpoElasticTensor<dim>::update_one_particle_property(const unsigned int data_position,
                                                          const Point<dim> &,
                                                          const Vector<double> &,
                                                          const std::vector<Tensor<1,dim>> &,
                                                          const ArrayView<double> &data) const
      {
        // Get a reference to the CPO particle property.
        const Particle::Property::CrystalPreferredOrientation<dim> &cpo_particle_property =
          this->get_particle_world().get_property_manager().template get_matching_property<Particle::Property::CrystalPreferredOrientation<dim>>();


        SymmetricTensor<2,6> S_average = voigt_average_elastic_tensor(cpo_particle_property,
                                                                      cpo_data_position,
                                                                      data);

        Particle::Property::CpoElasticTensor<dim>::set_elastic_tensor(data_position,
                                                                      data,
                                                                      S_average);


      }


      template <int dim>
      SymmetricTensor<2,6>
      CpoElasticTensor<dim>::get_elastic_tensor(unsigned int cpo_data_position,
                                                const ArrayView<double> &data)
      {
        SymmetricTensor<2,6> elastic_tensor;
        for (unsigned int i = 0; i < SymmetricTensor<2,6>::n_independent_components ; ++i)
          elastic_tensor[SymmetricTensor<2,6>::unrolled_to_component_indices(i)] = data[cpo_data_position + i];
        return elastic_tensor;
      }


      template <int dim>
      void
      CpoElasticTensor<dim>::set_elastic_tensor(unsigned int cpo_data_position,
                                                const ArrayView<double> &data,
                                                SymmetricTensor<2,6> &elastic_tensor)
      {
        for (unsigned int i = 0; i < SymmetricTensor<2,6>::n_independent_components ; ++i)
          data[cpo_data_position + i] = elastic_tensor[SymmetricTensor<2,6>::unrolled_to_component_indices(i)];
      }


      template<int dim>
      Tensor<4,3>
      CpoElasticTensor<dim>::rotate_4th_order_tensor(const Tensor<4,3> &input_tensor, const Tensor<2,3> &rotation_tensor)
      {
        Tensor<4,3> output;

        for (unsigned short int i1 = 0; i1 < 3; i1++)
          {
            for (unsigned short int i2 = 0; i2 < 3; i2++)
              {
                for (unsigned short int i3 = 0; i3 < 3; i3++)
                  {
                    for (unsigned short int i4 = 0; i4 < 3; i4++)
                      {
                        for (unsigned short int j1 = 0; j1 < 3; j1++)
                          {
                            for (unsigned short int j2 = 0; j2 < 3; j2++)
                              {
                                for (unsigned short int j3 = 0; j3 < 3; j3++)
                                  {
                                    for (unsigned short int j4 = 0; j4 < 3; j4++)
                                      {
                                        output[i1][i2][i3][i4] = output[i1][i2][i3][i4] + rotation_tensor[i1][j1]*rotation_tensor[i2][j2]*rotation_tensor[i3][j3]*rotation_tensor[i4][j4]*input_tensor[j1][j2][j3][j4];
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
          }

        return output;
      }

      template<int dim>
      SymmetricTensor<2,6>
      CpoElasticTensor<dim>::rotate_6x6_matrix(const Tensor<2,6> &input_tensor, const Tensor<2,3> &rotation_tensor)
      {
        // we can represent the roation of the 4th order tensor as a rotation in the voigt
        // notation by computing $C'=MCM^{-1}$. Because M is orhtogonal we can replace $M^{-1}$
        // with $M^T$ resutling in $C'=MCM^{T}$ (Carcione, J. M. (2007). Wave Fields in Real Media:
        // Wave Propagation in Anisotropic, Anelastic, Porous and Electromagnetic Media. Netherlands:
        // Elsevier Science. Pages 8-9).
        Tensor<2,6> rotation_matrix;
        // top left block
        rotation_matrix[0][0] = rotation_tensor[0][0] * rotation_tensor[0][0];
        rotation_matrix[1][0] = rotation_tensor[1][0] * rotation_tensor[1][0];
        rotation_matrix[2][0] = rotation_tensor[2][0] * rotation_tensor[2][0];
        rotation_matrix[0][1] = rotation_tensor[0][1] * rotation_tensor[0][1];
        rotation_matrix[1][1] = rotation_tensor[1][1] * rotation_tensor[1][1];
        rotation_matrix[2][1] = rotation_tensor[2][1] * rotation_tensor[2][1];
        rotation_matrix[0][2] = rotation_tensor[0][2] * rotation_tensor[0][2];
        rotation_matrix[1][2] = rotation_tensor[1][2] * rotation_tensor[1][2];
        rotation_matrix[2][2] = rotation_tensor[2][2] * rotation_tensor[2][2];

        // top right block
        rotation_matrix[0][3] = 2.0 * rotation_tensor[0][1] * rotation_tensor[0][2];
        rotation_matrix[1][3] = 2.0 * rotation_tensor[1][1] * rotation_tensor[1][2];
        rotation_matrix[2][3] = 2.0 * rotation_tensor[2][1] * rotation_tensor[2][2];
        rotation_matrix[0][4] = 2.0 * rotation_tensor[0][2] * rotation_tensor[0][0];
        rotation_matrix[1][4] = 2.0 * rotation_tensor[1][2] * rotation_tensor[1][0];
        rotation_matrix[2][4] = 2.0 * rotation_tensor[2][2] * rotation_tensor[2][0];
        rotation_matrix[0][5] = 2.0 * rotation_tensor[0][0] * rotation_tensor[0][1];
        rotation_matrix[1][5] = 2.0 * rotation_tensor[1][0] * rotation_tensor[1][1];
        rotation_matrix[2][5] = 2.0 * rotation_tensor[2][0] * rotation_tensor[2][1];

        // bottom left block
        rotation_matrix[3][0] = rotation_tensor[1][0] * rotation_tensor[2][0];
        rotation_matrix[4][0] = rotation_tensor[2][0] * rotation_tensor[0][0];
        rotation_matrix[5][0] = rotation_tensor[0][0] * rotation_tensor[1][0];
        rotation_matrix[3][1] = rotation_tensor[1][1] * rotation_tensor[2][1];
        rotation_matrix[4][1] = rotation_tensor[2][1] * rotation_tensor[0][1];
        rotation_matrix[5][1] = rotation_tensor[0][1] * rotation_tensor[1][1];
        rotation_matrix[3][2] = rotation_tensor[1][2] * rotation_tensor[2][2];
        rotation_matrix[4][2] = rotation_tensor[2][2] * rotation_tensor[0][2];
        rotation_matrix[5][2] = rotation_tensor[0][2] * rotation_tensor[1][2];

        // bottom right block
        rotation_matrix[3][3] = rotation_tensor[1][1] * rotation_tensor[2][2] + rotation_tensor[1][2] * rotation_tensor[2][1];
        rotation_matrix[4][3] = rotation_tensor[0][1] * rotation_tensor[2][2] + rotation_tensor[0][2] * rotation_tensor[2][1];
        rotation_matrix[5][3] = rotation_tensor[0][1] * rotation_tensor[1][2] + rotation_tensor[0][2] * rotation_tensor[1][1];
        rotation_matrix[3][4] = rotation_tensor[1][0] * rotation_tensor[2][2] + rotation_tensor[1][2] * rotation_tensor[2][0];
        rotation_matrix[4][4] = rotation_tensor[0][2] * rotation_tensor[2][0] + rotation_tensor[0][0] * rotation_tensor[2][2];
        rotation_matrix[5][4] = rotation_tensor[0][2] * rotation_tensor[1][0] + rotation_tensor[0][0] * rotation_tensor[1][2];
        rotation_matrix[3][5] = rotation_tensor[1][1] * rotation_tensor[2][0] + rotation_tensor[1][0] * rotation_tensor[2][1];
        rotation_matrix[4][5] = rotation_tensor[0][0] * rotation_tensor[2][1] + rotation_tensor[0][1] * rotation_tensor[2][0];
        rotation_matrix[5][5] = rotation_tensor[0][0] * rotation_tensor[1][1] + rotation_tensor[0][1] * rotation_tensor[1][0];

        Tensor<2,6> rotation_matrix_tranposed = transpose(rotation_matrix);

        return symmetrize((rotation_matrix*input_tensor)*rotation_matrix_tranposed);
      }



      template<int dim>
      SymmetricTensor<2,6>
      CpoElasticTensor<dim>::transform_4th_order_tensor_to_6x6_matrix(const Tensor<4,3> &input_tensor)
      {
        SymmetricTensor<2,6> output;

        for (unsigned short int i = 0; i < 3; i++)
          {
            output[i][i] = input_tensor[i][i][i][i];
          }

        for (unsigned short int i = 1; i < 3; i++)
          {
            output[0][i] = 0.5*(input_tensor[0][0][i][i] + input_tensor[i][i][0][0]);
            //symmetry: output[0][i] = output[i][0];
          }
        output[1][2]=0.5*(input_tensor[1][1][2][2]+input_tensor[2][2][1][1]);
        //symmetry: output[2][1]=output[1][2];

        for (unsigned short int i = 0; i < 3; i++)
          {
            output[i][3]=0.25*(input_tensor[i][i][1][2]+input_tensor[i][i][2][1]+ input_tensor[1][2][i][i]+input_tensor[2][1][i][i]);
            //symmetry: output[3][i]=output[i][3];
          }

        for (unsigned short int i = 0; i < 3; i++)
          {
            output[i][4]=0.25*(input_tensor[i][i][0][2]+input_tensor[i][i][2][0]+ input_tensor[0][2][i][i]+input_tensor[2][0][i][i]);
            //symmetry: output[4][i]=output[i][4];
          }

        for (unsigned short int i = 0; i < 3; i++)
          {
            output[i][5]=0.25*(input_tensor[i][i][0][1]+input_tensor[i][i][1][0]+input_tensor[0][1][i][i]+input_tensor[1][0][i][i]);
            //symmetry: output[5][i]=output[i][5];
          }

        output[3][3]=0.25*(input_tensor[1][2][1][2]+input_tensor[1][2][2][1]+input_tensor[2][1][1][2]+input_tensor[2][1][2][1]);
        output[4][4]=0.25*(input_tensor[0][2][0][2]+input_tensor[0][2][2][0]+input_tensor[2][0][0][2]+input_tensor[2][0][2][0]);
        output[5][5]=0.25*(input_tensor[1][0][1][0]+input_tensor[1][0][0][1]+input_tensor[0][1][1][0]+input_tensor[0][1][0][1]);

        output[3][4]=0.125*(input_tensor[1][2][0][2]+input_tensor[1][2][2][0]+input_tensor[2][1][0][2]+input_tensor[2][1][2][0]+input_tensor[0][2][1][2]+input_tensor[0][2][2][1]+input_tensor[2][0][1][2]+input_tensor[2][0][2][1]);
        //symmetry: output[4][3]=output[3][4];
        output[3][5]=0.125*(input_tensor[1][2][0][1]+input_tensor[1][2][1][0]+input_tensor[2][1][0][1]+input_tensor[2][1][1][0]+input_tensor[0][1][1][2]+input_tensor[0][1][2][1]+input_tensor[1][0][1][2]+input_tensor[1][0][2][1]);
        //symmetry: output[5][3]=output[3][5];
        output[4][5]=0.125*(input_tensor[0][2][0][1]+input_tensor[0][2][1][0]+input_tensor[2][0][0][1]+input_tensor[2][0][1][0]+input_tensor[0][1][0][2]+input_tensor[0][1][2][0]+input_tensor[1][0][0][2]+input_tensor[1][0][2][0]);
        //symmetry: output[5][4]=output[4][5];

        return output;
      }

      template<int dim>
      Tensor<4,3>
      CpoElasticTensor<dim>::transform_6x6_matrix_to_4th_order_tensor(const SymmetricTensor<2,6> &input_tensor)
      {
        Tensor<4,3> output;

        for (unsigned short int i = 0; i < 3; i++)
          for (unsigned short int j = 0; j < 3; j++)
            for (unsigned short int k = 0; k < 3; k++)
              for (unsigned short int l = 0; l < 3; l++)
                {
                  // The first part of the inline if statment gets the diagonal.
                  // The second part is never higher then 5 (which is the limit of the tensor index)
                  // because to reach this part the variables need to be different, which results in
                  // at least a minus 1.
                  const unsigned short int p = (i == j ? i : 6 - i - j);
                  const unsigned short int q = (k == l ? k : 6 - k - l);
                  output[i][j][k][l] = input_tensor[p][q];
                }
        return output;
      }

      template<int dim>
      Tensor<1,21>
      CpoElasticTensor<dim>::transform_6x6_matrix_to_21D_vector(const SymmetricTensor<2,6> &input)
      {
        return Tensor<1,21,double> (
        {
          input[0][0],           // 0  // 1
          input[1][1],           // 1  // 2
          input[2][2],           // 2  // 3
          sqrt(2)*input[1][2],   // 3  // 4
          sqrt(2)*input[0][2],   // 4  // 5
          sqrt(2)*input[0][1],   // 5  // 6
          2*input[3][3],         // 6  // 7
          2*input[4][4],         // 7  // 8
          2*input[5][5],         // 8  // 9
          2*input[0][3],         // 9  // 10
          2*input[1][4],         // 10 // 11
          2*input[2][5],         // 11 // 12
          2*input[2][3],         // 12 // 13
          2*input[0][4],         // 13 // 14
          2*input[1][5],         // 14 // 15
          2*input[1][3],         // 15 // 16
          2*input[2][4],         // 16 // 17
          2*input[0][5],         // 17 // 18
          2*sqrt(2)*input[4][5], // 18 // 19
          2*sqrt(2)*input[3][5], // 19 // 20
          2*sqrt(2)*input[3][4]  // 20 // 21
        });

      }


      template<int dim>
      SymmetricTensor<2,6>
      CpoElasticTensor<dim>::transform_21D_vector_to_6x6_matrix(const Tensor<1,21> &input)
      {
        SymmetricTensor<2,6> result;

        constexpr double sqrt_2_inv = 1/sqrt(2);

        result[0][0] = input[0];
        result[1][1] = input[1];
        result[2][2] = input[2];
        result[1][2] = sqrt_2_inv * input[3];
        result[0][2] = sqrt_2_inv * input[4];
        result[0][1] = sqrt_2_inv * input[5];
        result[3][3] = 0.5 * input[6];
        result[4][4] = 0.5 * input[7];
        result[5][5] = 0.5 * input[8];
        result[0][3] = 0.5 * input[9];
        result[1][4] = 0.5 * input[10];
        result[2][5] = 0.5 * input[11];
        result[2][3] = 0.5 * input[12];
        result[0][4] = 0.5 * input[13];
        result[1][5] = 0.5 * input[14];
        result[1][3] = 0.5 * input[15];
        result[2][4] = 0.5 * input[16];
        result[0][5] = 0.5 * input[17];
        result[4][5] = 0.5 * sqrt_2_inv * input[18];
        result[3][5] = 0.5 * sqrt_2_inv * input[19];
        result[3][4] = 0.5 * sqrt_2_inv * input[20];

        return result;

      }


      template<int dim>
      Tensor<1,21>
      CpoElasticTensor<dim>::transform_4th_order_tensor_to_21D_vector(const Tensor<4,3> &input_tensor)
      {
        return Tensor<1,21,double> (
        {
          input_tensor[0][0][0][0],           // 0  // 1
          input_tensor[1][1][1][1],           // 1  // 2
          input_tensor[2][2][2][2],           // 2  // 3
          sqrt(2)*0.5*(input_tensor[1][1][2][2] + input_tensor[2][2][1][1]),   // 3  // 4
          sqrt(2)*0.5*(input_tensor[0][0][2][2] + input_tensor[2][2][0][0]),   // 4  // 5
          sqrt(2)*0.5*(input_tensor[0][0][1][1] + input_tensor[1][1][0][0]),   // 5  // 6
          0.5*(input_tensor[1][2][1][2]+input_tensor[1][2][2][1]+input_tensor[2][1][1][2]+input_tensor[2][1][2][1]),         // 6  // 7
          0.5*(input_tensor[0][2][0][2]+input_tensor[0][2][2][0]+input_tensor[2][0][0][2]+input_tensor[2][0][2][0]),         // 7  // 8
          0.5*(input_tensor[1][0][1][0]+input_tensor[1][0][0][1]+input_tensor[0][1][1][0]+input_tensor[0][1][0][1]),         // 8  // 9
          0.5*(input_tensor[0][0][1][2]+input_tensor[0][0][2][1]+input_tensor[1][2][0][0]+input_tensor[2][1][0][0]),         // 9  // 10
          0.5*(input_tensor[1][1][0][2]+input_tensor[1][1][2][0]+input_tensor[0][2][1][1]+input_tensor[2][0][1][1]),         // 10 // 11
          0.5*(input_tensor[2][2][0][1]+input_tensor[2][2][1][0]+input_tensor[0][1][2][2]+input_tensor[1][0][2][2]),         // 11 // 12
          0.5*(input_tensor[2][2][1][2]+input_tensor[2][2][2][1]+input_tensor[1][2][2][2]+input_tensor[2][1][2][2]),         // 12 // 13
          0.5*(input_tensor[0][0][0][2]+input_tensor[0][0][2][0]+input_tensor[0][2][0][0]+input_tensor[2][0][0][0]),         // 13 // 14
          0.5*(input_tensor[1][1][0][1]+input_tensor[1][1][1][0]+input_tensor[0][1][1][1]+input_tensor[1][0][1][1]),         // 14 // 15
          0.5*(input_tensor[1][1][1][2]+input_tensor[1][1][2][1]+input_tensor[1][2][1][1]+input_tensor[2][1][1][1]),         // 15 // 16
          0.5*(input_tensor[2][2][0][2]+input_tensor[2][2][2][0]+input_tensor[0][2][2][2]+input_tensor[2][0][2][2]),         // 16 // 17
          0.5*(input_tensor[0][0][0][1]+input_tensor[0][0][1][0]+input_tensor[0][1][0][0]+input_tensor[1][0][0][0]),         // 17 // 18
          sqrt(2)*0.25*(input_tensor[0][2][0][1]+input_tensor[0][2][1][0]+input_tensor[2][0][0][1]+input_tensor[2][0][1][0]+input_tensor[0][1][0][2]+input_tensor[0][1][2][0]+input_tensor[1][0][0][2]+input_tensor[1][0][2][0]), // 18 // 19
          sqrt(2)*0.25*(input_tensor[1][2][0][1]+input_tensor[1][2][1][0]+input_tensor[2][1][0][1]+input_tensor[2][1][1][0]+input_tensor[0][1][1][2]+input_tensor[0][1][2][1]+input_tensor[1][0][1][2]+input_tensor[1][0][2][1]), // 19 // 20
          sqrt(2)*0.25*(input_tensor[1][2][0][2]+input_tensor[1][2][2][0]+input_tensor[2][1][0][2]+input_tensor[2][1][2][0]+input_tensor[0][2][1][2]+input_tensor[0][2][2][1]+input_tensor[2][0][1][2]+input_tensor[2][0][2][1])  // 20 // 21
        });

      }

      template <int dim>
      UpdateTimeFlags
      CpoElasticTensor<dim>::need_update() const
      {
        return update_output_step;
      }

      template <int dim>
      UpdateFlags
      CpoElasticTensor<dim>::get_needed_update_flags () const
      {
        return update_default;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      CpoElasticTensor<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int>> property_information;

        property_information.push_back(std::make_pair("cpo_elastic_tensor",SymmetricTensor<2,6>::n_independent_components));

        return property_information;
      }

      template <int dim>
      void
      CpoElasticTensor<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Crystal Preferred Orientation");
            {
              prm.declare_entry ("Number of grains per particle", "50",
                                 Patterns::Integer (1),
                                 "The number of grains of each different mineral "
                                 "each particle contains.");

              prm.enter_subsection("Initial grains");
              {
                prm.enter_subsection("Uniform grains and random uniform rotations");
                {
                  prm.declare_entry ("Minerals", "Olivine: Karato 2008, Enstatite",
                                     Patterns::List(Patterns::Anything()),
                                     "This determines what minerals and fabrics or fabric selectors are used used for the LPO calculation. "
                                     "The options are Olivine: Passive, A-fabric, Olivine: B-fabric, Olivine: C-fabric, Olivine: D-fabric, "
                                     "Olivine: E-fabric, Olivine: Karato 2008 or Enstatite. Passive sets all RRSS entries to the maximum. The "
                                     "Karato 2008 selector selects a fabric based on stress and water content as defined in "
                                     "figure 4 of the Karato 2008 review paper (doi: 10.1146/annurev.earth.36.031207.124120).");
                }
                prm.leave_subsection ();
              }
              prm.leave_subsection ();
            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }


      template <int dim>
      void
      CpoElasticTensor<dim>::parse_parameters (ParameterHandler &prm)
      {

        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Crystal Preferred Orientation");
            {
              n_grains = prm.get_integer("Number of grains per particle");
              prm.enter_subsection("Initial grains");
              {
                prm.enter_subsection("Uniform grains and random uniform rotations");
                {
                  n_minerals = dealii::Utilities::split_string_list(prm.get("Minerals")).size();
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();


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
                                        "defined as the deformation the Voigt average of the "
                                        "elastic tensors of the cpo grains of the minerals."
                                        "Currently only Olivine and Enstatite are supported.")
    }
  }
}
