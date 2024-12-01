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

#include <aspect/particle/property/elastic_tensor_decomposition.h>
#include <aspect/particle/property/cpo_elastic_tensor.h>
#include <aspect/particle/property/crystal_preferred_orientation.h>
#include <aspect/particle/manager.h>

#include <aspect/utilities.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      namespace Utilities
      {
        std::array<unsigned int, 3>
        indexed_even_permutation(const unsigned int index)
        {
          // there are 6 permutations, but only the odd or even are needed. We use the even
          // permutation here.
          switch (index)
            {
              case 0 :
                return {{0,1,2}};
              case 1 :
                return {{1,2,0}};
              case 2:
                return {{2,0,1}};
              /*case 3:
                return {0,2,1};
              case 4 :
                return {1,0,2};
              case 5:
                return {2,1,0};*/
              default:
                AssertThrow(false,ExcMessage("Provided index larger then 2 (" + std::to_string(index)+ ")."));
                return {{0,0,0}};
            }

        }



        SymmetricTensor<2,3>
        compute_voigt_stiffness_tensor(const SymmetricTensor<2,6> &elastic_matrix)
        {
          /**
           * the Voigt stiffness tensor (see Browaeys and Chevrot, 2004)
           * It defines the stress needed to cause an isotropic strain in the
           * material
           */
          SymmetricTensor<2,3> voigt_stiffness_tensor;
          voigt_stiffness_tensor[0][0]=elastic_matrix[0][0]+elastic_matrix[5][5]+elastic_matrix[4][4];
          voigt_stiffness_tensor[1][1]=elastic_matrix[5][5]+elastic_matrix[1][1]+elastic_matrix[3][3];
          voigt_stiffness_tensor[2][2]=elastic_matrix[4][4]+elastic_matrix[3][3]+elastic_matrix[2][2];
          voigt_stiffness_tensor[1][0]=elastic_matrix[0][5]+elastic_matrix[1][5]+elastic_matrix[3][4];
          voigt_stiffness_tensor[2][0]=elastic_matrix[0][4]+elastic_matrix[2][4]+elastic_matrix[3][5];
          voigt_stiffness_tensor[2][1]=elastic_matrix[1][3]+elastic_matrix[2][3]+elastic_matrix[4][5];

          return voigt_stiffness_tensor;
        }



        SymmetricTensor<2,3>
        compute_dilatation_stiffness_tensor(const SymmetricTensor<2,6> &elastic_matrix)
        {
          /**
           * The dilatational stiffness tensor (see Browaeys and Chevrot, 2004)
           * It defines the stress to cause isotropic dilatation in the material.
           */
          SymmetricTensor<2,3> dilatation_stiffness_tensor;
          for (size_t i = 0; i < 3; i++)
            {
              dilatation_stiffness_tensor[0][0]=elastic_matrix[0][i]+dilatation_stiffness_tensor[0][0];
              dilatation_stiffness_tensor[1][1]=elastic_matrix[1][i]+dilatation_stiffness_tensor[1][1];
              dilatation_stiffness_tensor[2][2]=elastic_matrix[2][i]+dilatation_stiffness_tensor[2][2];
              dilatation_stiffness_tensor[1][0]=elastic_matrix[5][i]+dilatation_stiffness_tensor[1][0];
              dilatation_stiffness_tensor[2][0]=elastic_matrix[4][i]+dilatation_stiffness_tensor[2][0];
              dilatation_stiffness_tensor[2][1]=elastic_matrix[3][i]+dilatation_stiffness_tensor[2][1];
            }
          return dilatation_stiffness_tensor;
        }



        Tensor<2,3>
        compute_unpermutated_SCCS(const SymmetricTensor<2,3> &dilatation_stiffness_tensor,
                                  const SymmetricTensor<2,3> &voigt_stiffness_tensor)
        {
          // computing the eigenvector of the dilation and Voigt stiffness matrices and then averaging them by bysection.
          const std::array<std::pair<double,Tensor<1,3,double>>, 3> voigt_eigenvectors_a = eigenvectors(voigt_stiffness_tensor, SymmetricTensorEigenvectorMethod::jacobi);
          const std::array<std::pair<double,Tensor<1,3,double>>, 3> dilatation_eigenvectors_a = eigenvectors(dilatation_stiffness_tensor, SymmetricTensorEigenvectorMethod::jacobi);


          std::array<Tensor<1,3,double>,3> unpermutated_SCCS;
          // Averaging dilatation and voigt eigenvectors
          // To do this we need to find the eigenvectors which are closest to each other and average those.
          // The next function looks for the smallest angle and returns the corresponding vector index for
          // that vector.
          int vector_index_signed = 0;
          for (unsigned int i1 = 0; i1 < 3; i1++)
            {
              vector_index_signed = 0;
              double smallest_angle = std::numeric_limits<double>::max(); // angle D VeCtor
              for (unsigned int i2 = 0; i2 < 3; i2++)
                {
                  double dv_dot_product = dilatation_eigenvectors_a[i1].second*voigt_eigenvectors_a[i2].second;
                  // limit the dot product between 1 and -1 so we can use the arccos function safely.
                  if (std::abs(dv_dot_product) >= 1.0)
                    dv_dot_product = std::copysign(1.0,dv_dot_product);
                  // Compute the angle between the vectors and account for that vector in the opposite
                  // direction are the same (0 == 180 degrees). So limit the direction of the vectors between
                  // 0 and 90 degrees such that it represents the minimum angle between the two lines.
                  const double angle = dv_dot_product < 0.0 ? std::acos(-1.0)-std::acos(dv_dot_product) : std::acos(dv_dot_product);
                  // store this if the angle is smaller
                  if (angle < smallest_angle)
                    {
                      vector_index_signed = (dv_dot_product < 0.0) ? -i2 : i2;
                      smallest_angle = angle;
                    }
                }

              // Adds to the dilatation eigenvectors to the Voigt eigenvectors with the smallest angle
              // Note that the voigt eigenvector is multiplied with a value which can be negative, which means it would be a subtraction.
              // Lastly we normalize the unpermutated_SCCS.
              unpermutated_SCCS[i1] = 0.5*(dilatation_eigenvectors_a[i1].second + static_cast<double>(vector_index_signed)*voigt_eigenvectors_a[std::abs(vector_index_signed)].second);
              unpermutated_SCCS[i1] /= unpermutated_SCCS[i1].norm();
            }

          return Tensor<2,3>(
          {
            {unpermutated_SCCS[0][0],unpermutated_SCCS[0][1],unpermutated_SCCS[0][2]},
            {unpermutated_SCCS[1][0],unpermutated_SCCS[1][1],unpermutated_SCCS[1][2]},
            {unpermutated_SCCS[2][0],unpermutated_SCCS[2][1],unpermutated_SCCS[2][2]}
          });
        }



        std::array<std::array<double,3>,7>
        compute_elastic_tensor_SCCS_decompositions(
          const Tensor<2,3> &unpermutated_SCCS,
          const SymmetricTensor<2,6> &elastic_matrix)
        {
          /**
           * Try the different permutations to determine what is the best hexagonal projection.
           * This is based on Browaeys and Chevrot (2004), GJI (doi: 10.1111/j.1365-246X.2004.024115.x),
           * which states at the end of paragraph 3.3 that "... an important property of an orthogonal projection
           * is that the distance between a vector $X$ and its orthogonal projection $X_H = p(X)$ on a given
           * subspace is minimum. These two features ensure that the decomposition is optimal once a 3-D Cartesian
           * coordinate system is chosen.". The other property they talk about is that "The space of elastic
           * vectors has a finite dimension [...], i.e. using a different norm from eq. (2.3 will change distances
           * but not the resulting decomposition.".
           */
          std::array<Tensor<2,3>,3> permutated;
          std::array<SymmetricTensor<2,6>,3> rotated_elastic_matrix;
          // The norms variable contains the square norms of the different symmetry approximations of the elastic tensor.
          std::array<std::array<double,3>,7> squared_norms_to_projections;

          for (unsigned int permutation_i = 0; permutation_i < 3; permutation_i++)
            {
              std::array<unsigned int, 3> permutation = indexed_even_permutation(permutation_i);

              for (size_t j = 0; j < 3; j++)
                {
                  permutated[permutation_i][j] = unpermutated_SCCS[permutation[j]];
                }

              rotated_elastic_matrix[permutation_i] = aspect::Utilities::Tensors::rotate_voigt_stiffness_matrix((permutated[permutation_i]),elastic_matrix);

              const Tensor<1,21> full_elastic_vector_rotated = aspect::Utilities::Tensors::to_voigt_stiffness_vector(rotated_elastic_matrix[permutation_i]);


              const double full_norm_square = full_elastic_vector_rotated.norm_square();
              squared_norms_to_projections[6][permutation_i] = full_norm_square;

              // Get the monoclinic and higher symmetry axes, which can be comptued by taking specific elements from the full vector.
              // This means that this vector contains all symmetry axes, but the isotropic part is removed.
              // The following line would do the same as the lines below, but is is very slow. It has therefore been
              // replaced by the lines below.
              //auto monoclinic_and_higher_vector = projection_matrix_triclinic_to_monoclinic*full_elastic_vector_rotated;
              dealii::Tensor<1,21> monoclinic_and_higher_vector;
              monoclinic_and_higher_vector[0] = full_elastic_vector_rotated[0];
              monoclinic_and_higher_vector[1] = full_elastic_vector_rotated[1];
              monoclinic_and_higher_vector[2] = full_elastic_vector_rotated[2];
              monoclinic_and_higher_vector[3] = full_elastic_vector_rotated[3];
              monoclinic_and_higher_vector[4] = full_elastic_vector_rotated[4];
              monoclinic_and_higher_vector[5] = full_elastic_vector_rotated[5];
              monoclinic_and_higher_vector[6] = full_elastic_vector_rotated[6];
              monoclinic_and_higher_vector[7] = full_elastic_vector_rotated[7];
              monoclinic_and_higher_vector[8] = full_elastic_vector_rotated[8];
              monoclinic_and_higher_vector[11] = full_elastic_vector_rotated[11];
              monoclinic_and_higher_vector[14] = full_elastic_vector_rotated[14];
              monoclinic_and_higher_vector[17] = full_elastic_vector_rotated[17];
              monoclinic_and_higher_vector[20] = full_elastic_vector_rotated[20];

              // The triclinic vector is the full elastic tensor minux the monoclinic and higher symmetry axes vector.
              auto triclinic_vector = full_elastic_vector_rotated-monoclinic_and_higher_vector;
              squared_norms_to_projections[0][permutation_i] = triclinic_vector.norm_square();

              // Only the first 9 elements are now non-zero, so crop the vector to the first 9 elements.
              // The following line would do the same as the lines below, but it is slow. It has therefore been
              // replaced by the lines below.
              //auto orthorhombic_and_higher_vector = projection_matrix_monoclinic_to_orthorhombic*monoclinic_and_higher_vector;
              dealii::Tensor<1,9>  monoclinic_and_higher_vector_cropped;
              monoclinic_and_higher_vector_cropped[0] = monoclinic_and_higher_vector[0];
              monoclinic_and_higher_vector_cropped[1] = monoclinic_and_higher_vector[1];
              monoclinic_and_higher_vector_cropped[2] = monoclinic_and_higher_vector[2];
              monoclinic_and_higher_vector_cropped[3] = monoclinic_and_higher_vector[3];
              monoclinic_and_higher_vector_cropped[4] = monoclinic_and_higher_vector[4];
              monoclinic_and_higher_vector_cropped[5] = monoclinic_and_higher_vector[5];
              monoclinic_and_higher_vector_cropped[6] = monoclinic_and_higher_vector[6];
              monoclinic_and_higher_vector_cropped[7] = monoclinic_and_higher_vector[7];
              monoclinic_and_higher_vector_cropped[8] = monoclinic_and_higher_vector[8];
              dealii::Tensor<1,9> orthorhombic_and_higher_vector;
              orthorhombic_and_higher_vector[0] = monoclinic_and_higher_vector[0];
              orthorhombic_and_higher_vector[1] = monoclinic_and_higher_vector[1];
              orthorhombic_and_higher_vector[2] = monoclinic_and_higher_vector[2];
              orthorhombic_and_higher_vector[3] = monoclinic_and_higher_vector[3];
              orthorhombic_and_higher_vector[4] = monoclinic_and_higher_vector[4];
              orthorhombic_and_higher_vector[5] = monoclinic_and_higher_vector[5];
              orthorhombic_and_higher_vector[6] = monoclinic_and_higher_vector[6];
              orthorhombic_and_higher_vector[7] = monoclinic_and_higher_vector[7];
              orthorhombic_and_higher_vector[8] = monoclinic_and_higher_vector[8];

              // The monoclinic vector is the monoclinic and higher symmetry axes vector minux the orthoclinic and higher symmetry axes vector.
              auto monoclinic_vector = monoclinic_and_higher_vector_cropped-orthorhombic_and_higher_vector;
              squared_norms_to_projections[1][permutation_i] = monoclinic_vector.norm_square();


              auto tetragonal_and_higher_vector = projection_matrix_orthorhombic_to_tetragonal*orthorhombic_and_higher_vector;
              auto orthorhombic_vector = orthorhombic_and_higher_vector-tetragonal_and_higher_vector;
              squared_norms_to_projections[2][permutation_i] = orthorhombic_vector.norm_square();

              auto hexagonal_and_higher_vector = projection_matrix_tetragonal_to_hexagonal*tetragonal_and_higher_vector;
              auto tetragonal_vector = tetragonal_and_higher_vector-hexagonal_and_higher_vector;
              squared_norms_to_projections[3][permutation_i] = tetragonal_vector.norm_square();

              auto isotropic_vector = projection_matrix_hexagonal_to_isotropic*hexagonal_and_higher_vector;
              auto hexagonal_vector = hexagonal_and_higher_vector-isotropic_vector;
              squared_norms_to_projections[4][permutation_i] = hexagonal_vector.norm_square();
              squared_norms_to_projections[5][permutation_i] = isotropic_vector.norm_square();

            }
          return squared_norms_to_projections;

        }
      }



      template <int dim>
      void
      ElasticTensorDecomposition<dim>::initialize ()
      {
        const Particle::Property::Manager<dim> &manager = this->get_particle_manager(this->get_particle_manager_index()).get_property_manager();
        AssertThrow(manager.plugin_name_exists("crystal preferred orientation"),
                    ExcMessage("No cpo property plugin found."));
        AssertThrow(manager.plugin_name_exists("cpo elastic tensor"),
                    ExcMessage("No cpo elastic tensor property plugin found."));

        AssertThrow(manager.check_plugin_order("crystal preferred orientation","elastic tensor decomposition"),
                    ExcMessage("To use the elastic tensor decomposition plugin, the cpo plugin needs to be defined before this plugin."));

        AssertThrow(manager.check_plugin_order("cpo elastic tensor","elastic tensor decomposition"),
                    ExcMessage("To use the elastic tensor decomposition plugin, the cpo elastic tensor plugin needs to be defined before this plugin."));

        cpo_elastic_tensor_data_position = manager.get_data_info().get_position_by_plugin_index(manager.get_plugin_index_by_name("cpo elastic tensor"));
      }



      template <int dim>
      void
      ElasticTensorDecomposition<dim>::initialize_one_particle_property(const Point<dim> &,
                                                                        std::vector<double> &data) const
      {
        const SymmetricTensor<2,6> elastic_matrix = Particle::Property::CpoElasticTensor<dim>::get_elastic_tensor(cpo_elastic_tensor_data_position,
                                                    data);

        const SymmetricTensor<2,3> dilatation_stiffness_tensor_full = Property::Utilities::compute_dilatation_stiffness_tensor(elastic_matrix);
        const SymmetricTensor<2,3> voigt_stiffness_tensor_full = Property::Utilities::compute_voigt_stiffness_tensor(elastic_matrix);
        const Tensor<2,3> SCCS_full = Property::Utilities::compute_unpermutated_SCCS(dilatation_stiffness_tensor_full, voigt_stiffness_tensor_full);

        const std::array<std::array<double,3>,7 > norms = Property::Utilities::compute_elastic_tensor_SCCS_decompositions(SCCS_full, elastic_matrix);

        // get max hexagonal element index, which is the same as the permutation index
        const size_t max_hexagonal_element_index = std::max_element(norms[4].begin(),norms[4].end())-norms[4].begin();
        std::array<unsigned int, 3> permutation = Property::Utilities::indexed_even_permutation(max_hexagonal_element_index);
        // reorder the SCCS be the SCCS permutation which yields the largest hexagonal vector (percentage of anisotropy)
        Tensor<2,3> hexagonal_permutated_SCCS;
        for (size_t index = 0; index < 3; ++index)
          {
            hexagonal_permutated_SCCS[index] = SCCS_full[permutation[index]];
          }

        data.push_back(SCCS_full[0][0]);
        data.push_back(SCCS_full[0][1]);
        data.push_back(SCCS_full[0][2]);
        data.push_back(SCCS_full[1][0]);
        data.push_back(SCCS_full[1][1]);
        data.push_back(SCCS_full[1][2]);
        data.push_back(SCCS_full[2][0]);
        data.push_back(SCCS_full[2][1]);
        data.push_back(SCCS_full[2][2]);
        data.push_back(hexagonal_permutated_SCCS[2][0]);
        data.push_back(hexagonal_permutated_SCCS[2][1]);
        data.push_back(hexagonal_permutated_SCCS[2][2]);
        data.push_back(norms[6][0]);
        data.push_back(norms[0][0]); // triclinic
        data.push_back(norms[0][1]); // triclinic
        data.push_back(norms[0][2]); // triclinic
        data.push_back(norms[1][0]); // monoclinic
        data.push_back(norms[1][1]); // monoclinic
        data.push_back(norms[1][2]); // monoclinic
        data.push_back(norms[2][0]); // orthorhomic
        data.push_back(norms[2][1]); // orthorhomic
        data.push_back(norms[2][2]); // orthorhomic
        data.push_back(norms[3][0]); // tetragonal
        data.push_back(norms[3][1]); // tetragonal
        data.push_back(norms[3][2]); // tetragonal
        data.push_back(norms[4][0]); // hexagonal
        data.push_back(norms[4][1]); // hexagonal
        data.push_back(norms[4][2]); // hexagonal
        data.push_back(norms[5][0]); // isotropic

      }



      template <int dim>
      void
      ElasticTensorDecomposition<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &/*inputs*/,
                                                                  typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        const unsigned int data_position = this->data_position;
        for (auto &particle: particles)
          {
            ArrayView<double> data = particle.get_properties();
            const SymmetricTensor<2,6> elastic_matrix = Particle::Property::CpoElasticTensor<dim>::get_elastic_tensor(cpo_elastic_tensor_data_position,
                                                        data);

            const SymmetricTensor<2,3> dilatation_stiffness_tensor_full = Utilities::compute_dilatation_stiffness_tensor(elastic_matrix);
            const SymmetricTensor<2,3> voigt_stiffness_tensor_full = Utilities::compute_voigt_stiffness_tensor(elastic_matrix);
            const Tensor<2,3> SCCS_full = Utilities::compute_unpermutated_SCCS(dilatation_stiffness_tensor_full, voigt_stiffness_tensor_full);

            const std::array<std::array<double,3>,7 > norms = Utilities::compute_elastic_tensor_SCCS_decompositions(SCCS_full, elastic_matrix);

            // get max hexagonal element index, which is the same as the permutation index
            const size_t max_hexagonal_element_index = std::max_element(norms[4].begin(),norms[4].end())-norms[4].begin();
            std::array<unsigned int, 3> permutation = Utilities::indexed_even_permutation(max_hexagonal_element_index);
            // reorder the SCCS by the SCCS permutation which yields the largest hexagonal vector (percentage of anisotropy)
            Tensor<2,3> hexagonal_permutated_SCCS;
            for (size_t index = 0; index < 3; ++index)
              {
                hexagonal_permutated_SCCS[index] = SCCS_full[permutation[index]];
              }

            data[data_position]    = SCCS_full[0][0];
            data[data_position+1]  = SCCS_full[0][1];
            data[data_position+2]  = SCCS_full[0][2];
            data[data_position+3]  = SCCS_full[1][0];
            data[data_position+4]  = SCCS_full[1][1];
            data[data_position+5]  = SCCS_full[1][2];
            data[data_position+6]  = SCCS_full[2][0];
            data[data_position+7]  = SCCS_full[2][1];
            data[data_position+8]  = SCCS_full[2][2];
            data[data_position+9]  = hexagonal_permutated_SCCS[2][0];
            data[data_position+10] = hexagonal_permutated_SCCS[2][1];
            data[data_position+11] = hexagonal_permutated_SCCS[2][2];
            data[data_position+12] = norms[6][0];
            data[data_position+13] = norms[0][0]; // triclinic
            data[data_position+14] = norms[0][1]; // triclinic
            data[data_position+15] = norms[0][2]; // triclinic
            data[data_position+16] = norms[1][0]; // monoclinic
            data[data_position+17] = norms[1][1]; // monoclinic
            data[data_position+18] = norms[1][2]; // monoclinic
            data[data_position+19] = norms[2][0]; // orthorhomic
            data[data_position+20] = norms[2][1]; // orthorhomic
            data[data_position+21] = norms[2][2]; // orthorhomic
            data[data_position+22] = norms[3][0]; // tetragonal
            data[data_position+23] = norms[3][1]; // tetragonal
            data[data_position+24] = norms[3][2]; // tetragonal
            data[data_position+25] = norms[4][0]; // hexagonal
            data[data_position+26] = norms[4][1]; // hexagonal
            data[data_position+27] = norms[4][2]; // hexagonal
            data[data_position+28] = norms[5][0]; // isotropic
          }
      }



      template <int dim>
      UpdateTimeFlags
      ElasticTensorDecomposition<dim>::need_update() const
      {
        return update_output_step;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      ElasticTensorDecomposition<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int>> property_information =
        {
          {"cpo elastic axis e1",3},
          {"cpo elastic axis e2",3},
          {"cpo elastic axis e3",3},
          {"cpo elastic hexagonal axis",3},
          {"cpo elastic vector norm square",1},
          {"cpo elastic triclinic vector norm square p1",1},
          {"cpo elastic triclinic vector norm square p2",1},
          {"cpo elastic triclinic vector norm square p3",1},
          {"cpo elastic monoclinic vector norm square p1",1},
          {"cpo elastic monoclinic vector norm square p2",1},
          {"cpo elastic monoclinic vector norm square p3",1},
          {"cpo elastic orthorhombic vector norm square p1",1},
          {"cpo elastic orthorhombic vector norm square p2",1},
          {"cpo elastic orthorhombic vector norm square p3",1},
          {"cpo elastic tetragonal vector norm square p1",1},
          {"cpo elastic tetragonal vector norm square p2",1},
          {"cpo elastic tetragonal vector norm square p3",1},
          {"cpo elastic hexagonal vector norm square p1",1},
          {"cpo elastic hexagonal vector norm square p2",1},
          {"cpo elastic hexagonal vector norm square p3",1},
          {"cpo elastic isotropic vector norm square",1}
        };

        return property_information;
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(ElasticTensorDecomposition,
                                        "elastic tensor decomposition",
                                        "A plugin which decomposes the elastic tensor into different approximations "
                                        "(Isotropic, Hexagonal, Tetragonal, Orthorhombic, Monoclinic and Triclinic) "
                                        "and provides the eigenvectors of the tensor.")
    }
  }
}
