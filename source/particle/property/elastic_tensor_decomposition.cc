/*
 Copyright (C) 2023 by the authors of the ASPECT code.
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
#include <aspect/particle/property/elastic_tensor_decomposition.h>
#include <aspect/particle/property/cpo_elastic_tensor.h>
#include <aspect/particle/property/crystal_preferred_orientation.h>
#include <aspect/particle/world.h>

#include <aspect/utilities.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {

      template <int dim> SymmetricTensor<2,21> ElasticTensorDecomposition<dim>::projection_matrix_tric_to_mono;
      template <int dim> SymmetricTensor<2,9> ElasticTensorDecomposition<dim>::projection_matrix_mono_to_ortho;
      template <int dim> SymmetricTensor<2,9> ElasticTensorDecomposition<dim>::projection_matrix_ortho_to_tetra;
      template <int dim> SymmetricTensor<2,9> ElasticTensorDecomposition<dim>::projection_matrix_tetra_to_hexa;
      template <int dim> SymmetricTensor<2,9> ElasticTensorDecomposition<dim>::projection_matrix_hexa_to_iso;



      template <int dim>
      ElasticTensorDecomposition<dim>::ElasticTensorDecomposition ()
      {
        // setup projection matrices
        // projection_matrix_tric_to_mono
        projection_matrix_tric_to_mono[0][0] = 1.0;
        projection_matrix_tric_to_mono[1][1] = 1.0;
        projection_matrix_tric_to_mono[2][2] = 1.0;
        projection_matrix_tric_to_mono[3][3] = 1.0;
        projection_matrix_tric_to_mono[4][4] = 1.0;
        projection_matrix_tric_to_mono[5][5] = 1.0;
        projection_matrix_tric_to_mono[6][6] = 1.0;
        projection_matrix_tric_to_mono[7][7] = 1.0;
        projection_matrix_tric_to_mono[8][8] = 1.0;
        projection_matrix_tric_to_mono[11][11] = 1.0;
        projection_matrix_tric_to_mono[14][14] = 1.0;
        projection_matrix_tric_to_mono[17][17] = 1.0;
        projection_matrix_tric_to_mono[20][20] = 1.0;


        // projection_matrix_mono_to_ortho;
        projection_matrix_mono_to_ortho[0][0] = 1.0;
        projection_matrix_mono_to_ortho[1][1] = 1.0;
        projection_matrix_mono_to_ortho[2][2] = 1.0;
        projection_matrix_mono_to_ortho[3][3] = 1.0;
        projection_matrix_mono_to_ortho[4][4] = 1.0;
        projection_matrix_mono_to_ortho[5][5] = 1.0;
        projection_matrix_mono_to_ortho[6][6] = 1.0;
        projection_matrix_mono_to_ortho[7][7] = 1.0;
        projection_matrix_mono_to_ortho[8][8] = 1.0;

        // projection_matrix_ortho_to_tetra;
        projection_matrix_ortho_to_tetra[0][0] = 0.5;
        projection_matrix_ortho_to_tetra[0][1] = 0.5;
        projection_matrix_ortho_to_tetra[1][1] = 0.5;
        projection_matrix_ortho_to_tetra[2][2] = 1.0;
        projection_matrix_ortho_to_tetra[3][3] = 0.5;
        projection_matrix_ortho_to_tetra[3][4] = 0.5;
        projection_matrix_ortho_to_tetra[4][4] = 0.5;
        projection_matrix_ortho_to_tetra[5][5] = 1.0;
        projection_matrix_ortho_to_tetra[6][6] = 0.5;
        projection_matrix_ortho_to_tetra[6][7] = 0.5;
        projection_matrix_ortho_to_tetra[7][7] = 0.5;
        projection_matrix_ortho_to_tetra[8][8] = 1.0;

        // projection_matrix_tetra_to_hexa;
        projection_matrix_tetra_to_hexa[0][0] = 3./8.;
        projection_matrix_tetra_to_hexa[0][1] = 3./8.;
        projection_matrix_tetra_to_hexa[1][1] = 3./8.;
        projection_matrix_tetra_to_hexa[2][2] = 1.0;
        projection_matrix_tetra_to_hexa[3][3] = 0.5;
        projection_matrix_tetra_to_hexa[3][4] = 0.5;
        projection_matrix_tetra_to_hexa[4][4] = 0.5;
        projection_matrix_tetra_to_hexa[5][5] = 3./4.;
        projection_matrix_tetra_to_hexa[6][6] = 0.5;
        projection_matrix_tetra_to_hexa[6][7] = 0.5;
        projection_matrix_tetra_to_hexa[7][7] = 0.5;
        projection_matrix_tetra_to_hexa[8][8] = 0.5;
        projection_matrix_tetra_to_hexa[5][0] = 1./(4.*std::sqrt(2.0));
        projection_matrix_tetra_to_hexa[5][1] = 1./(4.*std::sqrt(2.0));
        projection_matrix_tetra_to_hexa[8][0] = 0.25;
        projection_matrix_tetra_to_hexa[8][1] = 0.25;
        projection_matrix_tetra_to_hexa[8][5] = -1/(2*std::sqrt(2.0));


        // projection_matrix_hexa_to_iso;
        projection_matrix_hexa_to_iso[0][0] = 3./15.;
        projection_matrix_hexa_to_iso[0][1] = 3./15.;
        projection_matrix_hexa_to_iso[0][2] = 3./15.;
        projection_matrix_hexa_to_iso[1][1] = 3./15.;
        projection_matrix_hexa_to_iso[1][2] = 3./15.;
        projection_matrix_hexa_to_iso[2][2] = 3./15.;
        projection_matrix_hexa_to_iso[3][3] = 4./15.;
        projection_matrix_hexa_to_iso[3][4] = 4./15.;
        projection_matrix_hexa_to_iso[3][5] = 4./15.;
        projection_matrix_hexa_to_iso[4][4] = 4./15.;
        projection_matrix_hexa_to_iso[4][5] = 4./15.;
        projection_matrix_hexa_to_iso[5][5] = 4./15.;
        projection_matrix_hexa_to_iso[6][6] = 1./5.;
        projection_matrix_hexa_to_iso[6][7] = 1./5.;
        projection_matrix_hexa_to_iso[6][8] = 1./5.;
        projection_matrix_hexa_to_iso[7][7] = 1./5.;
        projection_matrix_hexa_to_iso[7][8] = 1./5.;
        projection_matrix_hexa_to_iso[8][8] = 1./5.;

        projection_matrix_hexa_to_iso[0][3] = std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[0][4] = std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[0][5] = std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[1][3] = std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[1][4] = std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[1][5] = std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[2][3] = std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[2][4] = std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[2][5] = std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[0][6] = 2./15.;
        projection_matrix_hexa_to_iso[0][7] = 2./15.;
        projection_matrix_hexa_to_iso[0][8] = 2./15.;
        projection_matrix_hexa_to_iso[1][6] = 2./15.;
        projection_matrix_hexa_to_iso[1][7] = 2./15.;
        projection_matrix_hexa_to_iso[1][8] = 2./15.;
        projection_matrix_hexa_to_iso[2][6] = 2./15.;
        projection_matrix_hexa_to_iso[2][7] = 2./15.;
        projection_matrix_hexa_to_iso[2][8] = 2./15.;
        projection_matrix_hexa_to_iso[3][6] = -std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[3][7] = -std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[3][8] = -std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[4][6] = -std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[4][7] = -std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[4][8] = -std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[5][6] = -std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[5][7] = -std::sqrt(2.0)/15.;
        projection_matrix_hexa_to_iso[5][8] = -std::sqrt(2.0)/15.;
      }



      template <int dim>
      void
      ElasticTensorDecomposition<dim>::initialize ()
      {
        const Particle::Property::Manager<dim> &manager = this->get_particle_world().get_property_manager();
        AssertThrow(manager.plugin_name_exists("crystal preferred orientation"),
                    ExcMessage("No cpo property plugin found."));
        AssertThrow(manager.plugin_name_exists("cpo elastic tensor"),
                    ExcMessage("No cpo elastic tensor property plugin found."));

        AssertThrow(manager.check_plugin_order("crystal preferred orientation","elastic tensor decomposition"),
                    ExcMessage("To use the elastic tensor decomposition plugin, the cpo plugin need to be defined before this plugin."));

        AssertThrow(manager.check_plugin_order("cpo elastic tensor","elastic tensor decomposition"),
                    ExcMessage("To use the elastic tensor decomposition plugin, the cpo elastic tensor plugin need to be defined before this plugin."));

        cpo_elastic_tensor_data_position = manager.get_data_info().get_position_by_plugin_index(manager.get_plugin_index_by_name("cpo elastic tensor"));
      }



      template <int dim>
      void
      ElasticTensorDecomposition<dim>::initialize_one_particle_property(const Point<dim> &,
                                                                        std::vector<double> &data) const
      {
        const SymmetricTensor<2,6> elastic_matrix = Particle::Property::CpoElasticTensor<dim>::get_elastic_tensor(cpo_elastic_tensor_data_position,
                                                    data);

        const SymmetricTensor<2,3> dilatation_stiffness_tensor_full = compute_dilatation_stiffness_tensor(elastic_matrix);
        const SymmetricTensor<2,3> voigt_stiffness_tensor_full = compute_voigt_stiffness_tensor(elastic_matrix);
        Tensor<2,3> SCCS_full = compute_unpermutated_SCCS(dilatation_stiffness_tensor_full, voigt_stiffness_tensor_full);

        std::array<std::array<double,3>,7 > norms = compute_elastic_tensor_SCCS_decompositions(SCCS_full, elastic_matrix);

        // get max hexagonal element index, which is the same as the permutation index
        const size_t max_hexagonal_element_index = std::max_element(norms[4].begin(),norms[4].end())-norms[4].begin();
        std::array<unsigned short int, 3> perumation = indexed_even_permutation(max_hexagonal_element_index);
        // reorder the SCCS be the SCCS permutation which yields the largest hexagonal vector (percentage of anisotropy)
        Tensor<2,3> hexa_permutated_SCCS;
        for (size_t index = 0; index < 3; index++)
          {
            hexa_permutated_SCCS[index] = SCCS_full[perumation[index]];
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
        data.push_back(hexa_permutated_SCCS[2][0]);
        data.push_back(hexa_permutated_SCCS[2][1]);
        data.push_back(hexa_permutated_SCCS[2][2]);
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
      ElasticTensorDecomposition<dim>::update_one_particle_property(const unsigned int data_position,
                                                                    const Point<dim> &,
                                                                    const Vector<double> &,
                                                                    const std::vector<Tensor<1,dim>> &,
                                                                    const ArrayView<double> &data) const
      {
        const SymmetricTensor<2,6> elastic_matrix = Particle::Property::CpoElasticTensor<dim>::get_elastic_tensor(cpo_elastic_tensor_data_position,
                                                    data);


        const SymmetricTensor<2,3> dilatation_stiffness_tensor_full = compute_dilatation_stiffness_tensor(elastic_matrix);
        const SymmetricTensor<2,3> voigt_stiffness_tensor_full = compute_voigt_stiffness_tensor(elastic_matrix);
        Tensor<2,3> SCCS_full = compute_unpermutated_SCCS(dilatation_stiffness_tensor_full, voigt_stiffness_tensor_full);

        std::array<std::array<double,3>,7 > norms = compute_elastic_tensor_SCCS_decompositions(SCCS_full, elastic_matrix);

        // get max hexagonal element index, which is the same as the permutation index
        const size_t max_hexagonal_element_index = std::max_element(norms[4].begin(),norms[4].end())-norms[4].begin();
        std::array<unsigned short int, 3> perumation = indexed_even_permutation(max_hexagonal_element_index);
        // reorder the SCCS be the SCCS permutation which yields the largest hexagonal vector (percentage of anisotropy)
        Tensor<2,3> hexa_permutated_SCCS;
        for (size_t index = 0; index < 3; index++)
          {
            hexa_permutated_SCCS[index] = SCCS_full[perumation[index]];
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
        data[data_position+9]  = hexa_permutated_SCCS[2][0];
        data[data_position+10] = hexa_permutated_SCCS[2][1];
        data[data_position+11] = hexa_permutated_SCCS[2][2];
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



      template<int dim>
      std::array<unsigned short int, 3>
      ElasticTensorDecomposition<dim>::indexed_even_permutation(const unsigned short int index)
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



      template<int dim>
      SymmetricTensor<2,3>
      ElasticTensorDecomposition<dim>::compute_voigt_stiffness_tensor(const SymmetricTensor<2,6> &elastic_matrix)
      {
        /**
         * the Voigt stiffness tensor (see Browaeys and chevrot, 2004)
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



      template<int dim>
      SymmetricTensor<2,3>
      ElasticTensorDecomposition<dim>::compute_dilatation_stiffness_tensor(const SymmetricTensor<2,6> &elastic_matrix)
      {
        /**
         * The dilatational stiffness tensor (see Browaeys and chevrot, 2004)
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



      template<int dim>
      Tensor<2,3>
      ElasticTensorDecomposition<dim>::compute_unpermutated_SCCS(const SymmetricTensor<2,3> &dilatation_stiffness_tensor,
                                                                 const SymmetricTensor<2,3> &voigt_stiffness_tensor)
      {
        // computing the eigenvector of the dilation and voigt stiffness matrices and then averaging them by bysection.
        const std::array<std::pair<double,Tensor<1,3,double>>, 3> voigt_eigenvectors_a = eigenvectors(voigt_stiffness_tensor, SymmetricTensorEigenvectorMethod::jacobi);
        const std::array<std::pair<double,Tensor<1,3,double>>, 3> dilatation_eigenvectors_a = eigenvectors(dilatation_stiffness_tensor, SymmetricTensorEigenvectorMethod::jacobi);


        std::vector<Tensor<1,3,double>> unpermutated_SCCS(3);
        // averaging eigenvectors
        // the next function looks for the smallest angle
        // and returns the corresponding vecvo index for that
        // vector.
        size_t NDVC = 0;
        for (size_t i1 = 0; i1 < 3; i1++)
          {
            NDVC = 0;
            double ADVC = 10.0;
            //double SCN = 0.0;
            for (size_t i2 = 0; i2 < 3; i2++)
              {
                double dv_dot_product = dilatation_eigenvectors_a[i1].second*voigt_eigenvectors_a[i2].second;
                // limit the dot product between 1 and -1 so we can use the arccos function safely.
                if (std::abs(dv_dot_product) >= 1.0)
                  dv_dot_product = std::copysign(1.0,dv_dot_product);
                // compute the angle bewteen the vectors and account for that vectors in the oposit
                // direction are the same. So limit them between 0 and 90 degrees such that it
                // represents the minimum angle between the two lines.
                double ADV = dv_dot_product < 0.0 ? std::acos(-1.0)-std::acos(dv_dot_product) : std::acos(dv_dot_product);
                // store this if the angle is smaller
                if (ADV < ADVC)
                  {
                    NDVC=std::copysign(1.0, dv_dot_product)*i2;
                    ADVC = ADV;
                  }
              }

            // Adds/substracting to vecdi the vecvo with the smallest
            // angle times the i2/j index with the sign of SVD to vecdi
            // (effectively turning the eigenvector),and then nomalizing it.
            unpermutated_SCCS[i1] = 0.5*(dilatation_eigenvectors_a[i1].second + (double)NDVC*voigt_eigenvectors_a[std::abs((int)NDVC)].second);
            unpermutated_SCCS[i1] = unpermutated_SCCS[i1]/unpermutated_SCCS[i1].norm();
          }

        return Tensor<2,3>(
        {
          {unpermutated_SCCS[0][0],unpermutated_SCCS[0][1],unpermutated_SCCS[0][2]},
          {unpermutated_SCCS[1][0],unpermutated_SCCS[1][1],unpermutated_SCCS[1][2]},
          {unpermutated_SCCS[2][0],unpermutated_SCCS[2][1],unpermutated_SCCS[2][2]}
        });
      }



      template<int dim>
      std::array<std::array<double,3>,7>
      ElasticTensorDecomposition<dim>::compute_elastic_tensor_SCCS_decompositions(
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
        Tensor<2,3> permutated[3];
        SymmetricTensor<2,6> rotated_elastic_matrix[3];
        std::array<std::array<double,3>,7> norms;

        for (unsigned short int permutation_i = 0; permutation_i < 3; permutation_i++)
          {
            std::array<unsigned short int, 3> perumation = indexed_even_permutation(permutation_i);

            for (size_t j = 0; j < 3; j++)
              {
                permutated[permutation_i][j] = unpermutated_SCCS[perumation[j]];
              }

            rotated_elastic_matrix[permutation_i] = Utilities::Tensors::rotate_voigt_stiffness_matrix((permutated[permutation_i]),elastic_matrix);

            const Tensor<1,21> full_elastic_vector_rotated = Utilities::Tensors::to_voigt_stiffness_vector(rotated_elastic_matrix[permutation_i]);


            const double full_norm_square = full_elastic_vector_rotated.norm_square();
            norms[6][permutation_i] = full_norm_square;

            // The following line would do the same as the lines below, but is is very slow. It has therefore been
            // replaced by the lines below.
            //auto mono_and_higher_vector = projection_matrix_tric_to_mono*full_elastic_vector_rotated;
            dealii::Tensor<1,21> mono_and_higher_vector;
            mono_and_higher_vector[0] = full_elastic_vector_rotated[0];
            mono_and_higher_vector[1] = full_elastic_vector_rotated[1];
            mono_and_higher_vector[2] = full_elastic_vector_rotated[2];
            mono_and_higher_vector[3] = full_elastic_vector_rotated[3];
            mono_and_higher_vector[4] = full_elastic_vector_rotated[4];
            mono_and_higher_vector[5] = full_elastic_vector_rotated[5];
            mono_and_higher_vector[6] = full_elastic_vector_rotated[6];
            mono_and_higher_vector[7] = full_elastic_vector_rotated[7];
            mono_and_higher_vector[8] = full_elastic_vector_rotated[8];
            mono_and_higher_vector[11] = full_elastic_vector_rotated[11];
            mono_and_higher_vector[14] = full_elastic_vector_rotated[14];
            mono_and_higher_vector[17] = full_elastic_vector_rotated[17];
            mono_and_higher_vector[20] = full_elastic_vector_rotated[20];

            auto tric_vector = full_elastic_vector_rotated-mono_and_higher_vector;
            norms[0][permutation_i] = tric_vector.norm_square();

            // The following line would do the same as the lines below, but it is slow. It has therefore been
            // replaced by the lines below.
            //auto ortho_and_higher_vector = projection_matrix_mono_to_ortho*mono_and_higher_vector;
            dealii::Tensor<1,9>  mono_and_higher_vector_cropped;
            mono_and_higher_vector_cropped[0] = mono_and_higher_vector[0];
            mono_and_higher_vector_cropped[1] = mono_and_higher_vector[1];
            mono_and_higher_vector_cropped[2] = mono_and_higher_vector[2];
            mono_and_higher_vector_cropped[3] = mono_and_higher_vector[3];
            mono_and_higher_vector_cropped[4] = mono_and_higher_vector[4];
            mono_and_higher_vector_cropped[5] = mono_and_higher_vector[5];
            mono_and_higher_vector_cropped[6] = mono_and_higher_vector[6];
            mono_and_higher_vector_cropped[7] = mono_and_higher_vector[7];
            mono_and_higher_vector_cropped[8] = mono_and_higher_vector[8];
            dealii::Tensor<1,9> ortho_and_higher_vector;
            ortho_and_higher_vector[0] = mono_and_higher_vector[0];
            ortho_and_higher_vector[1] = mono_and_higher_vector[1];
            ortho_and_higher_vector[2] = mono_and_higher_vector[2];
            ortho_and_higher_vector[3] = mono_and_higher_vector[3];
            ortho_and_higher_vector[4] = mono_and_higher_vector[4];
            ortho_and_higher_vector[5] = mono_and_higher_vector[5];
            ortho_and_higher_vector[6] = mono_and_higher_vector[6];
            ortho_and_higher_vector[7] = mono_and_higher_vector[7];
            ortho_and_higher_vector[8] = mono_and_higher_vector[8];
            auto mono_vector = mono_and_higher_vector_cropped-ortho_and_higher_vector;
            norms[1][permutation_i] = mono_vector.norm_square();


            auto tetra_and_higher_vector = projection_matrix_ortho_to_tetra*ortho_and_higher_vector;
            auto ortho_vector = ortho_and_higher_vector-tetra_and_higher_vector;
            norms[2][permutation_i] = ortho_vector.norm_square();

            auto hexa_and_higher_vector = projection_matrix_tetra_to_hexa*tetra_and_higher_vector;
            auto tetra_vector = tetra_and_higher_vector-hexa_and_higher_vector;
            norms[3][permutation_i] = tetra_vector.norm_square();

            auto iso_vector = projection_matrix_hexa_to_iso*hexa_and_higher_vector;
            auto hexa_vector = hexa_and_higher_vector-iso_vector;
            norms[4][permutation_i] = hexa_vector.norm_square();
            norms[5][permutation_i] = iso_vector.norm_square();

          }
        return norms;

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
        std::vector<std::pair<std::string,unsigned int>> property_information;

        property_information.push_back(std::make_pair("cpo elastic axis e1",3));
        property_information.push_back(std::make_pair("cpo elastic axis e2",3));
        property_information.push_back(std::make_pair("cpo elastic axis e3",3));
        property_information.push_back(std::make_pair("cpo elastic hexagonal axis",3));
        property_information.push_back(std::make_pair("cpo elastic vector norm square",1));
        property_information.push_back(std::make_pair("cpo elastic triclinic vector norm square p1",1));
        property_information.push_back(std::make_pair("cpo elastic triclinic vector norm square p2",1));
        property_information.push_back(std::make_pair("cpo elastic triclinic vector norm square p3",1));
        property_information.push_back(std::make_pair("cpo elastic monoclinic vector norm square p1",1));
        property_information.push_back(std::make_pair("cpo elastic monoclinic vector norm square p2",1));
        property_information.push_back(std::make_pair("cpo elastic monoclinic vector norm square p3",1));
        property_information.push_back(std::make_pair("cpo elastic orthorhombic vector norm square p1",1));
        property_information.push_back(std::make_pair("cpo elastic orthorhombic vector norm square p2",1));
        property_information.push_back(std::make_pair("cpo elastic orthorhombic vector norm square p3",1));
        property_information.push_back(std::make_pair("cpo elastic tetragonal vector norm square p1",1));
        property_information.push_back(std::make_pair("cpo elastic tetragonal vector norm square p2",1));
        property_information.push_back(std::make_pair("cpo elastic tetragonal vector norm square p3",1));
        property_information.push_back(std::make_pair("cpo elastic hexagonal vector norm square p1",1));
        property_information.push_back(std::make_pair("cpo elastic hexagonal vector norm square p2",1));
        property_information.push_back(std::make_pair("cpo elastic hexagonal vector norm square p3",1));
        property_information.push_back(std::make_pair("cpo elastic isotropic vector norm square",1));

        return property_information;
      }



      template <int dim>
      void
      ElasticTensorDecomposition<dim>::declare_parameters (ParameterHandler &)
      {}



      template <int dim>
      void
      ElasticTensorDecomposition<dim>::parse_parameters (ParameterHandler &)
      {}
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
                                        "A plugin in which decomposes the elastic tensor "
                                        "into different approximations and provides the "
                                        "eigenvectors of the tensor.")
    }
  }
}
