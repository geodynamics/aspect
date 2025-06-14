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

#include <aspect/particle/property/cpo_bingham_average_euler.h>
#include <aspect/particle/manager.h>

#include <aspect/utilities.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      CpoBinghamAverageEA<dim>::initialize ()
      {
        const unsigned int my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
        this->random_number_generator.seed(random_number_seed+my_rank);
        const auto &manager = this->get_particle_manager(this->get_particle_manager_index()).get_property_manager();
        AssertThrow(manager.plugin_name_exists("crystal preferred orientation"),
                    ExcMessage("No crystal preferred orientation property plugin found."));

        AssertThrow(manager.check_plugin_order("crystal preferred orientation","cpo bingham average"),
                    ExcMessage("To use the cpo bingham average plugin, the cpo plugin need to be defined before this plugin."));

        cpo_data_position = manager.get_data_info().get_position_by_plugin_index(manager.get_plugin_index_by_name("crystal preferred orientation"));

      }



      template <int dim>
      void
      CpoBinghamAverageEA<dim>::initialize_one_particle_property(const Point<dim> &,
                                                                 std::vector<double> &data) const
      {
        std::vector<double> volume_fractions_grains(n_grains);
        std::vector<Tensor<2,3>> rotation_matrices_grains(n_grains);
        for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
          {
            // create volume fractions and rotation matrix vectors in the order that it is stored in the data array
            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                volume_fractions_grains[grain_i] = cpo_particle_property->get_volume_fractions_grains(cpo_data_position,data,mineral_i,grain_i);
                rotation_matrices_grains[grain_i] = cpo_particle_property->get_rotation_matrix_grains(cpo_data_position,data,mineral_i,grain_i);
              }

            const std::vector<Tensor<2,3>> weighted_rotation_matrices =
              Utilities::rotation_matrices_random_draw_volume_weighting(volume_fractions_grains,
                                                                        rotation_matrices_grains,
                                                                        n_samples,
                                                                        this->random_number_generator);
            const std::array<std::array<double,4>,3> bingham_average = compute_bingham_average(weighted_rotation_matrices);

            for (unsigned int i = 0; i < 3; ++i)
              for (unsigned int j = 0; j < 4; ++j)
                data.emplace_back(bingham_average[i][j]);
          }
      }



      template <int dim>
      void
      CpoBinghamAverageEA<dim>::update_particle_properties(const ParticleUpdateInputs<dim> &/*inputs*/,
                                                           typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        std::vector<double> volume_fractions_grains(n_grains);
        std::vector<Tensor<2,3>> rotation_matrices_grains(n_grains);

        for (auto &particle: particles)
          {
            ArrayView<double> data = particle.get_properties();
            for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
              {
                // create volume fractions and rotation matrix vectors in the order that it is stored in the data array
                for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
                  {
                    volume_fractions_grains[grain_i] = cpo_particle_property->get_volume_fractions_grains(cpo_data_position,data,mineral_i,grain_i);
                    rotation_matrices_grains[grain_i] = cpo_particle_property->get_rotation_matrix_grains(cpo_data_position,data,mineral_i,grain_i);
                  }

                const std::vector<Tensor<2,3>> weighted_rotation_matrices = Utilities::rotation_matrices_random_draw_volume_weighting(volume_fractions_grains, rotation_matrices_grains, n_samples, this->random_number_generator);
                std::array<std::array<double,4>,3> bingham_average = compute_bingham_average(weighted_rotation_matrices);

                for (unsigned int i = 0; i < 3; ++i)
                  for (unsigned int j = 0; j < 4; ++j)
                    {
                      data[this->data_position + mineral_i*12 + i*4 + j] = bingham_average[i][j];
                    }
              }
          }
      }



      template <int dim>
      std::array<std::array<double,4>,3>
      CpoBinghamAverageEA<dim>::compute_bingham_average(std::vector<Tensor<2,3>> matrices) const
      {
        SymmetricTensor<2,3> sum_matrix_a;
        SymmetricTensor<2,3> sum_matrix_b;
        SymmetricTensor<2,3> sum_matrix_c;

        // extracting the a, b and c orientations from the olivine a matrix
        // see https://courses.eas.ualberta.ca/eas421/lecturepages/orientation.html
        const unsigned int n_matrices = matrices.size();
        for (unsigned int i_grain = 0; i_grain < n_matrices; ++i_grain)
          {
            sum_matrix_a[0][0] += matrices[i_grain][0][0] * matrices[i_grain][0][0]; // SUM(l^2)
            sum_matrix_a[1][1] += matrices[i_grain][0][1] * matrices[i_grain][0][1]; // SUM(m^2)
            sum_matrix_a[2][2] += matrices[i_grain][0][2] * matrices[i_grain][0][2]; // SUM(n^2)
            sum_matrix_a[0][1] += matrices[i_grain][0][0] * matrices[i_grain][0][1]; // SUM(l*m)
            sum_matrix_a[0][2] += matrices[i_grain][0][0] * matrices[i_grain][0][2]; // SUM(l*n)
            sum_matrix_a[1][2] += matrices[i_grain][0][1] * matrices[i_grain][0][2]; // SUM(m*n)


            sum_matrix_b[0][0] += matrices[i_grain][1][0] * matrices[i_grain][1][0]; // SUM(l^2)
            sum_matrix_b[1][1] += matrices[i_grain][1][1] * matrices[i_grain][1][1]; // SUM(m^2)
            sum_matrix_b[2][2] += matrices[i_grain][1][2] * matrices[i_grain][1][2]; // SUM(n^2)
            sum_matrix_b[0][1] += matrices[i_grain][1][0] * matrices[i_grain][1][1]; // SUM(l*m)
            sum_matrix_b[0][2] += matrices[i_grain][1][0] * matrices[i_grain][1][2]; // SUM(l*n)
            sum_matrix_b[1][2] += matrices[i_grain][1][1] * matrices[i_grain][1][2]; // SUM(m*n)


            sum_matrix_c[0][0] += matrices[i_grain][2][0] * matrices[i_grain][2][0]; // SUM(l^2)
            sum_matrix_c[1][1] += matrices[i_grain][2][1] * matrices[i_grain][2][1]; // SUM(m^2)
            sum_matrix_c[2][2] += matrices[i_grain][2][2] * matrices[i_grain][2][2]; // SUM(n^2)
            sum_matrix_c[0][1] += matrices[i_grain][2][0] * matrices[i_grain][2][1]; // SUM(l*m)
            sum_matrix_c[0][2] += matrices[i_grain][2][0] * matrices[i_grain][2][2]; // SUM(l*n)
            sum_matrix_c[1][2] += matrices[i_grain][2][1] * matrices[i_grain][2][2]; // SUM(m*n)

          }
        const std::array<std::pair<double,Tensor<1,3,double>>, 3> eigenvectors_a = eigenvectors(sum_matrix_a, SymmetricTensorEigenvectorMethod::jacobi);
        const std::array<std::pair<double,Tensor<1,3,double>>, 3> eigenvectors_b = eigenvectors(sum_matrix_b, SymmetricTensorEigenvectorMethod::jacobi);
        const std::array<std::pair<double,Tensor<1,3,double>>, 3> eigenvectors_c = eigenvectors(sum_matrix_c, SymmetricTensorEigenvectorMethod::jacobi);

        // eigenvalues of all axes, used in the anisotropic viscosity material model to compute Hill's coefficients
        const double eigenvalue_a1 = eigenvectors_a[0].first/matrices.size();
        const double eigenvalue_a2 = eigenvectors_a[1].first/matrices.size();
        const double eigenvalue_a3 = eigenvectors_a[2].first/matrices.size();
        const double eigenvalue_b1 = eigenvectors_b[0].first/matrices.size();
        const double eigenvalue_b2 = eigenvectors_b[1].first/matrices.size();
        const double eigenvalue_b3 = eigenvectors_b[2].first/matrices.size();
        const double eigenvalue_c1 = eigenvectors_c[0].first/matrices.size();
        const double eigenvalue_c2 = eigenvectors_c[1].first/matrices.size();
        const double eigenvalue_c3 = eigenvectors_c[2].first/matrices.size();

        const Tensor<1,3,double> eigvec_a = eigenvectors_a[0].second;
        const Tensor<1,3,double> eigvec_b = eigenvectors_b[0].second;
        const Tensor<1,3,double> eigvec_c = eigenvectors_c[0].second;

        // build rotation matrix from the eigen vectors
        Tensor<2,3> R_CPO;
        R_CPO[0][0] = eigvec_a[0];
        R_CPO[1][0] = eigvec_a[1];
        R_CPO[2][0] = eigvec_a[2];
        R_CPO[0][1] = eigvec_b[0];
        R_CPO[1][1] = eigvec_b[1];
        R_CPO[2][1] = eigvec_b[2];
        R_CPO[0][2] = eigvec_c[0];
        R_CPO[1][2] = eigvec_c[1];
        R_CPO[2][2] = eigvec_c[2];

        // convert rotation matrix to euler angles phi1, theta, phi2
        Tensor<2,3> Rot = transpose(R_CPO);
        Rot = dealii::project_onto_orthogonal_tensors(Rot);
        // save euler angles in radians
        std::array<double,3> EA = Utilities::zxz_euler_angles_from_rotation_matrix(Rot); // in degrees
        const double phi1 = EA[0]*constants::degree_to_radians;
        const double theta = EA[1]*constants::degree_to_radians;
        const double phi2 = EA[2]*constants::degree_to_radians;

        return
        {
          {
            {{phi1, eigenvalue_a1, eigenvalue_a2, eigenvalue_a3}},
            {{theta, eigenvalue_b1, eigenvalue_b2, eigenvalue_b3}},
            {{phi2, eigenvalue_c1, eigenvalue_c2, eigenvalue_c3}}
          }
        };

      }



      template <int dim>
      UpdateTimeFlags
      CpoBinghamAverageEA<dim>::need_update() const
      {
        return update_output_step;
      }



      template <int dim>
      UpdateFlags
      CpoBinghamAverageEA<dim>::get_update_flags (const unsigned int /*component*/) const
      {
        return update_default;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      CpoBinghamAverageEA<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int>> property_information;
        property_information.reserve(6*n_minerals);
        for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
          {
            property_information.emplace_back("cpo mineral ea " + std::to_string(mineral_i) + " phi1",1);
            property_information.emplace_back("cpo mineral ea " + std::to_string(mineral_i) + " eigenvalues a axis",3);
            property_information.emplace_back("cpo mineral ea " + std::to_string(mineral_i) + " theta",1);
            property_information.emplace_back("cpo mineral ea " + std::to_string(mineral_i) + " eigenvalues b axis",3);
            property_information.emplace_back("cpo mineral ea " + std::to_string(mineral_i) + " phi2",1);
            property_information.emplace_back("cpo mineral ea " + std::to_string(mineral_i) + " eigenvalues c axis",3);
          }

        return property_information;
      }



      template <int dim>
      void
      CpoBinghamAverageEA<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("CPO Bingham Average Euler angles");
        {
          prm.declare_entry ("Random number seed", "1",
                             Patterns::Integer (0),
                             "The seed used to generate random numbers. This will make sure that "
                             "results are reproducible as long as the problem is run with the "
                             "same amount of MPI processes. It is implemented as final seed = "
                             "Random number seed + MPI Rank. ");

          prm.declare_entry ("Number of samples", "0",
                             Patterns::Double(0),
                             "This determines how many samples are taken when using the random "
                             "draw volume averaging. Setting it to zero means that the number of "
                             "samples is set to be equal to the number of grains.");
        }
        prm.leave_subsection ();
      }



      template <int dim>
      void
      CpoBinghamAverageEA<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("CPO Bingham Average Euler angles");
        {
          // Get a pointer to the CPO particle property.
          cpo_particle_property = std::make_unique<const Particle::Property::CrystalPreferredOrientation<dim>> (
                                    this->get_particle_manager(this->get_particle_manager_index()).get_property_manager().template get_matching_active_plugin<Particle::Property::CrystalPreferredOrientation<dim>>());

          random_number_seed = prm.get_integer ("Random number seed");
          n_grains = cpo_particle_property->get_number_of_grains();
          n_minerals = cpo_particle_property->get_number_of_minerals();
          n_samples = prm.get_integer("Number of samples");
          if (n_samples == 0)
            n_samples = n_grains;
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(CpoBinghamAverageEA,
                                        "cpo bingham average euler angle",
                                        "This is a particle property plugin which computes the Bingham "
                                        "average for the Crystal Preferred Orientation particle property "
                                        "plugin so that it can be visualized. The orientations are saved "
                                        "as Euler angles.")
    }
  }
}
