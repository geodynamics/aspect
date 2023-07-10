/*
  Copyright (C) 2022 - 2023 by the authors of the ASPECT code.

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

#include <aspect/particle/property/crystal_preferred_orientation.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <world_builder/grains.h>
#include <world_builder/world.h>

#include <aspect/citation_info.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {

      template <int dim>
      CrystalPreferredOrientation<dim>::CrystalPreferredOrientation ()
      {
        permutation_operator_3d[0][1][2]  = 1;
        permutation_operator_3d[1][2][0]  = 1;
        permutation_operator_3d[2][0][1]  = 1;
        permutation_operator_3d[0][2][1]  = -1;
        permutation_operator_3d[1][0][2]  = -1;
        permutation_operator_3d[2][1][0]  = -1;
      }



      template <int dim>
      void
      CrystalPreferredOrientation<dim>::initialize ()
      {
        CitationInfo::add("CPO");
        const unsigned int my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
        this->random_number_generator.seed(random_number_seed+my_rank);

        // Don't assert when called by the unit tester.
        if (this->simulator_is_past_initialization())
          {
            AssertThrow(this->introspection().compositional_name_exists("water"),
                        ExcMessage("Particle property CPO only works if"
                                   "there is a compositional field called water."));
            water_index = this->introspection().compositional_index_for_name("water");
          }
      }



      template <int dim>
      void
      CrystalPreferredOrientation<dim>::compute_random_rotation_matrix(Tensor<2,3> &rotation_matrix) const
      {

        const double dt = this -> get_time();
        // This function is based on an article in Graphic Gems III, written by James Arvo, Cornell University (p 116-120).
        // The original code can be found on  http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
        // and is licenced according to this website with the following licence:
        //
        // "The Graphics Gems code is copyright-protected. In other words, you cannot claim the text of the code as your own and
        // resell it. Using the code is permitted in any program, product, or library, non-commercial or commercial. Giving credit
        // is not required, though is a nice gesture. The code comes as-is, and if there are any flaws or problems with any Gems
        // code, nobody involved with Gems - authors, editors, publishers, or webmasters - are to be held responsible. Basically,
        // don't be a jerk, and remember that anything free comes with no guarantee.""
        //
        // The book states in the preface the following: "As in the first two volumes, all of the C and C++ code in this book is in
        // the public domain, and is yours to study, modify, and use."

        // first generate three random numbers between 0 and 1 and multiply them with 2 PI or 2 for z. Note that these are not the same as phi_1, theta and phi_2.

        boost::random::uniform_real_distribution<double> uniform_distribution(0,1);
        double one = uniform_distribution(this->random_number_generator);
        double two = uniform_distribution(this->random_number_generator);
        double three = uniform_distribution(this->random_number_generator);

        double theta = 2.0 * numbers::PI * one; // Rotation about the pole (Z)
        double phi = 2.0 * numbers::PI * two; // For direction of pole deflection.
        double z = dt != 0.0 ? 0.1 * three : 2.0 * three; //For magnitude of pole deflection.

        // Compute a vector V used for distributing points over the sphere
        // via the reflection I - V Transpose(V).  This formulation of V
        // will guarantee that if x[1] and x[2] are uniformly distributed,
        // the reflected points will be uniform on the sphere.  Note that V
        // has length sqrt(2) to eliminate the 2 in the Householder matrix.

        double r  = std::sqrt( z );
        double Vx = std::sin( phi ) * r;
        double Vy = std::cos( phi ) * r;
        double Vz = std::sqrt( 2.f - z );

        // Compute the row vector S = Transpose(V) * R, where R is a simple
        // rotation by theta about the z-axis.  No need to compute Sz since
        // it's just Vz.

        double st = std::sin( theta );
        double ct = std::cos( theta );
        double Sx = Vx * ct - Vy * st;
        double Sy = Vx * st + Vy * ct;

        // Construct the rotation matrix  ( V Transpose(V) - I ) R, which
        // is equivalent to V S - R.

        rotation_matrix[0][0] = Vx * Sx - ct;
        rotation_matrix[0][1] = Vx * Sy - st;
        rotation_matrix[0][2] = Vx * Vz;
        rotation_matrix[1][0] = Vy * Sx + st;
        rotation_matrix[1][1] = Vy * Sy - ct;
        rotation_matrix[1][2] = Vy * Vz;
        rotation_matrix[2][0] = Vz * Sx;
        rotation_matrix[2][1] = Vz * Sy;
        rotation_matrix[2][2] = 1.0 - z;   // This equals Vz * Vz - 1.0
      }



      template <int dim>
      void
      CrystalPreferredOrientation<dim>::initialize_one_particle_property(const Point<dim> &,
                                                                         std::vector<double> &data) const
      {
        // the layout of the data vector per perticle is the following:
        // 1. M mineral times
        //    1.1  olivine deformation type   -> 1 double, at location
        //                                      => data_position + 0 + mineral_i * (n_grains * 10 + 2)
        //    2.1. Mineral volume fraction    -> 1 double, at location
        //                                      => data_position + 1 + mineral_i *(n_grains * 10 + 2)
        //    2.2. N grains times:
        //         2.1. volume fraction grain -> 1 double, at location:
        //                                      => data_position + 2 + i_grain * 10 + mineral_i *(n_grains * 10 + 2), or
        //                                      => data_position + 2 + i_grain * (2 * Tensor<2,3>::n_independent_components+ 2) + mineral_i * (n_grains * 10 + 2)
        //         2.2. rotation matrix grain -> 9 (Tensor<2,dim>::n_independent_components) doubles, starts at:
        //                                      => data_position + 3 + i_grain * 10 + mineral_i * (n_grains * 10 + 2), or
        //                                      => data_position + 3 + i_grain * (2 * Tensor<2,3>::n_independent_components+ 2) + mineral_i * (n_grains * 10 + 2)
        //
        // Note that we store exactly the same number of grains of all minerals (e.g. olivine and enstatite
        // grains), although their volume fractions may not be the same. We need a minimum amount
        // of grains per tracer to perform reliable statistics on it. This minimum is the same for all phases.
        // and enstatite.
        //
        // Furthermore, for this plugin the following dims are always 3. When using 2d an infinitely thin 3d domain is assumed.
        //
        // The rotation matrix is a direction cosine matrix, representing the orientation of the grain in the domain.
        // The fabric is determined later in the computations, so initialize it to -1.
        std::vector<double> deformation_type(n_minerals, -1.0);
        std::vector<std::vector<double >>volume_fractions_grains(n_minerals);
        std::vector<std::vector<Tensor<2,3>>> rotation_matrices_grains(n_minerals);

        for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
          {
            volume_fractions_grains[mineral_i].resize(n_grains);
            rotation_matrices_grains[mineral_i].resize(n_grains);

            // This will be set by the initial grain subsection.
            bool use_world_builder = false;
            if (use_world_builder)
              {
#ifdef ASPECT_WITH_WORLD_BUILDER
                AssertThrow(false,
                            ExcMessage("Not implemented."))
#else
                AssertThrow(false,
                            ExcMessage("The world builder was requested but not provided. Make sure that aspect is "
                                       "compiled with the World Builder and that you provide a world builder file in the input."))
#endif
              }
            else
              {
                // set volume fraction

                /*
                To use D-Rex uncomment the first statement, setting the initial volume fraction = 1/n_grains
                To use D-Rex ++ uncomment the other lines
                */

                //const double initial_volume_fraction = 1.0/n_grains;
                const double large_grains_size = 1e-4;
                const double small_grain_size = large_grains_size/1000.;
                const double initial_volume_fraction_large = (4./3.)*numbers::PI*pow(0.5*large_grains_size,3);
                const double initial_volume_fraction_small = (4./3.)*numbers::PI*pow(0.5*small_grain_size,3);

                for (unsigned int grain_i = 0; grain_i < n_grains ; ++grain_i)
                  {
                    //set volume fraction
                    if (grain_i < n_grains/10.)
                      {
                        volume_fractions_grains[mineral_i][grain_i] = initial_volume_fraction_large;
                      }
                    else
                      {
                        volume_fractions_grains[mineral_i][grain_i] = initial_volume_fraction_small;
                      }

                    // set a uniform random rotation_matrix per grain
                    // volume_fractions_grains[mineral_i][grain_i] = initial_volume_fraction;
                    this->compute_random_rotation_matrix(rotation_matrices_grains[mineral_i][grain_i]);
                  }
              }
          }

        for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
          {
            data.emplace_back(deformation_type[mineral_i]);
            data.emplace_back(volume_fractions_minerals[mineral_i]);
            for (unsigned int grain_i = 0; grain_i < n_grains ; ++grain_i)
              {
                data.emplace_back(volume_fractions_grains[mineral_i][grain_i]);
                for (unsigned int i = 0; i < Tensor<2,3>::n_independent_components ; ++i)
                  {
                    const dealii::TableIndices<2> index = Tensor<2,3>::unrolled_to_component_indices(i);
                    data.emplace_back(rotation_matrices_grains[mineral_i][grain_i][index]);
                  }
              }
          }
      }



      template <int dim>
      void
      CrystalPreferredOrientation<dim>::update_one_particle_property(const unsigned int data_position,
                                                                     const Point<dim> &position,
                                                                     const Vector<double> &solution,
                                                                     const std::vector<Tensor<1,dim>> &gradients,
                                                                     const ArrayView<double> &data) const
      {

        // STEP 1: Load data and preprocess it.

        // need access to the pressure, viscosity,
        // get velocity
        Tensor<1,dim> velocity;
        for (unsigned int i = 0; i < dim; ++i)
          velocity[i] = solution[this->introspection().component_indices.velocities[i]];

        // get velocity gradient tensor.
        Tensor<2,dim> velocity_gradient;
        for (unsigned int d=0; d<dim; ++d)
          velocity_gradient[d] = gradients[d];

        // Calculate strain rate from velocity gradients
        const SymmetricTensor<2,dim> strain_rate = symmetrize (velocity_gradient);
        const SymmetricTensor<2,dim> deviatoric_strain_rate
          = (this->get_material_model().is_compressible()
             ?
             strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
             :
             strain_rate);

        const double pressure = solution[this->introspection().component_indices.pressure];
        const double temperature = solution[this->introspection().component_indices.temperature];
        const double water_content = solution[this->introspection().component_indices.compositional_fields[water_index]];

        // get the composition of the particle
        std::vector<double> compositions;
        for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
          {
            const unsigned int solution_component = this->introspection().component_indices.compositional_fields[i];
            compositions.push_back(solution[solution_component]);
          }

        const double dt = this->get_timestep();

        // even in 2d we need 3d strain-rates and velocity gradient tensors. So we make them 3d by
        // adding an extra dimension which is zero.
        SymmetricTensor<2,3> strain_rate_3d;
        strain_rate_3d[0][0] = strain_rate[0][0];
        strain_rate_3d[0][1] = strain_rate[0][1];
        //sym: strain_rate_3d[1][0] = strain_rate[1][0];
        strain_rate_3d[1][1] = strain_rate[1][1];

        if (dim == 3)
          {
            strain_rate_3d[0][2] = strain_rate[0][2];
            strain_rate_3d[1][2] = strain_rate[1][2];
            //sym: strain_rate_3d[2][0] = strain_rate[0][2];
            //sym: strain_rate_3d[2][1] = strain_rate[1][2];
            strain_rate_3d[2][2] = strain_rate[2][2];
          }
        Tensor<2,3> velocity_gradient_3d;
        velocity_gradient_3d[0][0] = velocity_gradient[0][0];
        velocity_gradient_3d[0][1] = velocity_gradient[0][1];
        velocity_gradient_3d[1][0] = velocity_gradient[1][0];
        velocity_gradient_3d[1][1] = velocity_gradient[1][1];
        if (dim == 3)
          {
            velocity_gradient_3d[0][2] = velocity_gradient[0][2];
            velocity_gradient_3d[1][2] = velocity_gradient[1][2];
            velocity_gradient_3d[2][0] = velocity_gradient[2][0];
            velocity_gradient_3d[2][1] = velocity_gradient[2][1];
            velocity_gradient_3d[2][2] = velocity_gradient[2][2];
          }

        for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
          {

            /**
            * Now we have loaded all the data and can do the actual computation.
            * The computation consists of two parts. The first part is computing
            * the derivatives for the directions and grain sizes. Then those
            * derivatives are used to advect the particle properties.
            */
            double sum_volume_mineral = 0;
            std::pair<std::vector<double>, std::vector<Tensor<2,3>>>
            derivatives_grains = this->compute_derivatives(data_position,
                                                           data,
                                                           mineral_i,
                                                           strain_rate_3d,
                                                           velocity_gradient_3d,
                                                           position,
                                                           temperature,
                                                           pressure,
                                                           velocity,
                                                           compositions,
                                                           strain_rate,
                                                           deviatoric_strain_rate,
                                                           water_content);

            switch (advection_method)
              {
                case AdvectionMethod::forward_euler:

                  sum_volume_mineral = this->advect_forward_euler(data_position,
                                                                  data,
                                                                  mineral_i,
                                                                  dt,
                                                                  derivatives_grains);

                  break;

                case AdvectionMethod::backward_euler:
                  sum_volume_mineral = this->advect_backward_euler(data_position,
                                                                   data,
                                                                   mineral_i,
                                                                   dt,
                                                                   derivatives_grains);

                  break;
              }

            // normalize the volume fractions back to a total of 1 for each mineral
            const double inv_sum_volume_mineral = 1.0;///sum_volume_mineral;

            Assert(std::isfinite(inv_sum_volume_mineral),
                   ExcMessage("inv_sum_volume_mineral is not finite. sum_volume_enstatite = "
                              + std::to_string(sum_volume_mineral)));

            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                const double volume_fraction_grains = get_volume_fractions_grains(data_position,data,mineral_i,grain_i)*inv_sum_volume_mineral;
                set_volume_fractions_grains(data_position,data,mineral_i,grain_i,volume_fraction_grains);
                Assert(isfinite(get_volume_fractions_grains(data_position,data,mineral_i,grain_i)),
                       ExcMessage("volume_fractions_grains[mineral_i]" + std::to_string(grain_i) + "] is not finite: "
                                  + std::to_string(get_volume_fractions_grains(data_position,data,mineral_i,grain_i)) + ", inv_sum_volume_mineral = "
                                  + std::to_string(inv_sum_volume_mineral) + "."));

                /**
                 * Correct direction rotation matrices numerical error (orthnormality) after integration
                 * Follows same method as in matlab version from Thissen (see https://github.com/cthissen/Drex-MATLAB/)
                 * of finding the nearest orthonormal matrix using the SVD
                 */
                Tensor<2,3> rotation_matrix = get_rotation_matrix_grains(data_position,data,mineral_i,grain_i);
                for (size_t i = 0; i < 3; i++)
                  for (size_t j = 0; j < 3; j++)
                    Assert(abs(rotation_matrix[i][j]) <= 1.0,
                           ExcMessage("rotation_matrix[" + std::to_string(i) + "][" + std::to_string(j) +
                                      "] is larger than one: " + std::to_string(rotation_matrix[i][j]) + " (" + std::to_string(rotation_matrix[i][j]-1.0) + "). rotation_matrix = \n"
                                      + std::to_string(rotation_matrix[0][0]) + " " + std::to_string(rotation_matrix[0][1]) + " " + std::to_string(rotation_matrix[0][2]) + "\n"
                                      + std::to_string(rotation_matrix[1][0]) + " " + std::to_string(rotation_matrix[1][1]) + " " + std::to_string(rotation_matrix[1][2]) + "\n"
                                      + std::to_string(rotation_matrix[2][0]) + " " + std::to_string(rotation_matrix[2][1]) + " " + std::to_string(rotation_matrix[2][2])));

                for (size_t i = 0; i < 3; ++i)
                  {
                    for (size_t j = 0; j < 3; ++j)
                      {
                        Assert(!std::isnan(rotation_matrix[i][j]), ExcMessage("rotation_matrix is nan before orthogonalization."));
                      }
                  }

                rotation_matrix = dealii::project_onto_orthogonal_tensors(rotation_matrix);
                for (size_t i = 0; i < 3; ++i)
                  for (size_t j = 0; j < 3; ++j)
                    {
                      // I don't think this should happen with the projection, but D-Rex
                      // does not do the orthogonal projection, but just clamps the values
                      // to 1 and -1.
                      Assert(std::fabs(rotation_matrix[i][j]) <= 1.0,
                             ExcMessage("The rotation_matrix has a entry larger than 1."));

                      Assert(!std::isnan(rotation_matrix[i][j]),
                             ExcMessage("rotation_matrix is nan after orthoganalization: "
                                        + std::to_string(rotation_matrix[i][j])));

                      Assert(abs(rotation_matrix[i][j]) <= 1.0,
                             ExcMessage("3. rotation_matrix[" + std::to_string(i) + "][" + std::to_string(j) +
                                        "] is larger than one: "
                                        + std::to_string(rotation_matrix[i][j]) + " (" + std::to_string(rotation_matrix[i][j]-1.0) + "). rotation_matrix = \n"
                                        + std::to_string(rotation_matrix[0][0]) + " " + std::to_string(rotation_matrix[0][1]) + " " + std::to_string(rotation_matrix[0][2]) + "\n"
                                        + std::to_string(rotation_matrix[1][0]) + " " + std::to_string(rotation_matrix[1][1]) + " " + std::to_string(rotation_matrix[1][2]) + "\n"
                                        + std::to_string(rotation_matrix[2][0]) + " " + std::to_string(rotation_matrix[2][1]) + " " + std::to_string(rotation_matrix[2][2])));
                    }

                set_rotation_matrix_grains(data_position,data,mineral_i,grain_i,rotation_matrix);
              }
          }
      }



      template <int dim>
      UpdateTimeFlags
      CrystalPreferredOrientation<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      InitializationModeForLateParticles
      CrystalPreferredOrientation<dim>::late_initialization_mode () const
      {
        return InitializationModeForLateParticles::interpolate;
      }



      template <int dim>
      UpdateFlags
      CrystalPreferredOrientation<dim>::get_needed_update_flags () const
      {
        return update_values | update_gradients;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      CrystalPreferredOrientation<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int>> property_information;

        for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
          {
            property_information.push_back(std::make_pair("cpo mineral " + std::to_string(mineral_i) + " type",1));
            property_information.push_back(std::make_pair("cpo mineral " + std::to_string(mineral_i) + " volume fraction",1));

            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                property_information.push_back(std::make_pair("cpo mineral " + std::to_string(mineral_i) + " grain " + std::to_string(grain_i) + " volume fraction",1));

                for (unsigned int index = 0; index < Tensor<2,3>::n_independent_components; ++index)
                  {
                    property_information.push_back(std::make_pair("cpo mineral " + std::to_string(mineral_i) + " grain " + std::to_string(grain_i) + " rotation_matrix " + std::to_string(index),1));
                  }
              }
          }

        return property_information;
      }



      template <int dim>
      double
      CrystalPreferredOrientation<dim>::advect_forward_euler(const unsigned int cpo_index,
                                                             const ArrayView<double> &data,
                                                             const unsigned int mineral_i,
                                                             const double dt,
                                                             const std::pair<std::vector<double>, std::vector<Tensor<2,3>>> &derivatives) const
      {
        double sum_volume_fractions = 0;
        Tensor<2,3> rotation_matrix;
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            // Do the volume fraction of the grain
            Assert(std::isfinite(get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i)),ExcMessage("volume_fractions[grain_i] is not finite before it is set."));
            double volume_fraction_grains = get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i);
            volume_fraction_grains =  volume_fraction_grains + dt * volume_fraction_grains * derivatives.first[grain_i];
            set_volume_fractions_grains(cpo_index,data,mineral_i,grain_i, volume_fraction_grains);
            Assert(std::isfinite(get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i)),ExcMessage("volume_fractions[grain_i] is not finite. grain_i = "
                   + std::to_string(grain_i) + ", volume_fractions[grain_i] = " + std::to_string(get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i))
                   + ", derivatives.first[grain_i] = " + std::to_string(derivatives.first[grain_i])));

            sum_volume_fractions += get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i);

            // Do the rotation matrix for this grain
            rotation_matrix = get_rotation_matrix_grains(cpo_index,data,mineral_i,grain_i);
            rotation_matrix += dt * rotation_matrix * derivatives.second[grain_i];
            set_rotation_matrix_grains(cpo_index,data,mineral_i,grain_i,rotation_matrix);
          }

        Assert(sum_volume_fractions != 0, ExcMessage("The sum of all grain volume fractions of a mineral is equal to zero. This should not happen."));
        return sum_volume_fractions;
      }



      template <int dim>
      double
      CrystalPreferredOrientation<dim>::advect_backward_euler(const unsigned int cpo_index,
                                                              const ArrayView<double> &data,
                                                              const unsigned int mineral_i,
                                                              const double dt,
                                                              const std::pair<std::vector<double>, std::vector<Tensor<2,3>>> &derivatives) const
      {
        double sum_volume_fractions = 0;
        Tensor<2,3> cosine_ref;
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            // Do the volume fraction of the grain
            double vf_old = get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i);
            double vf_new = get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i);
            Assert(std::isfinite(vf_new),ExcMessage("vf_new is not finite before it is set."));
            for (size_t iteration = 0; iteration < property_advection_max_iterations; ++iteration)
              {
                Assert( std::isfinite(vf_new),ExcMessage("vf_new is not finite before it is set. grain_i = "
                                                         + std::to_string(grain_i) + ", volume_fractions[grain_i] = " + std::to_string(get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i))
                                                         + ", derivatives.first[grain_i] = " + std::to_string(derivatives.first[grain_i])));

                vf_new = get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i) + (dt * vf_new * derivatives.first[grain_i]);

                Assert(std::isfinite(get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i)),ExcMessage("volume_fractions[grain_i] is not finite. grain_i = "
                       + std::to_string(grain_i) + ", volume_fractions[grain_i] = " + std::to_string(get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i))
                       + ", derivatives.first[grain_i] = " + std::to_string(derivatives.first[grain_i])));
                if (std::fabs(vf_new-vf_old) < property_advection_tolerance)
                  {
                    break;
                  }
                vf_old = vf_new;
              }

            set_volume_fractions_grains(cpo_index,data,mineral_i,grain_i,vf_new);
            sum_volume_fractions += vf_new;

            // Do the rotation matrix for this grain
            cosine_ref = get_rotation_matrix_grains(cpo_index,data,mineral_i,grain_i);
            Tensor<2,3> cosine_old = cosine_ref;
            Tensor<2,3> cosine_new = cosine_ref;

            for (size_t iteration = 0; iteration < property_advection_max_iterations; ++iteration)
              {
                cosine_new = cosine_ref + (dt * cosine_new * derivatives.second[grain_i]);

                if ((cosine_new-cosine_old).norm() < property_advection_tolerance)
                  {
                    break;
                  }
                cosine_old = cosine_new;
              }

            set_rotation_matrix_grains(cpo_index,data,mineral_i,grain_i,cosine_new);

          }


        Assert(sum_volume_fractions != 0, ExcMessage("The sum of all grain volume fractions of a mineral is equal to zero. This should not happen."));
        return sum_volume_fractions;
      }



      template <int dim>
      std::pair<std::vector<double>, std::vector<Tensor<2,3>>>
      CrystalPreferredOrientation<dim>::compute_derivatives(const unsigned int cpo_index,
                                                            const ArrayView<double> &data,
                                                            const unsigned int mineral_i,
                                                            const SymmetricTensor<2,3> &strain_rate_3d,
                                                            const Tensor<2,3> &velocity_gradient_tensor,
                                                            const Point<dim> &position,
                                                            const double temperature,
                                                            const double pressure,
                                                            const Tensor<1,dim> &velocity,
                                                            const std::vector<double> &compositions,
                                                            const SymmetricTensor<2,dim> &strain_rate,
                                                            const SymmetricTensor<2,dim> &deviatoric_strain_rate,
                                                            const double water_content) const
      {
        std::pair<std::vector<double>, std::vector<Tensor<2,3>>> derivatives;
        switch (cpo_derivative_algorithm)
          {
            case CPODerivativeAlgorithm::spin_tensor:
            {
              return compute_derivatives_spin_tensor(velocity_gradient_tensor);
              break;
            }
            case CPODerivativeAlgorithm::drex_2004:
            {

              const DeformationType deformation_type = determine_deformation_type(deformation_type_selector[mineral_i],
                                                                                  position,
                                                                                  temperature,
                                                                                  pressure,
                                                                                  velocity,
                                                                                  compositions,
                                                                                  strain_rate,
                                                                                  deviatoric_strain_rate,
                                                                                  water_content);

              set_deformation_type(cpo_index,data,mineral_i,static_cast<unsigned int>(deformation_type));

              const std::array<double,4> ref_resolved_shear_stress = reference_resolved_shear_stress_from_deformation_type(deformation_type);

              return compute_derivatives_drex_2004(cpo_index,
                                                   data,
                                                   mineral_i,
                                                   strain_rate_3d,
                                                   velocity_gradient_tensor,
                                                   ref_resolved_shear_stress);
              break;
            }
            case CPODerivativeAlgorithm::drexpp:
            {
              const double pressure = 2.0e9;
              const DeformationType deformation_type = determine_deformation_type(deformation_type_selector[mineral_i],
                                                                                  position,
                                                                                  temperature,
                                                                                  pressure,
                                                                                  velocity,
                                                                                  compositions,
                                                                                  strain_rate,
                                                                                  deviatoric_strain_rate,
                                                                                  water_content);

              set_deformation_type(cpo_index,data,mineral_i,static_cast<unsigned int>(deformation_type));

              /*
                The following part of the algorithm describes the rate of recrystalization of new strain-free grains. This is done by using the
                Johnson- Avrami-Mehl_Kolmogodorov theory of transformation. The equations and the values for the variables are taken from Cross and Skemer
                (2019).

                Incremental recrystalization fraction = n * beta * (strain - critical strain)^(n - 1) * exp( -beta * (strain - critical strain)^ n )
                 where
                 n -> avrami exponent
                 beta -> rate of transformation
                  where beta -> C * exp(g * (T/T_melt))
                  C and g are constantsS
                 critical strain -> strain at which dynamic recrystalization start
              */

              const double strain = this->get_time() * std::sqrt(std::max(-second_invariant(deviatoric_strain_rate), 0.));
              const double const_C = exp(-10.0);
              const double const_g = 13.8;
              const double homologous_T = 1770+273.15; // Have to take this formulation from Katz et al (2003).
              const double avrami_exponet = 1.48;
              const double strain_critical = 0.25;
              const double rate_of_transformation = const_C*exp(const_g*temperature/homologous_T);
              double aggregate_recrystalization_increment;
              std::cout<<"strain at timestep = "<<strain<<std::endl;
              if (strain_critical < strain)
                aggregate_recrystalization_increment = avrami_exponet * rate_of_transformation * std::pow((strain - strain_critical),avrami_exponet - 1) * exp( -1 * rate_of_transformation * std::pow(( strain - strain_critical), avrami_exponet)); // TODO: compute actual value with equation 5
              else
                aggregate_recrystalization_increment = 0;
              const std::array<double,4> ref_resolved_shear_stress = reference_resolved_shear_stress_from_deformation_type(deformation_type);


              // Compute the diffusion strain with the diffusion viscosity
              const double gravity_norm = this->get_gravity_model().gravity_vector(position).norm();
              const double reference_density = (this->get_adiabatic_conditions().is_initialized())
                                               ?
                                               this->get_adiabatic_conditions().density(position)
                                               :
                                               3300.;

              // The phase index is set to invalid_unsigned_int, because it is only used internally
              // in phase_average_equation_of_state_outputs to loop over all existing phases
              MaterialModel::MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(temperature,
                                                                                      pressure,
                                                                                      this->get_geometry_model().depth(position),
                                                                                      gravity_norm*reference_density,
                                                                                      numbers::invalid_unsigned_int);
              std::vector<double> phase_function_values(phase_function.n_phase_transitions(), 0.0);
              // Compute value of phase functions
              for (unsigned int j=0; j < phase_function.n_phase_transitions(); j++)
                {
                  phase_inputs.phase_index = j;
                  phase_function_values[j] = phase_function.compute_value(phase_inputs);
                }

              const std::vector<double> volume_fractions = MaterialModel::MaterialUtilities::compute_composition_fractions(compositions);
              // the diffusion_pre_viscosities is the diffusion viscosity withouth the grainsize
              std::vector<double> diffusion_pre_viscosities(volume_fractions.size(),std::numeric_limits<double>::quiet_NaN());
              std::vector<double> diffusion_grain_size_exponent(volume_fractions.size(),std::numeric_limits<double>::quiet_NaN());
              std::vector<double> dislocation_viscosities(volume_fractions.size(),std::numeric_limits<double>::quiet_NaN());
              const double strain_rate_inv = std::max(std::sqrt(std::max(-second_invariant(deviatoric_strain_rate), 0.)),min_strain_rate);
              for (unsigned int composition = 0; composition < volume_fractions.size(); composition++)
                {
                  const MaterialModel::Rheology::DiffusionCreepParameters p_dif = rheology_diff->compute_creep_parameters(composition,
                                                                                  phase_function_values,
                                                                                  phase_function.n_phase_transitions_for_each_composition());
                  //const double grain_size = 1e-3;
                  // Power law creep equation
                  //    viscosity = 0.5 * A^(-1) * d^(m) * exp((E + P*V)/(RT))
                  // A: prefactor,
                  // d: grain size, m: grain size exponent, E: activation energy, P: pressure,
                  // V; activation volume, R: gas constant, T: temperature.
                  diffusion_pre_viscosities[composition] = 0.5 / p_dif.prefactor *
                                                           std::exp((p_dif.activation_energy +
                                                                     pressure*p_dif.activation_volume)/
                                                                    (constants::gas_constant*temperature));

                  diffusion_grain_size_exponent[composition] = p_dif.grain_size_exponent;

                  const MaterialModel::Rheology::DislocationCreepParameters p_dis = rheology_disl->compute_creep_parameters(composition,
                                                                                    phase_function_values,
                                                                                    phase_function.n_phase_transitions_for_each_composition());

                  // Power law creep equation
                  //    viscosity = 0.5 * A^(-1) * d^(m) * exp((E + P*V)/(RT))
                  // A: prefactor,
                  // d: grain size, m: grain size exponent, E: activation energy, P: pressure,
                  // V; activation volume, R: gas constant, T: temperature.

                  dislocation_viscosities[composition] = 0.5 * std::pow(p_dis.prefactor,-1/p_dis.stress_exponent) *
                                                         std::exp((p_dis.activation_energy + pressure*p_dis.activation_volume)/
                                                                  (constants::gas_constant*temperature*p_dis.stress_exponent)) *
                                                         std::pow(strain_rate_inv,((1. - p_dis.stress_exponent)/p_dis.stress_exponent));

                }

              // now compute the normal viscosity to be able to computes the stress
              // Create the material model inputs and outputs to
              // retrieve the current viscosity.

              MaterialModel::MaterialModelInputs<dim> in = MaterialModel::MaterialModelInputs<dim>(1,compositions.size());
              in.pressure[0] = pressure;
              in.temperature[0] = temperature;
              in.position[0] = position;
              in.strain_rate[0] = strain_rate;
              in.composition[0] = compositions;

              in.requested_properties = MaterialModel::MaterialProperties::viscosity;

              MaterialModel::MaterialModelOutputs<dim> out(1,
                                                           this->n_compositional_fields());

              this->get_material_model().evaluate(in, out);

              // Compressive stress is positive in geoscience applications.
              SymmetricTensor<2, dim> stress = pressure * unit_symmetric_tensor<dim>();

              // Add elastic stresses if existent.
              AssertThrow(this->get_parameters().enable_elasticity == false, ExcMessage("Elasticity not supported when computing the CPO stress"));

              const double eta = out.viscosities[0];

              stress += -2. * eta * deviatoric_strain_rate;


              // Compute the deviatoric stress tensor after elastic stresses were added.
              const SymmetricTensor<2, dim> deviatoric_stress = deviator(stress);

              // Compute the second moment invariant of the deviatoric stress
              // in the same way as the second moment invariant of the deviatoric
              // strain rate is computed in the viscoplastic material model.
              // TODO check that this is valid for the compressible case.
              //const double stress_invariant = (std::sqrt(std::max(-second_invariant(deviatoric_stress), 0.)));

              const std::array< double, dim > eigenvalues = dealii::eigenvalues(deviatoric_stress);
              double differential_stress = eigenvalues[0]-eigenvalues[dim-1];
              const double dislocation_creep_strain_rate = differential_stress/(2.0*dislocation_viscosities[0]);

              /*
               Calculation of the recrystalized grain gize using the paleowattmeter from Austin and Evans (2009)
               and the values are taken from Dannberg et al (2017)

               The equation is of the form
               d_recrystalized = [(pre_exponent_numerator)/(pre_exponent_denominator)]*exp(exponent_term)^(1/1+p_g)
               where
               pre_exponent_num = c * gamma * k_g
                where
                 c -> geometric constant
                 gamma -> average specific grain boundary energy
                 k_g -> experimentally derived prefactor

               pre_exponent_den = lambda * stress * dislocation_strain_rate * p_g
                where
                 lambda -> fraction of work
                 p_g -> grain growth exponent

               exponent term = -( E_g + P * V_g )/( R * T )
               E_g -> Activation energy for grain growth
               V_g -> activation volume for grain growth
               R -> boltzmann constant
               T -> temperature
              */

              long double recrystalized_grain_size;
              const double t = this -> get_time();

              // The terms related to the paleowattmeter are declared here
              // Note : These values prescribe to the piezometer for olivine alone

              const double C      = 3;
              const double gamma  = 1.0;
              const double k_g    = 1.92 * std::pow(10,-10);
              const double lambda = 0.1;
              const double p_g    = 3;
              const double E_g    = 400 * std::pow(10,3);
              const double V_g    = 0;
              const double R      = 8.314;

              if (t!=0)
                {
                  const long double pre_exponent_term_num = C * gamma * k_g;
                  const long double pre_exponent_term_den = lambda * differential_stress * dislocation_creep_strain_rate * p_g;
                  const long double pre_exponent_term = pre_exponent_term_num/pre_exponent_term_den;
                  const long double exponent_term = (E_g + pressure * V_g) / (R * temperature);
                  recrystalized_grain_size =  std::pow(pre_exponent_term*std::exp(-1 * exponent_term),1/(1+p_g));
                }
              else
                {
                  recrystalized_grain_size = 0.5;
                }

              const long double half_recrystalized_grain_size = 0.5 * recrystalized_grain_size;
              long double recrystalized_grain_volume = (4./3.)*numbers::PI*half_recrystalized_grain_size*half_recrystalized_grain_size*half_recrystalized_grain_size;

              return compute_derivatives_drexpp(cpo_index,
                                                data,
                                                mineral_i,
                                                strain_rate_3d,
                                                velocity_gradient_tensor,
                                                ref_resolved_shear_stress,
                                                recrystalized_grain_volume,
                                                aggregate_recrystalization_increment,
                                                volume_fractions,
                                                diffusion_pre_viscosities,
                                                diffusion_grain_size_exponent,
                                                dislocation_viscosities);
              break;
            }
            default:
              AssertThrow(false, ExcMessage("Internal error."));
              break;
          }
        AssertThrow(false, ExcMessage("Internal error."));
        return derivatives;
      }



      template <int dim>
      std::pair<std::vector<double>, std::vector<Tensor<2,3>>>
      CrystalPreferredOrientation<dim>::compute_derivatives_spin_tensor(const Tensor<2,3> &velocity_gradient_tensor) const
      {
        // dA/dt = W * A, where W is the spin tensor and A is the rotation matrix
        // The spin tensor is defined as W = 0.5 * ( L - L^T ), where L is the velocity gradient tensor.
        const Tensor<2,3> spin_tensor = -0.5 *(velocity_gradient_tensor - dealii::transpose(velocity_gradient_tensor));

        return std::pair<std::vector<double>, std::vector<Tensor<2,3>>>(std::vector<double>(n_grains,0.0), std::vector<Tensor<2,3>>(n_grains, spin_tensor));
      }


      template <int dim>
      std::pair<std::vector<double>, std::vector<Tensor<2,3>>>
      CrystalPreferredOrientation<dim>::compute_derivatives_drex_2004(const unsigned int cpo_index,
                                                                      const ArrayView<double> &data,
                                                                      const unsigned int mineral_i,
                                                                      const SymmetricTensor<2,3> &strain_rate_3d,
                                                                      const Tensor<2,3> &velocity_gradient_tensor,
                                                                      const std::array<double,4> ref_resolved_shear_stress,
                                                                      const bool prevent_nondimensionalization) const
      {
        // This if statement is only there for the unit test. In normal situations it should always be set to false,
        // because the nondimensionalization should always be done (in this exact way), unless you really know what
        // you are doing.
        double nondimensionalization_value = 1.0;
        if (!prevent_nondimensionalization)
          {
            const std::array< double, 3 > eigenvalues = dealii::eigenvalues(strain_rate_3d);
            nondimensionalization_value = std::max(std::abs(eigenvalues[0]),std::abs(eigenvalues[2]));

            Assert(!std::isnan(nondimensionalization_value), ExcMessage("The second invariant of the strain rate is not a number."));
          }


        // Make the strain-rate and velocity gradient tensor non-dimensional
        // by dividing it through the second invariant
        const Tensor<2,3> strain_rate_nondimensional = nondimensionalization_value != 0 ? strain_rate_3d/nondimensionalization_value : strain_rate_3d;
        const Tensor<2,3> velocity_gradient_tensor_nondimensional = nondimensionalization_value != 0 ? velocity_gradient_tensor/nondimensionalization_value : velocity_gradient_tensor;

        // create output variables
        std::vector<double> deriv_volume_fractions(n_grains);
        std::vector<Tensor<2,3>> deriv_a_cosine_matrices(n_grains);

        // create shortcuts
        const std::array<double, 4> &tau = ref_resolved_shear_stress;

        std::vector<double> strain_energy(n_grains);
        double mean_strain_energy = 0;

        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            // Compute the Schmidt tensor for this grain (nu), s is the slip system.
            // We first compute beta_s,nu (equation 5, Kaminski & Ribe, 2001)
            // Then we use the beta to calculate the Schmidt tensor G_{ij} (Eq. 5, Kaminski & Ribe, 2001)
            Tensor<2,3> G;
            Tensor<1,3> w;
            Tensor<1,4> beta({1.0, 1.0, 1.0, 1.0});
            std::array<Tensor<1,3>,4> slip_normal_reference {{Tensor<1,3>({0,1,0}),Tensor<1,3>({0,0,1}),Tensor<1,3>({0,1,0}),Tensor<1,3>({1,0,0})}};
            std::array<Tensor<1,3>,4> slip_direction_reference {{Tensor<1,3>({1,0,0}),Tensor<1,3>({1,0,0}),Tensor<1,3>({0,0,1}),Tensor<1,3>({0,0,1})}};

            // these are variables we only need for olivine, but we need them for both
            // within this if block and the next ones
            // Ordered vector where the first entry is the max/weakest and the last entry is the inactive slip system.
            std::array<unsigned int,4> indices {};

            // compute G and beta
            Tensor<1,4> bigI;
            const Tensor<2,3> rotation_matrix = get_rotation_matrix_grains(cpo_index,data,mineral_i,grain_i);
            const Tensor<2,3> rotation_matrix_transposed = transpose(rotation_matrix);
            for (unsigned int slip_system_i = 0; slip_system_i < 4; ++slip_system_i)
              {
                const Tensor<1,3> slip_normal_global = rotation_matrix_transposed*slip_normal_reference[slip_system_i];
                const Tensor<1,3> slip_direction_global = rotation_matrix_transposed*slip_direction_reference[slip_system_i];
                const Tensor<2,3> slip_cross_product = outer_product(slip_direction_global,slip_normal_global);
                bigI[slip_system_i] = scalar_product(slip_cross_product,strain_rate_nondimensional);
              }

            if (bigI.norm() < 1e-10)
              {
                // In this case there is no shear, only (possibly) a rotation. So \gamma_y and/or G should be zero.
                // Which is the default value, so do nothing.
              }
            else
              {
                // compute the element wise absolute value of the element wise
                // division of BigI by tau (tau = ref_resolved_shear_stress).
                std::array<double,4> q_abs;
                for (unsigned int i = 0; i < 4; ++i)
                  {
                    q_abs[i] = std::abs(bigI[i] / tau[i]);
                  }

                // here we find the indices starting at the largest value and ending at the smallest value
                // and assign them to special variables. Because all the variables are absolute values,
                // we can set them to a negative value to ignore them. This should be faster then deleting
                // the element, which would require allocation. (not tested)
                for (unsigned int slip_system_i = 0; slip_system_i < 4; ++slip_system_i)
                  {
                    indices[slip_system_i] = std::distance(q_abs.begin(),std::max_element(q_abs.begin(), q_abs.end()));
                    q_abs[indices[slip_system_i]] = -1;
                  }

                // compute the ordered beta vector, which is the relative slip rates of the active slip systems.
                // Test whether the max element is not equal to zero.
                Assert(bigI[indices[0]] != 0.0, ExcMessage("Internal error: bigI is zero."));
                beta[indices[0]] = 1.0; // max q_abs, weak system (most deformation) "s=1"

                const double ratio = tau[indices[0]]/bigI[indices[0]];
                for (unsigned int slip_system_i = 1; slip_system_i < 4-1; ++slip_system_i)
                  {
                    beta[indices[slip_system_i]] = std::pow(std::abs(ratio * (bigI[indices[slip_system_i]]/tau[indices[slip_system_i]])), stress_exponent);
                  }
                beta[indices.back()] = 0.0;
                for (unsigned int slip_system_i = 0; slip_system_i < 4; ++slip_system_i)

                  // Now compute the crystal rate of deformation tensor.
                  for (unsigned int i = 0; i < 3; ++i)
                    {
                      for (unsigned int j = 0; j < 3; j++)
                        {
                          G[i][j] = 2.0 * (beta[0] * rotation_matrix[0][i] * rotation_matrix[1][j]
                                           + beta[1] * rotation_matrix[0][i] * rotation_matrix[2][j]
                                           + beta[2] * rotation_matrix[2][i] * rotation_matrix[1][j]
                                           + beta[3] * rotation_matrix[2][i] * rotation_matrix[0][j]);
                        }
                    }
              }

            // Now calculate the analytic solution to the deformation minimization problem
            // compute gamma (equation 7, Kaminiski & Ribe, 2001)

            // Top is the numerator and bottom is the denominator in equation 7.
            double top = 0;
            double bottom = 0;
            for (unsigned int i = 0; i < 3; ++i)
              {
                // Following the actual Drex implementation we use i+2, which differs
                // from the EPSL paper, which says gamma_nu depends on i+1
                const unsigned int i_offset = (i==0) ? (i+2) : (i-1);

                top = top - (velocity_gradient_tensor_nondimensional[i][i_offset]-velocity_gradient_tensor_nondimensional[i_offset][i])*(G[i][i_offset]-G[i_offset][i]);
                bottom = bottom - (G[i][i_offset]-G[i_offset][i])*(G[i][i_offset]-G[i_offset][i]);

                for (unsigned int j = 0; j < 3; ++j)
                  {
                    top = top + 2.0 * G[i][j]*velocity_gradient_tensor_nondimensional[i][j];
                    bottom = bottom + 2.0* G[i][j] * G[i][j];
                  }
              }
            // see comment on if all BigI are zero. In that case gamma should be zero.
            const double gamma = (bottom != 0.0) ? top/bottom : 0.0;
            // compute w (equation 8, Kaminiski & Ribe, 2001)
            // w is the Rotation rate vector of the crystallographic axes of grain
            w[0] = 0.5*(velocity_gradient_tensor_nondimensional[2][1]-velocity_gradient_tensor_nondimensional[1][2]) - 0.5*(G[2][1]-G[1][2])*gamma;
            w[1] = 0.5*(velocity_gradient_tensor_nondimensional[0][2]-velocity_gradient_tensor_nondimensional[2][0]) - 0.5*(G[0][2]-G[2][0])*gamma;
            w[2] = 0.5*(velocity_gradient_tensor_nondimensional[1][0]-velocity_gradient_tensor_nondimensional[0][1]) - 0.5*(G[1][0]-G[0][1])*gamma;

            // Compute strain energy for this grain (abbreviated Estr)
            // For olivine: DREX only sums over 1-3. But Christopher Thissen's matlab
            // code (https://github.com/cthissen/Drex-MATLAB) corrected
            // this and writes each term using the indices created when calculating bigI.
            // Note tau = RRSS = (tau_m^s/tau_o), this why we get tau^(p-n)
            for (unsigned int slip_system_i = 0; slip_system_i < 4; ++slip_system_i)
              {
                const double rhos = std::pow(tau[indices[slip_system_i]],exponent_p-stress_exponent) *
                                    std::pow(std::abs(gamma*beta[indices[slip_system_i]]),exponent_p/stress_exponent);
                strain_energy[grain_i] += rhos * std::exp(-nucleation_efficiency * rhos * rhos);

                Assert(isfinite(strain_energy[grain_i]), ExcMessage("strain_energy[" + std::to_string(grain_i) + "] is not finite: " + std::to_string(strain_energy[grain_i])
                                                                    + ", rhos (" + std::to_string(slip_system_i) + ") = " + std::to_string(rhos)
                                                                    + ", nucleation_efficiency = " + std::to_string(nucleation_efficiency) + "."));
              }

            // compute the derivative of the rotation matrix: \frac{\partial a_{ij}}{\partial t}
            // (Eq. 9, Kaminski & Ribe 2001)
            const double volume_fraction_grain = get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i);
            if (volume_fraction_grain >= threshold_GBS/n_grains)
              {
                deriv_a_cosine_matrices[grain_i] = permutation_operator_3d * w  * nondimensionalization_value;
                // volume averaged strain energy
                mean_strain_energy += volume_fraction_grain * strain_energy[grain_i];

                Assert(isfinite(mean_strain_energy), ExcMessage("mean_strain_energy when adding grain " + std::to_string(grain_i) + " is not finite: " + std::to_string(mean_strain_energy)
                                                                + ", volume_fraction_grain = " + std::to_string(volume_fraction_grain) + "."));
              }
            else
              {
                strain_energy[grain_i] = 0;
              }
          }

        // Change of volume fraction of grains by grain boundary migration
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            // Different than D-Rex. Here we actually only compute the derivative and do not multiply it with the volume_fractions. We do that when we advect.
            deriv_volume_fractions[grain_i] = get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i) * mobility * (mean_strain_energy - strain_energy[grain_i]) * nondimensionalization_value;
            Assert(isfinite(deriv_volume_fractions[grain_i]),
                   ExcMessage("deriv_volume_fractions[" + std::to_string(grain_i) + "] is not finite: "
                              + std::to_string(deriv_volume_fractions[grain_i])));
          }

        return std::pair<std::vector<double>, std::vector<Tensor<2,3>>>(deriv_volume_fractions, deriv_a_cosine_matrices);
      }



      template <int dim>
      void
      CrystalPreferredOrientation<dim>::recrystalize_grains(const unsigned int cpo_index,
                                                            const ArrayView<double> &data,
                                                            const unsigned int mineral_i,
                                                            const double recrystalized_grain_volume,
                                                            const std::vector<double> &recrystalization_fractions,
                                                            std::vector<double> &strain_energy) const
      {
        // make a vector for which the first entry contains the index of the smallest vector, the second entry contains the
        // index of the second smallest vector, etc.
        std::vector<std::size_t> permutation_vector(n_grains);
        std::iota(permutation_vector.begin(), permutation_vector.end(), 0);
        std::sort(permutation_vector.begin(), permutation_vector.end(),
                  [&](std::size_t grain_i, std::size_t grain_j)
        {
          return get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i) < get_volume_fractions_grains(cpo_index,data,mineral_i,grain_j);
        });

        unsigned int permutation_vector_counter = 0;
        Tensor<2,3> main_rotation_matrix = get_rotation_matrix_grains(cpo_index,data,mineral_i,permutation_vector[permutation_vector_counter]);

        for (size_t i = 0; i < 3; i++)
          for (size_t j = 0; j < 3; j++)
            Assert(abs(main_rotation_matrix[i][j]) <= 1.0,
                   ExcMessage("rotation_matrix[" + std::to_string(i) + "][" + std::to_string(j) +
                              "] is larger than one: " + std::to_string(main_rotation_matrix[i][j]) + " (" + std::to_string(main_rotation_matrix[i][j]-1.0) + "). rotation_matrix = \n"
                              + std::to_string(main_rotation_matrix[0][0]) + " " + std::to_string(main_rotation_matrix[0][1]) + " " + std::to_string(main_rotation_matrix[0][2]) + "\n"
                              + std::to_string(main_rotation_matrix[1][0]) + " " + std::to_string(main_rotation_matrix[1][1]) + " " + std::to_string(main_rotation_matrix[1][2]) + "\n"
                              + std::to_string(main_rotation_matrix[2][0]) + " " + std::to_string(main_rotation_matrix[2][1]) + " " + std::to_string(main_rotation_matrix[2][2])));

        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            const double grain_volume = get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i);
            double temp_total_volume = 0.;
            for (unsigned int i = 0; i < permutation_vector.size(); ++i)
              {
                temp_total_volume += get_volume_fractions_grains(cpo_index,data,mineral_i,i);
              }

            size_t n_recrystalized_grains = std::floor((recrystalization_fractions[grain_i]*grain_volume)/recrystalized_grain_volume);

            if (n_recrystalized_grains > 0)
              {
                const double grain_volume_left = grain_volume-n_recrystalized_grains*recrystalized_grain_volume;
                set_volume_fractions_grains(cpo_index,data,mineral_i,grain_i,grain_volume_left);

                // compute the volume of n_recrystalized_grains+1 smallest grains
                double replaced_grain_volume = 0.0;
                for (unsigned int recrystalize_grain_i = 0; recrystalize_grain_i < n_recrystalized_grains+1; ++recrystalize_grain_i)
                  {
                    replaced_grain_volume += get_volume_fractions_grains(cpo_index,data,mineral_i,permutation_vector[permutation_vector_counter+recrystalize_grain_i]);
                  }

                set_volume_fractions_grains(cpo_index,data,mineral_i,permutation_vector[permutation_vector_counter],replaced_grain_volume);

                Tensor<2,3> random_rotation_matrix;
                this->compute_random_rotation_matrix(random_rotation_matrix);
                for (size_t i = 0; i < 3; i++)
                  for (size_t j = 0; j < 3; j++)
                    Assert(abs(random_rotation_matrix[i][j]) <= 1.0,
                           ExcMessage("rotation_matrix[" + std::to_string(i) + "][" + std::to_string(j) +
                                      "] is larger than one: " + std::to_string(random_rotation_matrix[i][j]) + " (" + std::to_string(random_rotation_matrix[i][j]-1.0) + "). rotation_matrix = \n"
                                      + std::to_string(random_rotation_matrix[0][0]) + " " + std::to_string(random_rotation_matrix[0][1]) + " " + std::to_string(random_rotation_matrix[0][2]) + "\n"
                                      + std::to_string(random_rotation_matrix[1][0]) + " " + std::to_string(random_rotation_matrix[1][1]) + " " + std::to_string(random_rotation_matrix[1][2]) + "\n"
                                      + std::to_string(random_rotation_matrix[2][0]) + " " + std::to_string(random_rotation_matrix[2][1]) + " " + std::to_string(random_rotation_matrix[2][2])));
                set_rotation_matrix_grains(cpo_index,data,mineral_i,permutation_vector[permutation_vector_counter],random_rotation_matrix);
                strain_energy[permutation_vector[permutation_vector_counter]] = 0.0;

                permutation_vector_counter++;

                Tensor<2,3> random_deflection;
                for (unsigned int recrystalize_grain_i = 0; recrystalize_grain_i < n_recrystalized_grains; ++recrystalize_grain_i)
                  {
                    set_volume_fractions_grains(cpo_index,data,mineral_i,permutation_vector[permutation_vector_counter],recrystalized_grain_volume);

                    // TODO: compute random rotation matrix between two degrees and apply it
                    // Have to test this to ensure that the misorientation between the two tensors is only 10 degrees

                    this->compute_random_rotation_matrix(random_deflection);
                    main_rotation_matrix = main_rotation_matrix*random_deflection;
                    for (size_t i = 0; i < 3; i++)
                      for (size_t j = 0; j < 3; j++)
                        Assert(abs(main_rotation_matrix[i][j]) <= 1.0,
                               ExcMessage("rotation_matrix[" + std::to_string(i) + "][" + std::to_string(j) +
                                          "] is larger than one: " + std::to_string(main_rotation_matrix[i][j]) + " (" + std::to_string(main_rotation_matrix[i][j]-1.0) + "). rotation_matrix = \n"
                                          + std::to_string(main_rotation_matrix[0][0]) + " " + std::to_string(main_rotation_matrix[0][1]) + " " + std::to_string(main_rotation_matrix[0][2]) + "\n"
                                          + std::to_string(main_rotation_matrix[1][0]) + " " + std::to_string(main_rotation_matrix[1][1]) + " " + std::to_string(main_rotation_matrix[1][2]) + "\n"
                                          + std::to_string(main_rotation_matrix[2][0]) + " " + std::to_string(main_rotation_matrix[2][1]) + " " + std::to_string(main_rotation_matrix[2][2])));

                    set_rotation_matrix_grains(cpo_index,data,mineral_i,permutation_vector[permutation_vector_counter],main_rotation_matrix*random_deflection);

                    strain_energy[permutation_vector[permutation_vector_counter]] = 0.0;
                    permutation_vector_counter++;
                  }
              }

            temp_total_volume = 0.;
            for (unsigned int i = 0; i < permutation_vector.size(); ++i)
              {
                temp_total_volume += get_volume_fractions_grains(cpo_index,data,mineral_i,i);
              }
          }
      }


      template <int dim>
      std::pair<std::vector<double>, std::vector<Tensor<2,3>>>
      CrystalPreferredOrientation<dim>::compute_derivatives_drexpp(const unsigned int cpo_index,
                                                                   const ArrayView<double> &data,
                                                                   const unsigned int mineral_i,
                                                                   const SymmetricTensor<2,3> &strain_rate,
                                                                   const Tensor<2,3> &velocity_gradient_tensor,
                                                                   const std::array<double,4> ref_resolved_shear_stress,
                                                                   const double recrystalized_grain_volume,
                                                                   const double aggregate_recrystalization_increment,
                                                                   const std::vector<double> &volume_fractions,
                                                                   const std::vector<double> &diffusion_pre_viscosities,
                                                                   const std::vector<double> &diffusion_grain_size_exponent,
                                                                   const std::vector<double> &dislocation_viscosities) const
      {
        // create output variables
        std::vector<double> deriv_volume_fractions(n_grains);
        std::vector<Tensor<2,3>> deriv_a_cosine_matrices(n_grains);

        // create shorcuts
        const std::array<double, 4> &tau = ref_resolved_shear_stress;

        std::vector<double> strain_energy(n_grains);
        std::vector<double> recrystalized_fractions(n_grains);
        std::vector<double> subgrain_rotation_fractions(n_grains);
        std::vector<double> grain_boundary_sliding_fractions(n_grains);
        std::vector<Tensor<1,3>> spin_vectors(n_grains);
        double total_subgrain_rotation_fraction = 0.0;
        double total_grain_boundary_sliding_fraction = 0.0;

        // first compute the strain energy and G for all grains
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            // Compute the Schmidt tensor for this grain (nu), s is the slip system.
            // We first compute beta_s,nu (equation 5, Kaminski & Ribe, 2001)
            // Then we use the beta to calculate the Schmidt tensor G_{ij} (Eq. 5, Kaminski & Ribe, 2001)
            //Tensor<2,3> G;
            Tensor<1,4> beta({1.0, 1.0, 1.0, 1.0});
            std::array<Tensor<1,3>,4> slip_normal_reference {{Tensor<1,3>({0,1,0}),Tensor<1,3>({0,0,1}),Tensor<1,3>({0,1,0}),Tensor<1,3>({1,0,0})}};
            std::array<Tensor<1,3>,4> slip_direction_reference {{Tensor<1,3>({1,0,0}),Tensor<1,3>({1,0,0}),Tensor<1,3>({0,0,1}),Tensor<1,3>({0,0,1})}};

            // these are variables we only need for olivine, but we need them for both
            // within this if block and the next ones
            // Ordered vector where the first entry is the max/weakest and the last entry is the inactive slip system.
            std::array<unsigned int,4> indices {};

            // compute G and beta
            Tensor<1,4> bigI;
            const Tensor<2,3> rotation_matrix = get_rotation_matrix_grains(cpo_index,data,mineral_i,grain_i);
            for (size_t i = 0; i < 3; i++)
              for (size_t j = 0; j < 3; j++)
                Assert(abs(rotation_matrix[i][j]) <= 1.0,
                       ExcMessage("rotation_matrix[" + std::to_string(i) + "][" + std::to_string(j) +
                                  "] is larger than one: " + std::to_string(rotation_matrix[i][j]) + " (" + std::to_string(rotation_matrix[i][j]-1.0) + "). rotation_matrix = \n"
                                  + std::to_string(rotation_matrix[0][0]) + " " + std::to_string(rotation_matrix[0][1]) + " " + std::to_string(rotation_matrix[0][2]) + "\n"
                                  + std::to_string(rotation_matrix[1][0]) + " " + std::to_string(rotation_matrix[1][1]) + " " + std::to_string(rotation_matrix[1][2]) + "\n"
                                  + std::to_string(rotation_matrix[2][0]) + " " + std::to_string(rotation_matrix[2][1]) + " " + std::to_string(rotation_matrix[2][2])));
            const Tensor<2,3> rotation_matrix_transposed = transpose(rotation_matrix);
            for (unsigned int slip_system_i = 0; slip_system_i < 4; ++slip_system_i)
              {
                const Tensor<1,3> slip_normal_global = rotation_matrix_transposed*slip_normal_reference[slip_system_i];
                const Tensor<1,3> slip_direction_global = rotation_matrix_transposed*slip_direction_reference[slip_system_i];
                const Tensor<2,3> slip_cross_product = outer_product(slip_direction_global,slip_normal_global);
                bigI[slip_system_i] = scalar_product(slip_cross_product,strain_rate);
              }

            Tensor<2,3> schmidt_tensor;

            // The value in this if statement is arbitrary. Has to be further checked out to for a more "realistic value".
            if (bigI.norm() < 1e-15)
              {
                // In this case there is no shear, only (possibly) a rotation. So \gamma_y and/or G should be zero.
                // Which is the default value, so do nothing.
              }
            else
              {
                // compute the element wise absolute value of the element wise
                // division of BigI by tau (tau = ref_resolved_shear_stress).
                std::array<double,4> q_abs;
                for (unsigned int i = 0; i < 4; ++i)
                  {
                    q_abs[i] = std::abs(bigI[i] / tau[i]);
                  }

                // here we find the indices starting at the largest value and ending at the smallest value
                // and assign them to special variables. Because all the variables are absolute values,
                // we can set them to a negative value to ignore them. This should be faster then deleting
                // the element, which would require allocation. (not tested)
                for (unsigned int slip_system_i = 0; slip_system_i < 4; ++slip_system_i)
                  {
                    indices[slip_system_i] = std::distance(q_abs.begin(),std::max_element(q_abs.begin(), q_abs.end()));
                    q_abs[indices[slip_system_i]] = -1;
                  }

                // compute the ordered beta vector, which is the relative slip rates of the active slip systems.
                // Test whether the max element is not equal to zero.
                Assert(bigI[indices[0]] != 0.0, ExcMessage("Internal error: bigI is zero."));
                beta[indices[0]] = 1.0; // max q_abs, weak system (most deformation) "s=1"

                const double ratio = tau[indices[0]]/bigI[indices[0]];
                for (unsigned int slip_system_i = 1; slip_system_i < 4-1; ++slip_system_i)
                  {
                    beta[indices[slip_system_i]] = std::pow(std::abs(ratio * (bigI[indices[slip_system_i]]/tau[indices[slip_system_i]])), stress_exponent);
                  }
                beta[indices.back()] = 0.0;


                // Now compute the crystal rate of deformation tensor.

                for (unsigned int i = 0; i < 3; ++i)
                  {
                    for (unsigned int j = 0; j < 3; j++)
                      {
                        schmidt_tensor[i][j] = 2.0 * (beta[0] * rotation_matrix[0][i] * rotation_matrix[1][j]
                                                      + beta[1] * rotation_matrix[0][i] * rotation_matrix[2][j]
                                                      + beta[2] * rotation_matrix[2][i] * rotation_matrix[1][j]
                                                      + beta[3] * rotation_matrix[2][i] * rotation_matrix[0][j]);

                      }

                  }
              }


            // Now calculate the analytic solution to the deformation minimization problem
            // compute gamma (equation 7, Kaminiski & Ribe, 2001)

            // Top is the numerator and bottom is the denominator in equation 7.
            double top = 0;
            double bottom = 0;
            for (unsigned int i = 0; i < 3; ++i)
              {
                // Following the actual Drex implementation we use i+2, which differs
                // from the EPSL paper, which says gamma_nu depends on i+1
                const unsigned int i_offset = (i==0) ? (i+2) : (i-1);

                top = top - (velocity_gradient_tensor[i][i_offset]-velocity_gradient_tensor[i_offset][i])*(schmidt_tensor[i][i_offset]-schmidt_tensor[i_offset][i]);
                bottom = bottom - (schmidt_tensor[i][i_offset]-schmidt_tensor[i_offset][i])*(schmidt_tensor[i][i_offset]-schmidt_tensor[i_offset][i]);

                for (unsigned int j = 0; j < 3; ++j)
                  {
                    top = top + 2.0 * schmidt_tensor[i][j]*velocity_gradient_tensor[i][j];
                    bottom = bottom + 2.0* schmidt_tensor[i][j] * schmidt_tensor[i][j];
                  }
              }
            // see comment on if all BigI are zero. In that case gamma should be zero.
            const double gamma = (bottom != 0.0) ? top/bottom : 0.0;

            // compute w (equation 8, Kaminiski & Ribe, 2001)
            // w is the Rotation rate vector of the crystallographic axes of grain
            spin_vectors[grain_i] = Tensor<1,3>
                                    (
            {
              0.5*(velocity_gradient_tensor[2][1]-velocity_gradient_tensor[1][2]) - 0.5*(schmidt_tensor[2][1]-schmidt_tensor[1][2]) *gamma,
              0.5*(velocity_gradient_tensor[0][2]-velocity_gradient_tensor[2][0]) - 0.5*(schmidt_tensor[0][2]-schmidt_tensor[2][0]) *gamma,
              0.5*(velocity_gradient_tensor[1][0]-velocity_gradient_tensor[0][1]) - 0.5*(schmidt_tensor[1][0]-schmidt_tensor[0][1]) *gamma
            });

            // Compute strain energy for this grain (abbreviated Estr)
            // For olivine: DREX only sums over 1-3. But Christopher Thissen's matlab
            // code (https://github.com/cthissen/Drex-MATLAB) corrected
            // this and writes each term using the indices created when calculating bigI.
            // Note tau = RRSS = (tau_m^s/tau_o), this why we get tau^(p-n)
            long double alpha = 0.0;
            subgrain_rotation_fractions[grain_i] = 4.0;
            for (unsigned int slip_system_i = 0; slip_system_i < 4; ++slip_system_i)
              {
                const double rhos = std::pow(tau[indices[slip_system_i]],exponent_p-stress_exponent) *
                                    std::pow(std::abs(gamma*beta[indices[slip_system_i]]),exponent_p/stress_exponent);
                strain_energy[grain_i] += rhos;
                alpha += std::exp(-nucleation_efficiency * rhos * rhos);
                Assert(isfinite(strain_energy[grain_i]), ExcMessage("strain_energy[" + std::to_string(grain_i) + "] is not finite: " + std::to_string(strain_energy[grain_i])
                                                                    + ", rhos (" + std::to_string(slip_system_i) + ") = " + std::to_string(rhos)
                                                                    + ", nucleation_efficiency = " + std::to_string(nucleation_efficiency) + "."));

                Assert(isfinite(alpha), ExcMessage("alpha is not finite: " + std::to_string(alpha)
                                                   + ", rhos (" + std::to_string(slip_system_i) + ") = " + std::to_string(rhos)
                                                   + ", nucleation_efficiency = " + std::to_string(nucleation_efficiency) + "."));
              }

            subgrain_rotation_fractions[grain_i] = 4.0-alpha;
            total_subgrain_rotation_fraction += subgrain_rotation_fractions[grain_i];
          }
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            // compute the diffusion viscosity
            // get grain_size
            const double grain_volume = get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i);
            const double grain_size = 2.0*std::cbrt(grain_volume * 3.0/(4.0*numbers::PI));

            // compute the viscosities
            std::vector<double> diffusion_viscosities(volume_fractions.size(),std::numeric_limits<double>::quiet_NaN());
            std::vector<double> composite_viscosities(volume_fractions.size(),std::numeric_limits<double>::quiet_NaN());
            for (unsigned int composition = 0; composition < volume_fractions.size(); ++composition)
              {
                diffusion_viscosities[composition] = diffusion_pre_viscosities[composition] * std::pow(grain_size, diffusion_grain_size_exponent[composition]);
                composite_viscosities[composition] = diffusion_viscosities[composition]+dislocation_viscosities[composition];
              }
            const double diffusion_viscosity = MaterialModel::MaterialUtilities::average_value(volume_fractions, diffusion_viscosities, MaterialModel::MaterialUtilities::harmonic);
            const double composite_viscosity = MaterialModel::MaterialUtilities::average_value(volume_fractions, composite_viscosities, MaterialModel::MaterialUtilities::harmonic);

            grain_boundary_sliding_fractions[grain_i] = 1 - (std::log(diffusion_viscosity)/std::log(composite_viscosity));
            AssertThrow(grain_boundary_sliding_fractions[grain_i]>0.0, ExcMessage("diffusion strain-rate larger than total strain-rate. "
                                                                                  "composite_viscosity = " + std::to_string(composite_viscosity)
                                                                                  + ", diffusion_viscosity = " + std::to_string(diffusion_viscosity)
                                                                                  + ", grain_boundary_sliding_fractions[grain_i] = " + std::to_string(grain_boundary_sliding_fractions[grain_i])
                                                                                 ));
            total_grain_boundary_sliding_fraction += grain_boundary_sliding_fractions[grain_i];
          }

        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            recrystalized_fractions[grain_i] = total_grain_boundary_sliding_fraction != 0 && total_subgrain_rotation_fraction !=0 ?
                                               (grain_boundary_sliding_fractions[grain_i]/total_grain_boundary_sliding_fraction)
                                               * (subgrain_rotation_fractions[grain_i]/total_subgrain_rotation_fraction)
                                               * aggregate_recrystalization_increment
                                               :
                                               0.0;
          }

        this->recrystalize_grains(cpo_index,
                                  data,
                                  mineral_i,
                                  recrystalized_grain_volume,
                                  recrystalized_fractions,
                                  strain_energy);

        double mean_strain_energy = 0.0;
        // Change of volume fraction of grains by grain boundary migration
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            // compute the derivative of the rotation matrix: \frac{\partial a_{ij}}{\partial t}
            // (Eq. 9, Kaminski & Ribe 2001)
            deriv_a_cosine_matrices[grain_i] = 0;
            const double volume_fraction_grain = get_volume_fractions_grains(cpo_index,data,mineral_i,grain_i);
            if (volume_fraction_grain >= (threshold_GBS*1e-12) && strain_energy[grain_i] > 0.0) // TODO: Change to actual grain size based on the rheology
              {
                deriv_a_cosine_matrices[grain_i] = permutation_operator_3d * spin_vectors[grain_i];

                // volume averaged strain energy
                mean_strain_energy += volume_fraction_grain * strain_energy[grain_i];

                Assert(isfinite(mean_strain_energy), ExcMessage("mean_strain_energy when adding grain " + std::to_string(grain_i) + " is not finite: " + std::to_string(mean_strain_energy)
                                                                + ", volume_fraction_grain = " + std::to_string(volume_fraction_grain) + "."));
              }
            else
              {
                strain_energy[grain_i] = 0;
              }

            // Different than D-Rex. Here we actually only compute the derivative and do not multiply it with the volume_fractions. We do that when we advect.
            deriv_volume_fractions[grain_i] = get_volume_fraction_mineral(cpo_index,data,mineral_i) * mobility * (mean_strain_energy - strain_energy[grain_i]);

            Assert(isfinite(deriv_volume_fractions[grain_i]),
                   ExcMessage("deriv_volume_fractions[" + std::to_string(grain_i) + "] is not finite: "
                              + std::to_string(deriv_volume_fractions[grain_i])));
          }

        return std::pair<std::vector<double>, std::vector<Tensor<2,3>>>(deriv_volume_fractions, deriv_a_cosine_matrices);
      }


      template<int dim>
      DeformationType
      CrystalPreferredOrientation<dim>::determine_deformation_type(const DeformationTypeSelector deformation_type_selector,
                                                                   const Point<dim> &position,
                                                                   const double temperature,
                                                                   const double pressure,
                                                                   const Tensor<1,dim> &velocity,
                                                                   const std::vector<double> &compositions,
                                                                   const SymmetricTensor<2,dim> &strain_rate,
                                                                   const SymmetricTensor<2,dim> &deviatoric_strain_rate,
                                                                   const double water_content) const
      {
        // Now compute what type of deformation takes place.
        switch (deformation_type_selector)
          {
            case DeformationTypeSelector::passive:
              return DeformationType::passive;
            case DeformationTypeSelector::olivine_a_fabric:
              return DeformationType::olivine_a_fabric;
            case DeformationTypeSelector::olivine_b_fabric:
              return DeformationType::olivine_b_fabric;
            case DeformationTypeSelector::olivine_c_fabric:
              return DeformationType::olivine_c_fabric;
            case DeformationTypeSelector::olivine_d_fabric:
              return DeformationType::olivine_d_fabric;
            case DeformationTypeSelector::olivine_e_fabric:
              return DeformationType::olivine_e_fabric;
            case DeformationTypeSelector::enstatite:
              return DeformationType::enstatite;
            case DeformationTypeSelector::olivine_karato_2008:
              // construct the material model inputs and outputs
              // Since this function is only evaluating one particle,
              // we use 1 for the amount of quadrature points.
              MaterialModel::MaterialModelInputs<dim> material_model_inputs(1,this->n_compositional_fields());
              material_model_inputs.position[0] = position;
              material_model_inputs.temperature[0] = temperature;
              material_model_inputs.pressure[0] = pressure;
              material_model_inputs.velocity[0] = velocity;
              material_model_inputs.composition[0] = compositions;
              material_model_inputs.strain_rate[0] = strain_rate;

              MaterialModel::MaterialModelOutputs<dim> material_model_outputs(1,this->n_compositional_fields());
              this->get_material_model().evaluate(material_model_inputs, material_model_outputs);
              double eta = material_model_outputs.viscosities[0];

              const SymmetricTensor<2,dim> stress = 2*eta*deviatoric_strain_rate +
                                                    pressure * unit_symmetric_tensor<dim>();
              const std::array< double, dim > eigenvalues = dealii::eigenvalues(stress);
              double differential_stress = eigenvalues[0]-eigenvalues[dim-1];
              return determine_deformation_type_karato_2008(differential_stress, water_content);

          }

        AssertThrow(false, ExcMessage("Internal error. Deformation type not implemented."));
        return DeformationType::passive;
      }


      template<int dim>
      DeformationType
      CrystalPreferredOrientation<dim>::determine_deformation_type_karato_2008(const double stress, const double water_content) const
      {
        constexpr double MPa = 1e6;
        constexpr double ec_line_slope = -500./1050.;
        if (stress > (380. - 0.05 * water_content)*MPa)
          {
            if (stress > (625. - 2.5 * water_content)*MPa)
              {
                return DeformationType::olivine_b_fabric;
              }
            else
              {
                return DeformationType::olivine_d_fabric;
              }
          }
        else
          {
            if (stress < (625.0 -2.5 * water_content)*MPa)
              {
                return DeformationType::olivine_a_fabric;
              }
            else
              {
                if (stress < (500.0 + ec_line_slope*-100. + ec_line_slope * water_content)*MPa)
                  {
                    return DeformationType::olivine_e_fabric;
                  }
                else
                  {
                    return DeformationType::olivine_c_fabric;
                  }
              }
          }
      }


      template<int dim>
      std::array<double,4>
      CrystalPreferredOrientation<dim>::reference_resolved_shear_stress_from_deformation_type(DeformationType deformation_type,
          double max_value) const
      {
        std::array<double,4> ref_resolved_shear_stress;
        switch (deformation_type)
          {
            // from Kaminski and Ribe, GJI 2004 and
            // Becker et al., 2007 (http://www-udc.ig.utexas.edu/external/becker/preprints/bke07.pdf)
            case DeformationType::olivine_a_fabric :
              ref_resolved_shear_stress[0] = 1;
              ref_resolved_shear_stress[1] = 2;
              ref_resolved_shear_stress[2] = 3;
              ref_resolved_shear_stress[3] = max_value;
              break;

            // from Kaminski and Ribe, GJI 2004 and
            // Becker et al., 2007 (http://www-udc.ig.utexas.edu/external/becker/preprints/bke07.pdf)
            case DeformationType::olivine_b_fabric :
              ref_resolved_shear_stress[0] = 3;
              ref_resolved_shear_stress[1] = 2;
              ref_resolved_shear_stress[2] = 1;
              ref_resolved_shear_stress[3] = max_value;
              break;

            // from Kaminski and Ribe, GJI 2004 and
            // Becker et al., 2007 (http://www-udc.ig.utexas.edu/external/becker/preprints/bke07.pdf)
            case DeformationType::olivine_c_fabric :
              ref_resolved_shear_stress[0] = 3;
              ref_resolved_shear_stress[1] = max_value;
              ref_resolved_shear_stress[2] = 2;
              ref_resolved_shear_stress[3] = 1;
              break;

            // from Kaminski and Ribe, GRL 2002 and
            // Becker et al., 2007 (http://www-udc.ig.utexas.edu/external/becker/preprints/bke07.pdf)
            case DeformationType::olivine_d_fabric :
              ref_resolved_shear_stress[0] = 1;
              ref_resolved_shear_stress[1] = 1;
              ref_resolved_shear_stress[2] = 3;
              ref_resolved_shear_stress[3] = max_value;
              break;

            // Kaminski, Ribe and Browaeys, GJI, 2004 (same as in the matlab code) and
            // Becker et al., 2007 (http://www-udc.ig.utexas.edu/external/becker/preprints/bke07.pdf)
            case DeformationType::olivine_e_fabric :
              ref_resolved_shear_stress[0] = 2;
              ref_resolved_shear_stress[1] = 1;
              ref_resolved_shear_stress[2] = max_value;
              ref_resolved_shear_stress[3] = 3;
              break;

            // from Kaminski and Ribe, GJI 2004.
            // Todo: this one is not used in practice, since there is an optimization in
            // the code. So maybe remove it in the future.
            case DeformationType::enstatite :
              ref_resolved_shear_stress[0] = max_value;
              ref_resolved_shear_stress[1] = max_value;
              ref_resolved_shear_stress[2] = max_value;
              ref_resolved_shear_stress[3] = 1;
              break;

            default:
              AssertThrow(false,
                          ExcMessage("Deformation type enum with number " + std::to_string(static_cast<unsigned int>(deformation_type))
                                     + " was not found."));
              break;
          }
        return ref_resolved_shear_stress;
      }

      template<int dim>
      unsigned int
      CrystalPreferredOrientation<dim>::get_number_of_grains() const
      {
        return n_grains;
      }



      template<int dim>
      unsigned int
      CrystalPreferredOrientation<dim>::get_number_of_minerals() const
      {
        return n_minerals;
      }



      template <int dim>
      void
      CrystalPreferredOrientation<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Crystal Preferred Orientation");
            {
              prm.declare_entry ("Random number seed", "1",
                                 Patterns::Integer (0),
                                 "The seed used to generate random numbers. This will make sure that "
                                 "results are reproducible as long as the problem is run with the "
                                 "same number of MPI processes. It is implemented as final seed = "
                                 "user seed + MPI Rank. ");

              prm.declare_entry ("Number of grains per particle", "50",
                                 Patterns::Integer (1),
                                 "The number of grains of each different mineral "
                                 "each particle contains.");

              prm.declare_entry ("Property advection method", "Backward Euler",
                                 Patterns::Anything(),
                                 "Options: Forward Euler, Backward Euler");

              prm.declare_entry ("Property advection tolerance", "1e-10",
                                 Patterns::Double(0),
                                 "The Backward Euler property advection method involve internal iterations. "
                                 "This option allows for setting a tolerance. When the norm of tensor new - tensor old is "
                                 "smaller than this tolerance, the iteration is stopped.");

              prm.declare_entry ("Property advection max iterations", "100",
                                 Patterns::Integer(0),
                                 "The Backward Euler property advection method involve internal iterations. "
                                 "This option allows for setting the maximum number of iterations. Note that when the iteration "
                                 "is ended by the max iteration amount an assert is thrown.");

              prm.declare_entry ("CPO derivatives algorithm", "Spin tensor",
                                 Patterns::List(Patterns::Anything()),
                                 "Options: Spin tensor, D-Rex 2004, D-Rex++");

              prm.enter_subsection("Initial grains");
              {
                prm.declare_entry("Model name","Uniform grains and random uniform rotations",
                                  Patterns::Anything(),
                                  "The model used to initialize the CPO for all particles. "
                                  "Currently 'Uniform grains and random uniform rotations' is the only valid option.");

                prm.declare_entry ("Minerals", "Olivine: Karato 2008, Enstatite",
                                   Patterns::List(Patterns::Anything()),
                                   "This determines what minerals and fabrics or fabric selectors are used used for the LPO/CPO calculation. "
                                   "The options are Olivine: Passive, A-fabric, Olivine: B-fabric, Olivine: C-fabric, Olivine: D-fabric, "
                                   "Olivine: E-fabric, Olivine: Karato 2008 or Enstatite. Passive sets all RRSS entries to the maximum. The "
                                   "Karato 2008 selector selects a fabric based on stress and water content as defined in "
                                   "figure 4 of the Karato 2008 review paper (doi: 10.1146/annurev.earth.36.031207.124120).");

                prm.declare_entry ("Volume fractions minerals", "0.7, 0.3",
                                   Patterns::List(Patterns::Double(0)),
                                   "The volume fractions for the different minerals. "
                                   "There need to be the same number of values as there are minerals."
                                   "Note that the currently implemented scheme is incompressible and "
                                   "does not allow chemical interaction or the formation of new phases");
              }
              prm.leave_subsection ();

              prm.enter_subsection("D-Rex 2004");
              {
                prm.declare_entry ("Mobility", "125",
                                   Patterns::Double(0),
                                   "The dimensionless intrinsic grain boundary mobility for both olivine and enstatite.");

                prm.declare_entry ("Volume fractions minerals", "0.5, 0.5",
                                   Patterns::List(Patterns::Double(0)),
                                   "The volume fraction for the different minerals. "
                                   "There need to be the same amount of values as there are minerals");

                prm.declare_entry ("Stress exponents", "3.5",
                                   Patterns::Double(0),
                                   "This is the power law exponent that characterizes the rheology of the "
                                   "slip systems. It is used in equation 11 of Kaminski et al., 2004.");

                prm.declare_entry ("Exponents p", "1.5",
                                   Patterns::Double(0),
                                   "This is exponent p as defined in equation 11 of Kaminski et al., 2004. ");

                prm.declare_entry ("Nucleation efficiency", "5",
                                   Patterns::Double(0),
                                   "This is the dimensionless nucleation rate as defined in equation 8 of "
                                   "Kaminski et al., 2004. ");

                prm.declare_entry ("Threshold GBS", "0.3",
                                   Patterns::Double(0),
                                   "The Dimensionless Grain Boundary Sliding (GBS) threshold. "
                                   "This is a grain size threshold below which grain deform by GBS and "
                                   "become strain-free grains.");
              }
              prm.leave_subsection();

              prm.enter_subsection("D-Rex++");
              {
                prm.declare_entry ("Mobility", "125",
                                   Patterns::Double(0),
                                   "The dimensionless intrinsic grain boundary mobility for both olivine and enstatite.");

                prm.declare_entry ("Volume fractions minerals", "0.5, 0.5",
                                   Patterns::List(Patterns::Double(0)),
                                   "The volume fraction for the different minerals. "
                                   "There need to be the same amount of values as there are minerals");

                prm.declare_entry ("Stress exponents", "3.5",
                                   Patterns::Double(0),
                                   "This is the power law exponent that characterizes the rheology of the "
                                   "slip systems. It is used in equation 11 of Kaminski et al., 2004.");

                prm.declare_entry ("Exponents p", "1.5",
                                   Patterns::Double(0),
                                   "This is exponent p as defined in equation 11 of Kaminski et al., 2004. ");

                prm.declare_entry ("Nucleation efficiency", "5",
                                   Patterns::Double(0),
                                   "This is the dimensionless nucleation rate as defined in equation 8 of "
                                   "Kaminski et al., 2004. ");

                prm.declare_entry ("Threshold GBS", "0.3",
                                   Patterns::Double(0),
                                   "The Dimensionless Grain Boundary Sliding (GBS) threshold. "
                                   "This is a grain size threshold below which grain deform by GBS and "
                                   "become strain-free grains.");
              }
              prm.leave_subsection();
            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }



      template <int dim>
      void
      CrystalPreferredOrientation<dim>::parse_parameters (ParameterHandler &prm)
      {
        AssertThrow(dim == 3, ExcMessage("CPO computations are currently only supported for 3d models. "
                                         "2d computations will work when this assert is removed, but you will need to make sure that the "
                                         "correct 3d strain-rate and velocity gradient tensors are provided to the algorithm."));

        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Crystal Preferred Orientation");
            {
              random_number_seed = prm.get_integer ("Random number seed");
              n_grains = prm.get_integer("Number of grains per particle");

              property_advection_tolerance = prm.get_double("Property advection tolerance");
              property_advection_max_iterations = prm.get_integer ("Property advection max iterations");

              const std::string temp_cpo_derivative_algorithm = prm.get("CPO derivatives algorithm");

              if (temp_cpo_derivative_algorithm == "Spin tensor")
                {
                  cpo_derivative_algorithm = CPODerivativeAlgorithm::spin_tensor;
                }
              else if (temp_cpo_derivative_algorithm ==  "D-Rex 2004")
                {
                  cpo_derivative_algorithm = CPODerivativeAlgorithm::drex_2004;
                }
              else if (temp_cpo_derivative_algorithm ==  "D-Rex++")
                {
                  cpo_derivative_algorithm = CPODerivativeAlgorithm::drexpp;
                }
              else
                {
                  AssertThrow(false,
                              ExcMessage("The CPO derivatives algorithm needs to be one of the following: "
                                         "Spin tensor, D-Rex 2004."))
                }

              const std::string temp_advection_method = prm.get("Property advection method");
              if (temp_advection_method == "Forward Euler")
                {
                  advection_method = AdvectionMethod::forward_euler;
                }
              else if (temp_advection_method == "Backward Euler")
                {
                  advection_method = AdvectionMethod::backward_euler;
                }
              else
                {
                  AssertThrow(false, ExcMessage("particle property advection method not found: \"" + temp_advection_method + "\""));
                }

              prm.enter_subsection("Initial grains");
              {
                const std::string model_name = prm.get("Model name");
                AssertThrow(model_name == "Uniform grains and random uniform rotations",
                            ExcMessage("No model named " + model_name + "for CPO particle property initialization. "
                                       + "Only the model \"Uniform grains and random uniform rotations\" is available."));
                const std::vector<std::string> temp_deformation_type_selector = dealii::Utilities::split_string_list(prm.get("Minerals"));
                //prm.enter_subsection("Uniform grains and random uniform rotations");
                //{
                n_minerals = temp_deformation_type_selector.size();
                deformation_type_selector.resize(n_minerals);

                for (size_t mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
                  {
                    if (temp_deformation_type_selector[mineral_i] == "Passive")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::passive;
                      }
                    else if (temp_deformation_type_selector[mineral_i] == "Olivine: Karato 2008")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::olivine_karato_2008;
                      }
                    else if (temp_deformation_type_selector[mineral_i] ==  "Olivine: A-fabric")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::olivine_a_fabric;
                      }
                    else if (temp_deformation_type_selector[mineral_i] ==  "Olivine: B-fabric")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::olivine_b_fabric;
                      }
                    else if (temp_deformation_type_selector[mineral_i] ==  "Olivine: C-fabric")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::olivine_c_fabric;
                      }
                    else if (temp_deformation_type_selector[mineral_i] ==  "Olivine: D-fabric")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::olivine_d_fabric;
                      }
                    else if (temp_deformation_type_selector[mineral_i] ==  "Olivine: E-fabric")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::olivine_e_fabric;
                      }
                    else if (temp_deformation_type_selector[mineral_i] ==  "Enstatite")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::enstatite;
                      }
                    else
                      {
                        AssertThrow(false,
                                    ExcMessage("The fabric needs to be assigned one of the following comma-delimited values: Olivine: Karato 2008, "
                                               "Olivine: A-fabric, Olivine: B-fabric, Olivine: C-fabric, Olivine: D-fabric,"
                                               "Olivine: E-fabric, Enstatite, Passive."))
                      }
                  }

                volume_fractions_minerals = Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Volume fractions minerals")));
                double volume_fractions_minerals_sum = 0;
                for (auto fraction : volume_fractions_minerals)
                  {
                    volume_fractions_minerals_sum += fraction;
                  }

                AssertThrow(abs(volume_fractions_minerals_sum-1.0) < 2.0 * std::numeric_limits<double>::epsilon(),
                            ExcMessage("The sum of the CPO volume fractions should be one."));
                //}
                //prm.leave_subsection();
              }
              prm.leave_subsection();

              prm.enter_subsection("D-Rex 2004");
              {
                mobility = prm.get_double("Mobility");
                volume_fractions_minerals = Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Volume fractions minerals")));
                stress_exponent = prm.get_double("Stress exponents");
                exponent_p = prm.get_double("Exponents p");
                nucleation_efficiency = prm.get_double("Nucleation efficiency");
                threshold_GBS = prm.get_double("Threshold GBS");
              }
              prm.leave_subsection();

              prm.enter_subsection("D-Rex++");
              {
                mobility = prm.get_double("Mobility");
                volume_fractions_minerals = Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Volume fractions minerals")));
                stress_exponent = prm.get_double("Stress exponents");
                exponent_p = prm.get_double("Exponents p");
                nucleation_efficiency = prm.get_double("Nucleation efficiency");
                threshold_GBS = prm.get_double("Threshold GBS");
              }
              prm.leave_subsection();
            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();


        prm.enter_subsection("Material model");
        {
          prm.enter_subsection ("Visco Plastic");
          {
            // Phase transition parameters
            phase_function.initialize_simulator (this->get_simulator());
            phase_function.parse_parameters (prm);

            // Retrieve the list of composition names
            const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

            // Establish that a background field is required here
            const bool has_background_field = true;

            thermal_diffusivities = Utilities::parse_map_to_double_array (prm.get("Thermal diffusivities"),
                                                                          list_of_composition_names,
                                                                          has_background_field,
                                                                          "Thermal diffusivities");

            define_conductivities = prm.get_bool ("Define thermal conductivities");

            thermal_conductivities = Utilities::parse_map_to_double_array (prm.get("Thermal conductivities"),
                                                                           list_of_composition_names,
                                                                           has_background_field,
                                                                           "Thermal conductivities");

            rheology_diff = std::make_unique<MaterialModel::Rheology::DiffusionCreep<dim>>();
            rheology_diff->initialize_simulator (this->get_simulator());
            rheology_diff->parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(phase_function.n_phases_for_each_composition()));

            rheology_disl = std::make_unique<MaterialModel::Rheology::DislocationCreep<dim>>();
            rheology_disl->initialize_simulator (this->get_simulator());
            rheology_disl->parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(phase_function.n_phases_for_each_composition()));

            rheology_vipl = std::make_unique<MaterialModel::Rheology::ViscoPlastic<dim>>();
            rheology_vipl->initialize_simulator (this->get_simulator());
            rheology_vipl->parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(phase_function.n_phases_for_each_composition()));
            min_strain_rate = rheology_vipl->min_strain_rate;
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(CrystalPreferredOrientation,
                                        "crystal preferred orientation",
                                        "The plugin manages and computes the evolution of Lattice/Crystal Preferred Orientations (LPO/CPO) "
                                        "on particles. Each ASPECT particle can be assigned many grains. Each grain is assigned a size and a orientation "
                                        "matrix. This allows for CPO evolution tracking with polycrystalline kinematic CrystalPreferredOrientation evolution models such "
                                        "as D-Rex (Kaminski and Ribe, 2001; Kaminski et al., 2004).")
    }
  }
}
