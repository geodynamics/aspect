/*
  Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/interpolator/quaternion_average.h>
#include <aspect/particle/property/crystal_preferred_orientation.h>
#include <aspect/particle/manager.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      template <int dim>
      std::vector<std::vector<double>>
      QuaternionAverage<dim>::properties_at_points(const ParticleHandler<dim> &particle_handler,
                                                   const std::vector<Point<dim>> &positions,
                                                   const ComponentMask &selected_properties,
                                                   const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        const typename ParticleHandler<dim>::particle_iterator_range particle_range =
          particle_handler.particles_in_cell(cell);

        const unsigned int n_particles = std::distance(particle_range.begin(),particle_range.end());
        const unsigned int n_particle_properties = particle_handler.n_properties_per_particle();
        AssertThrow(selected_properties.size() == not_quaternion_properties.size(), ExcMessage("component masks have unequal length") );
        ComponentMask selected_properties_not_quaternions = selected_properties & not_quaternion_properties;

        // average quaternion for each mineral
        std::vector<std::array<double,4>> mean_quat (n_minerals);

        // unit weights for average
        // can be adapted to distance weighting in the future
        const std::vector<double> unit_weights(n_particles, 1.0);

        // average mineral agregates separately
        for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
          {
            // collect all quaternions in the cell
            std::vector<std::array<double,4>> quat_array;
            quat_array.reserve(n_particles);
            for (const auto &particle : particle_range)
              {
                const ArrayView<const double> &particle_properties = particle.get_properties();
                std::array<double,4> quaternion;
                for (unsigned int j = 0; j < 4; ++j)
                  quaternion[j] = particle_properties[quaternion_data_pos[mineral_i] + j];
                quat_array.push_back(quaternion);
              }
            if (symmetry_group == "triclinic")
              mean_quat[mineral_i] = markley_average(quat_array, unit_weights);
            else
              mean_quat[mineral_i] = symmetry_average(quat_array, unit_weights);
          }


        std::vector<std::vector<double>> cell_properties(positions.size(),
                                                          std::vector<double>(n_particle_properties,
                                                                              numbers::signaling_nan<double>()));
        cell_properties = base_interpolator->properties_at_points(particle_handler,
                                                                  positions,
                                                                  selected_properties_not_quaternions,
                                                                  cell);

        for (unsigned int mineral_i = 0; mineral_i < n_minerals; mineral_i++)
          for (unsigned int pos_index=0; pos_index<positions.size(); pos_index++)
            for (unsigned int j = 0; j<4; j++)
              cell_properties[pos_index][quaternion_data_pos[mineral_i]+j] = mean_quat[mineral_i][j];

        return cell_properties;
      }



      template <int dim>
      void
      QuaternionAverage<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Interpolator");
        {
          prm.enter_subsection("Quaternion average");
          {
            prm.declare_entry("Base interpolation scheme", "cell average",
                              Patterns::Selection(Interpolator::get_valid_interpolator_names_pattern<dim>()),
                              "Scheme that is used to interpolate everything but quaternions");
            prm.declare_entry("Symmetry group", "triclinic",
                              Patterns::Selection("triclinic|orthorhombic"),
                              "symmetry group of underlying rotations that are supposed to be averaged"
                              "this determines which symmetry operators are used to project rotations"
                              "into the fundamental zone of the respective symmetry group");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

      }



      template <int dim>
      void
      QuaternionAverage<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Interpolator");
        {
          prm.enter_subsection("Quaternion average");
          {
            symmetry_group = prm.get("Symmetry group");
            AssertThrow( prm.get("Base interpolation scheme") != "quaternion average",
                         ExcMessage("You may not use the ``quaternion average'' for averaging other objects than quaternions") );

            // create the base model;
            // it will get a chance to read its parameters below
            // and initialize its SimulatorAccess base class after we leave the current section
            base_interpolator = create_particle_interpolator<dim>(prm.get("Base interpolation scheme"));
            if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_interpolator.get()))
              sim->initialize_simulator (this->get_simulator());

            base_interpolator -> set_particle_manager_index(this->get_particle_manager_index());
            base_interpolator -> initialize();

            // get quaternion data position and number of minerals
            const auto &particle_property_manager = this->get_particle_manager(this->get_particle_manager_index()).get_property_manager();
            const auto &particle_property_information = particle_property_manager.get_data_info();

            std::vector<std::string> active_plugin_names = particle_property_information.get_property_names();
            // remove this in case you want to average another rotation based particle property
            AssertThrow(particle_property_manager.plugin_name_exists("crystal preferred orientation"),
                        ExcMessage("No crystal preferred orientation property plugin found."));
            AssertThrow(particle_property_manager.plugin_name_exists("cpo bingham average"),
                        ExcMessage("No cpo bingham average property plugin found."));

            const auto &cpo_particle_property = particle_property_manager.template get_matching_active_plugin<Particle::Property::CrystalPreferredOrientation<dim>>();
            n_minerals = cpo_particle_property.get_number_of_minerals();
            quaternion_data_pos.resize(n_minerals);

            const unsigned int n_property_components = particle_property_information.n_components();

            std::vector<bool> components_not_quaternions (n_property_components, true);

            for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
              {
                AssertThrow(particle_property_manager.get_data_info().fieldname_exists("cpo mineral " + std::to_string(mineral_i) + " q.n"), ExcMessage("quaternions must be a particle property.\noutput rotation as in subsection Bingham average has to be set to quaternion."));

                quaternion_data_pos[mineral_i] = particle_property_manager.get_data_info().get_position_by_field_name("cpo mineral " + std::to_string(mineral_i) + " q.w");

                for (unsigned int j=0; j<4; j++)
                  components_not_quaternions[quaternion_data_pos[mineral_i]+j] = false;
              }
            not_quaternion_properties = ComponentMask(components_not_quaternions);
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        base_interpolator -> parse_parameters(prm);
      }



      template <int dim>
      std::array<double,4>
      QuaternionAverage<dim>::markley_average(const std::vector<std::array<double,4>> &quat_array, const std::vector<double> &weights) const
      {
        AssertThrow(quat_array.size()==weights.size(), ExcMessage("unequal length of quaternion array and weight array"));

        // use LAPACK to compute eigenvalues/vectors of a 4x4 real symmetric matrix
        LAPACKFullMatrix<double> K(4,4);
        double normalization = 0;

        for (unsigned int n=0; n < quat_array.size(); n++)
          {
            for (unsigned int i1=0; i1<4; i1++)
              {
                for (unsigned int i2 =0; i2<i1+1; i2++)
                  {
                    K(i1,i2) += quat_array[n][i1]*quat_array[n][i2]*weights[n];
                    if ((n == quat_array.size()-1) && (i1 !=i2)) K(i2,i1) = K(i1,i2);
                  }
              }
            normalization += weights[n];
          }

        // ensures matrix is traceless i.e. all eigenvalues add up to 0
        K *= 4.0;
        for (unsigned int j=0; j<4; j++) K(j,j) -= normalization;

        Vector<double> eigenvalues;
        FullMatrix<double> eigenvectors;

        // only pick eigenvalues larger than 0 -> reduce search space
        // maximum eiganvalue has to be > 0 as eigenvalues add up to one
        K.compute_eigenvalues_symmetric(0,std::numeric_limits<double>::infinity(), 0, eigenvalues, eigenvectors);

        unsigned int last_idx = eigenvalues.size()-1;

        // stay in positive quaternion half space
        if (eigenvectors(0, last_idx) < 0)
          return {{-eigenvectors(0,last_idx),-eigenvectors(1,last_idx),-eigenvectors(2,last_idx),-eigenvectors(3,last_idx)}};
        else
          return {{ eigenvectors(0,last_idx), eigenvectors(1,last_idx), eigenvectors(2,last_idx), eigenvectors(3,last_idx)}};
      }



      template <int dim>
      double
      QuaternionAverage<dim>::SO3_geodesic(const std::array<double,4> &quaternion1, const std::array<double,4> &quaternion2) const
      {
        double inner_product = 0;
        for (unsigned int j=0; j<4; j++) inner_product += quaternion1[j]*quaternion2[j];

        return 2*std::acos(std::min(std::abs(inner_product),1.0));
      }

      template <int dim>
      double
      QuaternionAverage<dim>::objective_function(const std::array<double,4> &mean_quat, const std::vector<std::array<double,4>> &quat_array, const std::vector<double> &weights) const
      {
        double objective_fct = 0.;
        double normalization = 0.;
        for (unsigned int j=0; j<quat_array.size(); j++)
          {
            objective_fct += weights[j]*Utilities::fixed_power<2>(SO3_geodesic(quat_array[j],mean_quat));
            normalization += weights[j];
          }

        return objective_fct/normalization;
      }

      template <int dim>
      std::array<double,4>
      QuaternionAverage<dim>::to_fundamental_zone(const std::array<double,4> &quaternion, const std::array<double,4> &reference_quaternion) const
      {
        std::vector<std::array<double,4>> possible_mutations;

        // q_i = q*g_i are the equivalent representations of one rotation
        // where g_i are the elements of the selected symmetry group G_cr
        if (symmetry_group == "orthorhombic")
          {
            // for orthorhombic symmetry gi = {+/-1, +/-i, +/-j, +/-k }
            // note that +/- are the same rotation and have the same geodesic
            possible_mutations.push_back(quaternion);
            possible_mutations.push_back({{-quaternion[1], quaternion[0], quaternion[3],-quaternion[2]}});
            possible_mutations.push_back({{-quaternion[2],-quaternion[3], quaternion[0], quaternion[1]}});
            possible_mutations.push_back({{-quaternion[3], quaternion[2],-quaternion[1], quaternion[0]}});
          }
        else
          AssertThrow(false, ExcMessage(symmetry_group + " is not a valid symmetry group or not yet implemented"));

        const unsigned int n_symmetry_elements = possible_mutations.size();

        // find minimal geodesic
        std::vector<double> geodesics (n_symmetry_elements);

        for (unsigned int j=0; j < n_symmetry_elements; j++)
          geodesics[j] = SO3_geodesic(possible_mutations[j], reference_quaternion);

        unsigned int best_idx = 0;
        for (unsigned int i = 1; i < n_symmetry_elements; ++i)
          if (geodesics[i] < geodesics[best_idx])
            best_idx = i;

        // corresponding to one optimal permutation
        std::array<double, 4> closest_quaternion = possible_mutations[best_idx];

        // remain in halfspace where q[0] > 0
        if (closest_quaternion[0] < 0)
          closest_quaternion = {{
              -closest_quaternion[0], -closest_quaternion[1],
              -closest_quaternion[2], -closest_quaternion[3]
            }
          };

        return closest_quaternion;
      }



      template <int dim>
      std::array<double,4>
      QuaternionAverage<dim>::symmetry_average(const std::vector<std::array<double,4>> &quat_array, const std::vector<double> &weights) const
      {
        const unsigned int n_particles = quat_array.size();
        // crafting an initial guess
        // avoids locking into a not optimal solution for few particles
        std::array<double,4> quat_mean = markley_average(quat_array, weights);

        // find maximal geodesic from initial mean
        // this corresponds to the rotation closest to the edge of the fundamental zone
        // for an antipodal distribution this is likely to be close to the actual mean
        std::vector<double> geodesics(n_particles);
        for (unsigned int j=0; j<n_particles; j++)
          geodesics[j] = SO3_geodesic(quat_array[j], quat_mean);

        unsigned int best_idx = 0;
        for (unsigned int i = 1; i < n_particles; ++i)
          if (geodesics[i] > geodesics[best_idx])
            best_idx = i;

        std::array<double,4> quat_0 = quat_array[best_idx];

        // initializing loop variables to be optimized
        double objective = objective_function(quat_mean, quat_array, weights);
        double new_objective = objective;
        std::vector<std::array<double,4>> symmetrized_quat_array(n_particles);

        // combinatoric optimization problem
        unsigned int max_iterations = 50;
        for (unsigned int n=0; n<max_iterations; n++)
          {
            // send all quaternions into fundamental zone of guess/mean
            for (unsigned int j=0; j<n_particles; j++)
              symmetrized_quat_array[j] = to_fundamental_zone(quat_array[j], quat_0);

            // standart average in SO(3)
            quat_mean = markley_average(symmetrized_quat_array, weights);

            // compute objective function
            new_objective = objective_function(quat_mean, symmetrized_quat_array, weights);

            // iterate as long as objective function decreases
            if (objective > new_objective)
              {
                objective = new_objective;
                quat_0 = quat_mean;
              }
            else return quat_mean;
            AssertThrow(n!=max_iterations-1, ExcMessage("maximum number of iterations reached in quaternion average aborting"));
          }
        // we should never reach this point, but just in case so nothing breaks
        return {{1.,0.,0.,0.}};
      }

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(QuaternionAverage,
                                            "quaternion average",
                                            "Return a quaternion average with taking into account the orthotropic symmetry group ")
    }
  }
}
