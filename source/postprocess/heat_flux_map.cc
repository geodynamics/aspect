/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


#include <aspect/postprocess/heat_flux_map.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/heating_model/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/boundary_heat_flux/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      /**
       * This function computes the consistent heat flux through each boundary face.
       * For reflecting boundaries the integrated heat flux is 0, for boundaries with prescribed heat flux
       * (Neumann boundary conditions) it is simply the integral of the prescribed heat flux over
       * the face, but for boundaries with prescribed temperature (Dirichlet boundary conditions) it
       * is computed using the consistent boundary flux method. The method is described in
       *
       * Gresho, P. M., Lee, R. L., Sani, R. L., Maslanik, M. K., & Eaton, B. E. (1987).
       * The consistent Galerkin FEM for computing derived boundary quantities in thermal and or fluids
       * problems. International Journal for Numerical Methods in Fluids, 7(4), 371-394.
       *
       * In summary, the method solves the temperature equation again on the boundary faces, with known
       * temperatures and solving for the boundary fluxes that satisfy the equation. Since the
       * equation is only formed on the faces and it can be solved using only diagonal matrices it is cheap,
       * and conceptually simpler methods like evaluating the temperature gradient on the face are
       * significantly less accurate.
       */
      template <int dim>
      std::vector<std::vector<std::pair<double, double> > >
      compute_heat_flux_through_boundary_faces (const SimulatorAccess<dim> &simulator_access)
      {
        std::vector<std::vector<std::pair<double, double> > > heat_flux_and_area(simulator_access.get_triangulation().n_active_cells(),
                                                                                 std::vector<std::pair<double, double> >(GeometryInfo<dim>::faces_per_cell,
                                                                                     std::pair<double,double>(0.0,0.0)));

        // Quadrature degree for assembling the consistent boundary flux equation, see Simulator::assemble_advection_system()
        // for a justification of the chosen quadrature degree.
        const unsigned int quadrature_degree = simulator_access.get_parameters().temperature_degree
                                               +
                                               (simulator_access.get_parameters().stokes_velocity_degree+1)/2;

        // Gauss quadrature in the interior for best accuracy.
        const QGauss<dim> quadrature_formula(quadrature_degree);
        // GLL quadrature on the surface to get a diagonal mass matrix.
        const QGaussLobatto<dim-1> quadrature_formula_face(quadrature_degree);

        // The CBF method involves both boundary and volume integrals on the
        // cells at the boundary. Construct FEValues objects for each of these integrations.
        FEValues<dim> fe_volume_values (simulator_access.get_mapping(),
                                        simulator_access.get_fe(),
                                        quadrature_formula,
                                        update_values |
                                        update_gradients |
                                        update_q_points |
                                        update_JxW_values);

        FEFaceValues<dim> fe_face_values (simulator_access.get_mapping(),
                                          simulator_access.get_fe(),
                                          quadrature_formula_face,
                                          update_JxW_values |
                                          update_values |
                                          update_gradients |
                                          update_normal_vectors |
                                          update_q_points);

        const unsigned int dofs_per_cell = simulator_access.get_fe().dofs_per_cell;
        const unsigned int n_q_points = quadrature_formula.size();
        const unsigned int n_face_q_points = quadrature_formula_face.size();

        // Vectors for solving CBF system.
        Vector<double> local_vector(dofs_per_cell);
        Vector<double> local_mass_matrix(dofs_per_cell);

        // The mass matrix may be stored in a vector as it is a diagonal matrix.
        LinearAlgebra::BlockVector mass_matrix(simulator_access.introspection().index_sets.system_partitioning,
                                               simulator_access.get_mpi_communicator());
        LinearAlgebra::BlockVector distributed_heat_flux_vector(simulator_access.introspection().index_sets.system_partitioning,
                                                                simulator_access.get_mpi_communicator());
        LinearAlgebra::BlockVector heat_flux_vector(simulator_access.introspection().index_sets.system_partitioning,
                                                    simulator_access.introspection().index_sets.system_relevant_partitioning,
                                                    simulator_access.get_mpi_communicator());
        LinearAlgebra::BlockVector rhs_vector(simulator_access.introspection().index_sets.system_partitioning,
                                              simulator_access.get_mpi_communicator());

        distributed_heat_flux_vector = 0.;
        heat_flux_vector = 0.;

        typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_volume_values.n_quadrature_points, simulator_access.n_compositional_fields());
        typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_volume_values.n_quadrature_points, simulator_access.n_compositional_fields());
        typename HeatingModel::HeatingModelOutputs heating_out(fe_volume_values.n_quadrature_points, simulator_access.n_compositional_fields());

        typename MaterialModel::Interface<dim>::MaterialModelInputs face_in(fe_face_values.n_quadrature_points, simulator_access.n_compositional_fields());
        typename MaterialModel::Interface<dim>::MaterialModelOutputs face_out(fe_face_values.n_quadrature_points, simulator_access.n_compositional_fields());

        std::vector<double> old_temperatures (n_q_points);
        std::vector<double> old_old_temperatures (n_q_points);
        std::vector<Tensor<1,dim> > temperature_gradients (n_q_points);

        const double time_step = simulator_access.get_timestep();
        const double old_time_step = simulator_access.get_old_timestep();

        const std::set<types::boundary_id> &fixed_temperature_boundaries =
          simulator_access.get_boundary_temperature_manager().get_fixed_temperature_boundary_indicators();

        const std::set<types::boundary_id> &fixed_heat_flux_boundaries =
          simulator_access.get_parameters().fixed_heat_flux_boundary_indicators;

        Vector<float> artificial_viscosity(simulator_access.get_triangulation().n_active_cells());
        simulator_access.get_artificial_viscosity(artificial_viscosity, true);

        // loop over all of the surface cells and evaluate the heat flux
        typename DoFHandler<dim>::active_cell_iterator
        cell = simulator_access.get_dof_handler().begin_active(),
        endc = simulator_access.get_dof_handler().end();

        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned() && cell->at_boundary())
            {
              // First check if we need to compute heat flux for this cell at all,
              // that is if it has a face on a boundary with prescribed temperatures
              // or prescribed heat fluxes
              bool compute_heat_flux = false;
              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                {
                  if (!cell->at_boundary(f))
                    continue;

                  const unsigned int boundary_id = cell->face(f)->boundary_id();
                  if (fixed_temperature_boundaries.find(boundary_id) != fixed_temperature_boundaries.end())
                    {
                      compute_heat_flux = true;
                      break;
                    }
                  if (fixed_heat_flux_boundaries.find(boundary_id) != fixed_heat_flux_boundaries.end())
                    {
                      compute_heat_flux = true;
                      break;
                    }
                }

              if (compute_heat_flux == false)
                continue;

              fe_volume_values.reinit (cell);
              in.reinit(fe_volume_values, cell, simulator_access.introspection(), simulator_access.get_solution(), true);
              simulator_access.get_material_model().evaluate(in, out);

              if (simulator_access.get_parameters().formulation_temperature_equation ==
                  Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
                {
                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      out.densities[q] = simulator_access.get_adiabatic_conditions().density(in.position[q]);
                    }
                }

              MaterialModel::MaterialAveraging::average (simulator_access.get_parameters().material_averaging,
                                                         cell,
                                                         fe_volume_values.get_quadrature(),
                                                         fe_volume_values.get_mapping(),
                                                         out);

              simulator_access.get_heating_model_manager().evaluate(in, out, heating_out);

              local_vector = 0.;
              local_mass_matrix = 0.;

              fe_volume_values[simulator_access.introspection().extractors.temperature].get_function_gradients (simulator_access.get_solution(), temperature_gradients);
              fe_volume_values[simulator_access.introspection().extractors.temperature].get_function_values (simulator_access.get_old_solution(), old_temperatures);
              fe_volume_values[simulator_access.introspection().extractors.temperature].get_function_values (simulator_access.get_old_old_solution(), old_old_temperatures);

              // Compute volume integrals on RHS of the CBF system
              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  double temperature_time_derivative;

                  if (simulator_access.get_timestep_number() > 1)
                    temperature_time_derivative = (1.0/time_step) *
                                                  (in.temperature[q] *
                                                   (2*time_step + old_time_step) / (time_step + old_time_step)
                                                   -
                                                   old_temperatures[q] *
                                                   (1 + time_step/old_time_step)
                                                   +
                                                   old_old_temperatures[q] *
                                                   (time_step * time_step) / (old_time_step * (time_step + old_time_step)));
                  else if (simulator_access.get_timestep_number() == 1)
                    temperature_time_derivative =
                      (in.temperature[q] - old_temperatures[q]) / time_step;
                  else
                    temperature_time_derivative = 0.0;

                  const double JxW = fe_volume_values.JxW(q);

                  const double density_c_P = out.densities[q] * out.specific_heat[q];
                  const double latent_heat_LHS = heating_out.lhs_latent_heat_terms[q];
                  const double material_prefactor = density_c_P + latent_heat_LHS;

                  const double artificial_viscosity_cell = static_cast<double>(artificial_viscosity(cell->active_cell_index()));
                  const double diffusion_constant = std::max(out.thermal_conductivities[q],
                                                             artificial_viscosity_cell);

                  for (unsigned int i = 0; i<dofs_per_cell; ++i)
                    {
                      local_vector(i) +=
                        // conduction term
                        (-diffusion_constant *
                         (fe_volume_values[simulator_access.introspection().extractors.temperature].gradient(i,q)
                          * temperature_gradients[q])
                         +
                         // advection term and time derivative
                         (- material_prefactor * (temperature_gradients[q] * in.velocity[q] + temperature_time_derivative)
                          // source terms
                          + heating_out.heating_source_terms[q])
                         * fe_volume_values[simulator_access.introspection().extractors.temperature].value(i,q))
                        * JxW;
                    }
                }

              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                {
                  if (!cell->at_boundary(f))
                    continue;

                  const unsigned int boundary_id = cell->face(f)->boundary_id();

                  if (fixed_temperature_boundaries.find(boundary_id) != fixed_temperature_boundaries.end())
                    {
                      fe_face_values.reinit (cell, f);

                      // Assemble the mass matrix for cell face. Since we are using GLL
                      // quadrature, the mass matrix will be diagonal, and we can just assemble it into a vector.
                      for (unsigned int q=0; q<n_face_q_points; ++q)
                        for (unsigned int i=0; i<dofs_per_cell; ++i)
                          local_mass_matrix(i) += fe_face_values[simulator_access.introspection().extractors.temperature].value(i,q) *
                                                  fe_face_values[simulator_access.introspection().extractors.temperature].value(i,q) *
                                                  fe_face_values.JxW(q);
                    }
                  else if (fixed_heat_flux_boundaries.find(boundary_id) != fixed_heat_flux_boundaries.end())
                    {
                      fe_face_values.reinit (cell, f);
                      face_in.reinit(fe_face_values, cell, simulator_access.introspection(), simulator_access.get_solution(), true);
                      simulator_access.get_material_model().evaluate(face_in, face_out);

                      if (simulator_access.get_parameters().formulation_temperature_equation ==
                          Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
                        {
                          for (unsigned int q=0; q<n_q_points; ++q)
                            {
                              face_out.densities[q] = simulator_access.get_adiabatic_conditions().density(face_in.position[q]);
                            }
                        }

                      std::vector<Tensor<1,dim> > heat_flux(n_face_q_points);
                      heat_flux = simulator_access.get_boundary_heat_flux().heat_flux(
                                    boundary_id,
                                    face_in,
                                    face_out,
#if DEAL_II_VERSION_GTE(9,0,0)
                                    fe_face_values.get_normal_vectors()
#else
                                    fe_face_values.get_all_normal_vectors()
#endif
                                  );

                      // For Neumann boundaries: Integrate the heat flux for the face, but still assemble the
                      // boundary terms for the Dirichlet boundary computations on this cell
                      for (unsigned int q=0; q < n_face_q_points; ++q)
                        {
                          const double normal_heat_flux = heat_flux[q] * fe_face_values.normal_vector(q);
                          const double JxW = fe_face_values.JxW(q);
                          heat_flux_and_area[cell->active_cell_index()][f].first += normal_heat_flux * JxW;
                          heat_flux_and_area[cell->active_cell_index()][f].second += JxW;

                          for (unsigned int i = 0; i<dofs_per_cell; ++i)
                            {
                              // heat flux boundary condition terms
                              local_vector(i) += - fe_face_values[simulator_access.introspection().extractors.temperature].value(i,q) *
                                                 normal_heat_flux * JxW;
                            }
                        }
                    }
                  else
                    continue;
                }

              cell->distribute_local_to_global(local_mass_matrix, mass_matrix);
              cell->distribute_local_to_global(local_vector, rhs_vector);
            }

        mass_matrix.compress(VectorOperation::add);
        rhs_vector.compress(VectorOperation::add);

        const IndexSet local_elements = mass_matrix.locally_owned_elements();
        for (unsigned int k=0; k<local_elements.n_elements(); ++k)
          {
            const unsigned int global_index = local_elements.nth_index_in_set(k);

            // Since the mass matrix is diagonal, we can just solve for the heat flux vector by dividing the
            // right-hand side by the mass matrix entry
            if (mass_matrix[global_index] > 1.e-15)
              distributed_heat_flux_vector[global_index] = rhs_vector[global_index] / mass_matrix[global_index];
          }

        distributed_heat_flux_vector.compress(VectorOperation::insert);
        heat_flux_vector = distributed_heat_flux_vector;

        std::vector<double> heat_flux_values(n_face_q_points);

        // Now integrate the heat flux for each face
        for (cell = simulator_access.get_dof_handler().begin_active(); cell!=endc; ++cell)
          if (cell->is_locally_owned() && cell->at_boundary())
            {
              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                if (cell->at_boundary(f))
                  {
                    const unsigned int boundary_id = cell->face(f)->boundary_id();

                    // Only fill output vector for non-Neumann boundaries (Neumann boundaries have been filled above)
                    if (fixed_heat_flux_boundaries.find(boundary_id) == fixed_heat_flux_boundaries.end())
                      {
                        fe_face_values.reinit (cell, f);

                        if (fixed_temperature_boundaries.find(boundary_id) != fixed_temperature_boundaries.end())
                          fe_face_values[simulator_access.introspection().extractors.temperature].get_function_values(heat_flux_vector, heat_flux_values);
                        else
                          for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                            heat_flux_values[q] = 0.0;

                        // Integrate the consistent heat flux and face area
                        for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                          {
                            heat_flux_and_area[cell->active_cell_index()][f].first += heat_flux_values[q] *
                                                                                      fe_face_values.JxW(q);
                            heat_flux_and_area[cell->active_cell_index()][f].second += fe_face_values.JxW(q);
                          }
                      }
                  }
            }
        return heat_flux_and_area;
      }
    }

    template <int dim>
    std::pair<std::string,std::string>
    HeatFluxMap<dim>::execute (TableHandler &)
    {
      std::vector<std::vector<std::pair<double, double> > > heat_flux_and_area =
        internal::compute_heat_flux_through_boundary_faces (*this);

      // have a stream into which we write the data. the text stream is then
      // later sent to processor 0
      std::ostringstream output;
      std::vector<std::pair<Point<dim>,double> > stored_values;

      // loop over all of the surface cells and evaluate the heat flux
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->at_boundary(f) &&
                (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top" ||
                 this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "bottom"))
              {
                // evaluate position of heat flow to write into output file
                const Point<dim> midpoint_at_surface = cell->face(f)->center();

                const double flux_density = heat_flux_and_area[cell->active_cell_index()][f].first /
                                            heat_flux_and_area[cell->active_cell_index()][f].second;

                // store final position and heat flow
                stored_values.push_back (std::make_pair(midpoint_at_surface, flux_density));
              }


      // Write the solution to an output stream
      for (unsigned int i=0; i<stored_values.size(); ++i)
        {
          output << stored_values[i].first
                 << ' '
                 << stored_values[i].second
                 << std::endl;
        }



      std::string filename = this->get_output_directory() +
                             "heat_flux." +
                             Utilities::int_to_string(this->get_timestep_number(), 5);
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        filename.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

      const unsigned int max_data_length = Utilities::MPI::max (output.str().size()+1,
                                                                this->get_mpi_communicator());
      const unsigned int mpi_tag = 567;

      // on processor 0, collect all of the data the individual processors send
      // and concatenate them into one file
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::ofstream file (filename.c_str());

          file << "# "
               << ((dim==2)? "x y" : "x y z")
               << " heat flux" << std::endl;

          // first write out the data we have created locally
          file << output.str();

          std::string tmp;
          tmp.resize (max_data_length, '\0');

          // then loop through all of the other processors and collect
          // data, then write it to the file
          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              MPI_Status status;
              // get the data. note that MPI says that an MPI_Recv may receive
              // less data than the length specified here. since we have already
              // determined the maximal message length, we use this feature here
              // rather than trying to find out the exact message length with
              // a call to MPI_Probe.
              MPI_Recv (&tmp[0], max_data_length, MPI_CHAR, p, mpi_tag,
                        this->get_mpi_communicator(), &status);

              // output the string. note that 'tmp' has length max_data_length,
              // but we only wrote a certain piece of it in the MPI_Recv, ended
              // by a \0 character. write only this part by outputting it as a
              // C string object, rather than as a std::string
              file << tmp.c_str();
            }
        }
      else
        // on other processors, send the data to processor zero. include the \0
        // character at the end of the string
        {
          MPI_Send (&output.str()[0], output.str().size()+1, MPI_CHAR, 0, mpi_tag,
                    this->get_mpi_communicator());
        }

      return std::pair<std::string,std::string>("Writing heat flux map:",
                                                filename);
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(HeatFluxMap,
                                  "heat flux map",
                                  "A postprocessor that computes the heat flux "
                                  "density across each boundary. The heat flux density is computed in "
                                  "outward direction, i.e., from the domain to the outside, "
                                  "using the consistent boundary flux method as described in "
                                  "Gresho, P. M., Lee, R. L., Sani, R. L., Maslanik, M. K., & Eaton, B. E. (1987). "
                                  "The consistent Galerkin FEM for computing derived boundary quantities in thermal and or fluids "
                                  "problems. International Journal for Numerical Methods in Fluids, 7(4), 371-394.")
  }
}
