//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/postprocess_visualization.h>
#include <aspect/simulator.h>
#include <aspect/equation_data.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      /**
       * A class derived from DataPostprocessor that takes an output vector and
       * computes a number of derived statistics from that.
       *
       * The member functions are all implementations of those declared in the base
       * class. See there for their meaning.
       */
      template <int dim>
      class Postprocessor : public DataPostprocessor<dim>
      {
        public:
          Postprocessor (const unsigned int partition,
                         const double       minimal_pressure);

          virtual
          void
          compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                             const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                             const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                             const std::vector<Point<dim> >                  &normals,
                                             const std::vector<Point<dim> >                  &evaluation_points,
                                             std::vector<Vector<double> >                    &computed_quantities) const;

          virtual std::vector<std::string> get_names () const;

          virtual unsigned int n_output_variables() const;

          virtual
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const;

          virtual UpdateFlags get_needed_update_flags () const;

        private:
          const unsigned int partition;
          const double       minimal_pressure;
      };


      template <int dim>
      Postprocessor<dim>::
      Postprocessor (const unsigned int partition,
                     const double       minimal_pressure)
        :
        partition (partition),
        minimal_pressure (minimal_pressure)
      {}


      template <int dim>
      std::vector<std::string>
      Postprocessor<dim>::get_names() const
      {
        std::vector<std::string> solution_names (dim, "velocity");
        solution_names.push_back ("p");
        solution_names.push_back ("T");
        solution_names.push_back ("friction_heating");
        solution_names.push_back ("partition");
        solution_names.push_back ("viscosity");
        solution_names.push_back ("non_adiabatic_pressure");
        solution_names.push_back ("non_adiabatic_temperature");
        solution_names.push_back ("density");
        return solution_names;
      }


      template <int dim>
      unsigned int
      Postprocessor<dim>::n_output_variables() const
      {
        // make our lives a bit easier here
        return get_names().size();
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      Postprocessor<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation>
        interpretation (dim,
                        DataComponentInterpretation::component_is_part_of_vector);

        interpretation.push_back (DataComponentInterpretation::component_is_scalar);
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);

        return interpretation;
      }


      template <int dim>
      UpdateFlags
      Postprocessor<dim>::get_needed_update_flags() const
      {
        return update_values | update_gradients | update_q_points;
      }


      template <int dim>
      void
      Postprocessor<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &/*dduh*/,
                                         const std::vector<Point<dim> >                  &/*normals*/,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (duh.size() == n_quadrature_points,                  ExcInternalError());
        Assert (computed_quantities.size() == n_quadrature_points,  ExcInternalError());
        Assert (uh[0].size() == dim+2,                              ExcInternalError());
        Assert (computed_quantities[0].size()==n_output_variables(),ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            // velocity; rescale in cm/year
            for (unsigned int d=0; d<dim; ++d)
              computed_quantities[q](d)
                = (uh[q](d) *  EquationData::year_in_seconds * 100);

            // pressure
            const double pressure_at_surface = 1e6;
            const double pressure = (uh[q](dim)-minimal_pressure) + pressure_at_surface;
            computed_quantities[q](dim) = pressure;

            // temperature
            const double temperature = uh[q](dim+1);
            computed_quantities[q](dim+1) = temperature;

            // friction heating
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = duh[q][d];
            const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);
            computed_quantities[q](dim+2) = 2 * EquationData::MaterialModel::eta(temperature, pressure, evaluation_points[q]) *
                                            strain_rate * strain_rate;

            computed_quantities[q](dim+3) = partition;

            computed_quantities[q](dim+4) = EquationData::MaterialModel::real_viscosity(temperature,
                                                                                        pressure,
                                                                                        evaluation_points[q],
                                                                                        strain_rate);

            computed_quantities[q](dim+5) = pressure - EquationData::adiabatic_pressure (evaluation_points[q]);

            computed_quantities[q](dim+6) = temperature -
                                            EquationData::adiabatic_temperature (evaluation_points[q]);

            computed_quantities[q](dim+7) = EquationData::MaterialModel::density(temperature, pressure, evaluation_points[q]);
          }
      }
    }


    template <int dim>
    Visualization<dim>::Visualization (const Simulator<dim> &simulator_object)
      :
      SimulatorAccess<dim> (simulator_object),
      // TODO: do something sensible here
      output_interval (50000),
      next_output_time (0),
      output_file_number (0)
    {}



    template <int dim>
    std::pair<std::string,std::string>
    Visualization<dim>::execute (TableHandler &)
    {
      const FESystem<dim> joint_fe (this->get_stokes_dof_handler().get_fe(), 1,
                                    this->get_temperature_dof_handler().get_fe(), 1);

      DoFHandler<dim> joint_dof_handler (this->get_triangulation());
      joint_dof_handler.distribute_dofs (joint_fe);
      Assert (joint_dof_handler.n_dofs() ==
              this->get_stokes_dof_handler().n_dofs() + this->get_temperature_dof_handler().n_dofs(),
              ExcInternalError());

      TrilinosWrappers::MPI::Vector joint_solution;
      joint_solution.reinit (joint_dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);

      {
        std::vector<unsigned int> local_joint_dof_indices (joint_fe.dofs_per_cell);
        std::vector<unsigned int> local_stokes_dof_indices (this->get_stokes_dof_handler().get_fe().dofs_per_cell);
        std::vector<unsigned int> local_temperature_dof_indices (this->get_temperature_dof_handler().get_fe().dofs_per_cell);

        typename DoFHandler<dim>::active_cell_iterator
        joint_cell       = joint_dof_handler.begin_active(),
        joint_endc       = joint_dof_handler.end(),
        stokes_cell      = this->get_stokes_dof_handler().begin_active(),
        temperature_cell = this->get_temperature_dof_handler().begin_active();
        for (; joint_cell!=joint_endc;
             ++joint_cell, ++stokes_cell, ++temperature_cell)
          if (joint_cell->is_locally_owned())
            {
              joint_cell->get_dof_indices (local_joint_dof_indices);
              stokes_cell->get_dof_indices (local_stokes_dof_indices);
              temperature_cell->get_dof_indices (local_temperature_dof_indices);

              for (unsigned int i=0; i<joint_fe.dofs_per_cell; ++i)
                if (joint_fe.system_to_base_index(i).first.first == 0)
                  {
                    Assert (joint_fe.system_to_base_index(i).second
                            <
                            local_stokes_dof_indices.size(),
                            ExcInternalError());

                    joint_solution(local_joint_dof_indices[i])
                      = this->get_stokes_solution()(local_stokes_dof_indices
                                                    [joint_fe.system_to_base_index(i).second]);
                  }
                else
                  {
                    Assert (joint_fe.system_to_base_index(i).first.first == 1,
                            ExcInternalError());
                    Assert (joint_fe.system_to_base_index(i).second
                            <
                            local_temperature_dof_indices.size(),
                            ExcInternalError());
                    joint_solution(local_joint_dof_indices[i])
                      = this->get_temperature_solution()(local_temperature_dof_indices
                                                         [joint_fe.system_to_base_index(i).second]);
                  }
            }
      }


      IndexSet locally_relevant_joint_dofs(joint_dof_handler.n_dofs());
      DoFTools::extract_locally_relevant_dofs (joint_dof_handler, locally_relevant_joint_dofs);
      TrilinosWrappers::MPI::Vector locally_relevant_joint_solution;
      locally_relevant_joint_solution.reinit (locally_relevant_joint_dofs, MPI_COMM_WORLD);
      locally_relevant_joint_solution = joint_solution;

      internal::Postprocessor<dim> postprocessor (this->get_triangulation().locally_owned_subdomain(),
                                                  this->get_stokes_solution().block(1).minimal_value());

      DataOut<dim> data_out;
      data_out.attach_dof_handler (joint_dof_handler);
      data_out.add_data_vector (locally_relevant_joint_solution, postprocessor);
      data_out.build_patches ();

      const std::string filename = ("bin/solution-" +
                                    Utilities::int_to_string (output_file_number, 5) +
                                    "." +
                                    Utilities::int_to_string
                                    (this->get_triangulation().locally_owned_subdomain(), 4) +
                                    ".vtu");

      //throttle output
      const unsigned int concurrent_writers = 10;
      unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      unsigned int nproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      for (unsigned int i=0; i<nproc; ++i)
        {
          if (i == myid)
            {
              std::ofstream output (filename.c_str());
              if (!output)
                std::cout << "ERROR: proc " << myid << " could not create " << filename << std::endl;
              data_out.write_vtu (output);
            }
          if (i%concurrent_writers == 0)
            {
              sleep(1);
              MPI_Barrier(MPI_COMM_WORLD);
            }

        }


      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          std::vector<std::string> filenames;
          for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); ++i)
            filenames.push_back (std::string("solution-") +
                                 Utilities::int_to_string (output_file_number, 5) +
                                 "." +
                                 Utilities::int_to_string(i, 4) +
                                 ".vtu");
          const std::string
          pvtu_master_filename = ("bin/solution-" +
                                  Utilities::int_to_string (output_file_number, 5) +
                                  ".pvtu");
          std::ofstream pvtu_master (pvtu_master_filename.c_str());
          data_out.write_pvtu_record (pvtu_master, filenames);

          const std::string
          visit_master_filename = ("bin/solution-" +
                                   Utilities::int_to_string (output_file_number, 5) +
                                   ".visit");
          std::ofstream visit_master (visit_master_filename.c_str());
          data_out.write_visit_record (visit_master, filenames);
        }

      // up the counter of the number of the file by one
      ++output_file_number;

      // return what should be printed to the screen
      return std::make_pair (std::string ("Writing graphical output:"),
                             std::string ("bin/solution-") + Utilities::int_to_string (output_file_number, 5));
    }


    template <int dim>
    void
    Visualization<dim>::save (std::map<std::string, std::string> &status_strings) const
    {}


    template <int dim>
    void
    Visualization<dim>::load (const std::map<std::string, std::string> &status_strings)
    {}
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    template class Visualization<deal_II_dimension>;
  }
}
