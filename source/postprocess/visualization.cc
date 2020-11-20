/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/mesh_deformation/interface.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include <algorithm>
#include <type_traits>

#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      /**
       * This Postprocessor will generate the output variables of velocity,
       * pressure, temperature, and compositional fields. They can not be added
       * directly if the velocity needs to be converted from m/s to m/year, so
       * this is what this class does.
       */
      template <int dim>
      class BaseVariablePostprocessor: public DataPostprocessor< dim >, public SimulatorAccess<dim>
      {
        public:

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const override
          {
            const double velocity_scaling_factor =
              this->convert_output_to_years() ? year_in_seconds : 1.0;
            const unsigned int n_q_points = input_data.solution_values.size();
            for (unsigned int q=0; q<n_q_points; ++q)
              for (unsigned int i=0; i<computed_quantities[q].size(); ++i)
                {
                  // scale velocities and fluid velocities by year_in_seconds if needed
                  if (this->introspection().component_masks.velocities[i] ||
                      (this->include_melt_transport()
                       && this->introspection().variable("fluid velocity").component_mask[i]))
                    computed_quantities[q][i] = input_data.solution_values[q][i] * velocity_scaling_factor;
                  else
                    computed_quantities[q][i] = input_data.solution_values[q][i];
                }
          }

          std::vector<std::string> get_names () const override
          {
            std::vector<std::string> solution_names (dim, "velocity");

            if (this->include_melt_transport())
              {
                solution_names.emplace_back("p_f");
                solution_names.emplace_back("p_c_bar");
                for (unsigned int i=0; i<dim; ++i)
                  solution_names.emplace_back("u_f");
              }
            solution_names.emplace_back("p");

            solution_names.emplace_back("T");
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              solution_names.push_back (this->introspection().name_for_compositional_index(c));

            return solution_names;
          }


          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const override
          {
            std::vector<DataComponentInterpretation::DataComponentInterpretation>
            interpretation (dim,
                            DataComponentInterpretation::component_is_part_of_vector);
            if (this->include_melt_transport())
              {
                interpretation.push_back (DataComponentInterpretation::component_is_scalar);
                interpretation.push_back (DataComponentInterpretation::component_is_scalar);
                for (unsigned int i=0; i<dim; ++i)
                  interpretation.push_back (DataComponentInterpretation::component_is_part_of_vector);
              }
            interpretation.push_back (DataComponentInterpretation::component_is_scalar); // p
            interpretation.push_back (DataComponentInterpretation::component_is_scalar); // T
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              interpretation.push_back (DataComponentInterpretation::component_is_scalar);

            return interpretation;
          }

          UpdateFlags get_needed_update_flags () const override
          {
            return update_values;
          }
      };

      /**
       * This Postprocessor will generate the output variables of mesh velocity
       * for when a deforming mesh is used.
       */
      template <int dim>
      class MeshDeformationPostprocessor: public DataPostprocessorVector< dim >, public SimulatorAccess<dim>
      {
        public:
          MeshDeformationPostprocessor ()
            : DataPostprocessorVector<dim>( "mesh_velocity", UpdateFlags(update_values) )
          {}


          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const override
          {
            // check that the first quadrature point has dim components
            Assert( computed_quantities[0].size() == dim,
                    ExcMessage("Unexpected dimension in mesh velocity postprocessor"));
            const double velocity_scaling_factor =
              this->convert_output_to_years() ? year_in_seconds : 1.0;
            const unsigned int n_q_points = input_data.solution_values.size();
            for (unsigned int q=0; q<n_q_points; ++q)
              for (unsigned int i=0; i<dim; ++i)
                computed_quantities[q][i] = input_data.solution_values[q][i] * velocity_scaling_factor;
          }
      };
    }


    namespace VisualizationPostprocessors
    {


      template <int dim>
      Interface<dim>::~Interface ()
      {}


      template <int dim>
      void
      Interface<dim>::initialize ()
      {}

      template <int dim>
      void
      Interface<dim>::update ()
      {}


      template <int dim>
      void
      Interface<dim>::declare_parameters (ParameterHandler &)
      {}



      template <int dim>
      void
      Interface<dim>::parse_parameters (ParameterHandler &)
      {}



      template <int dim>
      std::list<std::string>
      Interface<dim>::required_other_postprocessors () const
      {
        return std::list<std::string>();
      }



      template <int dim>
      void
      Interface<dim>::save (std::map<std::string,std::string> &) const
      {}


      template <int dim>
      void
      Interface<dim>::load (const std::map<std::string,std::string> &)
      {}
    }


    template <int dim>
    Visualization<dim>::OutputHistory::OutputHistory()
      :
      // Record that we have not yet written any meshes and so
      // there are none we can reference from future output operations.
      mesh_changed (true)
    {}



    template <int dim>
    Visualization<dim>::OutputHistory::~OutputHistory()
    {
      // Make sure that any thread that may still be running in the background,
      // writing data, finishes
      background_thread.join ();
    }



    template <int dim>
    template <class Archive>
    void Visualization<dim>::OutputHistory::serialize (Archive &ar, const unsigned int)
    {
      ar
      & last_mesh_file_name
      & times_and_pvtu_names
      & output_file_names_by_timestep
      & xdmf_entries
      ;

      // We do not serialize mesh_changed but use the default (true) from our
      // constructor. This will result in a new mesh file the first time we
      // create visualization output after resuming from a snapshot. Otherwise
      // we might get corrupted graphical output, because the ordering of
      // vertices can be different after resuming. (This is because when
      // deal.II builds a triangulation, it tears down and re-creates its mesh every time
      // p4est changes the partitioning; this process introduces a history where
      // cells and vertices may be numbered differently if the triangulation
      // had previous content than when -- as is done after a restart -- the
      // triangulation is empty.)
    }



    template <int dim>
    Visualization<dim>::Visualization ()
      :
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      last_output_time (std::numeric_limits<double>::quiet_NaN()),
      maximum_timesteps_between_outputs (std::numeric_limits<int>::max()),
      last_output_timestep (numbers::invalid_unsigned_int),
      output_file_number (numbers::invalid_unsigned_int)
    {}



    template <int dim>
    void Visualization<dim>::mesh_changed_signal()
    {
      cell_output_history.mesh_changed = true;
      face_output_history.mesh_changed = true;
    }



    template <int dim>
    template <typename DataOutType>
    void
    Visualization<dim>::write_master_files (const DataOutType &data_out,
                                            const std::string &solution_file_prefix,
                                            const std::vector<std::string> &filenames,
                                            OutputHistory                  &output_history) const
    {
      static_assert (std::is_same<DataOutType,DataOut<dim>>::value ||
                     std::is_same<DataOutType,DataOutFaces<dim>>::value,
                     "The only allowed template types of this function are "
                     "DataOut and DataOutFaces.");
      const bool is_cell_data_output = std::is_same<DataOutType,DataOut<dim>>::value;

      const double time_in_years_or_seconds = (this->convert_output_to_years() ?
                                               this->get_time() / year_in_seconds :
                                               this->get_time());
      const std::string pvtu_master_filename = (solution_file_prefix +
                                                ".pvtu");
      std::ofstream pvtu_master ((this->get_output_directory() + "solution/" +
                                  pvtu_master_filename).c_str());
      data_out.write_pvtu_record (pvtu_master, filenames);

      // now also generate a .pvd file that matches simulation
      // time and corresponding .pvtu record
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        {
          // in case we output all nonlinear iterations, we only want one
          // entry per time step, so replace the last line with the current iteration
          if (this->get_nonlinear_iteration() == 0)
            output_history.times_and_pvtu_names.emplace_back(time_in_years_or_seconds, "solution/"+pvtu_master_filename);
          else
            output_history.times_and_pvtu_names.back() = (std::make_pair
                                                          (time_in_years_or_seconds, "solution/"+pvtu_master_filename));
        }
      else
        output_history.times_and_pvtu_names.emplace_back(time_in_years_or_seconds, "solution/"+pvtu_master_filename);

      const std::string pvd_master_filename = (this->get_output_directory() +
                                               (is_cell_data_output ? "solution.pvd" : "solution_surface.pvd"));
      std::ofstream pvd_master (pvd_master_filename.c_str());

      DataOutBase::write_pvd_record (pvd_master, output_history.times_and_pvtu_names);

      // finally, do the same for Visit via the .visit file for this
      // time step, as well as for all time steps together
      const std::string visit_master_filename = (this->get_output_directory()
                                                 + "solution/"
                                                 + solution_file_prefix
                                                 + ".visit");
      std::ofstream visit_master (visit_master_filename.c_str());

      DataOutBase::write_visit_record (visit_master, filenames);

      {
        // the global .visit file needs the relative path because it sits a
        // directory above
        std::vector<std::string> filenames_with_path;
        filenames_with_path.reserve(filenames.size());
        for (const auto &filename : filenames)
          {
            filenames_with_path.emplace_back("solution/" + filename);
          }

        if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
          {
            // in case we output all nonlinear iterations, we only want one
            // entry per time step, so replace the last line with the current iteration
            if (this->get_nonlinear_iteration() == 0)
              output_history.output_file_names_by_timestep.push_back (filenames_with_path);
            else
              output_history.output_file_names_by_timestep.back() = filenames_with_path;
          }
        else
          output_history.output_file_names_by_timestep.push_back (filenames_with_path);
      }

      std::ofstream global_visit_master ((this->get_output_directory() +
                                          (is_cell_data_output ? "solution.visit" : "solution_surface.visit")).c_str());

      std::vector<std::pair<double, std::vector<std::string> > > times_and_output_file_names;
      for (unsigned int timestep=0; timestep<output_history.times_and_pvtu_names.size(); ++timestep)
        times_and_output_file_names.push_back(std::make_pair(output_history.times_and_pvtu_names[timestep].first,
                                                             output_history.output_file_names_by_timestep[timestep]));
      DataOutBase::write_visit_record (global_visit_master, times_and_output_file_names);
    }



    template <int dim>
    void
    Visualization<dim>::update ()
    {
      //Call the .update() method for each visualization postprocessor.
      for (const auto &p : postprocessors)
        p->update();
    }



    template <int dim>
    template <typename DataOutType>
    std::string
    Visualization<dim>::write_data_out_data(DataOutType   &data_out,
                                            OutputHistory &output_history) const
    {
      static_assert (std::is_same<DataOutType,DataOut<dim>>::value ||
                     std::is_same<DataOutType,DataOutFaces<dim>>::value,
                     "The only allowed template types of this function are "
                     "DataOut and DataOutFaces.");
      const bool is_cell_data_output = std::is_same<DataOutType,DataOut<dim>>::value;

      const double time_in_years_or_seconds = (this->convert_output_to_years() ?
                                               this->get_time() / year_in_seconds :
                                               this->get_time());

      std::string solution_file_prefix =
        (is_cell_data_output ? "solution-" : "solution_surface-")
        + Utilities::int_to_string (output_file_number, 5);
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        solution_file_prefix.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

      if (output_format == "hdf5")
        {
          XDMFEntry new_xdmf_entry;
          const std::string h5_solution_file_name = "solution/"
                                                    + solution_file_prefix + ".h5";
          const std::string xdmf_filename = "solution.xdmf";
          // Filter redundant values if requested in the input file
          DataOutBase::DataOutFilter data_filter(
            DataOutBase::DataOutFilterFlags(filter_output, true));
          // If the mesh changed since the last output, make a new mesh file
          const std::string mesh_file_prefix =
            (is_cell_data_output ? "mesh-" : "mesh_surface-")
            + Utilities::int_to_string(output_file_number, 5);
          if (output_history.mesh_changed)
            output_history.last_mesh_file_name = "solution/" + mesh_file_prefix + ".h5";

          data_out.write_filtered_data(data_filter);
          data_out.write_hdf5_parallel(data_filter,
                                       output_history.mesh_changed,
                                       this->get_output_directory() + output_history.last_mesh_file_name,
                                       this->get_output_directory() + h5_solution_file_name,
                                       this->get_mpi_communicator());
          new_xdmf_entry = data_out.create_xdmf_entry(data_filter,
                                                      output_history.last_mesh_file_name,
                                                      h5_solution_file_name,
                                                      time_in_years_or_seconds, this->get_mpi_communicator());
          output_history.xdmf_entries.push_back(new_xdmf_entry);
          data_out.write_xdmf_file(output_history.xdmf_entries,
                                   this->get_output_directory() + xdmf_filename,
                                   this->get_mpi_communicator());
          output_history.mesh_changed = false;
        }
      else if (output_format == "vtu")
        {
          // Write master files (.pvtu,.pvd,.visit) on the master process
          const int my_id = Utilities::MPI::this_mpi_process(
                              this->get_mpi_communicator());
          if (my_id == 0)
            {
              std::vector<std::string> filenames;
              const unsigned int n_processes = Utilities::MPI::n_mpi_processes(
                                                 this->get_mpi_communicator());
              const unsigned int n_files =
                (group_files == 0) ?
                n_processes : std::min(group_files, n_processes);
              for (unsigned int i = 0; i < n_files; ++i)
                filenames.push_back(
                  solution_file_prefix + "."
                  + Utilities::int_to_string(i, 4) + ".vtu");
              write_master_files(data_out, solution_file_prefix, filenames, output_history);
            }
          const unsigned int n_processes = Utilities::MPI::n_mpi_processes(
                                             this->get_mpi_communicator());
          const unsigned int my_file_id = (
                                            group_files == 0 ? my_id : my_id % group_files);
          const std::string filename = this->get_output_directory() + "solution/"
                                       + solution_file_prefix + "."
                                       + Utilities::int_to_string(my_file_id, 4) + ".vtu";
          // pass time step number and time as metadata into the output file
          DataOutBase::VtkFlags vtk_flags;
          vtk_flags.cycle = this->get_timestep_number();
          vtk_flags.time = time_in_years_or_seconds;

          vtk_flags.write_higher_order_cells = write_higher_order_output;

          data_out.set_flags(vtk_flags);
          // Write as many files as processes. For this case we support writing in a
          // background thread and to a temporary location, so we first write everything
          // into a string that is written to disk in a writer function
          if ((group_files == 0) || (group_files >= n_processes))
            {
              // Put the content we want to write into a string object that
              // we can then write in the background
              const std::string *file_contents;
              {
                std::ostringstream tmp;
                data_out.write(tmp,
                               DataOutBase::parse_output_format(output_format));
                file_contents = new std::string(tmp.str());
              }
              if (write_in_background_thread)
                {
                  // Wait for all previous write operations to finish, should
                  // any be still active,
                  output_history.background_thread.join();
                  // then continue with writing our own data.
                  output_history.background_thread = Threads::new_thread(&writer,
                                                                         filename, temporary_output_location, file_contents);
                }
              else
                writer(filename, temporary_output_location, file_contents);
            }
          else
            // Just write one data file in parallel
            if (group_files == 1)
              {
                data_out.write_vtu_in_parallel(filename.c_str(),
                                               this->get_mpi_communicator());
              }
            else               // Write as many output files as 'group_files' groups
              {
                int color = my_id % group_files;
                MPI_Comm comm;
                int ierr = MPI_Comm_split(this->get_mpi_communicator(), color, my_id, &comm);
                AssertThrowMPI(ierr);
                data_out.write_vtu_in_parallel(filename.c_str(), comm);
                ierr = MPI_Comm_free(&comm);
                AssertThrowMPI(ierr);
              }
        }
      else   // Write in a different format than hdf5 or vtu. This case is supported, but is not
        // optimized for parallel output in that every process will write one file directly
        // into the output directory. This may or may not affect performance depending on
        // the model setup and the network file system type.
        {
          const unsigned int myid = Utilities::MPI::this_mpi_process(
                                      this->get_mpi_communicator());
          const std::string filename = this->get_output_directory() + "solution/"
                                       + solution_file_prefix + "." + Utilities::int_to_string(myid, 4)
                                       + DataOutBase::default_suffix(
                                         DataOutBase::parse_output_format(output_format));
          std::ofstream out(filename.c_str());
          AssertThrow(out,
                      ExcMessage(
                        "Unable to open file for writing: " + filename + "."));
          data_out.write(out, DataOutBase::parse_output_format(output_format));
        }

      return solution_file_prefix;
    }



    template <int dim>
    std::pair<std::string,std::string>
    Visualization<dim>::execute (TableHandler &statistics)
    {
      // if this is the first time we get here, set the last output time
      // to the current time - output_interval. this makes sure we
      // always produce data during the first time step
      if (std::isnan(last_output_time))
        {
          last_output_time = this->get_time() - output_interval;
          last_output_timestep = this->get_timestep_number();
        }

      // Return if graphical output is not requested at this time. Do not
      // return in the first timestep, or if the last output was more than
      // output_interval in time ago, or maximum_timesteps_between_outputs in
      // number of timesteps ago.
      // The comparison in number of timesteps is safe from integer overflow for
      // at most 2 billion timesteps , which is not likely to
      // be ever reached (both values are unsigned int,
      // and the default value of maximum_timesteps_between_outputs is
      // set to numeric_limits<int>::max())
      if ((this->get_time() < last_output_time + output_interval)
          && (this->get_timestep_number() < last_output_timestep + maximum_timesteps_between_outputs)
          && (this->get_timestep_number() != 0))
        return std::pair<std::string,std::string>();

      // up the counter of the number of the file by one, but not in
      // the very first output step. if we run postprocessors on all
      // iterations, only increase file number in the first nonlinear iteration
      const bool increase_file_number = (this->get_nonlinear_iteration() == 0) || (!this->get_parameters().run_postprocessors_on_nonlinear_iterations);
      if (output_file_number == numbers::invalid_unsigned_int)
        output_file_number = 0;
      else if (increase_file_number)
        ++output_file_number;

      internal::BaseVariablePostprocessor<dim> base_variables;
      base_variables.initialize_simulator (this->get_simulator());

      std::unique_ptr<internal::MeshDeformationPostprocessor<dim> > mesh_deformation_variables;

      DataOut<dim> data_out;
      data_out.attach_dof_handler (this->get_dof_handler());
      data_out.add_data_vector (this->get_solution(),
                                base_variables);

      // Also create an object for outputting information that lives on
      // the faces of the mesh. If there are postprocessors derived from
      // the VisualizationPostprocessors::SurfaceOnlyVisualization class, then
      // we will use this object for viz purposes.
      DataOutFaces<dim> data_out_faces;
      data_out_faces.attach_dof_handler (this->get_dof_handler());
      const bool have_face_viz_postprocessors
        = (std::find_if (postprocessors.begin(),
                         postprocessors.end(),
                         [](const std::unique_ptr<VisualizationPostprocessors::Interface<dim> > &p)
      {
        return (dynamic_cast<const VisualizationPostprocessors::SurfaceOnlyVisualization<dim>*>
                (p.get()) != nullptr);
      })
      != postprocessors.end());

      // If there is a deforming mesh, also attach the mesh velocity object
      if ( this->get_parameters().mesh_deformation_enabled && output_mesh_velocity)
        {
          mesh_deformation_variables = std_cxx14::make_unique<internal::MeshDeformationPostprocessor<dim>>();
          mesh_deformation_variables->initialize_simulator(this->get_simulator());
          data_out.add_data_vector (this->get_mesh_velocity(),
                                    *mesh_deformation_variables);
        }

      // then for each additional selected output variable
      // add the computed quantity as well. keep a list of
      // pointers to data vectors created by cell data visualization
      // postprocessors that will later be deleted
      std::list<std::unique_ptr<Vector<float> > > cell_data_vectors;
      for (const auto &p : postprocessors)
        {
          try
            {
              // There are two ways of writing visualization postprocessors:
              // - deriving from DataPostprocessor
              // - deriving from DataVectorCreator
              // treat them in turn. In both cases, the information can
              // be output on all cells, or via the faces on the surface
              // only (if the class in question is derived from
              // SurfaceOnlyVisualization), so we will have to switch between
              // the two ways when we send things to the data_out or
              // data_out_faces objects.
              if (const DataPostprocessor<dim> *viz_postprocessor
                  = dynamic_cast<const DataPostprocessor<dim>*>(& *p))
                {
                  if (dynamic_cast<const VisualizationPostprocessors::SurfaceOnlyVisualization<dim>*>
                      (& *p) == nullptr)
                    data_out.add_data_vector (this->get_solution(),
                                              *viz_postprocessor);
                  else
                    data_out_faces.add_data_vector (this->get_solution(),
                                                    *viz_postprocessor);
                }
              else if (const VisualizationPostprocessors::CellDataVectorCreator<dim> *
                       cell_data_creator
                       = dynamic_cast<const VisualizationPostprocessors::CellDataVectorCreator<dim>*>
                         (& *p))
                {
                  // get the data produced here
                  const std::pair<std::string, Vector<float> *>
                  cell_data = cell_data_creator->execute();
                  Assert (cell_data.second->size() ==
                          this->get_triangulation().n_active_cells(),
                          ExcMessage ("Cell data visualization postprocessors must generate "
                                      "vectors that have as many entries as there are active cells "
                                      "on the current processor."));

                  // store the pointer, then attach the vector to the DataOut object
                  cell_data_vectors.push_back (std::unique_ptr<Vector<float> >
                                               (cell_data.second));

                  if (dynamic_cast<const VisualizationPostprocessors::SurfaceOnlyVisualization<dim>*>
                      (& *p) == nullptr)
                    data_out.add_data_vector (*cell_data.second,
                                              cell_data.first,
                                              DataOut<dim>::type_cell_data);
                  else
                    data_out_faces.add_data_vector (*cell_data.second,
                                                    cell_data.first,
                                                    DataOutFaces<dim>::type_cell_data);
                }
              else
                // A viz postprocessor not derived from either DataPostprocessor
                // or CellDataVectorCreator? We don't know what to do with that!
                Assert (false,
                        ExcMessage("The visualization system found a visualization "
                                   "postprocessor class that is not either derived from "
                                   "DataPostprocessor or "
                                   "VisualizationPostprocessors::CellDataVectorCreator. "
                                   "ASPECT does not know what to do with these kinds of "
                                   "classes."));
            }
          // viz postprocessors that throw exceptions usually do not result in
          // anything good because they result in an unwinding of the stack
          // and, if only one processor triggers an exception, the
          // destruction of objects often causes a deadlock. thus, if
          // an exception is generated, catch it, print an error message,
          // and abort the program
          catch (std::exception &exc)
            {
              std::cerr << std::endl << std::endl
                        << "----------------------------------------------------"
                        << std::endl;
              std::cerr << "An exception happened on MPI process <"
                        << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                        << "> while running the visualization postprocessor <"
                        << typeid(*p).name()
                        << ">: " << std::endl
                        << exc.what() << std::endl
                        << "Aborting!" << std::endl
                        << "----------------------------------------------------"
                        << std::endl;

              // terminate the program!
              MPI_Abort (MPI_COMM_WORLD, 1);
            }
          catch (...)
            {
              std::cerr << std::endl << std::endl
                        << "----------------------------------------------------"
                        << std::endl;
              std::cerr << "An exception happened on MPI process <"
                        << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                        << "> while running the visualization postprocessor <"
                        << typeid(*p).name()
                        << ">: " << std::endl;
              std::cerr << "Unknown exception!" << std::endl
                        << "Aborting!" << std::endl
                        << "----------------------------------------------------"
                        << std::endl;

              // terminate the program!
              MPI_Abort (MPI_COMM_WORLD, 1);
            }

        }

      // Now build the patches. If selected, increase the output resolution.
      // Giving the mapping ensures that the case with mesh deformation works correctly.
      const unsigned int subdivisions = interpolate_output
                                        ?
                                        this->get_stokes_velocity_degree()
                                        :
                                        0;

      // Now get everything written for the DataOut case, and record this
      // in the statistics file
      std::string solution_file_prefix;
      {
        data_out.build_patches (this->get_mapping(),
                                subdivisions,
                                this->get_geometry_model().has_curved_elements()
                                ?
                                DataOut<dim>::curved_inner_cells
                                :
                                DataOut<dim>::no_curved_cells);

        solution_file_prefix
          = write_data_out_data(data_out, cell_output_history);
        statistics.add_value ("Visualization file name",
                              this->get_output_directory()
                              + "solution/"
                              + solution_file_prefix);
      }

      // Then do the same again for the face data case. We won't print the
      // output file name to screen (too much clutter on the screen already)
      // but still put it into the statistics file
      if (have_face_viz_postprocessors)
        {
          data_out_faces.build_patches (this->get_mapping(),
                                        subdivisions);

          const std::string face_solution_file_prefix
            = write_data_out_data(data_out_faces, face_output_history);
          statistics.add_value ("Surface visualization file name",
                                this->get_output_directory()
                                + "solution_surface/"
                                + face_solution_file_prefix);
        }

      // Increment the next time we need output:
      set_last_output_time (this->get_time());
      last_output_timestep = this->get_timestep_number();

      // Return what should be printed to the screen. This is a bit
      // late (the output has already been written, and this probably took
      // a good long while), but it's still good to provide a status
      // update.
      return std::make_pair (std::string ("Writing graphical output:"),
                             this->get_output_directory()
                             + "solution/"
                             + solution_file_prefix);
    }


    template <int dim>
    // We need to pass the arguments by value, as this function can be called on a separate thread:
    void Visualization<dim>::writer (const std::string filename, //NOLINT(performance-unnecessary-value-param)
                                     const std::string temporary_output_location, //NOLINT(performance-unnecessary-value-param)
                                     const std::string *file_contents)
    {
      std::string tmp_filename = filename;
      if (temporary_output_location != "")
        {
          tmp_filename = temporary_output_location + "/aspect.tmp.XXXXXX";

          // Create the temporary file; get at the actual filename
          // by using a C-style string that mkstemp will then overwrite
          std::vector<char> tmp_filename_x (tmp_filename.size()+1);
          std::strcpy(tmp_filename_x.data(), tmp_filename.c_str());
          const int tmp_file_desc = mkstemp(tmp_filename_x.data());
          tmp_filename = tmp_filename_x.data();

          // If we failed to create the temp file, just write directly to the target file.
          // We also provide a warning about this fact. There are places where
          // this fails *on every node*, so we will get a lot of warning messages
          // into the output; in these cases, just writing multiple pieces to
          // std::cerr will produce an unreadable mass of text; rather, first
          // assemble the error message completely, and then output it atomically
          if (tmp_file_desc == -1)
            {
              const std::string x = ("***** WARNING: could not create temporary file <"
                                     +
                                     tmp_filename
                                     +
                                     ">, will output directly to final location. This may negatively "
                                     "affect performance. (On processor "
                                     + Utilities::int_to_string(Utilities::MPI::this_mpi_process (MPI_COMM_WORLD))
                                     + ".)\n");

              std::cerr << x << std::flush;

              tmp_filename = filename;
            }
          else
            close(tmp_file_desc);
        }

      std::ofstream out(tmp_filename.c_str());

      AssertThrow (out, ExcMessage(std::string("Trying to write to file <") +
                                   filename +
                                   ">, but the file can't be opened!"))

      // now write and then move the tmp file to its final destination
      // if necessary
      out << *file_contents;
      out.close ();

      if (tmp_filename != filename)
        {
          std::string command = std::string("mv ") + tmp_filename + " " + filename;
          int error = system(command.c_str());

          AssertThrow(error == 0,
                      ExcMessage("Could not move " + tmp_filename + " to "
                                 + filename + ". On processor "
                                 + Utilities::int_to_string(Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)) + "."));
        }

      // destroy the pointer to the data we needed to write
      delete file_contents;
    }


    namespace
    {
      std::tuple
      <void *,
      void *,
      aspect::internal::Plugins::PluginList<VisualizationPostprocessors::Interface<2> >,
      aspect::internal::Plugins::PluginList<VisualizationPostprocessors::Interface<3> > > registered_visualization_plugins;
    }


    template <int dim>
    void
    Visualization<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Visualization");
        {
          prm.declare_entry ("Time between graphical output", "1e8",
                             Patterns::Double (0.),
                             "The time interval between each generation of "
                             "graphical output files. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");

          prm.declare_entry ("Time steps between graphical output", boost::lexical_cast<std::string>(std::numeric_limits<int>::max()),
                             Patterns::Integer(0),
                             "The maximum number of time steps between each generation of "
                             "graphical output files.");

          // now also see about the file format we're supposed to write in
          prm.declare_entry ("Output format", "vtu",
                             Patterns::Selection (DataOutBase::get_output_format_names ()),
                             "The file format to be used for graphical output. The list "
                             "of possible output formats that can be given here is documented "
                             "in the appendix of the manual where the current parameter "
                             "is described.");

          prm.declare_entry ("Number of grouped files", "16",
                             Patterns::Integer(0),
                             "VTU file output supports grouping files from several CPUs "
                             "into a given number of files using MPI I/O when writing on a parallel "
                             "filesystem. Select 0 for no grouping. This will disable "
                             "parallel file output and instead write one file per processor. "
                             "A value of 1 will generate one big file containing the whole "
                             "solution, while a larger value will create that many files "
                             "(at most as many as there are MPI ranks).");

          prm.declare_entry ("Write in background thread", "false",
                             Patterns::Bool(),
                             "File operations can potentially take a long time, blocking the "
                             "progress of the rest of the model run. Setting this variable to "
                             "`true' moves this process into a background thread, while the "
                             "rest of the model continues.");

          prm.declare_entry ("Temporary output location", "",
                             Patterns::Anything(),
                             "On large clusters it can be advantageous to first write the "
                             "output to a temporary file on a local file system and later "
                             "move this file to a network file system. If this variable is "
                             "set to a non-empty string it will be interpreted as a "
                             "temporary storage location.");

          prm.declare_entry ("Interpolate output", "true",
                             Patterns::Bool(),
                             "deal.II offers the possibility to linearly interpolate "
                             "output fields of higher order elements to a finer resolution. "
                             "This somewhat compensates the fact that most visualization "
                             "software only offers linear interpolation between grid points "
                             "and therefore the output file is a very coarse representation "
                             "of the actual solution field. Activating this option increases "
                             "the spatial resolution in each dimension by a factor equal "
                             "to the polynomial degree used for the velocity finite element "
                             "(usually 2). In other words, instead of showing one quadrilateral "
                             "or hexahedron in the visualization per cell on which \\aspect{} "
                             "computes, it shows multiple (for quadratic elements, it will "
                             "describe each cell of the mesh on which we compute as "
                             "$2\\times 2$ or $2\\times 2\\times 2$ cells in 2d and 3d, "
                             "respectively; correspondingly more subdivisions are used if "
                             "you use cubic, quartic, or even higher order elements for the "
                             "velocity)."
                             "\n\n"
                             "The effect of using this option can be seen in the following "
                             "picture showing a variation of the output produced with the "
                             "input files from Section~\\ref{sec:shell-simple-2d}:"
                             "\n\n"
                             "\\begin{center}"
                             "  \\includegraphics[width=0.5\\textwidth]{viz/parameters/build-patches}"
                             "\\end{center}"
                             "Here, the left picture shows one visualization cell per "
                             "computational cell (i.e., the option is switched off), "
                             "and the right picture shows the same simulation with the "
                             "option switched on (which is the default). The images "
                             "show the same data, demonstrating "
                             "that interpolating the solution onto bilinear shape functions as is "
                             "commonly done in visualizing data loses information."
                             "\n\n"
                             "Of course, activating this option also greatly increases the amount of "
                             "data \\aspect{} will write to disk: approximately by a factor of 4 in 2d, "
                             "and a factor of 8 in 3d, when using quadratic elements for the velocity, "
                             "and correspondingly more for even higher order elements.");

          prm.declare_entry ("Point-wise stress and strain", "false",
                             Patterns::Bool(),
                             "If set to true, quantities related to stress and strain are computed "
                             "in each vertex. Otherwise, an average per cell is computed.");

          prm.declare_entry ("Write higher order output", "false",
                             Patterns::Bool(),
                             "deal.II offers the possibility to write vtu files with higher order "
                             "representations of the output data. This means each cell will correctly "
                             "show the higher order representation of the output data instead of the "
                             "linear interpolation between vertices that ParaView and Visit usually show. "
                             "Note that activating this option is safe and recommended, but requires that "
                             "(i) ``Output format'' is set to ``vtu'', (ii) ``Interpolate output'' is "
                             "set to true, (iii) you use a sufficiently new version of Paraview "
                             "or Visit to read the files (Paraview version 5.5 or newer, and Visit version "
                             "to be determined), and (iv) you use deal.II version 9.1.0 or newer. "
                             "\n"
                             "The effect of using this option can be seen in the following "
                             "picture:"
                             "\n\n"
                             "\\begin{center}"
                             "  \\includegraphics[width=0.5\\textwidth]{viz/parameters/higher-order-output}"
                             "\\end{center}"
                             "The top figure shows the plain output without interpolation or higher "
                             "order output. The middle figure shows output that was interpolated as "
                             "discussed for the ``Interpolate output'' option. The bottom panel "
                             "shows higher order output that achieves better accuracy than the "
                             "interpolated output at a lower memory cost.");

          prm.declare_entry ("Filter output", "false",
                             Patterns::Bool(),
                             "deal.II offers the possibility to filter duplicate vertices for HDF5 "
                             "output files. This merges the vertices of adjacent cells and "
                             "therefore saves disk space, but misrepresents discontinuous "
                             "output properties. Activating this function reduces the disk space "
                             "by about a factor of $2^{dim}$ for HDF5 output, and currently has no "
                             "effect on other output formats. "
                             "\\note{\\textbf{Warning:} Setting this flag to true will result in "
                             "visualization output that does not accurately represent discontinuous "
                             "fields. This may be because you are using a discontinuous finite "
                             "element for the pressure, temperature, or compositional variables, "
                             "or because you use a visualization postprocessor that outputs "
                             "quantities as discontinuous fields (e.g., the strain rate, viscosity, "
                             "etc.). These will then all be visualized as \\textit{continuous} "
                             "quantities even though, internally, \\aspect{} considers them as "
                             "discontinuous fields.}");

          prm.declare_entry ("Output mesh velocity", "false",
                             Patterns::Bool(),
                             "For computations with deforming meshes, ASPECT uses an Arbitrary-Lagrangian-"
                             "Eulerian formulation to handle deforming the domain, so the mesh "
                             "has its own velocity field.  This may be written as an output field "
                             "by setting this parameter to true.");

          // Finally also construct a string for Patterns::MultipleSelection that
          // contains the names of all registered visualization postprocessors.
          // Also add a number of removed plugins that are now combined in 'material properties'
          // to keep compatibility with input files. These will be filtered out in parse_parameters().
          const std::string pattern_of_names
            = std::get<dim>(registered_visualization_plugins).get_pattern_of_names ()
              + "|density|specific heat|thermal conductivity|thermal diffusivity|thermal expansivity|viscosity";
          prm.declare_entry("List of output variables",
                            "",
                            Patterns::MultipleSelection(pattern_of_names),
                            "A comma separated list of visualization objects that should be run "
                            "whenever writing graphical output. By default, the graphical "
                            "output files will always contain the primary variables velocity, "
                            "pressure, and temperature. However, one frequently wants to also "
                            "visualize derived quantities, such as the thermodynamic phase "
                            "that corresponds to a given temperature-pressure value, or the "
                            "corresponding seismic wave speeds. The visualization objects do "
                            "exactly this: they compute such derived quantities and place them "
                            "into the output file. The current parameter is the place where "
                            "you decide which of these additional output variables you want "
                            "to have in your output file.\n\n"
                            "The following postprocessors are available:\n\n"
                            +
                            std::get<dim>(registered_visualization_plugins).get_description_string());
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // now declare the parameters of each of the registered
      // visualization postprocessors in turn
      std::get<dim>(registered_visualization_plugins).declare_parameters (prm);
    }


    template <int dim>
    void
    Visualization<dim>::parse_parameters (ParameterHandler &prm)
    {
      Assert (std::get<dim>(registered_visualization_plugins).plugins != nullptr,
              ExcMessage ("No postprocessors registered!?"));
      std::vector<std::string> viz_names;

      std::string visualization_subdirectory = this->get_output_directory() + "solution/";
      Utilities::create_directory (visualization_subdirectory,
                                   this->get_mpi_communicator(),
                                   true);

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Visualization");
        {
          output_interval = prm.get_double ("Time between graphical output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;

          maximum_timesteps_between_outputs = prm.get_integer("Time steps between graphical output");

          if (output_interval > 0.0)
            {
              // since we increase the time indicating when to write the next graphical output
              // every time we execute the visualization postprocessor, there is no good way to
              // figure out when to write graphical output for the nonlinear iterations if we do
              // not want to output every time step
              AssertThrow(this->get_parameters().run_postprocessors_on_nonlinear_iterations == false,
                          ExcMessage("Postprocessing nonlinear iterations is only supported if every time "
                                     "step is visualized, or in other words, if the 'Time between graphical "
                                     "output' in the Visualization postprocessor is set to zero."));
            }

          output_format   = prm.get ("Output format");
          group_files     = prm.get_integer("Number of grouped files");
          write_in_background_thread = prm.get_bool("Write in background thread");
          temporary_output_location = prm.get("Temporary output location");

          if (temporary_output_location != "")
            {
              // Check if a command-processor is available by calling system() with a
              // null pointer. System is guaranteed to return non-zero if it finds
              // a terminal and zero if there is none (like on the compute nodes of
              // some cluster architectures, e.g. IBM BlueGene/Q)
              AssertThrow(system((char *)nullptr) != 0,
                          ExcMessage("Usage of a temporary storage location is only supported if "
                                     "there is a terminal available to move the files to their final location "
                                     "after writing. The system() command did not succeed in finding such a terminal."));
            }

          interpolate_output = prm.get_bool("Interpolate output");
          filter_output = prm.get_bool("Filter output");
          pointwise_stress_and_strain = prm.get_bool("Point-wise stress and strain");
          write_higher_order_output = prm.get_bool("Write higher order output");

          if (write_higher_order_output == true)
            {
              AssertThrow(output_format == "vtu",
                          ExcMessage("The option 'Postprocess/Visualization/Write higher order output' requires the "
                                     "data output format to be set to 'vtu'."));
              AssertThrow(interpolate_output == true,
                          ExcMessage("The input parameter 'Postprocess/Visualization/Write higher order output' "
                                     "requires the input parameter 'Interpolate output' to be set to 'true'."));
            }

          output_mesh_velocity = prm.get_bool("Output mesh velocity");

          // now also see which derived quantities we are to compute
          viz_names = Utilities::split_string_list(prm.get("List of output variables"));
          AssertThrow(Utilities::has_unique_entries(viz_names),
                      ExcMessage("The list of strings for the parameter "
                                 "'Postprocess/Visualization/List of output variables' contains entries more than once. "
                                 "This is not allowed. Please check your parameter file."));

          // see if 'all' was selected (or is part of the list). if so
          // simply replace the list with one that contains all names
          if (std::find (viz_names.begin(),
                         viz_names.end(),
                         "all") != viz_names.end())
            {
              viz_names.clear();
              for (typename std::list<typename aspect::internal::Plugins::PluginList<VisualizationPostprocessors::Interface<dim> >::PluginInfo>::const_iterator
                   p = std::get<dim>(registered_visualization_plugins).plugins->begin();
                   p != std::get<dim>(registered_visualization_plugins).plugins->end(); ++p)
                viz_names.push_back (std::get<0>(*p));
            }

          // Unify material property visualization plugins into the 'material properties'
          // plugin to avoid duplicated code and multiple calls to the material model
          prm.enter_subsection("Material properties");
          {
            bool material_properties_enabled = std::find(viz_names.begin(),
                                                         viz_names.end(),
                                                         "material properties") != viz_names.end() ;

            std::set<std::string> deprecated_postprocessors = {"density",
                                                               "specific heat",
                                                               "thermal conductivity",
                                                               "thermal diffusivity",
                                                               "thermal expansivity",
                                                               "viscosity"
                                                              };

            // For all selected visualization plugins
            auto plugin_name = viz_names.begin();
            while (plugin_name != viz_names.end())
              {
                // Check if the current name is in the set of the deprecated names
                if (deprecated_postprocessors.count(*plugin_name) != 0)
                  {
                    // If there is no 'material properties' yet
                    if (material_properties_enabled == false)
                      {
                        // Set the current property name as the parameter for 'material properties'
                        prm.set("List of material properties",*plugin_name);
                        // Then replace the currently selected plugin with 'material properties'
                        *plugin_name = "material properties";
                        material_properties_enabled = true;
                        ++plugin_name;
                      }
                    else
                      {
                        // Add the current property name to the parameter of 'material properties'
                        std::string new_property_names = prm.get("List of material properties") + ", " + *plugin_name;
                        prm.set("List of material properties",new_property_names);
                        // Then delete the current plugin
                        plugin_name = viz_names.erase(plugin_name);
                      }
                  }
                else
                  ++plugin_name;
              }
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // then go through the list, create objects and let them parse
      // their own parameters
      for (unsigned int name=0; name<viz_names.size(); ++name)
        {
          VisualizationPostprocessors::Interface<dim> *
          viz_postprocessor = std::get<dim>(registered_visualization_plugins)
                              .create_plugin (viz_names[name],
                                              "Visualization plugins");

          // make sure that the postprocessor is indeed of type
          // dealii::DataPostprocessor or of type
          // VisualizationPostprocessors::CellDataVectorCreator
          Assert ((dynamic_cast<DataPostprocessor<dim>*>(viz_postprocessor)
                   != nullptr)
                  ||
                  (dynamic_cast<VisualizationPostprocessors::CellDataVectorCreator<dim>*>(viz_postprocessor)
                   != nullptr)
                  ,
                  ExcMessage ("Can't convert visualization postprocessor to type "
                              "dealii::DataPostprocessor or "
                              "VisualizationPostprocessors::CellDataVectorCreator!?"));

          postprocessors.push_back (std::unique_ptr<VisualizationPostprocessors::Interface<dim> >
                                    (viz_postprocessor));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*postprocessors.back()))
            sim->initialize_simulator (this->get_simulator());

          postprocessors.back()->parse_parameters (prm);
          postprocessors.back()->initialize ();
        }

      // Finally also set up a listener to check when the mesh changes
      cell_output_history.mesh_changed = true;
      face_output_history.mesh_changed = true;
      this->get_triangulation().signals.post_refinement.connect(
        [&]()
      {
        this->mesh_changed_signal();
      });
    }



    template <int dim>
    template <class Archive>
    void Visualization<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &last_output_time
      & last_output_timestep
      & output_file_number
      & cell_output_history
      & face_output_history
      ;
    }


    template <int dim>
    void
    Visualization<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      std::ostringstream os;
      aspect::oarchive oa (os);
      oa << (*this);

      status_strings["Visualization"] = os.str();

//TODO: do something about the visualization postprocessor plugins
    }


    template <int dim>
    void
    Visualization<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("Visualization") != status_strings.end())
        {
          std::istringstream is (status_strings.find("Visualization")->second);
          aspect::iarchive ia (is);
          ia >> (*this);
        }

//TODO: do something about the visualization postprocessor plugins
    }


    template <int dim>
    void
    Visualization<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then update the last supposed output
      // time
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
        }
    }



    template <int dim>
    void
    Visualization<dim>::
    register_visualization_postprocessor (const std::string &name,
                                          const std::string &description,
                                          void (*declare_parameters_function) (ParameterHandler &),
                                          VisualizationPostprocessors::Interface<dim> *(*factory_function) ())
    {
      std::get<dim>(registered_visualization_plugins).register_plugin (name,
                                                                       description,
                                                                       declare_parameters_function,
                                                                       factory_function);
    }



    template <int dim>
    std::list<std::string>
    Visualization<dim>::required_other_postprocessors () const
    {
      std::list<std::string> requirements;

      // loop over all of the viz postprocessors and collect what
      // they want. don't worry about duplicates, the postprocessor
      // manager will filter them out
      for (const auto &p : postprocessors)
        {
          const std::list<std::string> this_requirements = p->required_other_postprocessors();
          requirements.insert (requirements.end(),
                               this_requirements.begin(), this_requirements.end());
        }

      return requirements;
    }



    template <int dim>
    bool
    Visualization<dim>::output_pointwise_stress_and_strain () const
    {
      return pointwise_stress_and_strain;
    }



    template <int dim>
    void
    Visualization<dim>::write_plugin_graph (std::ostream &out)
    {
      // in contrast to all other plugins, the visualization
      // postprocessors do not actually connect to the central
      // Simulator class, but they are a sub-system of
      // the Postprocessor plugin system. indicate this
      // through the last argument of the function call
      std::get<dim>(registered_visualization_plugins)
      .write_plugin_graph ("Visualization postprocessor interface",
                           out,
                           typeid(Visualization<dim>).name());
    }


  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
#define INSTANTIATE(dim) \
  template class Interface<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }

  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<2> >::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<3> >::plugins = nullptr;
    }
  }


  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Visualization,
                                  "visualization",
                                  "A postprocessor that takes the solution and writes "
                                  "it into files that can be read by a graphical "
                                  "visualization program. Additional run time parameters "
                                  "are read from the parameter subsection 'Visualization'.")
  }
}
