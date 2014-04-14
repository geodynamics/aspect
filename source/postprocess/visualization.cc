/*
  Copyright (C) 2011, 2012, 2013, 2014 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/
/*  $Id$  */


#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/lexical_cast.hpp>

#include <math.h>
#include <stdio.h>
#include <unistd.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Interface<dim>::~Interface ()
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
      void
      Interface<dim>::save (std::map<std::string,std::string> &) const
      {}


      template <int dim>
      void
      Interface<dim>::load (const std::map<std::string,std::string> &)
      {}
    }


    template <int dim>
    Visualization<dim>::Visualization ()
      :
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      next_output_time (std::numeric_limits<double>::quiet_NaN()),
      output_file_number (0)
    {}



    template <int dim>
    Visualization<dim>::~Visualization ()
    {
      // make sure a thread that may still be running in the background,
      // writing data, finishes
      background_thread.join ();
    }



    template <int dim>
    void Visualization<dim>::mesh_changed_signal()
    {
      mesh_changed = true;
    }



    template <int dim>
    std::pair<std::string,std::string>
    Visualization<dim>::execute (TableHandler &statistics)
    {
      // if this is the first time we get here, set the next output time
      // to the current time. this makes sure we always produce data during
      // the first time step
      if (std::isnan(next_output_time))
        {
          next_output_time = this->get_time();
        }

      // see if graphical output is requested at this time
      if (this->get_time() < next_output_time)
        return std::pair<std::string,std::string>();


      // create a DataOut object on the heap; ownership of this
      // object will later be transferred to a different thread
      // that will write data in the background. the other thread
      // will then also destroy the object
      DataOut<dim> data_out;
      data_out.attach_dof_handler (this->get_dof_handler());

      // add the primary variables
      std::vector<std::string> solution_names (dim, "velocity");
      solution_names.push_back ("p");
      solution_names.push_back ("T");
      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        solution_names.push_back ("C_" + boost::lexical_cast<std::string>(c+1));


      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      interpretation (dim,
                      DataComponentInterpretation::component_is_part_of_vector);
      interpretation.push_back (DataComponentInterpretation::component_is_scalar);
      interpretation.push_back (DataComponentInterpretation::component_is_scalar);
      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);

      data_out.add_data_vector (this->get_solution(),
                                solution_names,
                                DataOut<dim>::type_dof_data,
                                interpretation);


//      data_out.add_data_vector(free_surface_dof_handler, mesh_vertices, "meshvertices");

      // then for each additional selected output variable
      // add the computed quantity as well. keep a list of
      // pointers to data vectors created by cell data visualization
      // postprocessors that will later be deleted
      std::list<std_cxx1x::shared_ptr<Vector<float> > > cell_data_vectors;
      for (typename std::list<std_cxx1x::shared_ptr<VisualizationPostprocessors::Interface<dim> > >::const_iterator
           p = postprocessors.begin(); p!=postprocessors.end(); ++p)
        {
          try
            {
              // there are two ways of writing visualization postprocessors:
              // - deriving from DataPostprocessor
              // - deriving from DataVectorCreator
              // treat them in turn
              if (const DataPostprocessor<dim> *viz_postprocessor
                  = dynamic_cast<const DataPostprocessor<dim>*>(& **p))
                {
                  data_out.add_data_vector (this->get_solution(),
                                            *viz_postprocessor);
                }
              else if (const VisualizationPostprocessors::CellDataVectorCreator<dim> *
                       cell_data_creator
                       = dynamic_cast<const VisualizationPostprocessors::CellDataVectorCreator<dim>*>
                         (& **p))
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
                  cell_data_vectors.push_back (std_cxx1x::shared_ptr<Vector<float> >
                                               (cell_data.second));

                  data_out.add_data_vector (*cell_data.second,
                                            cell_data.first,
                                            DataOut<dim>::type_cell_data);
                }
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
              std::cerr << "Exception on MPI process <"
                        << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                        << "> while running visualization postprocessor <"
                        << typeid(**p).name()
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
              std::cerr << "Exception on MPI process <"
                        << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                        << "> while running visualization postprocessor <"
                        << typeid(**p).name()
                        << ">: " << std::endl;
              std::cerr << "Unknown exception!" << std::endl
                        << "Aborting!" << std::endl
                        << "----------------------------------------------------"
                        << std::endl;

              // terminate the program!
              MPI_Abort (MPI_COMM_WORLD, 1);
            }

        }

      // now build the patches and see how we can output these
      data_out.build_patches ();

      std::string solution_file_prefix = "solution-" + Utilities::int_to_string (output_file_number, 5);
      std::string mesh_file_prefix = "mesh-" + Utilities::int_to_string (output_file_number, 5);
      if (output_format=="hdf5")
        {
          XDMFEntry new_xdmf_entry;
          std::string     h5_solution_file_name = solution_file_prefix + ".h5";
          std::string     xdmf_filename = this->get_output_directory() + "solution.xdmf";

          // Filter redundant values if the functionality is available in the current
          // version of deal.II, otherwise use the old data format
#if DEAL_II_VERSION_MAJOR*100 + DEAL_II_VERSION_MINOR > 800
          DataOutBase::DataOutFilter   data_filter(DataOutBase::DataOutFilterFlags(true, true));

          // If the mesh changed since the last output, make a new mesh file
          if (mesh_changed) last_mesh_file_name = mesh_file_prefix + ".h5";
          data_out.write_filtered_data(data_filter);
          data_out.write_hdf5_parallel(data_filter,
                                       mesh_changed,
                                       (this->get_output_directory()+last_mesh_file_name).c_str(),
                                       (this->get_output_directory()+h5_solution_file_name).c_str(),
                                       this->get_mpi_communicator());
          new_xdmf_entry = data_out.create_xdmf_entry(data_filter,
                                                      last_mesh_file_name.c_str(),
                                                      h5_solution_file_name.c_str(),
                                                      this->get_time(),
                                                      this->get_mpi_communicator());
#else
          data_out.write_hdf5_parallel((this->get_output_directory()+h5_solution_file_name).c_str(),
                                       this->get_mpi_communicator());
          new_xdmf_entry = data_out.create_xdmf_entry(h5_solution_file_name.c_str(),
                                                      this->get_time(),
                                                      this->get_mpi_communicator());
#endif
          xdmf_entries.push_back(new_xdmf_entry);
          data_out.write_xdmf_file(xdmf_entries, xdmf_filename.c_str(),
                                   this->get_mpi_communicator());
          mesh_changed = false;
        }
      else if ((output_format=="vtu") && (group_files!=0))
        {
          //TODO: There is some code duplication between the following two
          //code blocks. unify!
          AssertThrow(group_files==1, ExcNotImplemented());
          data_out.write_vtu_in_parallel((this->get_output_directory() + solution_file_prefix +
                                          ".vtu").c_str(),
                                         this->get_mpi_communicator());

          if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::vector<std::string> filenames;
              filenames.push_back (solution_file_prefix + ".vtu");
              const std::string pvtu_master_filename = (solution_file_prefix + ".pvtu");
              std::ofstream pvtu_master ((this->get_output_directory() +
                                          pvtu_master_filename).c_str());
              data_out.write_pvtu_record (pvtu_master, filenames);

              // now also generate a .pvd file that matches simulation
              // time and corresponding .pvtu record
              times_and_pvtu_names.push_back(std::make_pair(this->get_time(),
                                                            pvtu_master_filename));
              const std::string
              pvd_master_filename = (this->get_output_directory() + "solution.pvd");
              std::ofstream pvd_master (pvd_master_filename.c_str());
              data_out.write_pvd_record (pvd_master, times_and_pvtu_names);

              // finally, do the same for Visit via the .visit file for this
              // time step, as well as for all time steps together
              const std::string
              visit_master_filename = (this->get_output_directory() +
                                       solution_file_prefix +
                                       ".visit");
              std::ofstream visit_master (visit_master_filename.c_str());
              data_out.write_visit_record (visit_master, filenames);

              output_file_names_by_timestep.push_back (filenames);
#if (DEAL_II_MAJOR*100 + DEAL_II_MINOR) > 800
              std::ofstream global_visit_master ((this->get_output_directory() +
                                                  "solution.visit").c_str());
              data_out.write_visit_record (global_visit_master, output_file_names_by_timestep);
#endif
            }
        }
      else
        {
          // put the stuff we want to write into a string object that
          // we can then write in the background
          const std::string *file_contents;
          {
            std::ostringstream tmp;

            // if deal.II supports it (after 7.3.x), pass time step number and time as
            // metadata into the output file
            DataOutBase::VtkFlags vtk_flags;
#if (DEAL_II_MAJOR*100 + DEAL_II_MINOR) >= 704
            vtk_flags.cycle = this->get_timestep_number();
            vtk_flags.time = this->get_time();
#endif
            data_out.set_flags (vtk_flags);

            data_out.write (tmp, DataOutBase::parse_output_format(output_format));
            file_contents = new std::string (tmp.str());
          }

          // let the master processor write the master record for all the distributed
          // files
          if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::vector<std::string> filenames;
              for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++i)
                filenames.push_back (solution_file_prefix +
                                     "." +
                                     Utilities::int_to_string(i, 4) +
                                     DataOutBase::default_suffix
                                     (DataOutBase::parse_output_format(output_format)));
              const std::string
              pvtu_master_filename = (solution_file_prefix +
                                      ".pvtu");
              std::ofstream pvtu_master ((this->get_output_directory() +
                                          pvtu_master_filename).c_str());
              data_out.write_pvtu_record (pvtu_master, filenames);

              // now also generate a .pvd file that matches simulation
              // time and corresponding .pvtu record
              times_and_pvtu_names.push_back(std::pair<double,std::string>
                                             (this->get_time(), pvtu_master_filename));
              const std::string
              pvd_master_filename = (this->get_output_directory() + "solution.pvd");
              std::ofstream pvd_master (pvd_master_filename.c_str());
              data_out.write_pvd_record (pvd_master, times_and_pvtu_names);

              // finally, do the same for Visit via the .visit file for this
              // time step, as well as for all time steps together
              const std::string
              visit_master_filename = (this->get_output_directory() +
                                       solution_file_prefix +
                                       ".visit");
              std::ofstream visit_master (visit_master_filename.c_str());
              data_out.write_visit_record (visit_master, filenames);

              output_file_names_by_timestep.push_back (filenames);
#if (DEAL_II_MAJOR*100 + DEAL_II_MINOR) > 800
              std::ofstream global_visit_master ((this->get_output_directory() +
                                                  "solution.visit").c_str());
              data_out.write_visit_record (global_visit_master, output_file_names_by_timestep);
#endif
            }

          const std::string *filename
            = new std::string (this->get_output_directory() +
                               solution_file_prefix +
                               "." +
                               Utilities::int_to_string
                               (this->get_triangulation().locally_owned_subdomain(), 4) +
                               DataOutBase::default_suffix
                               (DataOutBase::parse_output_format(output_format)));

          // wait for all previous write operations to finish, should
          // any be still active
          background_thread.join ();

          // then continue with writing our own stuff
          background_thread = Threads::new_thread (&background_writer,
                                                   filename,
                                                   file_contents);
        }

      // record the file base file name in the output file
      statistics.add_value ("Visualization file name",
                            this->get_output_directory() + solution_file_prefix);

      // up the counter of the number of the file by one; also
      // up the next time we need output
      ++output_file_number;
      set_next_output_time (this->get_time());

      // return what should be printed to the screen.
      return std::make_pair (std::string ("Writing graphical output:"),
                             this->get_output_directory() + solution_file_prefix);
    }


    template <int dim>
    void Visualization<dim>::background_writer (const std::string *filename,
                                                const std::string *file_contents)
    {
      // write stuff into a (hopefully local) tmp file first. to do so first
      // find out whether $TMPDIR is set and if so put the file in there
      std::string tmp_filename;

      {
        // Try getting the environment variable for the temporary directory
        const char *tmp_filedir = getenv("TMPDIR");
        // If we can't, default to /tmp
        if (tmp_filedir)
          tmp_filename = tmp_filedir;
        else
          tmp_filename = "/tmp";
        tmp_filename += "/aspect.tmp.XXXXXX";

        // Create the temporary file; get at the actual filename
        // by using a C-style string that mkstemp will then overwrite
        char *tmp_filename_x = new char[tmp_filename.size()+1];
        std::strcpy(tmp_filename_x, tmp_filename.c_str());
        const int tmp_file_desc = mkstemp(tmp_filename_x);
        tmp_filename = tmp_filename_x;
        delete []tmp_filename_x;

        // If we failed to create the temp file, just write directly to the target file.
        // We also provide a warning about this fact. There are places where
        // this fails *on every node*, so we will get a lot of warning messages
        // into the output; in these cases, just writing multiple pieces to
        // std::cerr will produce an unreadable mass of text; rather, first
        // assemble the error message completely, and then output it atomically
        if (tmp_file_desc == -1)
          {
            std::string x = std::string(
                              "***** WARNING: could not create temporary file, will "
                              "output directly to final location. This may negatively "
                              "affect performance. (On processor ")
                            + Utilities::int_to_string(
                              Utilities::MPI::this_mpi_process (MPI_COMM_WORLD))
                            + ".)\n";

            std::cerr << x << std::flush;

            tmp_filename = *filename;
          }
      }

      // open the file. if we can't open it, abort if this is the "real"
      // file. re-try with the "real" file if we had tried to write to
      // a temporary file
    re_try_with_non_tmp_file:
      std::ofstream out (tmp_filename.c_str());
      if (!out)
        {
          if (tmp_filename == *filename)
            AssertThrow (false, ExcMessage(std::string("Trying to write to file <") +
                                           *filename +
                                           " but the file can't be opened!"))
            else
              {
                tmp_filename = *filename;
                goto re_try_with_non_tmp_file;
              }
        }

      // now write and then move the tmp file to its final destination
      // if necessary
      out << *file_contents;

      if (tmp_filename != *filename)
        {
          std::string command = std::string("mv ") + tmp_filename + " " + *filename;

          bool first_attempt = true;

        re_try:
          int error = system(command.c_str());

          // if the move failed, and this is the first time, sleep for a second in
          // hopes that it was just an NFS timeout, then try again. if it fails the
          // second time around, try writing to the final file directly.
          if (error != 0)
            {
              if (first_attempt == true)
                {
                  first_attempt = false;
                  sleep (1);
                  goto re_try;
                }
              else
                {
                  std::cerr << "***** WARNING: could not move " << tmp_filename
                            << " to " << *filename << ". Trying again to write directly to "
                            << *filename
                            << ". (On processor "
                            << Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) << ".)"
                            << std::endl;
                  tmp_filename = *filename;
                  goto re_try_with_non_tmp_file;
                }
            }
        }

      // destroy the pointers to the data we needed to write
      delete file_contents;
      delete filename;
    }


    namespace
    {
      std_cxx1x::tuple
      <void *,
      void *,
      aspect::internal::Plugins::PluginList<VisualizationPostprocessors::Interface<2> >,
      aspect::internal::Plugins::PluginList<VisualizationPostprocessors::Interface<3> > > registered_plugins;
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
                             Patterns::Double (0),
                             "The time interval between each generation of "
                             "graphical output files. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");

          // now also see about the file format we're supposed to write in
          prm.declare_entry ("Output format", "vtu",
                             Patterns::Selection (DataOutBase::get_output_format_names ()),
                             "The file format to be used for graphical output.");

          prm.declare_entry ("Number of grouped files", "0",
                             Patterns::Integer(0),
                             "VTU file output supports grouping files from several CPUs "
                             "into one file using MPI I/O when writing on a parallel "
                             "filesystem. Select 0 for no grouping. This will disable "
                             "parallel file output and instead write one file per processor "
                             "in a background thread. "
                             "A value of 1 will generate one big file containing the whole "
                             "solution.");

          // finally also construct a string for Patterns::MultipleSelection that
          // contains the names of all registered visualization postprocessors
          const std::string pattern_of_names
            = std_cxx1x::get<dim>(registered_plugins).get_pattern_of_names (true);
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
                            std_cxx1x::get<dim>(registered_plugins).get_description_string());
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // now declare the parameters of each of the registered
      // visualization postprocessors in turn
      std_cxx1x::get<dim>(registered_plugins).declare_parameters (prm);
    }


    template <int dim>
    void
    Visualization<dim>::parse_parameters (ParameterHandler &prm)
    {
      Assert (std_cxx1x::get<dim>(registered_plugins).plugins != 0,
              ExcMessage ("No postprocessors registered!?"));
      std::vector<std::string> viz_names;

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Visualization");
        {
          output_interval = prm.get_double ("Time between graphical output");
          output_format   = prm.get ("Output format");
          group_files     = prm.get_integer("Number of grouped files");

          // now also see which derived quantities we are to compute
          viz_names = Utilities::split_string_list(prm.get("List of output variables"));

          // see if 'all' was selected (or is part of the list). if so
          // simply replace the list with one that contains all names
          if (std::find (viz_names.begin(),
                         viz_names.end(),
                         "all") != viz_names.end())
            {
              viz_names.clear();
              for (typename std::list<typename aspect::internal::Plugins::PluginList<VisualizationPostprocessors::Interface<dim> >::PluginInfo>::const_iterator
                   p = std_cxx1x::get<dim>(registered_plugins).plugins->begin();
                   p != std_cxx1x::get<dim>(registered_plugins).plugins->end(); ++p)
                viz_names.push_back (std_cxx1x::get<0>(*p));
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // then go through the list, create objects and let them parse
      // their own parameters
      for (unsigned int name=0; name<viz_names.size(); ++name)
        {
          VisualizationPostprocessors::Interface<dim> *
          viz_postprocessor = std_cxx1x::get<dim>(registered_plugins)
                              .create_plugin (viz_names[name],
                                              "Visualization plugins",
                                              prm);

          // make sure that the postprocessor is indeed of type
          // dealii::DataPostprocessor or of type
          // VisualizationPostprocessors::CellDataVectorCreator
          Assert ((dynamic_cast<DataPostprocessor<dim>*>(viz_postprocessor)
                   != 0)
                  ||
                  (dynamic_cast<VisualizationPostprocessors::CellDataVectorCreator<dim>*>(viz_postprocessor)
                   != 0)
                  ,
                  ExcMessage ("Can't convert visualization postprocessor to type "
                              "dealii::DataPostprocessor or "
                              "VisualizationPostprocessors::CellDataVectorCreator!?"));

          postprocessors.push_back (std_cxx1x::shared_ptr<VisualizationPostprocessors::Interface<dim> >
                                    (viz_postprocessor));
        }
    }


    template <int dim>
    template <class Archive>
    void Visualization<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &next_output_time
      & output_file_number
      & times_and_pvtu_names
      & output_file_names_by_timestep
      & mesh_changed
      & last_mesh_file_name
#if DEAL_II_VERSION_MAJOR*100 + DEAL_II_VERSION_MINOR > 800
      & xdmf_entries
#endif
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

      // set next output time to something useful
      set_next_output_time (this->get_time());
    }


    template <int dim>
    void
    Visualization<dim>::set_next_output_time (const double current_time)
    {
      // if output_interval is positive, then set the next output interval to
      // a positive multiple.
      if (output_interval > 0)
        {
          // the current time is always in seconds, so we need to convert the output_interval to the same unit
          double output_interval_in_s = (this->convert_output_to_years()) ? (output_interval*year_in_seconds) : output_interval;

          // we need to compute the smallest integer that is bigger than current_time/my_output_interval,
          // even if it is a whole number already (otherwise we output twice in a row)
          next_output_time = (std::floor(current_time/output_interval_in_s)+1.0) * output_interval_in_s;
        }
    }


    template <int dim>
    void
    Visualization<dim>::initialize (const Simulator<dim> &simulator)
    {
      // first call the respective function in the base class
      SimulatorAccess<dim>::initialize (simulator);

      // pass initialization through to the various visualization
      // objects if they so desire
      for (typename std::list<std_cxx1x::shared_ptr<VisualizationPostprocessors::Interface<dim> > >::iterator
           p = postprocessors.begin();
           p != postprocessors.end(); ++p)
        // see if a given visualization plugin is in fact derived
        // from the SimulatorAccess class, and if so initialize it.
        // note that viz plugins need not necessarily derive from
        // SimulatorAccess if they don't need anything beyond the
        // solution variables to compute what they compute
        if (SimulatorAccess<dim> *x = dynamic_cast<SimulatorAccess<dim>*>(& **p))
          x->initialize (simulator);

      // Also set up a listener to check when the mesh changes
      mesh_changed = true;
      this->get_triangulation().signals.post_refinement.connect(std_cxx1x::bind(&Visualization::mesh_changed_signal, std_cxx1x::ref(*this)));
    }


    template <int dim>
    void
    Visualization<dim>::
    register_visualization_postprocessor (const std::string &name,
                                          const std::string &description,
                                          void (*declare_parameters_function) (ParameterHandler &),
                                          VisualizationPostprocessors::Interface<dim> *(*factory_function) ())
    {
      std_cxx1x::get<dim>(registered_plugins).register_plugin (name,
                                                               description,
                                                               declare_parameters_function,
                                                               factory_function);
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
    }
  }

  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<Postprocess::VisualizationPostprocessors::Interface<3> >::plugins = 0;
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
