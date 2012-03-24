//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <math.h>
#include <stdio.h>

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
          Postprocessor (const unsigned int                   partition,
                         const double                         convert_output_to_years,
                         const MaterialModel::Interface<dim> &material_model,
                         const AdiabaticConditions<dim>      &adiabatic_conditions);

          virtual
          void
          compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                             const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                             const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                             const std::vector<Point<dim> >                  &normals,
                                             const std::vector<Point<dim> >                  &evaluation_points,
                                             std::vector<Vector<double> >                    &computed_quantities) const;

          virtual std::vector<std::string> get_names () const;

          virtual
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const;

          virtual UpdateFlags get_needed_update_flags () const;

        private:
          const unsigned int                   partition;
          const double                         convert_output_to_years;
          const MaterialModel::Interface<dim> &material_model;
          const AdiabaticConditions<dim>      &adiabatic_conditions;
      };


      template <int dim>
      Postprocessor<dim>::
      Postprocessor (const unsigned int                   partition,
                     const double                         convert_output_to_years,
                     const MaterialModel::Interface<dim> &material_model,
                     const AdiabaticConditions<dim>      &adiabatic_conditions)
        :
        partition (partition),
        convert_output_to_years (convert_output_to_years),
        material_model(material_model),
        adiabatic_conditions (adiabatic_conditions)
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
        solution_names.push_back ("strain_rate");
        solution_names.push_back ("Vp");
        solution_names.push_back ("Vs");
        solution_names.push_back ("thermodynamic_phase");
        solution_names.push_back ("Cp");
        return solution_names;
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

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            // velocity; rescale in m/year if so requested
            for (unsigned int d=0; d<dim; ++d)
              computed_quantities[q](d)
                = (uh[q](d) *
                   (convert_output_to_years == true ?
                    year_in_seconds :
                    1));

            // pressure
            const double pressure = uh[q](dim);
            computed_quantities[q](dim) = pressure;

            // temperature
            const double temperature = uh[q](dim+1);
            computed_quantities[q](dim+1) = temperature;

            // friction heating
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = duh[q][d];

            const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);
            const SymmetricTensor<2,dim> compressible_strain_rate
              = (material_model.is_compressible()
                 ?
                 strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                 :
                 strain_rate);
            computed_quantities[q](dim+2) = 2 * material_model.viscosity(temperature,
                                                                         pressure,
                                                                         strain_rate,
                                                                         evaluation_points[q]) *
                                            compressible_strain_rate * compressible_strain_rate;

            computed_quantities[q](dim+3) = partition;

            computed_quantities[q](dim+4) = material_model.viscosity(temperature,
                                                                     pressure,
                                                                     strain_rate,
                                                                     evaluation_points[q]);

            computed_quantities[q](dim+5) = pressure - adiabatic_conditions.pressure (evaluation_points[q]);

            computed_quantities[q](dim+6) = temperature -
                                            adiabatic_conditions.temperature (evaluation_points[q]);

            computed_quantities[q](dim+7) = material_model.density(temperature, pressure, evaluation_points[q]);
            computed_quantities[q](dim+8) = std::sqrt(strain_rate*strain_rate);
            computed_quantities[q](dim+9) = material_model.seismic_Vp(temperature, pressure);
            computed_quantities[q](dim+10) = material_model.seismic_Vs(temperature, pressure);
            computed_quantities[q](dim+11) = material_model.thermodynamic_phase(temperature, pressure);
            computed_quantities[q](dim+11) = material_model.specific_heat(temperature, pressure, evaluation_points[q]);
          }
      }
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
    std::pair<std::string,std::string>
    Visualization<dim>::execute (TableHandler &statistics)
    {
      // if this is the first time we get here, set the next output time
      // to the current time. this makes sure we always produce data during
      // the first time step
      if (isnan(next_output_time))
        next_output_time = this->get_time();

      // see if graphical output is requested at this time
      if (this->get_time() < next_output_time)
        return std::pair<std::string,std::string>();


      internal::Postprocessor<dim> postprocessor (this->get_triangulation().locally_owned_subdomain(),
                                                  this->convert_output_to_years(),
                                                  this->get_material_model(),
                                                  this->get_adiabatic_conditions());

      Vector<float> estimated_error_per_cell(this->get_triangulation().n_active_cells());
      this->get_refinement_criteria(estimated_error_per_cell);

      Vector<float> Vs_anomaly(this->get_triangulation().n_active_cells());
      this->get_Vs_anomaly(Vs_anomaly);
      Vector<float> Vp_anomaly(this->get_triangulation().n_active_cells());
      this->get_Vp_anomaly(Vp_anomaly);

      // create a DataOut object on the heap; ownership of this
      // object will later be transferred to a different thread
      // that will write data in the background. the other thread
      // will then also destroy the object
      DataOut<dim> data_out;
      data_out.attach_dof_handler (this->get_dof_handler());
      data_out.add_data_vector (this->get_solution(), postprocessor);
      data_out.add_data_vector (estimated_error_per_cell, "error_indicator");
      data_out.add_data_vector (Vs_anomaly, "Vs_anomaly");
      data_out.add_data_vector (Vp_anomaly, "Vp_anomaly");
      data_out.build_patches ();

      if (output_format=="vtu" && group_files!=0)
        {
          AssertThrow(group_files==1, ExcNotImplemented());
          data_out.write_vtu_in_parallel((this->get_output_directory() + std::string("solution-") +
                                          Utilities::int_to_string (output_file_number, 5) +
                                          ".vtu").c_str(), MPI_COMM_WORLD);

          if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::vector<std::string> filenames;
              filenames.push_back (std::string("solution-") +
                                   Utilities::int_to_string (output_file_number, 5) +
                                   ".vtu");
              const std::string
              pvtu_master_filename = (this->get_output_directory() +
                                      "solution-" +
                                      Utilities::int_to_string (output_file_number, 5) +
                                      ".pvtu");
              std::ofstream pvtu_master (pvtu_master_filename.c_str());
              data_out.write_pvtu_record (pvtu_master, filenames);

              const std::string
              visit_master_filename = (this->get_output_directory() +
                                       "solution-" +
                                       Utilities::int_to_string (output_file_number, 5) +
                                       ".visit");
              std::ofstream visit_master (visit_master_filename.c_str());
              data_out.write_visit_record (visit_master, filenames);
            }
        }
      else
        {

          // put the stuff we want to write into a string object that
          // we can then write in the background
          const std::string *file_contents;
          {
            std::ostringstream tmp;
            data_out.write (tmp, DataOutBase::parse_output_format(output_format));

            file_contents = new std::string (tmp.str());
          }

          // let the master processor write the master record for all the distributed
          // files
          if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::vector<std::string> filenames;
              for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++i)
                filenames.push_back (std::string("solution-") +
                                     Utilities::int_to_string (output_file_number, 5) +
                                     "." +
                                     Utilities::int_to_string(i, 4) +
                                     DataOutBase::default_suffix
                                     (DataOutBase::parse_output_format(output_format)));
              const std::string
              pvtu_master_filename = (this->get_output_directory() +
                                      "solution-" +
                                      Utilities::int_to_string (output_file_number, 5) +
                                      ".pvtu");
              std::ofstream pvtu_master (pvtu_master_filename.c_str());
              data_out.write_pvtu_record (pvtu_master, filenames);

              const std::string
              visit_master_filename = (this->get_output_directory() +
                                       "solution-" +
                                       Utilities::int_to_string (output_file_number, 5) +
                                       ".visit");
              std::ofstream visit_master (visit_master_filename.c_str());
              data_out.write_visit_record (visit_master, filenames);
            }

          const std::string *filename
            = new std::string (this->get_output_directory() +
                               "solution-" +
                               Utilities::int_to_string (output_file_number, 5) +
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
                            this->get_output_directory() + "solution-" +
                            Utilities::int_to_string (output_file_number, 5));

      // up the counter of the number of the file by one; also
      // up the next time we need output
      ++output_file_number;
      set_next_output_time (this->get_time());

      // return what should be printed to the screen. note that we had
      // just incremented the number, so use the previous value
      return std::make_pair (std::string ("Writing graphical output:"),
                             this->get_output_directory() +"solution-" +
                             Utilities::int_to_string (output_file_number-1, 5));
    }


    template <int dim>
    void Visualization<dim>::background_writer (const std::string *filename,
                                                const std::string *file_contents)
    {
      // write stuff into a (hopefully local) tmp file first. to do so first
      // find out whether $TMPDIR is set and if so put the file in there
      char tmp_filename[1025];

      {
        FILE *fp = popen("mktemp ${TMPDIR:-/tmp}/tmp.XXXXXXXXXX", "r");
        AssertThrow (fp, ExcMessage("Couldn't call mktemp!"));

        char *s = fgets (tmp_filename, sizeof(tmp_filename)-1, fp);
        AssertThrow (s!=0, ExcMessage("Couldn't create temporary file!"));

        pclose (fp);

        // tmp_filename also contains the \n from running mktemp
        // strip it
        tmp_filename[strlen(tmp_filename)-1] = '\0';
      }

      // write the data into the file
      {
        std::ofstream output (tmp_filename);
        if (!output)
          std::cout << "***** ERROR: could not create " << tmp_filename
                    << " *****"
                    << std::endl;
        output << *file_contents;
      }

      // now move the file to its final destination on the global file system
      std::string command = std::string("mv ") + tmp_filename + " " + *filename;
      int error = system(command.c_str());
      AssertThrow (error == 0,
                   ExcMessage ("Could not move temporary file to its final location: "
                               +
                               command));

      // destroy the pointers to the data we needed to write
      delete file_contents;
      delete filename;
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
                             Patterns::Selection (DataOutInterface<dim>::get_output_format_names ()),
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
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Visualization<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Visualization");
        {
          output_interval = prm.get_double ("Time between graphical output");
          output_format   = prm.get ("Output format");
          group_files     = prm.get_integer("Number of grouped files");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    template <class Archive>
    void Visualization<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &next_output_time
      & output_file_number;
    }


    template <int dim>
    void
    Visualization<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      std::ostringstream os;
      boost::archive::text_oarchive oa (os);
      oa << (*this);

      status_strings["Visualization"] = os.str();
    }


    template <int dim>
    void
    Visualization<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("Visualization") != status_strings.end())
        {
          std::istringstream is (status_strings.find("Visualization")->second);
          boost::archive::text_iarchive ia (is);
          ia >> (*this);
        }

      // set next output time to something useful
      set_next_output_time (this->get_time());
    }


    template <int dim>
    void
    Visualization<dim>::set_next_output_time (const double current_time)
    {
      // if output_interval is positive, then set the next output interval to
      // a positive multiple; we need to interpret output_interval either
      // as years or as seconds
      if (output_interval > 0)
        {
          if (this->convert_output_to_years() == true)
            next_output_time = std::ceil(current_time / (output_interval * year_in_seconds)) *
                               (output_interval * year_in_seconds);
          else
            next_output_time = std::ceil(current_time / (output_interval )) *
                               output_interval;
        }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    template class Visualization<deal_II_dimension>;

    ASPECT_REGISTER_POSTPROCESSOR(Visualization,
                                  "visualization",
                                  "A postprocessor that takes the solution and writes "
                                  "it into files that can be read by a graphical "
                                  "visualization program. Additional run time parameters "
                                  "are read from the parameter subsection 'Visualization'.")
  }
}
