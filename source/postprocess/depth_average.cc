/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/postprocess/depth_average.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out_stack.h>


#include <math.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    template <class Archive>
    void DepthAverage<dim>::DataPoint::serialize (Archive &ar, const unsigned int version)
    {
      ar &time &values;
    }


    template <int dim>
    DepthAverage<dim>::DepthAverage ()
      :
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      next_output_time (std::numeric_limits<double>::quiet_NaN()),
      n_depth_zones (100)
    {}



    template <int dim>
    std::pair<std::string,std::string>
    DepthAverage<dim>::execute (TableHandler &statistics)
    {
      // if this is the first time we get here, set the next output time
      // to the current time. this makes sure we always produce data during
      // the first time step
      if (std::isnan(next_output_time))
        next_output_time = this->get_time();

      // see if output is requested at this time
      if (this->get_time() < next_output_time)
        return std::pair<std::string,std::string>();

      const unsigned int n_statistics = 7+this->n_compositional_fields();

      DataPoint data_point;
      data_point.time       = this->get_time();
      data_point.values.resize(n_statistics, std::vector<double> (n_depth_zones));

      // add temperature and the compositional fields that follow
      // it immediately
      {
        this->get_depth_average_temperature(data_point.values[0]);
        for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
          this->get_depth_average_composition(c, data_point.values[1+c]);
        this->get_adiabatic_conditions().get_adiabatic_temperature_profile(data_point.values[1+this->n_compositional_fields()]);
        this->get_depth_average_velocity_magnitude(data_point.values[2+this->n_compositional_fields()]);
        this->get_depth_average_sinking_velocity(data_point.values[3+this->n_compositional_fields()]);
        this->get_depth_average_Vs(data_point.values[4+this->n_compositional_fields()]);
        this->get_depth_average_Vp(data_point.values[5+this->n_compositional_fields()]);
        this->get_depth_average_viscosity(data_point.values[6+this->n_compositional_fields()]);
      }
      entries.push_back (data_point);

      const double max_depth = this->get_geometry_model().maximal_depth();

      // On the root process, write out the file. do this using the DataOutStack
      // class on a piecewise constant finite element space on
      // a 1d mesh with the correct subdivisions
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          Triangulation<1> mesh;
          GridGenerator::subdivided_hyper_cube(mesh, n_depth_zones, 0, max_depth);

          FE_DGQ<1> fe(0);
          DoFHandler<1> dof_handler (mesh);
          dof_handler.distribute_dofs(fe);
          Assert (dof_handler.n_dofs() == n_depth_zones, ExcInternalError());

          DataOutStack<1> data_out_stack;
          std::vector<std::string> variables;
          variables.push_back ("temperature");
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            variables.push_back(std::string("C_") + Utilities::int_to_string(c));
          variables.push_back ("adiabatic_temperature");
          variables.push_back ("velocity_magnitude");
          variables.push_back ("sinking_velocity");
          variables.push_back ("Vs");
          variables.push_back ("Vp");
          variables.push_back ("viscosity");
          Assert (variables.size() == n_statistics, ExcInternalError());

          for (unsigned int j=0; j<n_statistics; ++j)
            data_out_stack.declare_data_vector (variables[j],
                                                DataOutStack<1>::cell_vector);

          for (unsigned int i=0; i<entries.size(); ++i)
            {
              data_out_stack.new_parameter_value ((this->convert_output_to_years()
                                                   ?
                                                   entries[i].time / year_in_seconds
                                                   :
                                                   entries[i].time),
                                                  // declare the time step, which here is the difference
                                                  // between successive output times. we don't have anything
                                                  // for the first time step, however. we could do a zero
                                                  // delta, but that leads to invisible output. rather, we
                                                  // use an artificial value of one tenth of the first interval,
                                                  // if available
                                                  (i == 0 ?
                                                   (entries.size() > 1 ? (entries[1].time - entries[0].time)/10 : 0) :
                                                   entries[i].time - entries[i-1].time) /
                                                  (this->convert_output_to_years()
                                                   ?
                                                   year_in_seconds
                                                   :
                                                   1));

              data_out_stack.attach_dof_handler (dof_handler);

              Vector<double> tmp(n_depth_zones);
              for (unsigned int j=0; j<n_statistics; ++j)
                {
                  std::copy (entries[i].values[j].begin(),
                             entries[i].values[j].end(),
                             tmp.begin());
                  data_out_stack.add_data_vector (tmp,
                                                  variables[j]);
                }
              data_out_stack.build_patches ();
              data_out_stack.finish_parameter_value ();
            }


          const std::string filename = (this->get_output_directory() +
                                        "depth_average" +
                                        DataOutBase::default_suffix(output_format));
          std::ofstream f (filename.c_str());
          data_out_stack.write (f, output_format);
        }


      set_next_output_time (this->get_time());

      // return what should be printed to the screen. note that we had
      // just incremented the number, so use the previous value
      return std::make_pair (std::string ("Writing depth average"),
                             this->get_output_directory() +
                             "depth_average" +
                             DataOutBase::default_suffix(output_format));
    }


    template <int dim>
    void
    DepthAverage<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Depth average");
        {
          prm.declare_entry ("Time between graphical output", "1e8",
                             Patterns::Double (0),
                             "The time interval between each generation of "
                             "graphical output files. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Number of zones", "10",
                             Patterns::Integer (1),
                             "The number of zones in depth direction within which we "
                             "are to compute averages. By default, we subdivide the entire "
                             "domain into 10 depth zones and compute temperature and other "
                             "averages within each of these zones. However, if you have a "
                             "very coarse mesh, it may not make much sense to subdivide "
                             "the domain into so many zones and you may wish to choose "
                             "less than this default. It may also make computations slightly "
                             "faster. On the other hand, if you have an extremely highly "
                             "resolved mesh, choosing more zones might also make sense.");
          prm.declare_entry ("Output format", "gnuplot",
                             Patterns::Selection(DataOutBase::get_output_format_names()),
                             "The format in which the output shall be produced. The "
                             "format in which the output is generated also determiens "
                             "the extension of the file into which data is written.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    DepthAverage<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Depth average");
        {
          output_interval = prm.get_double ("Time between graphical output");
          n_depth_zones = prm.get_integer ("Number of zones");
          output_format = DataOutBase::parse_output_format(prm.get("Output format"));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    template <class Archive>
    void DepthAverage<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &next_output_time
      & entries;
    }


    template <int dim>
    void
    DepthAverage<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      std::ostringstream os;
      aspect::oarchive oa (os);
      oa << (*this);

      status_strings["DepthAverage"] = os.str();
    }


    template <int dim>
    void
    DepthAverage<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("DepthAverage") != status_strings.end())
        {
          std::istringstream is (status_strings.find("DepthAverage")->second);
          aspect::iarchive ia (is);
          ia >> (*this);
        }

      // set next output time to something useful
      set_next_output_time (this->get_time());
    }


    template <int dim>
    void
    DepthAverage<dim>::set_next_output_time (const double current_time)
    {
      // if output_interval is positive, then set the next output interval to
      // a positive multiple.
      if (output_interval > 0)
        {
          // the current time is always in seconds, so we need to convert the output_interval to the same unit
          const double output_interval_in_s = (this->convert_output_to_years() ?
                                               (output_interval*year_in_seconds) :
                                               output_interval);

          // we need to compute the smallest integer that is bigger than current_time/my_output_interval,
          // even if it is a whole number already (otherwise we output twice in a row)
          next_output_time = (std::floor(current_time/output_interval_in_s)+1.0) * output_interval_in_s;
        }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(DepthAverage,
                                  "depth average",
                                  "A postprocessor that computes depth averaged "
                                  "quantities and writes them into a file <depth_average.ext> "
                                  "in the output directory, where the extension of the file "
                                  "is determined by the output format you select. In addition "
                                  "to the output format, a number of other parameters also influence "
                                  "this postprocessor, and they can be set in the section "
                                  "\\texttt{Postprocess/Depth average} in the input file. ")
  }
}
