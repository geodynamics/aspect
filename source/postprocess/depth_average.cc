//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/postprocess/depth_average.h>
#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <math.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    template <class Archive>
    void DepthAverage<dim>::DataPoint::serialize (Archive &ar, const unsigned int version)
    {
      ar &time &depth &value;
    }


    template <int dim>
    DepthAverage<dim>::DepthAverage ()
      :
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      next_output_time (std::numeric_limits<double>::quiet_NaN())
    {}



    template <int dim>
    std::pair<std::string,std::string>
    DepthAverage<dim>::execute (TableHandler &statistics)
    {
      // if this is the first time we get here, set the next output time
      // to the current time. this makes sure we always produce data during
      // the first time step
      if (isnan(next_output_time))
        next_output_time = this->get_time();

      // see if output is requested at this time
      if (this->get_time() < next_output_time)
        return std::pair<std::string,std::string>();

      std::vector<double> temp;
      this->get_depth_average_temperature(temp);

      const double max_depth = this->get_geometry_model().maximal_depth();

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          // store data from the current step
          for (unsigned int i=0; i<temp.size(); ++i)
            {
              DataPoint data_point;
              data_point.time  = this->get_time();
              data_point.depth = max_depth*i/temp.size();
              data_point.value = temp[i];
              entries.push_back(data_point);
            }

          // leave a marker for the end of this time step. we'll later
          // use this to leave a blank line in the output file at this
          // position
          DataPoint data_point;
          data_point.depth = -1.0;
          entries.push_back(data_point);

          // write out the file
          const std::string filename = this->get_output_directory() + "depthaverageT.plt";
          std::ofstream f1 (filename.c_str());
          for (unsigned int i=0; i<entries.size(); ++i)
            {
              if (entries[i].depth==-1.0)
                f1 << std::endl;
              else
                f1 << entries[i].time << " "
                   << entries[i].depth << " "
                   << entries[i].value << std::endl;
            }
          entries.clear();
        }

      temp.clear();
      this->get_depth_average_velocity_magnitude(temp);

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          // store data from the current step
          for (unsigned int i=0; i<temp.size(); ++i)
            {
              DataPoint data_point;
              data_point.time  = this->get_time();
              data_point.depth = max_depth*i/temp.size();
              data_point.value = temp[i];
              entries.push_back(data_point);
            }

          // leave a marker for the end of this time step. we'll later
          // use this to leave a blank line in the output file at this
          // position
          DataPoint data_point;
          data_point.depth = -1.0;
          entries.push_back(data_point);

          // write out the file
          const std::string filename = this->get_output_directory() + "depthaveragevelocitymagnitude.plt";
          std::ofstream f2 (filename.c_str());
          for (unsigned int i=0; i<entries.size(); ++i)
            {
              if (entries[i].depth==-1)
                f2 << std::endl;
              else
                f2 << entries[i].time << " "
                   << entries[i].depth << " "
                   << entries[i].value << std::endl;
            }
          entries.clear();
        }

      temp.clear();
      this->get_depth_average_sinking_velocity(temp);

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          // store data from the current step
          for (unsigned int i=0; i<temp.size(); ++i)
            {
              DataPoint data_point;
              data_point.time  = this->get_time();
              data_point.depth = max_depth*i/temp.size();
              data_point.value = temp[i];
              entries.push_back(data_point);
            }

          // leave a marker for the end of this time step. we'll later
          // use this to leave a blank line in the output file at this
          // position
          DataPoint data_point;
          data_point.depth = -1.0;
          entries.push_back(data_point);

          // write out the file
          const std::string filename = this->get_output_directory() + "depthaveragesinkingvelocity.plt";
          std::ofstream f3 (filename.c_str());
          for (unsigned int i=0; i<entries.size(); ++i)
            {
              if (entries[i].depth == -1.0)
                f3 << std::endl;
              else
                f3 << entries[i].time << " "
                   << entries[i].depth << " "
                   << entries[i].value << std::endl;
            }
          entries.clear();
        }

      temp.clear();
      this->get_depth_average_Vs(temp);

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // store all data from the current step
          for (unsigned int i=0; i<temp.size(); ++i)
            {
              DataPoint data_point;
              data_point.time  = this->get_time();
              data_point.depth = max_depth*i/temp.size();
              data_point.value = temp[i];
              entries.push_back(data_point);
            }

          // leave a marker for the end of this time step. we'll later
          // use this to leave a blank line in the output file at this
          // position
          DataPoint data_point;
          data_point.depth = -1.0;
          entries.push_back(data_point);

          // write out the file
          const std::string filename = this->get_output_directory() + "depthaverageVs.plt";
          std::ofstream f4 (filename.c_str());
          for (unsigned int i=0; i<entries.size(); ++i)
            {
              if (entries[i].depth == -1.0)
                f4 << std::endl;
              else
                f4 << entries[i].time << " "
                   << entries[i].depth << " "
                   << entries[i].value << std::endl;
            }
          entries.clear();
        }

      temp.clear();
      this->get_depth_average_Vp(temp);

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          // store data from the current step
          for (unsigned int i=0; i<temp.size(); ++i)
            {
              DataPoint data_point;
              data_point.time  = this->get_time();
              data_point.depth = max_depth*i/temp.size();
              data_point.value = temp[i];
              entries.push_back(data_point);
            }

          // leave a marker for the end of this time step. we'll later
          // use this to leave a blank line in the output file at this
          // position
          DataPoint data_point;
          data_point.depth = -1.0;
          entries.push_back(data_point);

          // write out the file
          const std::string filename = this->get_output_directory() + "depthaverageVp.plt";
          std::ofstream f5 (filename.c_str());
          for (unsigned int i=0; i<entries.size(); ++i)
            {
              if (entries[i].depth == -1.0)
                f5 << std::endl;
              else
                f5 << entries[i].time << " "
                   << entries[i].depth << " "
                   << entries[i].value << std::endl;
            }
        }
      temp.clear();
      this->get_depth_average_viscosity(temp);

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          // store data from the current step
          for (unsigned int i=0; i<temp.size(); ++i)
            {
              DataPoint data_point;
              data_point.time  = this->get_time();
              data_point.depth = max_depth*i/temp.size();
              data_point.value = temp[i];
              entries.push_back(data_point);
            }

          // leave a marker for the end of this time step. we'll later
          // use this to leave a blank line in the output file at this
          // position
          DataPoint data_point;
          data_point.depth = -1.0;
          entries.push_back(data_point);

          // write out the file
          const std::string filename = this->get_output_directory() + "depthaverageViscosity.plt";
          std::ofstream f6 (filename.c_str());
          for (unsigned int i=0; i<entries.size(); ++i)
            {
              if (entries[i].depth == -1.0)
                f6 << std::endl;
              else
                f6 << entries[i].time << " "
                << entries[i].depth << " "
                << entries[i].value << std::endl;
            }
        }

      set_next_output_time (this->get_time());

      // return what should be printed to the screen. note that we had
      // just incremented the number, so use the previous value
      return std::make_pair (std::string ("Writing depth average"),
                             this->get_output_directory() + "depthaverageT.plt");
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
      boost::archive::text_oarchive oa (os);
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
          boost::archive::text_iarchive ia (is);
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
      // a positive multiple; we need to interpret output_interval either
      // as years or as seconds
      if (output_interval > 0)
        {
          if (this->convert_output_to_years() == true)
            next_output_time = std::ceil(current_time / (output_interval * year_in_seconds)) *
                               (output_interval * year_in_seconds);
          else
            next_output_time = std::ceil(current_time / (output_interval )) *
                               (output_interval * year_in_seconds);
        }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    template class DepthAverage<deal_II_dimension>;

    ASPECT_REGISTER_POSTPROCESSOR(DepthAverage,
                                  "depth average",
                                  "A postprocessor that computes depth averaged quantities and writes them out.")
  }
}
