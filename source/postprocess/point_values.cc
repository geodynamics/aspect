/*
  Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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


#include <aspect/postprocess/point_values.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/global.h>
#include <deal.II/numerics/vector_tools.h>

#include <math.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    PointValues<dim>::PointValues ()
      :
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      last_output_time (std::numeric_limits<double>::quiet_NaN()),
      evaluation_points_cartesian (std::vector<Point<dim> >() ),
      point_values (std::vector<std::pair<double, std::vector<Vector<double> > > >() ),
      use_natural_coordinates (false)
    {}

    template <int dim>
    std::pair<std::string,std::string>
    PointValues<dim>::execute (TableHandler &)
    {
      // if this is the first time we get here, set the next output time
      // to the current time. this makes sure we always produce data during
      // the first time step
      if (std::isnan(last_output_time))
        last_output_time = this->get_time() - output_interval;

      // see if output is requested at this time
      if (this->get_time() < last_output_time + output_interval)
        return std::pair<std::string,std::string>();

      // evaluate the solution at all of our evaluation points
      std::vector<Vector<double> >
      current_point_values (evaluation_points_cartesian.size(),
                            Vector<double> (this->introspection().n_components));

      for (unsigned int p=0; p<evaluation_points_cartesian.size(); ++p)
        {
          // try to evaluate the solution at this point. in parallel, the point
          // will be on only one processor's owned cells, so the others are
          // going to throw an exception. make sure at least one processor
          // finds the given point
          bool point_found = false;

          try
            {
              VectorTools::point_value(this->get_mapping(),
                                       this->get_dof_handler(),
                                       this->get_solution(),
                                       evaluation_points_cartesian[p],
                                       current_point_values[p]);
              point_found = true;
            }
          catch (const VectorTools::ExcPointNotAvailableHere &)
            {
              // ignore
            }

          // ensure that at least one processor found things
          const int n_procs = Utilities::MPI::sum (point_found ? 1 : 0, this->get_mpi_communicator());
          AssertThrow (n_procs > 0,
                       ExcMessage ("While trying to evaluate the solution at point " +
                                   Utilities::to_string(evaluation_points_cartesian[p][0]) + ", " +
                                   Utilities::to_string(evaluation_points_cartesian[p][1]) +
                                   (dim == 3
                                    ?
                                    ", " + Utilities::to_string(evaluation_points_cartesian[p][2])
                                    :
                                    "") + "), " +
                                   "no processors reported that the point lies inside the " +
                                   "set of cells they own. Are you trying to evaluate the " +
                                   "solution at a point that lies outside of the domain?"
                                  ));

          // Reduce all collected values into local Vector
          Utilities::MPI::sum (current_point_values[p], this->get_mpi_communicator(),
                               current_point_values[p]);

          // Normalize in cases where points are claimed by multiple processors
          if (n_procs > 1)
            current_point_values[p] /= n_procs;
        }

      // finally push these point values all onto the list we keep
      point_values.emplace_back (this->get_time(), current_point_values);

      // now write all of the data to the file of choice. start with a pre-amble that
      // explains the meaning of the various fields
      const std::string filename = (this->get_output_directory() +
                                    "point_values.txt");

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {

          std::ofstream f (filename.c_str());
          f << ("# <time> "
                "<evaluation_point_x> "
                "<evaluation_point_y> ")
            << (dim == 3 ? "<evaluation_point_z> " : "")
            << ("<velocity_x> "
                "<velocity_y> ")
            << (dim == 3 ? "<velocity_z> " : "")
            << "<pressure> <temperature>";
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            f << " <" << this->introspection().name_for_compositional_index(c) << ">";
          f << '\n';

          for (const auto &time_point : point_values)
            {
              Assert (time_point.second.size() == evaluation_points_cartesian.size(),
                      ExcInternalError());
              for (unsigned int i=0; i<evaluation_points_cartesian.size(); ++i)
                {
                  f << /* time = */ time_point.first / (this->convert_output_to_years() ? year_in_seconds : 1.)
                    << ' '
                    << /* location = */ evaluation_points_cartesian[i] << ' ';

                  for (unsigned int c=0; c<time_point.second[i].size(); ++c)
                    {
                      // output a data element. internally, we store all point
                      // values in the same format in which they were computed,
                      // but we convert velocities to meters per year if so
                      // requested
                      if ((this->introspection().component_masks.velocities[c] == false)
                          ||
                          (this->convert_output_to_years() == false))
                        f << time_point.second[i][c];
                      else
                        f << time_point.second[i][c] * year_in_seconds;

                      f << (c != time_point.second[i].size()-1 ? ' ' : '\n');
                    }
                }

              // have an empty line between time steps
              f << '\n';
            }

          AssertThrow (f, ExcMessage("Writing data to <" + filename +
                                     "> did not succeed in the `point values' "
                                     "postprocessor."));
        }


      // Update time
      set_last_output_time (this->get_time());

      // return what should be printed to the screen. note that we had
      // just incremented the number, so use the previous value
      return std::make_pair (std::string ("Writing point values:"),
                             filename);
    }


    template <int dim>
    void
    PointValues<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Point values");
        {
          prm.declare_entry ("Time between point values output", "0.",
                             Patterns::Double (0.),
                             "The time interval between each generation of "
                             "point values output. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry("Evaluation points", "",
                            // a list of points, separated by semicolons; each point has
                            // exactly 'dim' components/coordinates, separated by commas
                            Patterns::List (Patterns::List (Patterns::Double(), dim, dim, ","),
                                            0, Patterns::List::max_int_value, ";"),
                            "The list of points at which the solution should be evaluated. "
                            "Points need to be separated by semicolons, and coordinates of "
                            "each point need to be separated by commas.");
          prm.declare_entry("Use natural coordinates", "false",
                            Patterns::Bool (),
                            "Whether or not the Evaluation points are specified in "
                            "the natural coordinates of the geometry model, e.g. "
                            "radius, lon, lat for the chunk model. "
                            "Currently, natural coordinates for the spherical shell "
                            "and sphere geometries are not supported. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    PointValues<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Point values");
        {
          output_interval = prm.get_double ("Time between point values output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;

          const std::vector<std::string> point_list
            = Utilities::split_string_list(prm.get("Evaluation points"), ';');

          std::vector<std::array<double,dim> > evaluation_points;

          for (unsigned int p=0; p<point_list.size(); ++p)
            {
              const std::vector<std::string> coordinates
                = Utilities::split_string_list(point_list[p], ',');
              AssertThrow (coordinates.size() == dim,
                           ExcMessage ("In setting up the list of evaluation points for the <Point values> "
                                       "postprocessor, one of the evaluation points reads <"
                                       + point_list[p] +
                                       ">, but this does not correspond to a list of numbers with "
                                       "as many coordinates as you run your simulation in."));

              std::array<double,dim> point;
              for (unsigned int d=0; d<dim; ++d)
                point[d] = Utilities::string_to_double (coordinates[d]);
              evaluation_points.push_back (point);
            }

          use_natural_coordinates = prm.get_bool("Use natural coordinates");

          // Convert the vector of coordinate arrays in Cartesian or natural
          // coordinates to a vector of Point<dim> of Cartesian coordinates.
          evaluation_points_cartesian.resize(evaluation_points.size());
          for (unsigned int p=0; p<evaluation_points.size(); ++p)
            {
              if (use_natural_coordinates)
                evaluation_points_cartesian[p] = this->get_geometry_model().natural_to_cartesian_coordinates(evaluation_points[p]);
              else
                for (unsigned int i = 0; i < dim; i++)
                  evaluation_points_cartesian[p][i] = evaluation_points[p][i];
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    template <class Archive>
    void PointValues<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &evaluation_points_cartesian
      & point_values
      & last_output_time;
    }


    template <int dim>
    void
    PointValues<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      std::ostringstream os;
      aspect::oarchive oa (os);
      oa << (*this);

      status_strings["PointValues"] = os.str();
    }


    template <int dim>
    void
    PointValues<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("PointValues") != status_strings.end())
        {
          std::istringstream is (status_strings.find("PointValues")->second);
          aspect::iarchive ia (is);
          ia >> (*this);
        }
    }


    template <int dim>
    void
    PointValues<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then set the next output interval to
      // a positive multiple.
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
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(PointValues,
                                  "point values",
                                  "A postprocessor that evaluates the solution (i.e., velocity, pressure, "
                                  "temperature, and compositional fields along with other fields that "
                                  "are treated as primary variables) at the end of every time step or "
                                  "after a user-specified time interval "
                                  "at a given set of points and then writes this data into the file "
                                  "<point\\_values.txt> in the output directory. The points at which "
                                  "the solution should be evaluated are specified in the section "
                                  "\\texttt{Postprocess/Point values} in the input file."
                                  "\n\n"
                                  "In the output file, data is organized as (i) time, (ii) the 2 or 3 "
                                  "coordinates of the evaluation points, and (iii) followed by the "
                                  "values of the solution vector at this point. The time is provided "
                                  "in seconds or, if the "
                                  "global ``Use years in output instead of seconds'' parameter is "
                                  "set, in years. In the latter case, the velocity is also converted "
                                  "to meters/year, instead of meters/second."
                                  "\n\n"
                                  "\\note{Evaluating the solution of a finite element field at "
                                  "arbitrarily chosen points is an expensive process. Using this "
                                  "postprocessor will only be efficient if the number of evaluation "
                                  "points or output times is relatively small. If you need a very large number of "
                                  "evaluation points, you should consider extracting this "
                                  "information from the visualization program you use to display "
                                  "the output of the `visualization' postprocessor.}")
  }
}
