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


#include <aspect/postprocess/point_values.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <deal.II/numerics/vector_tools.h>

#include <math.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    PointValues<dim>::execute (TableHandler &)
    {
      // evaluate the solution at all of our evaluation points
      std::vector<Vector<double> >
      current_point_values (evaluation_points.size(),
                            Vector<double> (this->introspection().n_components));

      for (unsigned int p=0; p<evaluation_points.size(); ++p)
        {
          // try to evaluate the solution at this point. in parallel, the point
          // will be on only one processor's owned cells, so the others are
          // going to throw an exception. make sure at least one processor
          // finds the given point
          bool point_found = false;
          try
            {
              VectorTools::point_value(this->get_dof_handler(),
                                       this->get_solution(),
                                       evaluation_points[p],
                                       current_point_values[p]);
              point_found = true;
            }
          catch (const VectorTools::ExcPointNotAvailableHere &)
            {
              // ignore
            }

          // ensure that exactly one processor found things
          Assert (Utilities::MPI::sum (point_found ? 1 : 0, this->get_mpi_communicator()) > 0,
                  ExcMessage ("While trying to evaluate the values of the solution at "
                              "evaluation point " + Utilities::int_to_string(p) +
                              ", no processor reported that that point lies inside the "
                              "set of cells it owns. Are you trying to evaluate the "
                              "solution at a point that lies outside the domain?"));

          // now exchange things. because we have exactly one processor that found
          // the point, we can just add up that value plus all of the zero
          // vectors from the other processors
          //
          // at the time of writing this (where we require deal.II 8.2) the
          // Utilities::MPI::sum function can't sum Vector<double> arguments, so
          // convert everything into a std::vector and back
          std::vector<double> v (current_point_values[p].begin(),
                                 current_point_values[p].end());
          Utilities::MPI::sum (v, this->get_mpi_communicator(),
                               v);
          std::copy (v.begin(), v.end(), current_point_values[p].begin());
        }

      // finally push these point values all onto the list we keep
      point_values.push_back (std::make_pair (this->get_time(),
                                              current_point_values));

      // now write all of the data to the file of choice
      const std::string filename = (this->get_output_directory() +
                                    "point_values.txt");
      std::ofstream f (filename.c_str());
      for (std::vector<std::pair<double, std::vector<Vector<double> > > >::iterator
           time_point = point_values.begin();
           time_point != point_values.end();
           ++time_point)
        {
          Assert (time_point->second.size() == evaluation_points.size(),
                  ExcInternalError());
          for (unsigned int i=0; i<evaluation_points.size(); ++i)
            {
              f << /* time = */ time_point->first / (this->convert_output_to_years() ? year_in_seconds : 1.)
                << ' '
                << /* location = */ evaluation_points[i] << ' ';

              for (unsigned int c=0; c<time_point->second[i].size(); ++c)
                {
                  // output a data element. internally, we store all point
                  // values in the same format in which they were computed,
                  // but we convert velocities to meters per year if so
                  // requested
                  if ((this->introspection().component_masks.velocities[c] == false)
                      ||
                      (this->convert_output_to_years() == false))
                    f << time_point->second[i][c];
                  else
                    f << time_point->second[i][c] * year_in_seconds;

                  f << (c != time_point->second[i].size()-1 ? ' ' : '\n');
                }
            }

          // have an empty line between time steps
          f << '\n';
        }

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
          prm.declare_entry("Evaluation points", "",
                            // a list of points, separated by commas; each point has
                            // exactly 'dim' components/coordinates, separated by spaces
                            Patterns::List (Patterns::List (Patterns::Double(), dim, dim, ","),
                                            0, Patterns::List::max_int_value, ";"),
                            "The list of points at which the solution should be evaluated. "
                            "Points need to be separated by semicolons, and coordinates of "
                            "each point need to be separated by commas.");
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
          const std::vector<std::string> point_list
            = Utilities::split_string_list(prm.get("Evaluation points"), ';');
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

              Point<dim> point;
              for (unsigned int d=0; d<dim; ++d)
                point[d] = Utilities::string_to_double (coordinates[d]);
              evaluation_points.push_back (point);
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
      ar &evaluation_points
      & point_values;
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
                                  "are treated as primary variables) at the end of every time step "
                                  "at a given set of points and then writes this data into the file "
                                  "<point_values.txt> in the output directory. The points at which "
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
                                  "points is relatively small. If you need a very large number of "
                                  "evaluation points, you should consider extracting this "
                                  "information from the visualization program you use to display "
                                  "the output of the 'visualization' postprocessor.}")
  }
}
