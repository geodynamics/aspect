/*
  Copyright (C) 2021 - 2022 by the authors of the ASPECT code.

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

#include <aspect/termination_criteria/steady_heat_flux.h>
#include <aspect/postprocess/heat_flux_map.h>

namespace aspect
{
  namespace TerminationCriteria
  {
    namespace
    {
      /**
       * A function that trims the handed over list and removes all entries from the front that are
       * further back in time measured from the last entry than given by the first argument.
       * Additionally it makes sure to always keep two entries in the list, if the list had
       * two or more entries. Otherwise the function does not change the list.
       */
      void trim_time_heat_flux_list (const double necessary_time_in_steady_state,
                                     std::list<std::pair<double, double>> &time_heat_flux_list)
      {
        // Remove old times until we're at the correct time period
        // but ensure at least two entries remain in the list (one old, one current timestep)
        auto it = time_heat_flux_list.begin();
        while (time_heat_flux_list.back().first - (*it).first > necessary_time_in_steady_state &&
               std::distance(it,time_heat_flux_list.end()) > 2)
          ++it;

        time_heat_flux_list.erase(time_heat_flux_list.begin(), it);
      }
    }



    template <int dim>
    bool
    SteadyHeatFlux<dim>::execute()
    {
      const std::vector<std::vector<std::pair<double, double>>> heat_flux_and_area =
        Postprocess::internal::compute_heat_flux_through_boundary_faces (*this);

      double local_boundary_fluxes = 0.0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->at_boundary(f))
              {
                const types::boundary_id boundary_indicator
                  = cell->face(f)->boundary_id();
                if (boundary_indicators.find(boundary_indicator) != boundary_indicators.end())
                  local_boundary_fluxes += heat_flux_and_area[cell->active_cell_index()][f].first;
              }

      const double global_heat_flux_integral
        = Utilities::MPI::sum (local_boundary_fluxes, this->get_mpi_communicator());

      // Keep a list of times and heat fluxes at those times
      time_heat_flux.emplace_back(this->get_time(), global_heat_flux_integral);

      // If the length of the simulation time covered in the list is shorter than the
      // specified parameter, we must continue the simulation
      if ((time_heat_flux.size() <= 2)
          ||
          (time_heat_flux.back().first - time_heat_flux.front().first < necessary_time_in_steady_state))
        return false;

      // Remove old entries outside of current time window
      trim_time_heat_flux_list(necessary_time_in_steady_state,time_heat_flux);

      // Scan through the list and calculate the min, mean and max heat flux
      // We assume a linear change of heat flux between times
      double flux_min, flux_max, flux_prev, time_prev, flux_sum=0, flux_mean, deviation_max;
      flux_min = flux_max = flux_prev = time_heat_flux.front().second;
      time_prev = time_heat_flux.front().first;
      for (const auto &it : time_heat_flux)
        {
          flux_min = std::min(flux_min, it.second);
          flux_max = std::max(flux_max, it.second);
          flux_sum += ((it.second + flux_prev)/2.0)*(it.first-time_prev);
          time_prev = it.first;
          flux_prev = it.second;
        }

      flux_mean = flux_sum/(time_heat_flux.back().first-time_heat_flux.front().first);

      // If the min and max are within the acceptable deviation of the mean,
      // we are in steady state and return true, otherwise return false
      deviation_max = std::max(flux_mean - flux_min, flux_max - flux_mean);

      AssertThrow(std::abs(flux_mean) > std::numeric_limits<double>::min(),
                  ExcMessage("The average heat flux is close to 0.0. The "
                             "'steady state heat flux' plugin can not compute a "
                             "relative deviation of heat flux in this case."));

      if (deviation_max/flux_mean > allowed_relative_deviation)
        return false;

      return true;
    }


    template <int dim>
    void
    SteadyHeatFlux<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        prm.enter_subsection("Steady state heat flux");
        {
          prm.declare_entry ("Maximum relative deviation", "0.05",
                             Patterns::Double (0.),
                             "The maximum relative deviation of the heat flux in recent "
                             "simulation time for the system to be considered in "
                             "steady state. If the actual deviation is smaller "
                             "than this number, then the simulation will be terminated.");
          prm.declare_entry ("Time in steady state", "1e7",
                             Patterns::Double (0.),
                             "The minimum length of simulation time that the system "
                             "should be in steady state before termination. "
                             "Note that if the time step size is similar to or larger than "
                             "this value, the termination criterion will only have very few "
                             "(in the most extreme case, just two) heat flux values to check. "
                             "To ensure that a larger number of time steps are included in "
                             "the check for steady state, this value should be much larger "
                             "than the time step size. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Boundary indicators", "",
                             Patterns::List (Patterns::Anything()),
                             "A comma separated list of names denoting those boundaries "
                             "that should be taken into account for integrating the heat "
                             "flux. Note that the plugin will compute the integrated heat "
                             "flux over these boundaries (instead of taking them into "
                             "account individually)."
                             "\n\n"
                             "The names of the boundaries listed here can either be "
                             "numbers (in which case they correspond to the numerical "
                             "boundary indicators assigned by the geometry object), or they "
                             "can correspond to any of the symbolic names the geometry object "
                             "may have provided for each part of the boundary. You may want "
                             "to compare this with the documentation of the geometry model you "
                             "use in your model.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    SteadyHeatFlux<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        prm.enter_subsection("Steady state heat flux");
        {
          allowed_relative_deviation = prm.get_double ("Maximum relative deviation");
          necessary_time_in_steady_state = prm.get_double ("Time in steady state");
          necessary_time_in_steady_state *= this->convert_output_to_years() ? year_in_seconds : 1.0;

          try
            {
              const std::vector<types::boundary_id> x_boundary_indicators
                = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                                      (prm.get ("Boundary indicators")));
              boundary_indicators
                = std::set<types::boundary_id> (x_boundary_indicators.begin(),
                                                x_boundary_indicators.end());

              AssertThrow (boundary_indicators.size() != 0,
                           ExcMessage ("You have specified the steady state heat flux as one of "
                                       "the termination criteria, but the list of boundary indicators "
                                       "that should be used for computing the heat flux is empty."));
            }
          catch (const std::string &error)
            {
              AssertThrow (false, ExcMessage ("While parsing the entry <Termination criteria/Steady "
                                              "state heat flux>, there was an error. Specifically, "
                                              "the conversion function complained as follows: "
                                              + error));
            }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      AssertThrow (allowed_relative_deviation >= 0,
                   ExcMessage("Relative deviation must be greater than or equal to 0."));
      AssertThrow (necessary_time_in_steady_state > 0,
                   ExcMessage("Steady state minimum time period must be greater than 0."));
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace TerminationCriteria
  {
    ASPECT_REGISTER_TERMINATION_CRITERION(SteadyHeatFlux,
                                          "steady state heat flux",
                                          "A criterion that terminates the simulation when the integrated "
                                          "heat flux over a given list of boundaries stays within a certain "
                                          "range for a specified period of time."
                                          "\n\n"
                                          "The criterion considers the total heat flux over all boundaries "
                                          "listed by their boundary indicators, rather than each boundary "
                                          "separately. As a consequence, if the \\textit{sum} of heat fluxes "
                                          "over individual parts of the boundary no longer changes, then this "
                                          "criterion recommends termination, even if the heat flux over "
                                          "individual parts of the boundary continues to change.")
  }
}
