// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/ODE_statistics.h>

#include <aspect/simulator_signals.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    void
    ODEStatistics<dim>::initialize()
    {
      this->get_signals().post_ARKode_solve.connect(
        [&](const SimulatorAccess<dim> &/*simulator_access*/,
            const unsigned int iteration_count)
      {
        this->store_ODE_solver_history(iteration_count);
      });

      this->clear_data();
    }



    template <int dim>
    void
    ODEStatistics<dim>::clear_data()
    {
      total_iteration_count = 0;
      number_of_solves = 0;
    }



    template <int dim>
    void
    ODEStatistics<dim>::store_ODE_solver_history(const unsigned int iteration_count)
    {
      total_iteration_count += iteration_count;
      number_of_solves += 1;
    }



    template <int dim>
    std::pair<std::string,std::string>
    ODEStatistics<dim>::execute (TableHandler &statistics)
    {
      const double average_iteration_count = number_of_solves > 0
                                             ?
                                             total_iteration_count / number_of_solves
                                             :
                                             total_iteration_count;
      statistics.add_value("Average iterations for ODE solver",
                           average_iteration_count);

      clear_data();

      return std::make_pair (std::string(),std::string());
    }



    template <int dim>
    template <class Archive>
    void ODEStatistics<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &total_iteration_count;
      ar &number_of_solves;
    }



    template <int dim>
    void ODEStatistics<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      // Serialize into a stringstream. Put the following into a code
      // block of its own to ensure the destruction of the 'oa'
      // archive triggers a flush() on the stringstream so we can
      // query the completed string below.
      std::ostringstream os;
      {
        aspect::oarchive oa (os);
        oa << (*this);
      }

      status_strings["ODEStatistics"] = os.str();
    }



    template <int dim>
    void ODEStatistics<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("ODEStatistics") != status_strings.end())
        {
          std::istringstream is (status_strings.find("ODEStatistics")->second);
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
    ASPECT_REGISTER_POSTPROCESSOR(ODEStatistics,
                                  "ODE statistics",
                                  "A postprocessor that computes some statistics about "
                                  "ODEs solved during the model evolution, specifically, "
                                  "how many iterations are needed to solve these ODEs "
                                  "on average. ")
  }
}
