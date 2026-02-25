/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/postprocess/ODE_statistics.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/simulator.h>



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
#if DEAL_II_VERSION_GTE(9,6,0)
      const double average_iteration_count = number_of_solves > 0
                                             ?
                                             total_iteration_count / number_of_solves
                                             :
                                             total_iteration_count;
#else
      // The computation is not correct for dealii versions older than
      // June 2024.
      const double average_iteration_count = std::numeric_limits<double>::quiet_NaN();
#endif
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
