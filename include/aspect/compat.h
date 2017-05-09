/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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

#ifndef _aspect_compat_h
#define _aspect_compat_h

#include <aspect/global.h>

// C++11 related includes. Can be removed when we require C++11.
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/std_cxx11/bind.h>
#include <deal.II/base/std_cxx11/function.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/base/std_cxx11/unique_ptr.h>


// We would like to use a function from SolverControl that was introduced after
// deal.II 8.5. For older versions use this derived class instead that implements
// the function and the typedef guarantees it is found before deal.II's class.
#if !DEAL_II_VERSION_GTE(9,0,0)

#include <deal.II/lac/solver_control.h>

namespace aspect
{
  using namespace dealii;

  class SolverControl : public dealii::SolverControl
  {
    public:
      SolverControl(const unsigned int n           = 100,
                    const double       tol         = 1.e-10,
                    const bool         log_history = false,
                    const bool         log_result  = true)
        :
        dealii::SolverControl (n, tol, log_history, log_result)
      {}

      dealii::SolverControl::State
      check (const unsigned int step,
             const double check_value)
      {
        dealii::SolverControl::State return_value = dealii::SolverControl::check(step, check_value);

        if (step == 0)
          history_data.resize(history_data.size()+1);
        return return_value;
      }


      const std::vector<double> &get_history_data() const
      {
        Assert (history_data_enabled, ExcHistoryDataRequired());
        Assert (history_data.size() > 0,
                ExcMessage("The SolverControl object was asked for the solver history "
                           "data, but there is no data. Possibly you requested the data before the "
                           "solver was run."));

        return history_data;
      }
  };
}
#endif


#endif
