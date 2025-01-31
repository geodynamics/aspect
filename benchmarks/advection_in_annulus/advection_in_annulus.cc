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
#include <aspect/simulator.h>
#include <aspect/material_model/simple.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/postprocess/dynamic_topography.h>

#include <aspect/prescribed_stokes_solution/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  /**
   * This is a variation of the "annulus" benchmark in benchmark/annulus
   * that tests the influence of the stabilization method on the temperature
   * equation.
   *
   */
  namespace AdvectionInAnnulus
  {
    namespace AnalyticSolutions
    {
      const double A=2.0, B=-3.0/std::log(2.0), C=-1;
      const double outer_radius = 2.;
      const double rho_0 = 1000.;
      const double gravity = 1.;

      Tensor<1,2>
      Annulus_velocity (const Point<2> &pos,
                        const double k)
      {
        const double x = pos[0];
        const double y = pos[1];
        const double r = std::sqrt(x*x+y*y);
        const double theta = std::atan2(y,x);
        const double f_r = A*r + B/r;
        const double g_r = A*r/2 + B*std::log(r)/r + C/r;
        const double v_r = g_r*k*std::sin(k*theta);
        const double v_theta = f_r*std::cos(k*theta);
        const double v_x = std::cos(theta)*v_r - std::sin(theta)*v_theta;
        const double v_y = std::sin(theta)*v_r + std::cos(theta)*v_theta;
        return Point<2> (v_x,v_y);
      }

      double
      Annulus_pressure (const Point<2> &pos,
                        const double k)
      {
        const double x = pos[0];
        const double y = pos[1];
        const double r=std::sqrt(x*x+y*y);
        const double theta = std::atan2(y,x);
        const double f_r = 2*r + B/r;
        const double g_r = A*r/2 + B*std::log(r)/r + C/r;
        const double h_r=(2*g_r-f_r)/r;
        return k*h_r*std::sin(k*theta)+rho_0*gravity*(outer_radius-r);
      }

      template <int dim>
      double
      Annulus_normal_traction (const Point<dim> &pos,
                               const double k)
      {
        Assert (dim == 2, ExcNotImplemented());
        const double x = pos[0];
        const double y = pos[1];
        const double r=std::sqrt(x*x+y*y);
        const double theta = std::atan2(y,x);
        const double f_r = 2*r + B/r;
        const double g_r = A*r/2 + B*std::log(r)/r + C/r;
        return k * 3.*f_r/r * std::sin(k*theta) - rho_0 * g_r * (outer_radius - r);
      }

    }
  }
}


namespace aspect
{
  namespace PrescribedStokesSolution
  {
    /**
     * A class that implements the flow field of the annulus benchmark.
     *
     * @ingroup PrescribedStokesSolution
     */
    template <int dim>
    class AdvectionInAnnulus : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void stokes_solution (const Point<dim> &p, Vector<double> &value) const override;
    };
  }
}

namespace aspect
{
  namespace PrescribedStokesSolution
  {
    template <int dim>
    void AdvectionInAnnulus<dim>::stokes_solution (const Point<dim> &, Vector<double> &) const
    {
      return;
    }


    template <>
    void AdvectionInAnnulus<2>::stokes_solution (const Point<2> &p, Vector<double> &value) const
    {
      Tensor<1,2> velocity = aspect::AdvectionInAnnulus::AnalyticSolutions::Annulus_velocity (p, 4);
      value(0) = velocity[0];
      value(1) = velocity[1];

      double pressure = aspect::AdvectionInAnnulus::AnalyticSolutions::Annulus_pressure (p, 4);
      value(2) = pressure;       // pressure
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace PrescribedStokesSolution
  {
    ASPECT_REGISTER_PRESCRIBED_STOKES_SOLUTION(AdvectionInAnnulus,
                                               "advection in annulus",
                                               "This plugin prescribes the Stokes solution of "
                                               "the annulus benchmark.")
  }
}
