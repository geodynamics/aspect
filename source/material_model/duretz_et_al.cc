/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/material_model/duretz_et_al.h>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    namespace DuretzEtAl
    {
      // ------------------ implementation of the SolCx benchmark ----------------------------

      template <int dim>
      double
      SolCx<dim>::
      viscosity (const double,
                 const double,
                 const SymmetricTensor<2,dim> &,
                 const Point<dim> &p) const
      {
        // defined as given in the Duretz et al. paper
        return (p[0] < 0.5 ? 1 : eta_B);
      }


      template <int dim>
      double
      SolCx<dim>::
      reference_viscosity () const
      {
        return 1;
      }

      template <int dim>
      double
      SolCx<dim>::
      reference_density () const
      {
        return 1;
      }

      template <int dim>
      double
      SolCx<dim>::
      reference_thermal_expansion_coefficient () const
      {
        return 0;
      }

      template <int dim>
      double
      SolCx<dim>::
      specific_heat (const double,
                     const double,
                     const Point<dim> &) const
      {
        return 0;
      }

      template <int dim>
      double
      SolCx<dim>::
      reference_cp () const
      {
        return 0;
      }

      template <int dim>
      double
      SolCx<dim>::
      thermal_conductivity (const double,
                            const double,
                            const Point<dim> &) const
      {
        return 0;
      }

      template <int dim>
      double
      SolCx<dim>::
      reference_thermal_diffusivity () const
      {
        return 0;
      }

      template <int dim>
      double
      SolCx<dim>::
      density (const double,
               const double,
               const Point<dim> &p) const
      {
        // defined as given in the paper
        return -std::sin(numbers::PI*p[1])*std::cos(numbers::PI*p[0]);
      }


      template <int dim>
      double
      SolCx<dim>::
      thermal_expansion_coefficient (const double temperature,
                                     const double,
                                     const Point<dim> &) const
      {
        return 0;
      }


      template <int dim>
      double
      SolCx<dim>::
      compressibility (const double,
                       const double,
                       const Point<dim> &) const
      {
        return 0.0;
      }



      template <int dim>
      bool
      SolCx<dim>::
      viscosity_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }


      template <int dim>
      bool
      SolCx<dim>::
      density_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }

      template <int dim>
      bool
      SolCx<dim>::
      compressibility_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }

      template <int dim>
      bool
      SolCx<dim>::
      specific_heat_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }

      template <int dim>
      bool
      SolCx<dim>::
      thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      template <int dim>
      bool
      SolCx<dim>::
      is_compressible () const
      {
        return false;
      }

      template <int dim>
      void
      SolCx<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("SolCx");
          {
            prm.declare_entry ("Viscosity jump", "1e6",
                               Patterns::Double (0),
                               "Viscosity in the right half of the domain.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      SolCx<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("SolCx");
          {
            eta_B = prm.get_double ("Viscosity jump");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      double
      SolCx<dim>::get_eta_B() const
      {
        return eta_B;
      }

// ------------------ implementation of the SolKz benchmark ----------------------------


      template <int dim>
      double
      SolKz<dim>::
      viscosity (const double,
                 const double,
                 const SymmetricTensor<2,dim> &,
                 const Point<dim> &p) const
      {
        // defined as given in the Duretz et al. paper
        static const double B = 0.5 * std::log(1e6);
        return std::exp(2*B*p[1]);
      }


      template <int dim>
      double
      SolKz<dim>::
      reference_viscosity () const
      {
        return 1;
      }

      template <int dim>
      double
      SolKz<dim>::
      reference_density () const
      {
        return 1;
      }

      template <int dim>
      double
      SolKz<dim>::
      reference_thermal_expansion_coefficient () const
      {
        return 0;
      }

      template <int dim>
      double
      SolKz<dim>::
      specific_heat (const double,
                     const double,
                     const Point<dim> &) const
      {
        return 0;
      }

      template <int dim>
      double
      SolKz<dim>::
      reference_cp () const
      {
        return 0;
      }

      template <int dim>
      double
      SolKz<dim>::
      thermal_conductivity (const double,
                            const double,
                            const Point<dim> &) const
      {
        return 0;
      }

      template <int dim>
      double
      SolKz<dim>::
      reference_thermal_diffusivity () const
      {
        return 0;
      }

      template <int dim>
      double
      SolKz<dim>::
      density (const double,
               const double,
               const Point<dim> &p) const
      {
        // defined as given in the paper
        return -std::sin(2*p[1])*std::cos(3*numbers::PI*p[0]);
      }


      template <int dim>
      double
      SolKz<dim>::
      thermal_expansion_coefficient (const double temperature,
                                     const double,
                                     const Point<dim> &) const
      {
        return 0;
      }


      template <int dim>
      double
      SolKz<dim>::
      compressibility (const double,
                       const double,
                       const Point<dim> &) const
      {
        return 0.0;
      }



      template <int dim>
      bool
      SolKz<dim>::
      viscosity_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }


      template <int dim>
      bool
      SolKz<dim>::
      density_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }

      template <int dim>
      bool
      SolKz<dim>::
      compressibility_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }

      template <int dim>
      bool
      SolKz<dim>::
      specific_heat_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }

      template <int dim>
      bool
      SolKz<dim>::
      thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      template <int dim>
      bool
      SolKz<dim>::
      is_compressible () const
      {
        return false;
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace DuretzEtAl
    {
      ASPECT_REGISTER_MATERIAL_MODEL(SolCx,
                                     "SolCx",
                                     "A material model that corresponds to the 'SolCx' benchmark "
                                     "defined in Duretz et al., G-Cubed, 2011.")
      ASPECT_REGISTER_MATERIAL_MODEL(SolKz,
                                     "SolKz",
                                     "A material model that corresponds to the 'SolKz' benchmark "
                                     "defined in Duretz et al., G-Cubed, 2011.")
    }
  }
}
