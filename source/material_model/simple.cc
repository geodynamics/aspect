//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/material_model/simple.h>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    double
    Simple<dim>::
    viscosity (const double,
               const double,
               const Point<dim> &) const
    {
      return reference_eta;
    }


    template <int dim>
    double
    Simple<dim>::
    reference_viscosity () const
    {
      return reference_eta;
    }



    template <int dim>
    double
    Simple<dim>::
    specific_heat (const double,
                   const double,
                   const Point<dim> &) const
    {
      return 1250.0;
    }



    template <int dim>
    double
    Simple<dim>::
    thermal_conductivity (const double,
                          const double,
                          const Point<dim> &) const
    {
      return 4.7;
    }



    template <int dim>
    double
    Simple<dim>::
    density (const double temperature,
             const double,
             const Point<dim> &) const
    {
      const double thermal_expansion_coefficient_ = 2e-5;
      return (reference_density *
              (1 - thermal_expansion_coefficient_ * (temperature -
                                                     reference_temperature)));
    }


    template <int dim>
    double
    Simple<dim>::
    compressibility (const double,
                     const double,
                     const Point<dim> &) const
    {
      return 0.0;
    }



    template <int dim>
    bool
    Simple<dim>::
    is_compressible () const
    {
      return false;
    }



    template <int dim>
    void
    Simple<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple model");
        {
          prm.declare_entry ("reference_density", "3300",
                             Patterns::Double (),
                             "rho0 in kg / m^3");
          prm.declare_entry ("reference_temperature", "293",
                             Patterns::Double (),
                             "T0 in K");
          prm.declare_entry ("reference_eta", "5e24",
                             Patterns::Double (),
                             "eta0");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Simple<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple model");
        {
          reference_density = prm.get_double ("reference_density");
          reference_temperature = prm.get_double ("reference_temperature");
          reference_eta = prm.get_double ("reference_eta");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    template class Simple<deal_II_dimension>;

    ASPECT_REGISTER_MATERIAL_MODEL(Simple,
                                   "simple",
                                   "A simple material model that has constant values "
                                   "throughout the domain. Additional parameters are "
                                   "read from the parameter file in subsection "
                                   "'Simple model'.")
  }
}
