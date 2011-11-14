//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/model_simple.h>

using namespace dealii;

namespace aspect
{

  template <int dim>
  MaterialModel_Simple<dim>::~MaterialModel_Simple ()
  {}

  template <int dim>
  double
  MaterialModel_Simple<dim>::
  eta (const double temperature, const double pressure, const dealii::Point<dim> &position) const
  {
    return reference_eta;
  }

  template <int dim>
  double
  MaterialModel_Simple<dim>::
  real_viscosity (const double                 temperature,
                  const double                  pressure,
                  const dealii::Point<dim> &position,
                  const dealii::SymmetricTensor<2,dim> &strain_rate) const
  {
    // this is currently only used
    // in generating graphical
    // output
    return eta (temperature, pressure, position);

// we could use something more complicated, like the following:

    //                                       // This is the function used
    //                                       // for the lower mantle in the
    //                                       // Science paper by Stadler et
    //                                       // al., see
    //       // http://www.sciencemag.org/content/329/5995/1033/suppl/\DC1
    // const double d    = 9.3e4;
    // const double p    = 3;
    // const double A    = 1;
    // const double C_OH = 1e3;         /* ppm */
    // const double r    = 1;
    // const double E_a  = 335e3;       /* J/mol */
    // const double V_a  = 1.25e-6;     /* m^3/mol */
    // const double R    = 8.31447215;  /* J/mol/K */
    // const double n    = 1;

    // return (std::pow(std::pow(d,p)/(A*std::pow(C_OH,r)), 1./n)
    //              *
    //              std::pow(second_invariant(strain_rate), (1.-n)/n)
    //              *
    //              std::exp((E_a+pressure*V_a)/(n*R*temperature)));
  }

  template <int dim>
  double
  MaterialModel_Simple<dim>::
  specific_heat (const double temperature,
                 const double pressure) const
  {
    return 1250.0;

  }

  template <int dim>
  double
  MaterialModel_Simple<dim>::
  density (const double temperature,
           const double pressure,
           const dealii::Point<dim> &position) const
  {
    const double expansion_coefficient_ = expansion_coefficient (temperature, pressure, position);
    return (reference_density *
            (1 - expansion_coefficient_ * (temperature -
                                           reference_temperature)));
  }

  template <int dim>
  double
  MaterialModel_Simple<dim>::
  compressibility (const double temperature,
                   const double pressure,
                   const dealii::Point<dim> &position) const
  {
    return 0.0;
  }

  template <int dim>
  double
  MaterialModel_Simple<dim>::
  expansion_coefficient (const double temperature,
                         const double pressure,
                         const dealii::Point<dim> &position) const
  {
    return 2e-5;
  }

  template <int dim>
  void
  MaterialModel_Simple<dim>::declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection("ModelParameters");
    {
      prm.enter_subsection("MaterialModel_Simple");
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
  MaterialModel_Simple<dim>::parse_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection("ModelParameters");
    {
      prm.enter_subsection("MaterialModel_Simple");
      {
        reference_density = prm.get_double ("Time between graphical output");
        reference_temperature = prm.get_double ("Time between graphical output");
        reference_eta = prm.get_double ("Time between graphical output");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

}

// explicit instantiations
namespace aspect
{

  template class MaterialModel_Simple<deal_II_dimension>;

}
