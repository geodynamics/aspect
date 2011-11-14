//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/model_simple.h>

namespace aspect
{
  namespace
  {
    const double reference_density     = 3300;    /* kg / m^3   */
    const double reference_temperature = 293;     /* K          */
    const double reference_eta    = 5e24;
  }

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
}

// explicit instantiations
namespace aspect
{

  template class MaterialModel_Simple<deal_II_dimension>;

}
