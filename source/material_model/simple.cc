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
               const SymmetricTensor<2,dim> &,
               const Point<dim> &) const
    {
      return eta;
    }


    template <int dim>
    double
    Simple<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    double
    Simple<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    Simple<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_alpha;
    }

    template <int dim>
    double
    Simple<dim>::
    specific_heat (const double,
                   const double,
                   const Point<dim> &) const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    Simple<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    Simple<dim>::
    thermal_conductivity (const double,
                          const double,
                          const Point<dim> &) const
    {
      return k_value;
    }

    template <int dim>
    double
    Simple<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }

    template <int dim>
    double
    Simple<dim>::
    density (const double temperature,
             const double,
             const Point<dim> &) const
    {
      return (reference_rho *
              (1 - thermal_alpha * (temperature - reference_T)));
    }


    template <int dim>
    double
    Simple<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const Point<dim> &) const
    {
      return thermal_alpha;
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
    viscosity_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    Simple<dim>::
    density_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Simple<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Simple<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Simple<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
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
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the constant viscosity. Units: $kg/m/s$.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $JG/kgK$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
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
          reference_rho     = prm.get_double ("Reference density");
          reference_T = prm.get_double ("Reference temperature");
          eta                   = prm.get_double ("Viscosity");
          k_value               = prm.get_double ("Thermal conductivity");
          reference_specific_heat = prm.get_double ("Reference specific heat");
          thermal_alpha = prm.get_double ("Thermal expansion coefficient");
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
                                   "for all coefficients but the density. This model uses "
                                   "the formulation that assumes an incompressible medium "
                                   "despite the fact that the density follows the law "
                                   "$\\rho(T)=\\rho_0(1-\\beta(T-T_{\\text{ref}})$. The value for "
                                   "the components of this formula and additional "
                                   "parameters are read from the parameter file in subsection "
                                   "'Simple model'.")
  }
}
