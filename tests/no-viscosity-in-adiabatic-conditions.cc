#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class ReallySimple : public MaterialModel::Simple<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;

      /**
        * Return true if the viscosity() function returns something that
        * may depend on the variable identifies by the argument.
        */
        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    double
    ReallySimple<dim>::
    viscosity (const double temperature,
               const double,
               const std::vector<double> &,
               const SymmetricTensor<2,dim> &strain_rate,
               const Point<dim> &) const
    {
      // the initial conditions for the viscosity are 1.0. If we find a zero,
      // then we are getting here erroneously. at the time of writing this test,
      // this happened when computing the adiabatic conditions, but we don't
      // need the viscosity there, so why compute it?
      Assert (temperature != 0, ExcMessage ("We can't evaluate the viscosity if the "
					    "temperature has not been initialized!"));
	
      return 1;
    }

    template <int dim>
    bool
    ReallySimple<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return ((dependence & NonlinearDependence::strain_rate) != NonlinearDependence::none);
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ReallySimple,
                                   "really simple",
                                   "A simple material model that is like the "
				   "'Simple' model, but aborts when the temperature has not been set.")
  }
}

