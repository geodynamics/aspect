// use the same postprocessing facilities as for the 'compressibility'
// testcase
#include "compressibility.cc"

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class CompressibilityIteratedStokes : public MaterialModel::Simple<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double density (const double                  temperature,
                                const double                  pressure,
                                const std::vector<double>    &compositional_fields,
                                const Point<dim>             &position) const;

        virtual double compressibility (const double                  temperature,
                                        const double                  pressure,
                                        const std::vector<double>    &compositional_fields,
                                        const Point<dim>             &position) const;

      /**
        * Return true if the compressibility() function returns something that
        * is not zero.
        */
        virtual bool
        is_compressible () const;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    double
    CompressibilityIteratedStokes<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &composition,
             const Point<dim> &position) const
    {
      return 10.0/11.0*exp(pressure/100.0);
    }

    template <int dim>
    double
    CompressibilityIteratedStokes<dim>::
    compressibility (const double ,
                     const double pressure,
                     const std::vector<double> &,
                     const Point<dim> &) const
    {
      // compressibility = 1/rho drho/dp
      return  0.01;
    }

    template <int dim>
    bool
    CompressibilityIteratedStokes<dim>::
    is_compressible () const
    {
      return true;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(CompressibilityIteratedStokes,
                                   "compressibility iterated stokes",
                                   "A simple material model that is like the "
				   "'Simple' model, but has a non-zero compressibility.")
  }
}
