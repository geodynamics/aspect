// use the same postprocessing facilities as for the
// 'compressibility_iterated_stokes' testcase
#include "compressibility_iterated_stokes.cc"

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class NoAdiabaticHeating : public CompressibilityIteratedStokes<dim>
    {
      public:
        virtual double thermal_expansion_coefficient (const double,
                                                      const double,
                                                      const std::vector<double> &,
                                                      const Point<dim> &) const
        {
          return 0;
        }

        virtual double viscosity (const double,
                                  const double,
                                  const std::vector<double> &,
                                  const Point<dim> &) const
        {
          return 0;
        }
    };

  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(NoAdiabaticHeating,
                                   "no adiabatic heating",
                                   "As described in the .prm file.")
  }
}
