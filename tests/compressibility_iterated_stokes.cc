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

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

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
    void
    CompressibilityIteratedStokes<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      Simple<dim>::evaluate(in, out);
      for (unsigned int i=0; i < in.position.size(); ++i)
        {
          out.densities[i] = 10.0/11.0*exp(in.pressure[i]/100.0);
          out.compressibilities[i] = 0.01;
        }
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
