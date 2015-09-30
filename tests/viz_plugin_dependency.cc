// create a visualization postprocessor that requires something else
//
// this works just like the plugin_dependency* tests, just that it's a
// viz postprocessor, not a regular postprocessor that has
// dependencies

#include "solcx.cc"
#include <aspect/postprocess/visualization.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      class MyPostprocessor
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          MyPostprocessor ()
            :
            DataPostprocessorScalar<dim> ("my_postprocessor",
                                          update_default)
          {}

          virtual
          void
          compute_derived_quantities_vector (const std::vector<Vector<double> > &,
                                             const std::vector<std::vector<Tensor<1,dim> > > &,
                                             const std::vector<std::vector<Tensor<2,dim> > > &,
                                             const std::vector<Point<dim> > &,
                                             const std::vector<Point<dim> > &,
                                             std::vector<Vector<double> >                    &computed_quantities) const
          {
            Assert (computed_quantities[0].size() == 1, ExcInternalError());

            for (unsigned int q=0; q<computed_quantities.size(); ++q)
              {
                // not important what we do here :-)
                computed_quantities[q](0) = 0;
              }
          }

          virtual
          std::list<std::string>
          required_other_postprocessors () const
          {
            // select a postprocessor that is not selected in the .prm file
            std::list<std::string> deps;
            deps.push_back ("velocity statistics");
            return deps;
          }
      };
    }
  }
}


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MyPostprocessor,
                                                  "my postprocessor",
                                                  ".")
    }
  }
}

