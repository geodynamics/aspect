
#ifndef _aspect_material_model_thermal_conductivity_lookup_h
#define _aspect_material_model_thermal_conductivity_lookup_h

#include <aspect/material_model/thermal_conductivity/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      using namespace dealii;

      /**
       * A class that implements a pressure- and temperature-dependent thermal conductivity
       * based on density and specific heat capacity read from thermodynamic lookup tables
       * and user-defined thermal diffusivity.
       *
       * @ingroup MaterialModels
       */
      template <int dim>
      class Lookup : public Interface<dim>, public aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Function to compute the thermal conductivities in @p out given the
           * inputs in @p in.
           */
          void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                         MaterialModel::MaterialModelOutputs<dim> &out) const override;

          /**
           * Declare the parameters this plugin takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * Parameter to define thermal diffusivity.
           */
          std::vector<double> thermal_diffusivities;
      };
    }
  }
}

#endif