#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>
#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>
#include <aspect/assembly.h>

#include <deal.II/fe/fe_dgq.h>
#include <iostream>

using namespace dealii;

namespace aspect
{

  template <int dim>
  void my_signal(std::vector<VariableDeclaration<dim> > &variables)
  {
    std::cout << "* signals.edit_finite_element_variables:" << std::endl;

    VariableDeclaration<dim> dummy("dummy",
                                   std_cxx11::shared_ptr<FiniteElement<dim> > (
                                     new dealii::FE_DGQ<dim>(4)),
                                   1,
                                   1);
    variables.insert(variables.begin()+6, dummy);

    for (unsigned int i=0; i<variables.size(); ++i)
      {
        std::cout << " name=" << variables[i].name
                  << " fe=" << variables[i].fe->get_name()
                  << " multiplicity=" << variables[i].multiplicity
                  << " n_blocks=" << variables[i].n_blocks
                  << " n_components=" << variables[i].n_components()
                  << std::endl;
      }

    std::cout << std::endl;
  }


  template <int dim>
  class MeltMaterial:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_viscosity () const
      {
        return 1.0;
      }

      virtual double reference_darcy_coefficient () const
      {
        const double porosity = 0.01;
        const double permeability = 1.0 * std::pow(porosity, 3) * std::pow(1.0-porosity, 2);
        return permeability / 0.1;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {

        for (unsigned int i=0; i<in.position.size(); ++i)
          {
            out.viscosities[i] = 1.0;
            out.densities[i] = 1.0 + (in.composition[i][0]);
            out.thermal_expansion_coefficients[i] = 1.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
            for (unsigned int c=0; c<in.composition[i].size(); ++c)
              out.reaction_terms[i][c] = 0.0;
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim> >();

        if (melt_out != NULL)
          {
            const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
            for (unsigned int i=0; i<in.position.size(); ++i)
              {
                const double porosity = in.composition[i][porosity_idx];

                melt_out->compaction_viscosities[i] = 100.0;
                melt_out->fluid_viscosities[i]=0.1;
                melt_out->permeabilities[i]=1.0 * std::pow(porosity,3) * std::pow(1.0-porosity,2);
                melt_out->fluid_densities[i]=.1;
                melt_out->fluid_density_gradients[i] = Tensor<1,dim>();
              }
          }
      }

  };


}




template <int dim>
void signal_connector (aspect::SimulatorSignals<dim> &signals)
{
  std::cout << "* Connecting signals" << std::endl;
  signals.edit_finite_element_variables.connect(&aspect::my_signal<dim>);
}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)


// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(MeltMaterial,
                                 "MeltMaterial",
                                 "")
}
