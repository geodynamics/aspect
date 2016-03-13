#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>
#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>
#include <aspect/assembly.h>

#include <iostream>
#include <typeinfo>

unsigned int counter_without = 0, counter_with = 0;
bool quiet = true;


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;
    template <int dim>
    class AdditionalOutputs1 : public AdditionalMaterialOutputs<dim>
    {
      public:
        AdditionalOutputs1 (const unsigned int n_points)
        {
          additional_material_output1.resize(n_points);
        }

        std::vector<double> additional_material_output1;
    };


    template <int dim>
    class Material1 : public MaterialModel::Simple<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          MaterialModel::Simple<dim>::evaluate(in, out);

          AdditionalOutputs1<dim> *additional;

          additional = out.template get_additional_output<AdditionalOutputs1<dim> >();
          if (additional)
            ++counter_with;
          else
            ++counter_without;

          if (additional)
            additional->additional_material_output1[0] = 42.0;

          if (quiet)
            return;

          if (additional)
            std::cout << "* evaluate called with additional outputs!" << std::endl;
          else
            std::cout << "* evaluate called without additional outputs!" << std::endl;

        }


    };

  }


  template <int dim>
  class TestAssembler :
    public aspect::internal::Assembly::Assemblers::AssemblerBase<dim>
  {
    public:

      virtual void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &out)
      {
        std::cout << "* create_additional_material_model_outputs() called" << std::endl;

        if (out.template get_additional_output<MaterialModel::AdditionalOutputs1<dim> >() != NULL)
          return;

        std::cout << "   creating additional output!" << std::endl;
        out.additional_outputs.push_back(std::make_shared<MaterialModel::AdditionalOutputs1<dim> > (1));

      }

      virtual void assemble_stokes(internal::Assembly::Scratch::StokesSystem<dim>        &scratch,
                                   internal::Assembly::CopyData::StokesSystem<dim>       &copy)
      {
        MaterialModel::AdditionalOutputs1<dim> *additional
          = scratch.material_model_outputs.template get_additional_output<MaterialModel::AdditionalOutputs1<dim> >();

        std::cout << "* local_assemble_stokes call, have additional? " << (additional!=NULL) << std::endl;
        if (additional!=NULL)
          std::cout << "   value = " << additional->additional_material_output1[0] << std::endl;


      }


  };


  template <int dim>
  void set_assemblers1(const SimulatorAccess<dim> &,
                       internal::Assembly::AssemblerLists<dim> &assemblers,
                       std::vector<dealii::std_cxx11::shared_ptr<internal::Assembly::Assemblers::AssemblerBase<dim> > > &assembler_objects)
  {
    std::cout << "* set_assemblers()" << std::endl;

    std::cout << "called without: " << counter_without << " with: " << counter_with << std::endl;
    counter_without = 0;
    counter_with = 0;
    quiet = false;

    TestAssembler<dim> *test_assembler = new TestAssembler<dim>();
    assembler_objects.push_back(std_cxx11::shared_ptr<internal::Assembly::Assemblers::AssemblerBase<dim> >(test_assembler));

    assemblers.local_assemble_stokes_system
    .connect (std_cxx11::bind(&TestAssembler<dim>::assemble_stokes,
                              std_cxx11::ref (*test_assembler),
                              // discard cell,
                              // discard pressure_scaling,
                              // discard bool,

                              std_cxx11::_4,
                              std_cxx11::_5));


  }
}


template <int dim>
void signal_connector (aspect::SimulatorSignals<dim> &signals)
{
  std::cout << "* Connecting signals" << std::endl;
  signals.set_assemblers.connect (&aspect::set_assemblers1<dim>);

}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)


class f
{
  public:
    ~f()
    {
      std::cout << "called without: " << counter_without << " with: " << counter_with << std::endl;
    }
} instance;

namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Material1,
                                   "Material1",
                                   "")
  }
}
