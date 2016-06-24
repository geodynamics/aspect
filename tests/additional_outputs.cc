// make sure we use the additional material model outputs
#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>
#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

#include <iostream>
#include <typeinfo>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;
    template <int dim>
    class AdditionalOutputs1 : public AdditionalMaterialOutputs<dim>
    {
      public:
        AdditionalOutputs1 (const unsigned int n_points,
                            const unsigned int n_comp)
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
          AdditionalOutputs1<dim> *additional;

          additional = out.template get_additional_output<AdditionalOutputs1<dim> >();
          additional->additional_material_output1[0] = 42.0;

        }


    };

  }
}


int f()
{
  const int dim=2;

  using namespace aspect::MaterialModel;
  MaterialModelInputs<dim> in(1,1);
  MaterialModelOutputs<dim> out(1,1);


  if (out.get_additional_output<AdditionalOutputs1<dim> >() != NULL)
    throw "error";

  out.additional_outputs.push_back(std::make_shared<AdditionalOutputs1<dim> > (1, 1));

  struct empty {};

  if (out.get_additional_output<empty>() != NULL)
    throw "error";

  Material1<dim> mat;
  mat.evaluate(in, out);

  std::cout << out.get_additional_output<AdditionalOutputs1<dim> >()->additional_material_output1[0] << std::endl;

  // test const version of get_additional_output:
  {
    const MaterialModelOutputs<dim> &const_out = out;
    if (const_out.get_additional_output<empty>() != NULL)
      throw "error";
    const AdditionalOutputs1<dim> *a = const_out.get_additional_output<AdditionalOutputs1<dim> >();
    if (a == NULL)
      throw "error";
  }

  std::cout << "OK" << std::endl;
  exit(0);
  return 42;
}

// run this function by initializing a global variable by it
int i = f();
