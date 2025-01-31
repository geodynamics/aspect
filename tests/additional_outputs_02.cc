/*
  Copyright (C) 2022 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>
#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator/assemblers/interface.h>

#include <iostream>

unsigned int counter_without = 0, counter_with = 0;
bool quiet = true;


namespace aspect
{
  namespace MaterialModel
  {
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

          additional = out.template get_additional_output<AdditionalOutputs1<dim>>();
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
    public aspect::Assemblers::Interface<dim>
  {
    public:

      virtual void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        std::cout << "* create_additional_material_model_outputs() called" << std::endl;

        if (out.template get_additional_output<MaterialModel::AdditionalOutputs1<dim>>() != nullptr)
          return;

        std::cout << "   creating additional output!" << std::endl;
        out.additional_outputs.push_back(std::make_unique<MaterialModel::AdditionalOutputs1<dim>> (1));

      }

      virtual void execute(internal::Assembly::Scratch::ScratchBase<dim>        &scratch_base,
                           internal::Assembly::CopyData::CopyDataBase<dim>       &/*data_base*/) const
      {
        internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);

        MaterialModel::AdditionalOutputs1<dim> *additional
          = scratch.material_model_outputs.template get_additional_output<MaterialModel::AdditionalOutputs1<dim>>();

        std::cout << "* local_assemble_stokes call, have additional? " << (additional!=nullptr) << std::endl;
        if (additional!=nullptr)
          std::cout << "   value = " << additional->additional_material_output1[0] << std::endl;


      }


  };


  template <int dim>
  void set_assemblers1(const SimulatorAccess<dim> &,
                       Assemblers::Manager<dim> &assemblers)
  {
    std::cout << "* set_assemblers()" << std::endl;

    std::cout << "called without: " << counter_without << " with: " << counter_with << std::endl;
    counter_without = 0;
    counter_with = 0;
    quiet = false;

    TestAssembler<dim> *test_assembler = new TestAssembler<dim>();
    assemblers.stokes_system.push_back(std::unique_ptr<Assemblers::Interface<dim>>(test_assembler));
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
