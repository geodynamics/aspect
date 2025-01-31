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

#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class ShearThinning : public MaterialModel::Simple<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        void
        parse_parameters (ParameterHandler &prm);
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    ShearThinning<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      Simple<dim>::evaluate(in, out);
      if (in.requests_property(MaterialProperties::viscosity))
        for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
          out.viscosities[i] = 1./(1+in.strain_rate[i].norm());
    }

    template <int dim>
    void
    ShearThinning<dim>::parse_parameters (ParameterHandler &prm)
    {
      Simple<dim>::parse_parameters(prm);

      // Declare dependencies on solution variables that are different from simple
      this->model_dependence.viscosity = MaterialModel::NonlinearDependence::strain_rate;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ShearThinning,
                                   "shear thinning",
                                   "A simple material model that is like the "
                                   "'Simple' model, but has a viscosity equal to 1/|strain rate|.")
  }
}



#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    class ShearThinning : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &);
    };
  }
}


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ShearThinning<dim>::execute (TableHandler &)
    {
      // create a quadrature formula based on the temperature element alone.
      // be defensive about determining that what we think is the temperature
      // element is indeed the temperature element
      Assert (this->get_fe().n_base_elements() == 3+(this->n_compositional_fields()>0 ? 1 : 0),
              ExcNotImplemented());
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(2).degree+1);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_gradients         | update_values |
                               update_quadrature_points | update_JxW_values);

      std::vector<std::vector<double>> composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      typename MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());
      in.requested_properties = MaterialModel::MaterialProperties::viscosity;

      // compute the integral of the viscosity. since we're on a unit box,
      // this also is the average value
      double viscosity_integral = 0;
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                         in.temperature);
            fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                      in.pressure);
            fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(),
                in.strain_rate);
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                  composition_values[c]);

            in.position = fe_values.get_quadrature_points();
            for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
              {
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  in.composition[i][c] = composition_values[c][i];
              }

            this->get_material_model().evaluate(in, out);


            for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
              viscosity_integral += out.viscosities[q] * fe_values.JxW(q);
          }

      std::ostringstream screen_text;
      screen_text.precision(4);
      screen_text << Utilities::MPI::sum(viscosity_integral, this->get_mpi_communicator());

      return std::pair<std::string, std::string> ("Average viscosity:",
                                                  screen_text.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ShearThinning,
                                  "shear thinning",
                                  "A postprocessor that computes some statistics about "
                                  "the viscosity.")
  }
}
