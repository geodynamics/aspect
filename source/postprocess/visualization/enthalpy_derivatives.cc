/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/visualization/enthalpy_derivatives.h>
#include <aspect/simulator_access.h>

#include <aspect/material_model/damage_rheology.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/grid/grid_tools.h>
#include <algorithm>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      EnthalpyDerivatives<dim>::
      EnthalpyDerivatives ()
        :
        DataPostprocessor<dim> ()
      {}

      template <int dim>
      std::vector<std::string>
      EnthalpyDerivatives<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;

        solution_names.push_back ("enthalpy_derivative_temperature");
        solution_names.push_back ("enthalpy_derivative_pressure");
        solution_names.push_back ("enthalpy_points_temperature");
        solution_names.push_back ("enthalpy_points_pressure");

        return solution_names;
      }

      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      EnthalpyDerivatives<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;

        interpretation.push_back (DataComponentInterpretation::component_is_scalar);
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);

        return interpretation;
      }

      template <int dim>
      UpdateFlags
      EnthalpyDerivatives<dim>::
      get_needed_update_flags () const
      {
        return update_gradients | update_values  | update_q_points;
      }

      template <int dim>
      void
      EnthalpyDerivatives<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                         const std::vector<Point<dim> >                  &normals,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());

        typename MaterialModel::Interface<dim>::MaterialModelInputs in(n_quadrature_points,
                                                                       this->n_compositional_fields());
        typename MaterialModel::Interface<dim>::MaterialModelOutputs out(n_quadrature_points,
                                                                         this->n_compositional_fields());

        in.position = evaluation_points;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = duh[q][d];
            in.strain_rate[q] = symmetrize (grad_u);

            in.pressure[q]=uh[q][this->introspection().component_indices.pressure];
            in.temperature[q]=uh[q][this->introspection().component_indices.temperature];
            for (unsigned int i = 0; i < dim; ++i)
              in.velocity[q][i]=uh[q][this->introspection().component_indices.velocities[i]];

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
          }

        Point<dim> average_position;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            average_position += in.position[q];
          }
        average_position /= n_quadrature_points;
        typename DoFHandler<dim>::active_cell_iterator cell = GridTools::find_active_cell_around_point(this->get_dof_handler(),average_position);

        in.cell = &cell;

        const MaterialModel::DamageRheology<dim>* material_model =
            dynamic_cast<const MaterialModel::DamageRheology<dim> * > (&this->get_material_model());

        AssertThrow(&material_model != 0,
                    ExcMessage("The enthalpy postprocessor currently only works with the damage rheology "
                        "material model"));

        const std_cxx1x::array<std::pair<double, unsigned int>,2> dH =
            material_model->enthalpy_derivative(in);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            computed_quantities[q][0] = dH[0].first;
            computed_quantities[q][1] = dH[1].first;
            computed_quantities[q][2] = dH[0].second;
            computed_quantities[q][3] = dH[1].second;
          }
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(EnthalpyDerivatives,
                                                  "enthalpy derivatives",
                                                  "A visualization output object that generates output "
                                                  "for the material properties given by the material model."
                                                  "There are a number of other visualization postprocessors "
                                                  "that offer to write individual material properties. However, "
                                                  "they all individually have to evaluate the material model. "
                                                  "This is inefficient if one wants to output more than just "
                                                  "one or two of the fields provided by the material model. "
                                                  "The current postprocessor allows to output a (potentially "
                                                  "large) subsets of all of the information provided by "
                                                  "material models at once, with just a single material model "
                                                  "evaluation per output point.")
    }
  }
}
