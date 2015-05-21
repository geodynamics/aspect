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



#include <aspect/postprocess/visualization/internal_heating.h>
#include <aspect/simulator_access.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      InternalHeating<dim>::
      InternalHeating ()
        :
        DataPostprocessorScalar<dim> ("Internal heating rate",
                                      update_values | update_gradients | update_q_points)
      {}



      template <int dim>
      void
      InternalHeating<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const HeatingModel::Interface<dim> &heating_model=this->get_heating_model();
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components, ExcInternalError());

        typename MaterialModel::Interface<dim>::MaterialModelInputs in(n_quadrature_points, this->n_compositional_fields());
        typename MaterialModel::Interface<dim>::MaterialModelOutputs out(n_quadrature_points, this->n_compositional_fields());

        in.strain_rate.resize(0); // we do not need the viscosity
        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (n_quadrature_points));
        std::vector<double> heating_model_outputs(n_quadrature_points);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            in.temperature[q] = uh[q][this->introspection().component_indices.temperature];
            in.pressure[q]    = uh[q][this->introspection().component_indices.pressure];

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
          }

        in.position = evaluation_points;

        this->get_material_model().evaluate(in, out);
        heating_model.evaluate(in, out, heating_model_outputs);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          computed_quantities[q](0) = heating_model_outputs[q];
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(InternalHeating,
                                                  "internal heating",
                                                  "A visualization output object that generates output "
                                                  "for the internal heating rate. Units: W/kg")
    }
  }
}
