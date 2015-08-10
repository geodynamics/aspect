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


#include <aspect/postprocess/visualization/vertical_heat_flux.h>
#include <aspect/simulator_access.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      VerticalHeatFlux<dim>::
      VerticalHeatFlux ()
        :
        DataPostprocessorScalar<dim> ("vertical_heat_flux",
                                      update_values | update_q_points | update_gradients)
      {}



      template <int dim>
      void
      VerticalHeatFlux<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points,
                                                   this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        //Create vector for the temperature gradients.  All the other things
        //we need are in MaterialModelInputs/Outputs
        std::vector<Tensor<1,dim> > temperature_gradient(n_quadrature_points);

        in.position = evaluation_points;
        in.strain_rate.resize(0); // we do not need the viscosity
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            in.pressure[q]=uh[q][this->introspection().component_indices.pressure];
            in.temperature[q]=uh[q][this->introspection().component_indices.temperature];
            for (unsigned int d = 0; d < dim; ++d)
              {
                in.velocity[q][d]=uh[q][this->introspection().component_indices.velocities[d]];
                in.pressure_gradient[q][d] = duh[q][this->introspection().component_indices.pressure][d];
                temperature_gradient[q][d] = duh[q][this->introspection().component_indices.temperature][d];
              }

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
          }

        this->get_material_model().evaluate(in, out);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(in.position[q]);
            const Tensor<1,dim> vertical = -gravity/( gravity.norm() != 0.0 ?
                                                      gravity.norm() : 1.0 );
            const double advective_flux = (in.velocity[q] * vertical) * in.temperature[q] *
                                          out.densities[q]*out.specific_heat[q];
            const double conductive_flux = -(temperature_gradient[q]*vertical) *
                                           out.thermal_conductivities[q];
            computed_quantities[q](0) = advective_flux + conductive_flux;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(VerticalHeatFlux,
                                                  "vertical heat flux",
                                                  "A visualization output object that generates output "
                                                  "for the heat flux in the vertical direction, which is "
                                                  "the sum of the advective and the conductive heat flux, "
                                                  "with the sign convention of positive flux upwards.")
    }
  }
}
