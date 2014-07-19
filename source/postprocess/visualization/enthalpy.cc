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


#include <aspect/postprocess/visualization/enthalpy.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/damage_rheology.h>

#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Enthalpy<dim>::
      Enthalpy ()
        :
        DataPostprocessorScalar<dim> ("enthalpy",
                                      update_values | update_q_points)
      {}



      template <int dim>
      void
      Enthalpy<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());

        const MaterialModel::DamageRheology<dim>* material_model =
            dynamic_cast<const MaterialModel::DamageRheology<dim> * > (&this->get_material_model());

        AssertThrow(&material_model != 0,
                    ExcMessage("The enthalpy postprocessor currently only works with the damage rheology "
                        "material model"));

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const double pressure    = uh[q][this->introspection().component_indices.pressure];
            const double temperature = uh[q][this->introspection().component_indices.temperature];
            std::vector<double> composition(this->n_compositional_fields());
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              composition[c] = uh[q][this->introspection().component_indices.compositional_fields[c]];

            computed_quantities[q](0) = material_model->enthalpy(temperature,
                                                                 pressure,
                                                                 composition,
                                                                 evaluation_points[q]);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Enthalpy,
                                                  "enthalpy",
                                                  "A visualization output object that generates output "
                                                  "for the enthalpy. Only works with DamageRheology "
                                                  "material model at the moment.")
    }
  }
}
