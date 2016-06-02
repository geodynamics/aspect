/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/melt.h>
#include <aspect/melt.h>
#include <aspect/material_model/interface.h>

#include <deal.II/numerics/data_out.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {

      template <int dim>
      MeltDensity<dim>::
      MeltDensity ()
        :
        DataPostprocessorScalar<dim> ("melt_density",
                                      update_values | update_q_points)
      {}

      template <int dim>
      void
      MeltDensity<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &/*duh*/,
                                         const std::vector<std::vector<Tensor<2,dim> > > &/*dduh*/,
                                         const std::vector<Point<dim> >                  &/*normals*/,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        AssertThrow(this->include_melt_transport()==true, ExcMessage("include_melt_transport has to be on when using melt transport postprocessors."));

        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points, this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points, this->n_compositional_fields());
        create_melt_material_outputs(out);

        in.position = evaluation_points;
        in.strain_rate.resize(0); // we do not need the viscosity
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            in.pressure[q]=uh[q][this->introspection().component_indices.pressure];
            in.temperature[q]=uh[q][this->introspection().component_indices.temperature];

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
          }

        this->get_material_model().evaluate(in, out);
        MaterialModel::MeltOutputs<dim> *melt_outputs = out.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
        AssertThrow(melt_outputs != NULL, ExcMessage("Need MeltOutputs from the material model for computing the melt density."));

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          computed_quantities[q](0) = melt_outputs->fluid_densities[q];
      }

      template <int dim>
      CompactionViscosity<dim>::
      CompactionViscosity ()
        :
        DataPostprocessorScalar<dim> ("compaction_viscosity",
                                      update_values | update_gradients | update_q_points)
      {}

      template <int dim>
      void
      CompactionViscosity<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &/*dduh*/,
                                         const std::vector<Point<dim> >                  &/*normals*/,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        AssertThrow(this->include_melt_transport()==true, ExcMessage("include_melt_transport has to be on when using melt transport postprocessors."));

        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());
        Assert (duh[0].size() == this->introspection().n_components,          ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points, this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points, this->n_compositional_fields());
        create_melt_material_outputs(out);

        in.position = evaluation_points;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = duh[q][d];
            in.strain_rate[q] = symmetrize (grad_u);

            in.pressure[q]=uh[q][this->introspection().component_indices.pressure];
            in.temperature[q]=uh[q][this->introspection().component_indices.temperature];

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
          }

        this->get_material_model().evaluate(in, out);
        MaterialModel::MeltOutputs<dim> *melt_outputs = out.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
        AssertThrow(melt_outputs != NULL, ExcMessage("Need MeltOutputs from the material model for computing the compaction viscosity."));

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          computed_quantities[q](0) = melt_outputs->compaction_viscosities[q];
      }

      template <int dim>
      Permeability<dim>::
      Permeability ()
        :
        DataPostprocessorScalar<dim> ("permeability",
                                      update_values | update_q_points)
      {}

      template <int dim>
      void
      Permeability<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &/*duh*/,
                                         const std::vector<std::vector<Tensor<2,dim> > > &/*dduh*/,
                                         const std::vector<Point<dim> >                  &/*normals*/,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        AssertThrow(this->include_melt_transport()==true, ExcMessage("include_melt_transport has to be on when using melt transport postprocessors."));

        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points, this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points, this->n_compositional_fields());
        create_melt_material_outputs(out);

        in.position = evaluation_points;
        in.strain_rate.resize(0); // we do not need the viscosity
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            in.pressure[q]=uh[q][this->introspection().component_indices.pressure];
            in.temperature[q]=uh[q][this->introspection().component_indices.temperature];

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
          }

        this->get_material_model().evaluate(in, out);
        MaterialModel::MeltOutputs<dim> *melt_outputs = out.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
        AssertThrow(melt_outputs != NULL, ExcMessage("Need MeltOutputs from the material model for computing the permeabilities."));

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          computed_quantities[q](0) = melt_outputs->permeabilities[q];
      }

      template <int dim>
      MeltViscosity<dim>::
      MeltViscosity ()
        :
        DataPostprocessorScalar<dim> ("melt_viscosity",
                                      update_values | update_gradients | update_q_points)
      {}

      template <int dim>
      void
      MeltViscosity<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > & /*dduh*/,
                                         const std::vector<Point<dim> >                  &/*normals*/,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        AssertThrow(this->include_melt_transport()==true, ExcMessage("include_melt_transport has to be on when using melt transport postprocessors."));

        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,           ExcInternalError());
        Assert (duh[0].size() == this->introspection().n_components,          ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points, this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points, this->n_compositional_fields());
        create_melt_material_outputs(out);

        in.position = evaluation_points;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = duh[q][d];
            in.strain_rate[q] = symmetrize (grad_u);

            in.pressure[q]=uh[q][this->introspection().component_indices.pressure];
            in.temperature[q]=uh[q][this->introspection().component_indices.temperature];

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
          }

        this->get_material_model().evaluate(in, out);
        MaterialModel::MeltOutputs<dim> *melt_outputs = out.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
        AssertThrow(melt_outputs != NULL, ExcMessage("Need MeltOutputs from the material model for computing the melt viscosity."));

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          computed_quantities[q](0) = melt_outputs->fluid_viscosities[q];
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MeltDensity,
                                                  "melt density",
                                                  "A visualization output object that generates output "
                                                  "for melt density.")
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(CompactionViscosity,
                                                  "compaction viscosity",
                                                  "A visualization output object that generates output "
                                                  "for compaction viscosity.")
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Permeability,
                                                  "permeability",
                                                  "A visualization output object that generates output "
                                                  "for permeability.")
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MeltViscosity,
                                                  "melt viscosity",
                                                  "A visualization output object that generates output "
                                                  "for melt viscosity.")
    }
  }
}
