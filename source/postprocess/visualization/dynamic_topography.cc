/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/postprocess/visualization/dynamic_topography.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      DynamicTopography<dim>::
      DynamicTopography ()
        :
        DataPostprocessorScalar<dim> ("dynamic_topography",
                                      update_values | update_gradients | update_q_points)
      {}



      template <int dim>
      void
      DynamicTopography<dim>::
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
        Assert (uh[0].size() == dim+2+this->n_compositional_fields(), ExcInternalError());
        Assert (duh[0].size() == dim+2+this->n_compositional_fields(),ExcInternalError());

        typename MaterialModel::Interface<dim>::MaterialModelInputs in(n_quadrature_points,
                                                                       this->n_compositional_fields());
        typename MaterialModel::Interface<dim>::MaterialModelOutputs out(n_quadrature_points,
                                                                         this->n_compositional_fields());

        // fill the various fields necessary to evaluate the material
        // properties
        in.position = evaluation_points;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = duh[q][d];
            in.strain_rate[q] = symmetrize (grad_u);

            in.pressure[q]=uh[q][dim];
            in.temperature[q]=uh[q][dim+1];

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][dim+2+c];
          }

        // evaluate the material model
        this->get_material_model().evaluate(in, out);

        // for each of the evaluation points, compute the dynamic
        // topography and put it into the output vector
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const Point<dim> location = evaluation_points[q];
            const double viscosity = out.viscosities[q];
            const double density   = out.densities[q];

//TODO: We need to subtract 2/3*div(u) from the stress here in the compressible case
            const SymmetricTensor<2,dim> stress = 2 * viscosity * in.strain_rate[q];

            const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(location);
            const Tensor<1,dim> gravity_direction = gravity/gravity.norm();

            const double sigma_rr           = gravity_direction * (stress * gravity_direction);
            const double dynamic_topography = -sigma_rr / gravity.norm() / density;

            computed_quantities[q](0) = dynamic_topography;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(DynamicTopography,
                                                  "dynamic topography",
                                                  "A visualization output object that generates output "
                                                  "for the dynamic topography. The approach to determine the "
                                                  "dynamic topography requires us to compute the stress tensor and "
                                                  "evaluate the component of it in the direction in which "
                                                  "gravity acts. In other words, we compute "
                                                  "$\\sigma_{rr}={\\hat g}^T(2 * \\eta \\varepsilon(\\mathbf u))\\hat g$ "
                                                  "where $\\hat g = \\mathbf g/\\|\\mathbf g\\|$ is the direction of "
                                                  "the gravity vector $\\mathbf g$. From this, the dynamic "
                                                  "topography is computed using the formula "
                                                  "$h=\\frac{\\sigma_{rr}}{\\|\\mathbf g\\| \\rho}$ where $\\rho$ "
                                                  "is the density at the cell center."
                                                  "\n\n"
                                                  "Strictly speaking, the dynamic topography is of course a "
                                                  "quantity that is only of interest at the surface. However, "
                                                  "we compute it everywhere to make things fit into the framework "
                                                  "within which we produce data for visualization. You probably "
                                                  "only want to visualize whatever data this postprocessor generates "
                                                  "at the surface of your domain and simply ignore the rest of the "
                                                  "data generated.")
    }
  }
}
