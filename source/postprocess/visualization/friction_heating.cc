/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization.h>
#include <aspect/simulator.h>

#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessor that takes an output vector and
       * computes a variable that represents the friction heating term
       * $\tau:\varepsilon$.
       *
       * The member functions are all implementations of those declared in the base
       * class. See there for their meaning.
       */
      template <int dim>
      class FrictionHeating
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          FrictionHeating ();

          virtual
          void
          compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                             const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                             const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                             const std::vector<Point<dim> >                  &normals,
                                             const std::vector<Point<dim> >                  &evaluation_points,
                                             std::vector<Vector<double> >                    &computed_quantities) const;
      };


      template <int dim>
      FrictionHeating<dim>::
      FrictionHeating ()
        :
        DataPostprocessorScalar<dim> ("friction_heating",
                                      update_values | update_gradients | update_q_points)
      {}



      template <int dim>
      void
      FrictionHeating<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                         const std::vector<Point<dim> >                  &normals,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,  ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                 ExcInternalError());
        Assert (uh[0].size() == dim+2,                              ExcInternalError());
        Assert (duh[0].size() == dim+2,                             ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            // extract the primal variables
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = duh[q][d];

            const double pressure    = uh[q][dim];
            const double temperature = uh[q][dim+1];

            const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);
            const SymmetricTensor<2,dim> compressible_strain_rate
              = (this->get_material_model().is_compressible()
                 ?
                 strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                 :
                 strain_rate);
            computed_quantities[q](0) = 2 * this->get_material_model().viscosity(temperature,
                                                                                 pressure,
                                                                                 strain_rate,
                                                                                 evaluation_points[q]) *
                                        compressible_strain_rate * compressible_strain_rate;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(FrictionHeating,
                                                  "friction heating",
                                                  "A visualization output object that generates output "
                                                  "for the amount of friction heating often referred "
                                                  "to as $\\tau:\\epsilon$. More concisely, in the "
                                                  "incompressible case, the quantity that is output "
                                                  "is defined as "
                                                  "$\\eta \\varepsilon(\\mathbf u):\\varepsilon(\\mathbf u)$ "
                                                  "where $\\eta$ is itself a function of temperature, "
                                                  "pressure and strain rate. In the compressible case, "
                                                  "the quantity that's computed is "
                                                  "$\\eta [\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u)\\mathbf I]:"
                                                  "[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u)\\mathbf I]$.")
    }
  }
}
