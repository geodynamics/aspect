/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/maximum_horizontal_compressive_stress.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      void
      MaximumHorizontalCompressiveStress<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert ((computed_quantities[0].size() == dim),
                ExcInternalError());
        Assert (uh[0].size() == this->introspection().n_components,   ExcInternalError());
        Assert (duh[0].size() == this->introspection().n_components,  ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(n_quadrature_points,
                                                   this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        // collect input information to compute the viscosity at every evaluation point
        in.position = evaluation_points;
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = duh[q][d];
            in.strain_rate[q] = symmetrize (grad_u);

            in.pressure[q]=uh[q][this->introspection().component_indices.pressure];
            in.temperature[q]=uh[q][this->introspection().component_indices.temperature];

            for (unsigned int d = 0; d < dim; ++d)
              {
                in.velocity[q][d]=uh[q][this->introspection().component_indices.velocities[d]];
                in.pressure_gradient[q][d] = duh[q][this->introspection().component_indices.pressure][d];
              }

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[q][c] = uh[q][this->introspection().component_indices.compositional_fields[c]];
          }

        // then do compute the viscosity...
        this->get_material_model().evaluate(in, out);

        // ...and use it to compute the stresses and from that the
        // maximum compressive stress direction
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];
            const SymmetricTensor<2,dim> compressible_strain_rate
              = (this->get_material_model().is_compressible()
                 ?
                 strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                 :
                 strain_rate);

            const double eta = out.viscosities[q];

            // first compute the stress tensor, ignoring the pressure
            // for the moment (the pressure has no effect on the
            // direction since it just adds a multiple of the identity
            // matrix to the stress, but because it is large, it may
            // lead to numerical instabilities)
            //
            // note that the *compressive* stress is simply the
            // negative stress
            const SymmetricTensor<2,dim> compressive_stress = -2*eta*compressible_strain_rate;

            // then find a set of (dim-1) horizontal, unit-length, mutually orthogonal vectors
            const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector (in.position[q]);
            const Tensor<1,dim> vertical_direction = gravity/gravity.norm();
            std_cxx11::array<Tensor<1,dim>,dim-1 > orthogonal_directions
              = Utilities::orthogonal_vectors(vertical_direction);
            for (unsigned int i=0; i<orthogonal_directions.size(); ++i)
              orthogonal_directions[i] /= orthogonal_directions[i].norm();

            Tensor<1,dim> maximum_horizontal_compressive_stress;
            switch (dim)
              {
                // in 2d, there is only one horizontal direction, and
                // we have already computed it above. give it the
                // length of the compressive_stress (now taking into account the
                // pressure) in this direction
                case 2:
                {
                  maximum_horizontal_compressive_stress = orthogonal_directions[0] *
                                                          (orthogonal_directions[0] *
                                                           (compressive_stress +
                                                            in.pressure[q] * unit_symmetric_tensor<dim>()) *
                                                           orthogonal_directions[0]);
                  break;
                }

                // in 3d, use the formulas discussed in the
                // documentation of the plugin below
                case 3:
                {
                  const double a = orthogonal_directions[0] *
                                   compressive_stress *
                                   orthogonal_directions[0];
                  const double b = orthogonal_directions[1] *
                                   compressive_stress *
                                   orthogonal_directions[1];
                  const double c = orthogonal_directions[0] *
                                   compressive_stress *
                                   orthogonal_directions[1];

                  // compute the two stationary points of f(alpha)
                  const double alpha_1 = 1./2 * std::atan2 (c, a-b);
                  const double alpha_2 = alpha_1 + numbers::PI/2;

                  // then check the sign of f''(alpha) to determine
                  // which of the two stationary points is the maximum
                  const double f_double_prime_1 = 2*(b-a)*std::cos(2*alpha_1)
                                                  - 2*c*sin(2*alpha_1);
                  double alpha;
                  if (f_double_prime_1 < 0)
                    alpha = alpha_1;
                  else
                    {
                      Assert (/* f_double_prime_2 = */
                        2*(b-a)*std::cos(2*alpha_2) - 2*c*sin(2*alpha_2) <= 0,
                        ExcInternalError());
                      alpha = alpha_2;
                    }

                  // then re-assemble the maximum horizontal compressive_stress
                  // direction from alpha and the two horizontal
                  // directions
                  const Tensor<1,dim> n = std::cos(alpha) * orthogonal_directions[0] +
                                          std::sin(alpha) * orthogonal_directions[1];

                  // finally compute the actual direction * magnitude,
                  // now taking into account the pressure (with the
                  // correct sign for the *compressive* stress)
                  maximum_horizontal_compressive_stress
                    = n * (n * (compressive_stress
                                -
                                in.pressure[q] * unit_symmetric_tensor<dim>()) * n);

                  break;
                }

                default:
                  Assert (false, ExcNotImplemented());
              }

            for (unsigned int i=0; i<dim; ++i)
              computed_quantities[q](i) = maximum_horizontal_compressive_stress[i];
          }
      }


      template <int dim>
      std::vector<std::string>
      MaximumHorizontalCompressiveStress<dim>::get_names () const
      {
        return std::vector<std::string> (dim, "maximum_horizontal_compressive_stress");
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      MaximumHorizontalCompressiveStress<dim>::get_data_component_interpretation () const
      {
        return
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          (dim,
           DataComponentInterpretation::component_is_part_of_vector);
      }



      template <int dim>
      UpdateFlags
      MaximumHorizontalCompressiveStress<dim>::get_needed_update_flags () const
      {
        return update_gradients | update_values | update_q_points;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MaximumHorizontalCompressiveStress,
                                                  "maximum horizontal compressive stress",
                                                  "A plugin that computes the direction and magnitude of the "
                                                  "maximum horizontal component of the compressive stress as a vector "
                                                  "field. The direction of this vector can often be used to "
                                                  "visualize the principal mode of deformation (e.g., at "
                                                  "normal faults or extensional margins) and can be "
                                                  "correlated with seismic anisotropies. Recall that the "
                                                  "\\textit{compressive} stress is simply the negative stress, "
                                                  "$\\sigma_c=-\\sigma=-\\left["
                                                  "     2\\eta (\\varepsilon(\\mathbf u)"
                                                  "             - \\frac 13 \\text{trace}\\varepsilon(\\mathbf u) I)"
                                                  "     + pI\\right]$."
                                                  "\n\n"
                                                  "Following \\cite{LundTownend07}, we define the maximum horizontal "
                                                  "stress direction as that \\textit{horizontal} direction "
                                                  "$\\mathbf n$ that maximizes $\\mathbf n^T \\sigma_c \\mathbf n$. We "
                                                  "call a vector \\textit{horizontal} if it is perpendicular to the "
                                                  "gravity vector $\\mathbf g$."
                                                  "\n\n"
                                                  "In two space dimensions, $\\mathbf n$ is simply a vector that "
                                                  "is horizontal (we choose one of the two possible choices)."
                                                  "\n\n"
                                                  "In three space dimensions, given two horizontal, perpendicular, "
                                                  "unit length, but otherwise arbitrarily chosen vectors "
                                                  "$\\mathbf u,\\mathbf v$, we can express "
                                                  "$\\mathbf n = (\\cos \\alpha)\\mathbf u + (\\sin\\alpha)\\mathbf v$ "
                                                  "where $\\alpha$ maximizes the expression "
                                                  "\\begin{align*}"
                                                  "  f(\\alpha) = \\mathbf n^T \\sigma_c \\mathbf n"
                                                  "  = (\\mathbf u^T \\sigma_c \\mathbf u)(\\cos\\alpha)^2"
                                                  "    +2(\\mathbf u^T \\sigma_c \\mathbf v)(\\cos\\alpha)(\\sin\\alpha)"
                                                  "    +(\\mathbf v^T \\sigma_c \\mathbf v)(\\sin\\alpha)^2."
                                                  "\\end{align*}"
                                                  "\n\n"
                                                  "The maximum of $f(\\alpha)$ is attained where $f'(\\alpha)=0$. "
                                                  "Evaluating the derivative and using trigonometric identities, "
                                                  "one finds that $\\alpha$ has to satisfy the equation "
                                                  "\\begin{align*}"
                                                  "  \\tan(2\\alpha) = \\frac{\\mathbf u^T \\sigma_c \\mathbf v}"
                                                  "                          {\\mathbf u^T \\sigma_c \\mathbf u "
                                                  "                           - \\mathbf v^T \\sigma_c \\mathbf v}."
                                                  "\\end{align*}"
                                                  "Since the transform $\\alpha\\mapsto\\alpha+\\pi$ flips the "
                                                  "direction of $\\mathbf n$, we only need to seek a solution "
                                                  "to this equation in the interval $\\alpha\\in[0,\\pi)$. "
                                                  "These are given by "
                                                  "$\\alpha_1=\\frac 12 \\arctan \\frac{\\mathbf u^T \\sigma_c "
                                                  "\\mathbf v}{\\mathbf u^T \\sigma_c \\mathbf u - "
                                                  "\\mathbf v^T \\sigma_c \\mathbf v}$ and "
                                                  "$\\alpha_2=\\alpha_1+\\frac{\\pi}{2}$, one of which will "
                                                  "correspond to a minimum and the other to a maximum of "
                                                  "$f(\\alpha)$. One checks the sign of "
                                                  "$f''(\\alpha)=-2(\\mathbf u^T \\sigma_c \\mathbf u - "
                                                  "\\mathbf v^T \\sigma_c \\mathbf v)\\cos(2\\alpha) "
                                                  "- 2 (\\mathbf u^T \\sigma_c \\mathbf v) \\sin(2\\alpha)$ for "
                                                  "each of these to determine the $\\alpha$ that maximizes "
                                                  "$f(\\alpha)$, and from this immediately arrives at the correct "
                                                  "form for the maximum horizontal stress $\\mathbf n$."
                                                  "\n\n"
                                                  "The description above computes a \\textit{direction} vector "
                                                  "$\\mathbf n$. The quantity this plugin then outputs is "
                                                  "this direction scaled by the size of the horizontal stress in "
                                                  "this direction, i.e., the vector "
                                                  "$\\mathbf w = (\\mathbf n^T \\sigma_c \\mathbf n) \\; \\mathbf n$.")
    }
  }
}
