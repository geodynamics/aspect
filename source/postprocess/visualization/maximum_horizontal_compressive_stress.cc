/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/maximum_horizontal_compressive_stress.h>
#include <aspect/gravity_model/interface.h>
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
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert ((computed_quantities[0].size() == dim),
                ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,  ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        // Compute the viscosity...
        this->get_material_model().evaluate(in, out);

        // ...and use it to compute the stresses and from that the
        // maximum compressive stress direction
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];
            const SymmetricTensor<2,dim> deviatoric_strain_rate
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
            const SymmetricTensor<2,dim> compressive_stress = -2*eta*deviatoric_strain_rate;

            // then find a set of (dim-1) horizontal, unit-length, mutually orthogonal vectors
            const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector (in.position[q]);
            const Tensor<1,dim> vertical_direction = gravity/gravity.norm();
            std::array<Tensor<1,dim>,dim-1 > orthogonal_directions
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
                                                           ((compressive_stress
                                                             -
                                                             in.pressure[q] * unit_symmetric_tensor<dim>()) *
                                                            orthogonal_directions[0]));
                  break;
                }

                // in 3d, use the formulas discussed in the
                // documentation of the plugin below
                case 3:
                {
                  const double a = orthogonal_directions[0] *
                                   (compressive_stress *
                                    orthogonal_directions[0]);
                  const double b = orthogonal_directions[1] *
                                   (compressive_stress *
                                    orthogonal_directions[1]);
                  const double c = 2.0*orthogonal_directions[0] *
                                   (compressive_stress *
                                    orthogonal_directions[1]);

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
                  // correct sign in front of the pressure for the
                  // *compressive* stress)
                  //
                  // the magnitude is computed as discussed in the
                  // description of the plugin below
                  const double maximum_horizontal_compressive_stress_magnitude
                    = (n * ((compressive_stress
                             -
                             in.pressure[q] * unit_symmetric_tensor<dim>()) * n));
                  const Tensor<1,dim> n_perp = std::sin(alpha) * orthogonal_directions[0] -
                                               std::cos(alpha) * orthogonal_directions[1];

                  const double minimum_horizontal_compressive_stress_magnitude
                    = (n_perp * ((compressive_stress
                                  -
                                  in.pressure[q] * unit_symmetric_tensor<dim>()) * n_perp));

                  maximum_horizontal_compressive_stress
                    = n * (maximum_horizontal_compressive_stress_magnitude -
                           minimum_horizontal_compressive_stress_magnitude);

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
        return update_gradients | update_values | update_quadrature_points;
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
                                                  "correlated with seismic anisotropy. Recall that the "
                                                  "\\textit{compressive} stress is simply the negative stress, "
                                                  "$\\sigma_c=-\\sigma=-\\left["
                                                  "     2\\eta (\\varepsilon(\\mathbf u)"
                                                  "             - \\frac 13 (\\nabla \\cdot \\mathbf u) I)"
                                                  "     + pI\\right]$."
                                                  "\n\n"
                                                  "Following \\cite{LundTownend07}, we define the maximum horizontal "
                                                  "stress direction as that \\textit{horizontal} direction "
                                                  "$\\mathbf n$ that maximizes $\\mathbf n^T \\sigma_c \\mathbf n$. We "
                                                  "call a vector \\textit{horizontal} if it is perpendicular to the "
                                                  "gravity vector $\\mathbf g$."
                                                  "\n\n"
                                                  "In two space dimensions, $\\mathbf n$ is simply a vector that "
                                                  "is horizontal (we choose one of the two possible choices). "
                                                  "This direction is then scaled by the size of the horizontal stress "
                                                  "in this direction, i.e., the plugin outputs the vector "
                                                  "$\\mathbf w = (\\mathbf n^T \\sigma_c \\mathbf n) \\; \\mathbf n$."
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
                                                  "  \\tan(2\\alpha) = \\frac{2.0\\mathbf u^T \\sigma_c \\mathbf v}"
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
                                                  "The description above computes a 3d \\textit{direction} vector "
                                                  "$\\mathbf n$. If one were to scale this vector the same way "
                                                  "as done in 2d, i.e., with the magnitude of the stress in "
                                                  "this direction, one will typically get vectors whose length "
                                                  "is principally determined by the hydrostatic pressure at "
                                                  "a given location simply because the hydrostatic pressure "
                                                  "is the largest component of the overall stress. On the other "
                                                  "hand, the hydrostatic pressure does not determine any "
                                                  "principal direction because it is an isotropic, anti-compressive "
                                                  "force. As a consequence, there are often points in simulations "
                                                  "(e.g., at the center of convection rolls) where the stress has "
                                                  "no dominant horizontal direction, and the algorithm above will "
                                                  "then in essence choose a random direction because the stress "
                                                  "is approximately equal in all horizontal directions. If one "
                                                  "scaled the output by the magnitude of the stress in this "
                                                  "direction (i.e., approximately equal to the hydrostatic "
                                                  "pressure at this point), one would get randomly oriented vectors "
                                                  "at these locations with significant lengths."
                                                  "\n\n"
                                                  "To avoid this problem, we scale the maximal horizontal "
                                                  "compressive stress direction $\\mathbf n$ by the \\textit{difference} "
                                                  "between the stress in the maximal and minimal horizontal stress "
                                                  "directions. In other words, let "
                                                  "$\\mathbf n_\\perp=(\\sin \\alpha)\\mathbf u - (\\cos\\alpha)\\mathbf v$ "
                                                  "be the horizontal direction perpendicular to "
                                                  "$\\mathbf n$, then this plugin outputs the vector quantity "
                                                  "$\\mathbf w = (\\mathbf n^T \\sigma_c \\mathbf n "
                                                  "               -\\mathbf n^T_\\perp \\sigma_c \\mathbf n_\\perp) "
                                                  "              \\; \\mathbf n$. "
                                                  "In other words, the length of the vector produced indicates "
                                                  "\\textit{how dominant} the direction of maximal horizontal "
                                                  "compressive strength is."
                                                  "\n\n"
                                                  "Fig.~\\ref{fig:max-horizontal-compressive-stress} shows a "
                                                  "simple example for this kind of visualization in 3d."
                                                  "\n\n"
                                                  "\\begin{figure}"
                                                  "  \\includegraphics[width=0.3\\textwidth]"
                                                  "    {viz/plugins/maximum_horizontal_compressive_stress/temperature.png}"
                                                  "  \\hfill"
                                                  "  \\includegraphics[width=0.3\\textwidth]"
                                                  "    {viz/plugins/maximum_horizontal_compressive_stress/velocity.png}"
                                                  "  \\hfill"
                                                  "  \\includegraphics[width=0.3\\textwidth]"
                                                  "    {viz/plugins/maximum_horizontal_compressive_stress/horizontal-stress.png}"
                                                  "  \\caption{\\it Illustration of the `maximum horizontal "
                                                  "    compressive stress' visualization plugin. The left "
                                                  "    figure shows a ridge-like temperature anomaly. Together "
                                                  "    with no-slip boundary along all six boundaries, this "
                                                  "    results in two convection rolls (center). The maximal "
                                                  "    horizontal compressive strength at the bottom center "
                                                  "    of the domain is perpendicular to the ridge because "
                                                  "    the flow comes together there from the left and right, "
                                                  "    yielding a compressive force in left-right direction. "
                                                  "    At the top of the model, the flow separates outward, "
                                                  "    leading to a \\textit{negative} compressive stress "
                                                  "    in left-right direction; because there is no flow "
                                                  "    in front-back direction, the compressive strength "
                                                  "    in front-back direction is zero, making the along-ridge "
                                                  "    direction the dominant one. At the center of the "
                                                  "    convection rolls, both horizontal directions yield "
                                                  "    the same stress; the plugin therefore chooses an "
                                                  "    essentially arbitrary horizontal vector, but then "
                                                  "    uses a zero magnitude given that the difference "
                                                  "    between the maximal and minimal horizontal stress "
                                                  "    is zero at these points.}"
                                                  "  \\label{fig:max-horizontal-compressive-stress}"
                                                  "\\end{figure}"
                                                 )
    }
  }
}
