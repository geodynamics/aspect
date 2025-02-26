/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/shear_stress_eigenvectors.h>
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
      ShearStressEigenvectors<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert ((computed_quantities[0].size() == dim*4),
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
            const SymmetricTensor<2,dim> shear_stress = 2*eta*compressible_strain_rate;
            const std::array<std::pair<double,Tensor<1,dim,double> >, dim> eigenvector_array = eigenvectors(shear_stress, SymmetricTensorEigenvectorMethod::jacobi);

            // Check whether sigma1 > sigma2 > simga3
            Assert(eigenvector_array[0].first > eigenvector_array[1].first, ExcMessage("the eigenvalue in array 0 should be larger then in array 1."));
            Assert(eigenvector_array[1].first > eigenvector_array[2].first, ExcMessage("the eigenvalue in array 1 should be larger then in array 2."));

            // Check whether the vectors are orthogonal
            Assert(std::fabs(eigenvector_array[0].second * eigenvector_array[1].second) < 1e-8, ExcMessage("Eigenvector in 0 is not orthogonal to eigenvector 1:"
                   "eigenvector_array[1].second * eigenvector_array[1].second = " + std::to_string(eigenvector_array[0].second * eigenvector_array[1].second)));
            Assert(std::fabs(eigenvector_array[1].second * eigenvector_array[2].second) < 1e-8, ExcMessage("Eigenvector in 1 is not orthogonal to eigenvector 2:"
                   "eigenvector_array[1].second * eigenvector_array[2].second = " + std::to_string(eigenvector_array[1].second * eigenvector_array[2].second)));
            Assert(std::fabs(eigenvector_array[0].second * eigenvector_array[2].second) < 1e-8, ExcMessage("Eigenvector in 0 is not orthogonal to eigenvector 2:"
                   "eigenvector_array[0].second * eigenvector_array[2].second = " + std::to_string(eigenvector_array[0].second * eigenvector_array[2].second)));

            const Tensor<1,dim> tau_direction = eigenvector_array[2].second - eigenvector_array[0].second;

            const double tau_magnitude = 0.5 * (eigenvector_array[0].first - eigenvector_array[2].first);


            Point<dim> p = in.position[q];
            double x = p[0];
            double y = p[1];
            double z = p[2];
            double e = 1.1e-1;
            Tensor<1,dim> ev0 = eigenvector_array[0].second;
            Tensor<1,dim> ev1 = eigenvector_array[1].second;
            Tensor<1,dim> ev2 = eigenvector_array[2].second;
            //if(fabs(x-0.5) < e && fabs(y-0.25) < e && fabs(z-0.375) < e || fabs(x-0.625) < e && fabs(y-0.375) < e && fabs(z-0.5) < e)
            {
              //std::cout << x << "," << y << "," << z << ", strain rate = " << compressible_strain_rate << std::endl;
              //std::cout << x << "," << y << "," << z << ": tau_magnitude = " << tau_magnitude <<  ", eigenvector 0 = " << eigenvector_array[0].second << ", 1 = " << eigenvector_array[1].second<< ", 2 = " << eigenvector_array[2].second << " = " << tau_direction << ", norm = " << tau_direction/tau_direction.norm() << std::endl;
              std::cout << x << " " << y << " " << z << " "
                        << eigenvector_array[0].first << " "<< eigenvector_array[1].first << " "<< eigenvector_array[2].first << " "
                        << ev0 << " " << ev1 << " " << ev2
                        << tau_direction << std::endl;
            }

            for (unsigned int i=0; i<dim; ++i)
              {
                // the eigenvectors are unit
                AssertThrow(std::fabs(eigenvector_array[0].second.norm()-1) < 1e-8,ExcMessage("eigenvector array 0 norm is not one: " + std::to_string(eigenvector_array[0].second.norm())));
                AssertThrow(std::fabs(eigenvector_array[1].second.norm()-1) < 1e-8,ExcMessage("eigenvector array 1 norm is not one: " + std::to_string(eigenvector_array[1].second.norm())));
                AssertThrow(std::fabs(eigenvector_array[2].second.norm()-1) < 1e-8,ExcMessage("eigenvector array 2 norm is not one: " + std::to_string(eigenvector_array[2].second.norm())));
                Tensor<1,dim> scaled_eigenvector = eigenvector_array[i].second * eigenvector_array[i].first;
                for (unsigned int j=0; j<dim; ++j)
                  {
                    computed_quantities[q](i*dim+j) = scaled_eigenvector[j];
                  }
              }

            for (unsigned int i=0; i<dim; ++i)
              computed_quantities[q](dim*dim+i) = eigenvector_array[i].first;

          }
      }


      template <int dim>
      std::vector<std::string>
      ShearStressEigenvectors<dim>::get_names () const
      {
        std::vector<std::string> names;
        switch (dim)
          {
            case 2:
              names.emplace_back("shear_stress_eigenvectors_sigma_1");
              names.emplace_back("shear_stress_eigenvectors_sigma_1");
              names.emplace_back("shear_stress_eigenvectors_sigma_1");
              names.emplace_back("shear_stress_eigenvectors_sigma_2");
              names.emplace_back("shear_stress_eigenvectors_sigma_2");
              names.emplace_back("shear_stress_eigenvectors_sigma_2");
              names.emplace_back("shear_stress_eigenvalue_sigma_1");
              names.emplace_back("shear_stress_eigenvalue_sigma_2");
              break;

            case 3:
              names.emplace_back("shear_stress_eigenvectors_sigma_1");
              names.emplace_back("shear_stress_eigenvectors_sigma_1");
              names.emplace_back("shear_stress_eigenvectors_sigma_1");
              names.emplace_back("shear_stress_eigenvectors_sigma_2");
              names.emplace_back("shear_stress_eigenvectors_sigma_2");
              names.emplace_back("shear_stress_eigenvectors_sigma_2");
              names.emplace_back("shear_stress_eigenvectors_sigma_3");
              names.emplace_back("shear_stress_eigenvectors_sigma_3");
              names.emplace_back("shear_stress_eigenvectors_sigma_3");
              names.emplace_back("shear_stress_eigenvalue_sigma_1");
              names.emplace_back("shear_stress_eigenvalue_sigma_2");
              names.emplace_back("shear_stress_eigenvalue_sigma_3");
              break;

            default:
              Assert (false, ExcNotImplemented());
          }
        return names;
      }

      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      ShearStressEigenvectors<dim>::get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> components
        (dim*dim,DataComponentInterpretation::component_is_part_of_vector);
        components.emplace_back(DataComponentInterpretation::component_is_scalar);
        components.emplace_back(DataComponentInterpretation::component_is_scalar);
        components.emplace_back(DataComponentInterpretation::component_is_scalar);
        return components;
      }



      template <int dim>
      UpdateFlags
      ShearStressEigenvectors<dim>::get_needed_update_flags () const
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(ShearStressEigenvectors,
                                                  "shear stress eigenvectors",
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
                                                  "The description above computes a 3d \\textit{direction} vector "
                                                  "$\\mathbf n$. If one were to scale this vector the same way "
                                                  "as done in 2d, i.e., with the magnitude of the stress in "
                                                  "this direction, one will typically get vectors whose length "
                                                  "is principally determined by the hydrostatic pressure at "
                                                  "a given location simply because the hydrostatic pressure "
                                                  "is the largest component of the overall stress. On the other "
                                                  "hand, the hydrostatic pressure does not determine any "
                                                  "principle direction because it is an isotropic, anti-compressive "
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

