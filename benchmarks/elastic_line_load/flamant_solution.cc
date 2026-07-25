/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>
#include <aspect/geometry_model/box.h>
#include <aspect/postprocess/visualization.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {

      template <int dim>
      class FlamantSolution : public DataPostprocessorVector<dim>, public SimulatorAccess<dim>, public Interface<dim>
      {
        public:

          /**
           * Default constructor.
           */
          FlamantSolution();

          /**
           * Initialization function. Take references to the material model and
           * initial conditions model to get the parameters necessary for computing
           * the analytical solution for the shape of the solitary wave and store them.
           */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          std::vector<std::string>
          get_names () const override;

          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const override;

          double
          flamant_sigma_xx (const double x, const double y, const double p0, const double a) const;

          double
          flamant_sigma_xy (const double x, const double y, const double p0, const double a) const;

          double
          flamant_sigma_yy (const double x, const double y, const double p0, const double a) const;

        private:

      } ;

      template <int dim>
      FlamantSolution<dim>::
      FlamantSolution ()
        :
        DataPostprocessorVector<dim> ("flamant_analytic_stress",
                                      update_values | update_quadrature_points | update_gradients),
        Interface<dim>()
      {}

      template <int dim>
      void
      FlamantSolution<dim>::evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                                  std::vector<Vector<double>> &computed_quantities) const
      {
        AssertThrow(dim==2,
                    ExcMessage("The Flamant solution postprocessor is only implemented for 2D problems."));
        AssertThrow(Plugins::plugin_type_matches<GeometryModel::Box<dim>>(this->get_geometry_model()),
                    ExcMessage("The Flamant solution postprocessor requires a Box geometry model to be set."));


        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(in.n_evaluation_points(),
                                                     this->n_compositional_fields());
        this->get_material_model().evaluate(in, out);

        const double load_width = 20e3;
        const double p0         = 2.7e6;

        for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
          {
            const double x = input_data.evaluation_points[q][0];
            const double y = input_data.evaluation_points[q][1];

            computed_quantities[q][0] = -flamant_sigma_xx(x, y, p0, load_width);
            computed_quantities[q][1] =  flamant_sigma_xy(x, y, p0, load_width);
            computed_quantities[q][2] = -flamant_sigma_yy(x, y, p0, load_width);
          }
      }

      template <int dim>
      std::vector<std::string>
      FlamantSolution<dim>::get_names () const
      {
        std::vector<std::string> names;
        names.emplace_back("flamant_sigma_xx");
        names.emplace_back("flamant_sigma_xy");
        names.emplace_back("flamant_sigma_yy");
        return names;
      }

      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      FlamantSolution<dim>::get_data_component_interpretation () const
      {
        return
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          (3, DataComponentInterpretation::component_is_scalar);
      }

      template <int dim>
      double
      FlamantSolution<dim>::flamant_sigma_xx (const double x, const double y, const double p0, const double a) const
      {
        const auto &box = dynamic_cast<const GeometryModel::Box<dim> &>(this->get_geometry_model());

        const double x_extent = box.get_extents()[0];
        const double y_extent = box.get_extents()[1];
        const double xR = x_extent/2+a;
        const double xL = x_extent/2-a;
        const double theta1 = std::atan((x-xR)/(y_extent-y));
        const double theta2 = std::atan((x-xL)/(y_extent-y));
        const double val = p0 / M_PI * (theta2 - theta1 - 0.5 * (std::sin(2 * theta2) -
                                                                 std::sin(2 * theta1)));
        return -val;
      }

      template <int dim>
      double
      FlamantSolution<dim>::flamant_sigma_yy (const double x, const double y, const double p0, const double a) const
      {
        const auto &box = dynamic_cast<const GeometryModel::Box<dim> &>(this->get_geometry_model());
        const double x_extent = box.get_extents()[0];
        const double y_extent = box.get_extents()[1];
        const double xR = x_extent/2+a;
        const double xL = x_extent/2-a;
        const double theta1 = std::atan((x-xR)/(y_extent-y));
        const double theta2 = std::atan((x-xL)/(y_extent-y));
        const double val = p0 / M_PI * (theta2 - theta1 + 0.5 * (std::sin(2 * theta2) -
                                                                 std::sin(2 * theta1)));
        return -val;
      }

      template <int dim>
      double
      FlamantSolution<dim>::flamant_sigma_xy (const double x, const double y, const double p0, const double a) const
      {
        const auto &box = dynamic_cast<const GeometryModel::Box<dim> &>(this->get_geometry_model());
        const double x_extent = box.get_extents()[0];
        const double y_extent = box.get_extents()[1];
        const double xR = x_extent/2+a;
        const double xL = x_extent/2-a;
        const double theta1 = std::atan((x-xR)/(y_extent-y));
        const double theta2 = std::atan((x-xL)/(y_extent-y));
        const double val = p0 / M_PI * ((std::sin(theta2) * std::sin(theta2) -
                                         std::sin(theta1) * std::sin(theta1)));
        return -val;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(FlamantSolution,
                                                  "flamant analytic stress",
                                                  "A simple material model that is like the "
                                                  "'Simple' model, but tracks the finite strain as compositional "
                                                  "fields. More precisely, the model assumes that the first 4 (in 2D) "
                                                  "or 9 (in 3D) compositional fields contain the components "
                                                  "of the deformation gradient tensor, $\\mathbf F$, which can "
                                                  "be polar-decomposed into the left stretching tensor "
                                                  "$\\mathbf L$ (the finite strain we are interested in), and the "
                                                  "rotation tensor $\\mathbf Q$. See the corresponding cookbook in "
                                                  "the manual for more detailed information.")
    }
  }
}
