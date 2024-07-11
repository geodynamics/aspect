/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

#include <aspect/material_model/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>
#include <aspect/simulator.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

const double c = 1.0;


namespace aspect
{
  template <int dim>
  class TestMeltMaterial:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_darcy_coefficient () const
      {
        const double porosity = 0.01;
        const double permeability = porosity * porosity;
        return permeability / 1.0;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            const double porosity = in.composition[i][porosity_idx];
            const double x = in.position[i](0);
            const double z = in.position[i](1);
            out.viscosities[i] = std::exp(c * porosity);
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
            out.densities[i] = 1.0;
            // This is the RHS we need to use for the manufactured solution.
            // We calculated it by subtracting the first term of the RHS (\Nabla (K_D \rho_f g)) from the LHS
            // we computed using our analytical solution.
            for (unsigned int c=0; c<in.composition[i].size(); ++c)
              out.reaction_terms[i][c] = 0.0;
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim>>();

        if (melt_out != nullptr)
          {
            const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                double porosity = in.composition[i][porosity_idx];
                melt_out->compaction_viscosities[i] = std::exp(c * porosity);
                melt_out->fluid_viscosities[i] = 1.0;
                melt_out->permeabilities[i] = porosity * porosity;
                melt_out->fluid_density_gradients[i] = Tensor<1,dim>();
                melt_out->fluid_densities[i] = 1.0;
              }
          }
      }
  };




  template <int dim>
  class RefFunction : public Function<dim>
  {
    public:
      RefFunction () : Function<dim>(2*dim+5) {}
      virtual void vector_value (const Point<dim>   &p,
                                 Vector<double>   &values) const
      {
        double x = p(0);
        double z = p(1);
        const double porosity = 0.15 + 1.0/20.0 * (numbers::PI/2.0 + std::atan(x + 2*z));

        values[0]= 1.0;       // x vel
        values[1]= 1.0;    // z vel
        values[2]= std::exp(z/10.0);  // p_f
        values[3]= 0;  // p_c
        values[4]= 1.0;       // x melt vel
        values[5]= 1.0;    // y melt vel
        values[6]= std::exp(z/10.0);  // p_s
        values[7]= 0; // T
        values[8]= porosity;

        // We have to scale the compaction pressure solution to p_c_bar using sqrt(K_D / ref_K_D).
        // K_D is equal to the porosity (as defined in the material model).
        const double K_D = values[8] * values[8];
        const double ref_K_D = 0.01 * 0.01;
        const double p_c_scale = std::sqrt(K_D / ref_K_D);

        if (p_c_scale > 0)
          values[3] /= p_c_scale;
      }
  };

  /**
    * A postprocessor that evaluates the accuracy of the solution
    * by using the L2 norm.
    */
  template <int dim>
  class ConvergenceMeltPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      /**
       * Generate graphical output from the current solution.
       */
      virtual
      std::pair<std::string,std::string>
      execute (TableHandler &statistics);

  };

  template <int dim>
  std::pair<std::string,std::string>
  ConvergenceMeltPostprocessor<dim>::execute (TableHandler &statistics)
  {
    RefFunction<dim> ref_func;
    const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities +2);

    const unsigned int n_total_comp = this->introspection().n_components;

    Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_f (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_c (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_c_bar (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_u_f (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_porosity (this->get_triangulation().n_active_cells());

    ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                        n_total_comp);
    ComponentSelectFunction<dim> comp_p_f(dim, n_total_comp);
    ComponentSelectFunction<dim> comp_p_c(dim+1, n_total_comp);
    ComponentSelectFunction<dim> comp_u_f(std::pair<unsigned int, unsigned int>(dim+2,dim+2+dim),
                                          n_total_comp);
    ComponentSelectFunction<dim> comp_p(dim+2+dim, n_total_comp);
    ComponentSelectFunction<dim> comp_porosity(dim+2+dim+2, n_total_comp);

    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_u,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_u);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p_f,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p_f);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p_c_bar,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p_c);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_porosity,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_porosity);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_u_f,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_u_f);


    // Loop over all cells to compute the error for p_c from p_c_bar
    const QGauss<dim> quadrature(this->get_parameters().stokes_velocity_degree+1);
    FEValues<dim> fe_values (this->get_mapping(),
                             this->get_fe(),
                             quadrature,
                             update_quadrature_points | update_values | update_gradients | update_JxW_values);

    MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
    MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());

    MeltHandler<dim>::create_material_model_outputs(out);

    typename DoFHandler<dim>::active_cell_iterator
    cell = this->get_dof_handler().begin_active(),
    endc = this->get_dof_handler().end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          in.reinit(fe_values, cell, this->introspection(), this->get_solution());

          this->get_material_model().evaluate(in, out);

          const double p_c_scale = dynamic_cast<const MaterialModel::MeltInterface<dim>*>(&this->get_material_model())->p_c_scale(in, out, this->get_melt_handler(), true);

          const unsigned int i = cell->active_cell_index();
          cellwise_errors_p_c[i] = cellwise_errors_p_c_bar[i] * p_c_scale;
        }

    const double u_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_u, VectorTools::L2_norm);
    const double p_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p, VectorTools::L2_norm);
    const double p_f_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p_f, VectorTools::L2_norm);
    const double p_c_bar_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p_c_bar, VectorTools::L2_norm);
    const double p_c_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p_c, VectorTools::L2_norm);
    const double poro_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_porosity, VectorTools::L2_norm);
    const double u_f_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_u_f, VectorTools::L2_norm);

    std::ostringstream os;
    os << std::scientific << u_l2
       << ", " << p_l2
       << ", " << p_f_l2
       << ", " << p_c_bar_l2
       << ", " << p_c_l2
       << ", " << poro_l2
       << ", " << u_f_l2;

    return std::make_pair("Errors u_L2, p_L2, p_f_L2, p_c_bar_L2, p_c_L2, porosity_L2, u_f_L2:", os.str());
  }


  template <int dim>
  class PressureBdry:

    public BoundaryFluidPressure::Interface<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const types::boundary_id boundary_indicator,
        const typename MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
        const typename MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
        const std::vector<Tensor<1,dim>> &normal_vectors,
        std::vector<double> &output
      ) const
      {
        for (unsigned int q=0; q<output.size(); ++q)
          {
            const double z = material_model_inputs.position[q][1];
            Tensor<1,dim> gravity;
            gravity[0] = 0.0;
            gravity[dim-1] = 1.0;
            output[q] = 0.1 * std::exp(z/10.0) * gravity * normal_vectors[q];
          }
      }



  };

}

// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(TestMeltMaterial,
                                 "test melt material",
                                 "")


  ASPECT_REGISTER_POSTPROCESSOR(ConvergenceMeltPostprocessor,
                                "melt error calculation",
                                "A postprocessor that compares the numerical solution to the analytical "
                                "solution derived for incompressible melt transport in a 2D box as described "
                                "in the manuscript and reports the error.")

  ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(PressureBdry,
                                                "PressureBdry",
                                                "A fluid pressure boundary condition that prescribes the "
                                                "gradient of the fluid pressure at the boundaries as "
                                                "calculated in the analytical solution. ")

}
