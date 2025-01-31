/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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
#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>
#include <aspect/boundary_fluid_pressure/interface.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  template <int dim>
  class MeltMaterial:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_darcy_coefficient () const
      {
        const double porosity = 0.01;
        const double permeability = 1.0 + 0.001 / (1.0 - porosity);
        return permeability / 2.0;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {

        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            const double y = in.position[i][dim-1] + 0.1;
            out.viscosities[i] = 3.0/4.0;
            out.densities[i] = 2 + Utilities::fixed_power<2>(y)/0.001;
            out.thermal_expansion_coefficients[i] = 1.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
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
                melt_out->compaction_viscosities[i] = 1.0;
                melt_out->fluid_viscosities[i]= 2.0;
                melt_out->permeabilities[i]= 1.0 + 0.001 / (1.0 - porosity);
                melt_out->fluid_densities[i]= 2.0;
                melt_out->fluid_density_gradients[i] = Tensor<1,dim>();
              }
          }
      }

  };





  template <int dim>
  class RefFunction : public Function<dim>
  {
    public:
      RefFunction () : Function<dim>(9) {}
      virtual void vector_value (const Point<dim>   &p,
                                 Vector<double>   &values) const
      {
        double x = p(0);
        double y = p(1);

        const double porosity = 1.0 - 0.001/Utilities::fixed_power<2>(y+0.1);
        const double p_f = 1.0 - (y+0.1);
        const double p_s = 0.0;

        values[0]=0;  // x vel
        values[1]=0.5 * Utilities::fixed_power<2>(y+0.1);  // y vel
        values[2]=p_f;  // p_f
        values[3]=(1.0-porosity)*(p_s-p_f); // p_c
        values[4]=-1; // u_f_x
        values[5]=-1; // u_f_y
        values[6]=p_s;  // p_s
        values[7]=0; // T
        values[8]=porosity; // porosity
      }
  };



  template <int dim>
  class MMPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      virtual
      std::pair<std::string,std::string>
      execute (TableHandler &statistics);
  };


  template <int dim>
  std::pair<std::string,std::string>
  MMPostprocessor<dim>::execute (TableHandler &statistics)
  {
    RefFunction<dim> ref_func;
    const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities +2);

    Vector<float> cellwise_errors_porosity (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());

    const unsigned int n_total_comp = this->introspection().n_components;
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
                                       cellwise_errors_p,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p_f);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_porosity,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_porosity);

    const double u_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_u.norm_sqr(),this->get_mpi_communicator()));
    const double p_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_p.norm_sqr(),this->get_mpi_communicator()));
    const double poro_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_porosity.norm_sqr(),this->get_mpi_communicator()));

    std::ostringstream os;
    os << std::scientific << u_l2
       << ", " << p_l2
       << ", " << poro_l2;

    return std::make_pair("Errors u_L2, p_f_L2, porosity_L2:", os.str());

  }

  template <int dim>
  class PressureBdry:

    public BoundaryFluidPressure::Interface<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const types::boundary_id boundary_indicator,
        const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
        const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
        const std::vector<Tensor<1,dim>> &normal_vectors,
        std::vector<double> &output
      ) const
      {
        for (unsigned int q=0; q<output.size(); ++q)
          {
            Tensor<1,dim> gradient;
            gradient[0] = 0.0;
            gradient[1] = -1.0;
            output[q] = gradient * normal_vectors[q];
          }
      }
  };

}


// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(MeltMaterial,
                                 "MeltMaterial4",
                                 "")

  ASPECT_REGISTER_POSTPROCESSOR(MMPostprocessor,
                                "MMPostprocessor",
                                "")

  ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(PressureBdry,
                                                "PressureBdry",
                                                "A fluid pressure boundary condition that prescribes the "
                                                "gradient of the fluid pressure at the boundaries as "
                                                "calculated in the analytical solution. ")

}
