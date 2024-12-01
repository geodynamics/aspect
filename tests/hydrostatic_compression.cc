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
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/simulator.h>
#include <aspect/simulator/assemblers/interface.h>
#include <aspect/gravity_model/interface.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  template <int dim>
  class TestMaterial:
    public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      virtual bool is_compressible () const
      {
        return true;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        const double alpha = -0.2;     // T0 = 100

        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            const double x = in.position[i](0);
            const double z = in.position[i](1);
            const double p = in.pressure[i];
            const double T = in.temperature[i];
            out.viscosities[i] = 1.0;
            out.densities[i] = std::exp(alpha*T)*p;
            out.thermal_expansion_coefficients[i] = -alpha;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 1.0/p;
          }
      }
  };


  template <int dim>
  class Gravity : public aspect::GravityModel::Interface<dim>
  {
    public:
      virtual Tensor<1,dim> gravity_vector (const Point<dim> &pos) const
      {
        const double x=pos[0];
        const double z=pos[1];
        const double pi = numbers::PI;

        const double T0 = 100.0;
        const double alpha = -0.2;

        Tensor<1,dim> gravity;
        gravity[0] = -(0.1e1 + 0.2e1 * alpha * (double) pi * std::cos((double) (2 * pi * x)) * std::cos((double) (pi * z))) / std::exp(alpha * (T0 + std::sin((double) (2 * pi * x)) * std::cos((double) (pi * z))));
        gravity[1] = -(0.1e1 - alpha * std::sin((double) (2 * pi * x)) * (double) pi * std::sin((double) (pi * z))) / std::exp(alpha * (T0 + std::sin((double) (2 * pi * x)) * std::cos((double) (pi * z))));


        return gravity;
      }
  };


  template <int dim>
  class RefFunction : public Function<dim>
  {
    public:
      RefFunction () : Function<dim>(dim+2) {}
      virtual void vector_value (const Point<dim>   &p,
                                 Vector<double>   &values) const
      {
        double x = p(0);
        double z = p(1);
        const double pi = numbers::PI;
        const double u1 = 1.0;
        const double w1 = 0.5;

        values[0] = -2.0*u1*pi*std::sin((double) (2*pi*x)) + u1*std::cos((double) (2*pi*x)) - w1*(pi*pi+1.0)*std::cos((double) (pi*z));
        values[1] = -w1*pi*std::sin(pi*z)+w1*std::cos(pi*z)-u1*(4*pi*pi+1)*std::cos(2*pi*x);
        values[2] = 0.5-z;
      }
  };

  /**
    * A postprocessor that evaluates the accuracy of the solution
    * by using the L2 norm.
    */
  template <int dim>
  class ConvergencePostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
  ConvergencePostprocessor<dim>::execute (TableHandler &statistics)
  {
    RefFunction<dim> ref_func;
    const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

    const unsigned int n_total_comp = this->introspection().n_components;

    Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());

    ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                        n_total_comp);
    ComponentSelectFunction<dim> comp_p(dim, n_total_comp);

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
                                       &comp_p);

    const double u_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_u, VectorTools::L2_norm);
    const double p_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p, VectorTools::L2_norm);

    std::ostringstream os;
    os << std::scientific
       << "ndofs= " << this->get_solution().size()
       << " u_L2= " << u_l2
       << " p_L2= " << p_l2
       ;

    return std::make_pair("Errors", os.str());
  }

  template <int dim>
  class ForceAssembler :
    public aspect::Assemblers::Interface<dim>, public SimulatorAccess<dim>
  {

    public:

      virtual
      void
      execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
               internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
      {
        internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
        internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

        const Introspection<dim> &introspection = this->introspection();
        const FiniteElement<dim> &fe = this->get_fe();
        const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
        const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
        const double pressure_scaling = this->get_pressure_scaling();

        for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
          {
            const unsigned int component_index_i = fe.system_to_component_index(i).first;

            if (component_index_i == introspection.component_indices.pressure
                ||
                (
                  component_index_i >= introspection.component_indices.velocities[0]
                  &&
                  component_index_i <= introspection.component_indices.velocities[dim-1]
                ))
              {
                data.local_dof_indices[i_stokes] = scratch.local_dof_indices[i];
                ++i_stokes;
              }
            ++i;
          }

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            const double x = scratch.finite_element_values.quadrature_point(q)[0];
            const double z = scratch.finite_element_values.quadrature_point(q)[1];

            const double pi = numbers::PI;
            const double u1 = 1.0;
            const double w1 = 0.5;
            const double alpha = -0.2;

            Tensor<1,dim> force_u;
            const double force_p = 0.0;

            force_u[0]= -0.32e2 / 0.3e1 * u1 * std::pow(pi, 0.3e1) * std::sin(0.2e1 * pi * x)
                        + 0.16e2 / 0.3e1 * u1 * pi * pi * std::cos(0.2e1 * pi * x)
                        - w1 * (pi * pi + 0.1e1) * pi * pi * std::cos(pi * z)
                        + (0.5e0 - z) * (0.1e1 + 0.2e1 * alpha * pi * std::cos(0.2e1 * pi * x) * std::cos(pi * z));

            force_u[1]= -0.4e1 / 0.3e1 * w1 * std::pow(pi, 0.3e1) * std::sin(pi * z)
                        + 0.4e1 / 0.3e1 * w1 * pi * pi * std::cos(pi * z)
                        - 0.4e1 * u1 * (0.4e1 * pi * pi + 0.1e1) * pi * pi * std::cos(0.2e1 * pi * x)
                        - 0.1e1 + (0.5e0 - z) * (0.1e1 - alpha * std::sin(0.2e1 * pi * x) * pi * std::sin(pi * z));


            for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
              {
                const unsigned int component_index_i = fe.system_to_component_index(i).first;

                if (component_index_i == introspection.component_indices.pressure
                    ||
                    (
                      component_index_i >= introspection.component_indices.velocities[0]
                      &&
                      component_index_i <= introspection.component_indices.velocities[dim-1]
                    ))
                  {
                    scratch.phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].value (i,q);
                    scratch.phi_p[i_stokes]   = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);

                    ++i_stokes;
                  }
                ++i;
              }

            const double JxW = scratch.finite_element_values.JxW(q);

            for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
              {
                data.local_rhs(i) += (
                                       (force_u * scratch.phi_u[i]
                                        + pressure_scaling * force_p * scratch.phi_p[i])
                                     )
                                     * JxW;
              }
          }
      }

  };


  template <int dim>
  void add_force_assembler(const SimulatorAccess<dim> &simulator_access,
                           Assemblers::Manager<dim> &assemblers)
  {

    Assemblers::Interface<dim> *assembler = new ForceAssembler<dim>();
    assemblers.stokes_system.push_back(std::unique_ptr<Assemblers::Interface<dim>>(assembler));
  }
}




template <int dim>
void signal_connector (aspect::SimulatorSignals<dim> &signals)
{
  signals.set_assemblers.connect (&aspect::add_force_assembler<dim>);
}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)





// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_GRAVITY_MODEL(Gravity,
                                "MyGravity",
                                "")

  ASPECT_REGISTER_MATERIAL_MODEL(TestMaterial,
                                 "test material",
                                 "")

  ASPECT_REGISTER_POSTPROCESSOR(ConvergencePostprocessor,
                                "error calculation",
                                "")
}
