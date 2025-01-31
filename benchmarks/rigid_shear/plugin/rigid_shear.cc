/*
  Copyright (C) 2019 - 2024 by the authors of the ASPECT code.

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

#include <aspect/gravity_model/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/postprocess/particles.h>
#include <aspect/particle/property/interface.h>
#include <aspect/particle/manager.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>



namespace aspect
{
  namespace RigidShearBenchmark
  {
    /**
     * This is the "Rigid shear" benchmark based on a suggestion in the following paper:
     * @code
     *  @Article{KKSCND97,
     *  author =       {P. E. van Keken and S. D. King and H. Schmeling and U. R. Christensen and D. Neumeister and M.-P. Doin},
     *  title =        {A comparison of methods for the modeling of thermochemical convection},
     *  journal =      {J. Geoph. Res.},
     *  year =         1997,
     *  volume =       102,
     *  pages =        {22477--22495}}
     * @endcode
     *
     * The results of the modification are published in Gassmoeller,
     * Lokavarapu, Bangerth, Puckett, 2019, "Evaluating the Accuracy of Hybrid
     * Finite Element/Particle-In-Cell Methods for Modeling Incompressible
     * Stokes Flow".
     */
    namespace AnalyticSolutions
    {
      double
      phase_function(const double t)
      {
        return std::exp(t)-1.;
      }

      template <int dim>
      double
      density(const Point<dim> &p,
              const double t,
              const bool transient)
      {
        const double pi = numbers::PI;
        const double phase = (transient == true) ? phase_function(t) : 0.0;

        return std::sin(pi*(p[0]-phase)) * std::sin(pi*p[1]) + 2.0;
      }

      /**
       * The exact solution for the Rigid Shear benchmark.
       */
      template <int dim>
      class FunctionRigidShear : public Function<dim>
      {
        public:
          FunctionRigidShear(const unsigned int n_components,
                             const bool transient)
            : Function<dim>(n_components),
              transient(transient) {}

          void vector_value(const Point<dim> &p,
                            Vector<double> &values) const override
          {
            const double pi = numbers::PI;
            const double t = this->get_time();
            const double phase = (transient == true) ? phase_function(t) : 0.0;

            values[0] = std::sin(pi*(p[0]-phase)) * std::cos(pi*p[1]);

            if (transient == true)
              values[0] += std::exp(t);

            values[1] = -std::cos(pi*(p[0]-phase)) * std::sin(pi*p[1]);
            values[2] = 2.0 * pi * std::cos(pi*(p[0]-phase)) * std::cos(pi*p[1]);
            values[4] = density(p,t,transient);
            return;
          }

        private:
          bool transient;
      };
    }



    /**
     * A material model for the stationary form of the rigid shear benchmark. All properties
     * are defined in dependence of position.
     */
    template <int dim>
    class RigidShearMaterial : public MaterialModel::Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override
        {
          const double t = (this->simulator_is_past_initialization()) ? this->get_time() : 0.0;

          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              if (use_analytical_density == true)
                out.densities[i] = AnalyticSolutions::density(in.position[i],t,use_transient_flow_solution);
              else
                out.densities[i] = in.composition[i][0];

              out.viscosities[i] = 1.0;
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0.0;
              out.thermal_conductivities[i] = 0.0;

              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                out.reaction_terms[i][c] = 0.0;
            }
        }



        bool is_compressible() const override
        {
          return false;
        }



        bool use_transient_solution() const
        {
          return use_transient_flow_solution;
        }



        static
        void
        declare_parameters (ParameterHandler &prm)
        {
          //create a global section in the parameter file for parameters
          //that describe this benchmark. note that we declare them here
          //in the material model, but other kinds of plugins (e.g., the gravity
          //model below) may also read these parameters even though they do not
          //declare them
          prm.enter_subsection("Rigid shear benchmark");
          {
            prm.declare_entry("Use analytical density", "true",
                              Patterns::Bool(),
                              "Whether to use the analytical density solution, or to look for a "
                              "compositional field named <density_field> to use as density.");

            prm.declare_entry("Use transient solution", "false",
                              Patterns::Bool(),
                              "Whether to use the transient flow solution, or to use the "
                              "default steady-state solution.");
          }
          prm.leave_subsection();
        }

        void
        parse_parameters(ParameterHandler &prm) override
        {
          prm.enter_subsection("Rigid shear benchmark");
          {
            use_analytical_density = prm.get_bool ("Use analytical density");
            use_transient_flow_solution = prm.get_bool ("Use transient solution");
          }
          prm.leave_subsection();

          this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
          this->model_dependence.density = MaterialModel::NonlinearDependence::none;
          this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
          this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
          this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
        }

        bool use_analytical_density;
        bool use_transient_flow_solution;
    };




    /**
     * Gravity model for the Rigid shear benchmark
     */
    template <int dim>
    class RigidShearGravity : public aspect::GravityModel::Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        Tensor<1,dim> gravity_vector (const Point<dim> &pos) const override
        {
          const double pi = numbers::PI;
          const double t = (this->simulator_is_past_initialization()) ? this->get_time() : 0.0;

          const RigidShearMaterial<dim> &
          material_model
            = Plugins::get_plugin_as_type<const RigidShearMaterial<dim>>(this->get_material_model());
          const bool transient = material_model.use_transient_solution();
          const double phase = (transient == true) ? AnalyticSolutions::phase_function(t) : 0.0;

          const double forcing_term = - 4.0 * pi * pi * std::cos(pi*(pos[0] - phase)) * std::sin(pi*pos[1]);

          Tensor<1,dim> gravity;
          gravity[0] =  0.0;
          gravity[1] =  forcing_term / AnalyticSolutions::density(pos,t,transient);

          return gravity;
        }
    };

    /**
     * A postprocessor that evaluates the accuracy of the solution.
     *
     * The implementation of error evaluators that correspond to the
     * benchmarks defined in the paper Gassmoeller et al. referenced above.
     */
    template <int dim>
    class RigidShearPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate error output for velocity, pressure, and density.
         */
        std::pair<std::string, std::string>
        execute(TableHandler &statistics) override
        {
          const RigidShearMaterial<dim> &
          material_model
            = Plugins::get_plugin_as_type<const RigidShearMaterial<dim>>(this->get_material_model());

          AnalyticSolutions::FunctionRigidShear<dim> ref_func(this->introspection().n_components,
                                                              material_model.use_transient_solution());
          ref_func.set_time(this->get_time());

          AssertThrow(Plugins::plugin_type_matches<const RigidShearMaterial<dim>>(this->get_material_model()),
                      ExcMessage(
                        "Postprocessor RigidShearPostprocessor only works with the material model RigidShearMaterial."));

          const QGauss<dim> quadrature_formula(this->introspection().polynomial_degree.velocities + 2);

          Vector<double> cellwise_errors_u(this->get_triangulation().n_active_cells());
          Vector<double> cellwise_errors_p(this->get_triangulation().n_active_cells());
          Vector<double> cellwise_errors_ul2(this->get_triangulation().n_active_cells());
          Vector<double> cellwise_errors_pl2(this->get_triangulation().n_active_cells());
          Vector<double> cellwise_errors_rhol2(this->get_triangulation().n_active_cells());

          ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0, dim),
                                              this->get_fe().n_components());
          ComponentSelectFunction<dim> comp_p(dim,
                                              this->get_fe().n_components());
          ComponentSelectFunction<dim> comp_rho(dim+2,
                                                this->get_fe().n_components());

          VectorTools::integrate_difference(this->get_mapping(), this->get_dof_handler(),
                                            this->get_solution(),
                                            ref_func,
                                            cellwise_errors_u,
                                            quadrature_formula,
                                            VectorTools::L1_norm,
                                            &comp_u);
          VectorTools::integrate_difference(this->get_mapping(), this->get_dof_handler(),
                                            this->get_solution(),
                                            ref_func,
                                            cellwise_errors_p,
                                            quadrature_formula,
                                            VectorTools::L1_norm,
                                            &comp_p);
          VectorTools::integrate_difference(this->get_mapping(), this->get_dof_handler(),
                                            this->get_solution(),
                                            ref_func,
                                            cellwise_errors_ul2,
                                            quadrature_formula,
                                            VectorTools::L2_norm,
                                            &comp_u);
          VectorTools::integrate_difference(this->get_mapping(), this->get_dof_handler(),
                                            this->get_solution(),
                                            ref_func,
                                            cellwise_errors_pl2,
                                            quadrature_formula,
                                            VectorTools::L2_norm,
                                            &comp_p);
          VectorTools::integrate_difference(this->get_mapping(), this->get_dof_handler(),
                                            this->get_solution(),
                                            ref_func,
                                            cellwise_errors_rhol2,
                                            quadrature_formula,
                                            VectorTools::L2_norm,
                                            &comp_rho);

          const double u_l1 = Utilities::MPI::sum(cellwise_errors_u.l1_norm(), this->get_mpi_communicator());
          const double p_l1 = Utilities::MPI::sum(cellwise_errors_p.l1_norm(), this->get_mpi_communicator());
          const double u_l2 = std::sqrt(
                                Utilities::MPI::sum(cellwise_errors_ul2.norm_sqr(), this->get_mpi_communicator()));
          const double p_l2 = std::sqrt(
                                Utilities::MPI::sum(cellwise_errors_pl2.norm_sqr(), this->get_mpi_communicator()));
          const double rho_l2 = std::sqrt(
                                  Utilities::MPI::sum(cellwise_errors_rhol2.norm_sqr(), this->get_mpi_communicator()));

          statistics.add_value ("u_L2",
                                u_l2);
          statistics.set_precision ("u_L2", 14);
          statistics.set_scientific ("u_L2", true);

          statistics.add_value ("p_L2",
                                p_l2);
          statistics.set_precision ("p_L2", 14);
          statistics.set_scientific ("p_L2", true);

          statistics.add_value ("rho_L2",
                                rho_l2);
          statistics.set_precision ("rho_L2", 14);
          statistics.set_scientific ("rho_L2", true);

          std::ostringstream os;
          os << std::scientific << u_l1
             << ", " << p_l1
             << ", " << u_l2
             << ", " << p_l2
             << ", " << rho_l2;

          return std::make_pair("Errors u_L1, p_L1, u_L2, p_L2, rho_L2:", os.str());
        }
    };
  }
}

// explicit instantiations
namespace aspect
{
  namespace RigidShearBenchmark
  {
    ASPECT_REGISTER_GRAVITY_MODEL(RigidShearGravity,
                                  "rigid shear",
                                  "A gravity model.")

    ASPECT_REGISTER_MATERIAL_MODEL(RigidShearMaterial,
                                   "rigid shear",
                                   "A material model for the stationary form of the rigid shear benchmark. All properties "
                                   "are defined in dependence of position.")

    ASPECT_REGISTER_POSTPROCESSOR(RigidShearPostprocessor,
                                  "rigid shear",
                                  "The implementation of error evaluators that correspond to the "
                                  "benchmarks defined in the paper Gassmoeller et al. referenced above.")
  }
}
