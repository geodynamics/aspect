/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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

#include <aspect/initial_composition/function.h>
#include <deal.II/base/parsed_function.h>


namespace aspect
{
  /**
   * This is the analytic solution, postprocessor, and material model setup
   * for the time dependent annular flow benchmark from Section 5 of
   *
   * Gassmoeller, Lokavarapu, Bangerth, Puckett (2019):
   * Evaluating the Accuracy of Hybrid Finite Element/Particle-In-Cell
   * Methods for Modeling Incompressible Stokes Flow. Geophys. J. Int.
   * submitted.
   *
   * It features a spherical annulus with a circular flow and a
   * time-indepedent analytical, but a time-dependent numerical solution,
   * which allows to quantify how errors in the advection equation
   * influence the accuracy of the Stokes equation. In this example
   * the benchmark uses particles to carry density, but it can be
   * used for other advection methods as well.
   */
  namespace AnalyticSolutions
  {
    template <int dim>
    void analytic_solution(
      double pos[],
      double vel[],
      double *pressure,
      double *density,
      std::shared_ptr<Functions::ParsedFunction<dim> > pressure_function,
      std::shared_ptr<Functions::ParsedFunction<dim> > velocity_function,
      std::shared_ptr<Functions::ParsedFunction<dim> > density_function)
    {
      /****************************************************************************************/
      /****************************************************************************************/
      /* Output */
      for (unsigned int i=0; i < dim; i++)
        vel[i] = velocity_function->value(Point<dim>(pos[0],pos[1]), i);

      (*pressure) = pressure_function->value(Point<dim>(pos[0], pos[1]));

      (*density) = density_function->value(Point<dim>(pos[0], pos[1]));
    }

    /**
     * The exact solution for the benchmark.
     */
    template <int dim>
    class FunctionStreamline : public Function<dim>
    {
      public:
        FunctionStreamline (std::shared_ptr<Functions::ParsedFunction<dim> > pressure,
                            std::shared_ptr<Functions::ParsedFunction<dim> > velocity,
                            std::shared_ptr<Functions::ParsedFunction<dim> > density,
                            const unsigned int n_components)
          :
          Function<dim>(n_components),
          pressure_function (pressure),
          velocity_function (velocity),
          density_function (density)
        {}

        virtual void vector_value (const Point< dim > &p,
                                   Vector< double >   &values) const
        {
          double pos[2]= {p(0),p(1)};

          AnalyticSolutions::analytic_solution<dim>
          (pos,
           &values[0], &values[2], &values[4], pressure_function, velocity_function, density_function);
        }
      private:
        std::shared_ptr<Functions::ParsedFunction<dim> > pressure_function;
        std::shared_ptr<Functions::ParsedFunction<dim> > velocity_function;
        std::shared_ptr<Functions::ParsedFunction<dim> > density_function;
    };
  }

  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class TimeDependentAnnulus : public MaterialModel::Interface<dim>
    {
      private:
        std::shared_ptr<Functions::ParsedFunction<dim> > density_function, pressure_function, velocity_function;

      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              if (use_analytic_density)
                out.densities[i] = density_function->value(in.position[i]);
              else
                out.densities[i] = in.composition[i][0];

              out.viscosities[i] = 1;
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.thermal_conductivities[i] = 0.0;
            }
        }

        virtual bool is_compressible () const
        {
          return false;
        }

        static
        void
        declare_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Material model");
          {
            prm.enter_subsection("Time dependent annulus");
            {
              prm.declare_entry ("Use analytic density", "false",
                                 Patterns::Bool(),
                                 "Whether to use the analytic density for the computations, or whatever"
                                 "density is stored in the first compositional field.");

              prm.enter_subsection("Analytical density");
              {
                prm.enter_subsection("Function");
                {
                  Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();

              prm.enter_subsection("Analytical pressure");
              {
                prm.enter_subsection("Function");
                {
                  Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();

              prm.enter_subsection("Analytical velocity");
              {
                prm.enter_subsection("Function");
                {
                  Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();

            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Material model");
          {
            prm.enter_subsection("Time dependent annulus");
            {
              use_analytic_density = prm.get_bool ("Use analytic density");

              prm.enter_subsection("Analytical density");
              {
                prm.enter_subsection("Function");
                {
                  density_function.reset(new Functions::ParsedFunction<dim>(1));
                  density_function->parse_parameters(prm);
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();

              prm.enter_subsection("Analytical pressure");
              {
                prm.enter_subsection("Function");
                {
                  pressure_function.reset(new Functions::ParsedFunction<dim>(1));
                  pressure_function->parse_parameters(prm);
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();

              prm.enter_subsection("Analytical velocity");
              {
                prm.enter_subsection("Function");
                {
                  try
                    {
                      velocity_function.reset (new Functions::ParsedFunction<dim>(dim));
                      velocity_function->parse_parameters (prm);
                    }
                  catch (...)
                    {
                      std::cerr << "ERROR: FunctionParser failed to parse\n"
                                << "\t'Analytical velocity.Function'\n"
                                << "with expression\n"
                                << "\t'" << prm.get("Function expression") << "'\n"
                                << "More information about the cause of the parse error \n"
                                << "is shown below.\n";
                      throw;
                    }
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();

          // Declare dependencies on solution variables
          this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
          this->model_dependence.density = MaterialModel::NonlinearDependence::none;
          this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
          this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
          this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
        }

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const
        {
          return 1;
        }

        /**
         * Returns the analytic solutions of this model. See the
         * corresponding member variable of this class for more information.
         */
        std::shared_ptr<Functions::ParsedFunction<dim>>  get_pressure() const
        {
          return pressure_function;
        }

        std::shared_ptr<Functions::ParsedFunction<dim>>  get_velocity() const
        {
          return velocity_function;
        }

        std::shared_ptr<Functions::ParsedFunction<dim>>  get_density() const
        {
          return density_function;
        }

      private:
        /**
         * Whether to use the analytic density for the benchmark or whatever density is
         * stored in the first compositional field.
         */
        bool use_analytic_density;
    };
  }

  namespace Postprocess
  {
    /**
     * A postprocessor that evaluates the accuracy of the solution.
     *
     */
    template <int dim>
    class TimeDependentAnnulus : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      private:
        std::shared_ptr<Function<dim> > ref_func;

      public:
        virtual
        void
        initialize ()
        {
          const MaterialModel::TimeDependentAnnulus<dim> &material_model =
            Plugins::get_plugin_as_type<const MaterialModel::TimeDependentAnnulus<dim> >(this->get_material_model());

          ref_func.reset (new AnalyticSolutions::FunctionStreamline<dim>(material_model.get_pressure(),
                                                                         material_model.get_velocity(),
                                                                         material_model.get_density(),
                                                                         1 + dim + 1 + this->n_compositional_fields()));
        }

        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics)
        {
          const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

          Vector<float> cellwise_errors_ul2 (this->get_triangulation().n_active_cells());
          Vector<float> cellwise_errors_pl2 (this->get_triangulation().n_active_cells());
          Vector<float> cellwise_errors_rhol2 (this->get_triangulation().n_active_cells());

          ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                              this->get_fe().n_components());
          ComponentSelectFunction<dim> comp_p(dim, this->get_fe().n_components());
          ComponentSelectFunction<dim> comp_rho(dim+2, this->get_fe().n_components());

          VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                             this->get_solution(),
                                             *ref_func,
                                             cellwise_errors_ul2,
                                             quadrature_formula,
                                             VectorTools::L2_norm,
                                             &comp_u);
          VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                             this->get_solution(),
                                             *ref_func,
                                             cellwise_errors_pl2,
                                             quadrature_formula,
                                             VectorTools::L2_norm,
                                             &comp_p);
          VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                             this->get_solution(),
                                             *ref_func,
                                             cellwise_errors_rhol2,
                                             quadrature_formula,
                                             VectorTools::L2_norm,
                                             &comp_rho);

          const double u_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_ul2, VectorTools::L2_norm);
          const double p_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_pl2, VectorTools::L2_norm);
          const double rho_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_rhol2, VectorTools::L2_norm);

          statistics.add_value ("Error_u_l2", u_l2);
          statistics.add_value ("Error_p_l2", p_l2);
          statistics.add_value ("Error_rho_l2", rho_l2);

          statistics.set_scientific ("Error_u_l2", true);
          statistics.set_scientific ("Error_p_l2", true);
          statistics.set_scientific ("Error_rho_l2", true);

          std::ostringstream os;
          os << std::scientific
             << u_l2
             << ", " << p_l2
             << ", " << rho_l2;

          return std::make_pair("Errors u_L2, p_L2, rho_L2", os.str());
        }
    };
  }
}
