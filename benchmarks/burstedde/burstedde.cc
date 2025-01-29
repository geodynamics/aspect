/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/postprocess/interface.h>


#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  /**
   * This is the "Burstedde" benchmark defined in the following paper:
   * @code
   *  @Article{busa13,
   *    author =       {Burstedde, Carsten and Stadler, Georg and Alisic, Laura and Wilcox, Lucas C and Tan, Eh and Gurnis, Michael and Ghattas, Omar},
   *    title =        {Large-scale adaptive mantle convection simulation},
   *    journal =      {Geophysical Journal International},
   *    year =         2013,
   *    volume =       192,
   *    number =       {3},
   *    publisher =    {Oxford University Press},
   *    pages =        {889--906}}
   * @endcode
   *
   */
  namespace BursteddeBenchmark
  {
    namespace AnalyticSolutions
    {
      Tensor<1,3>
      burstedde_velocity (const Point<3> &pos,
                          const double /*eta*/)
      {
        const double x = pos[0];
        const double y = pos[1];
        const double z = pos[2];

        // create a Point<3> (because it has a constructor that takes
        // three doubles) and return it (it automatically converts to
        // the necessary Tensor<1,3>).
        return Point<3> (x+x*x+x*y+x*x*x*y,
                         y+x*y+y*y+x*x*y*y,
                         -2.*z-3.*x*z-3.*y*z-5.*x*x*y*z);
      }

      double
      burstedde_pressure (const Point<3> &pos,
                          const double /*eta*/)
      {
        const double x = pos[0];
        const double y = pos[1];
        const double z = pos[2];

        return x*y*z+x*x*x*y*y*y*z-5./32.;
      }


      /**
       * The exact solution for the Burstedde benchmark.
       */
      template <int dim>
      class FunctionBurstedde : public Function<dim>
      {
        public:
          FunctionBurstedde (const double beta)
            :
            Function<dim>(dim+2),
            beta_(beta)
          {}

          void vector_value (const Point<dim>   &pos,
                             Vector<double>   &values) const override
          {
            Assert (dim == 3, ExcNotImplemented());
            Assert (values.size() >= 4, ExcInternalError());

            const Point<3> p (pos[0], pos[1], pos[2]);

            const Tensor<1,3> v = AnalyticSolutions::burstedde_velocity (p, beta_);
            values[0] = v[0];
            values[1] = v[1];
            values[2] = v[2];

            values[3] = AnalyticSolutions::burstedde_pressure (p, beta_);
          }

        private:
          const double beta_;
      };
    }



    template <int dim>
    class BursteddeBoundary : public BoundaryVelocity::Interface<dim>
    {
      public:
        /**
         * Constructor.
         */
        BursteddeBoundary();

        /**
         * Return the boundary velocity as a function of position.
         */
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const override;

      private:
        const double beta;
    };

    template <int dim>
    BursteddeBoundary<dim>::BursteddeBoundary ()
      :
      beta (0)
    {}



    template <>
    Tensor<1,2>
    BursteddeBoundary<2>::
    boundary_velocity (const types::boundary_id ,
                       const Point<2> &/*p*/) const
    {
      Assert (false, ExcNotImplemented());
      return Tensor<1,2>();
    }



    template <>
    Tensor<1,3>
    BursteddeBoundary<3>::
    boundary_velocity (const types::boundary_id ,
                       const Point<3> &p) const
    {
      return AnalyticSolutions::burstedde_velocity (p, beta);
    }




    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class BursteddeMaterial : public MaterialModel::Interface<dim>
    {
      public:
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              const Point<dim> &p = in.position[i];

              const double x=p[0];
              const double y=p[1];
              const double z=p[2];
              const double mu=std::exp(1. - beta * (x*(1.-x)+y*(1.-y) + z*(1.-z)));

              out.viscosities[i] = mu;
              out.thermal_conductivities[i] = 0.0;
              out.densities[i] = 1.0;

              out.thermal_expansion_coefficients[i] = 0.0;
              out.compressibilities[i] = 0.0;

              out.specific_heat[i] = 0.0;

              // Pressure derivative of entropy at the given positions.
              out.entropy_derivative_pressure[i] = 0.0;
              // Temperature derivative of entropy at the given positions.
              out.entropy_derivative_temperature[i] = 0.0;
              // Change in composition due to chemical reactions at the
              // given positions. The term reaction_terms[i][c] is the
              // change in compositional field c at point i.
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                out.reaction_terms[i][c] = 0.0;

            }

        }







        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.
         * Incompressibility does not necessarily imply that the density is
         * constant; rather, it may still depend on temperature or pressure.
         * In the current context, compressibility means whether we should
         * solve the continuity equation as $\nabla \cdot (\rho \mathbf u)=0$
         * (compressible Stokes) or as $\nabla \cdot \mathbf{u}=0$
         * (incompressible Stokes).
         */
        bool is_compressible () const override;
        /**
         * @}
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Returns the viscosity value in the inclusion
         */
        double get_beta() const;

        /**
         * viscosity value in the inclusion
         */
        double beta;
    };



    template <int dim>
    bool
    BursteddeMaterial<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    BursteddeMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      //create a global section in the parameter file for parameters
      //that describe this benchmark. note that we declare them here
      //in the material model, but other kinds of plugins (e.g., the gravity
      //model below) may also read these parameters even though they do not
      //declare them
      prm.enter_subsection("Burstedde benchmark");
      {
        prm.declare_entry("Viscosity parameter", "20",
                          Patterns::Double (0),
                          "Viscosity in the Burstedde benchmark.");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    BursteddeMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Burstedde benchmark");
      {
        beta = prm.get_double ("Viscosity parameter");
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
      this->model_dependence.density = MaterialModel::NonlinearDependence::none;
      this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
      this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
    }


    template <int dim>
    double
    BursteddeMaterial<dim>::get_beta() const
    {
      return beta;
    }

    /**
     * Gravity model for the Burstedde benchmark
     */

    template <int dim>
    class BursteddeGravity : public aspect::GravityModel::Interface<dim>
    {
      public:
        Tensor<1,dim> gravity_vector (const Point<dim> &pos) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
        */
        void
        parse_parameters (ParameterHandler &prm) override;

        double beta;
    };

    template <int dim>
    Tensor<1,dim>
    BursteddeGravity<dim>::
    gravity_vector(const Point<dim> &pos) const
    {

      const double x=pos[0];
      const double y=pos[1];
      const double z=pos[2];
      const double mu=std::exp(1. - beta * (x*(1.-x)+y*(1.-y) + z*(1.-z)));

      const double dmudx=-beta*(1.-2.*x)*mu;
      const double dmudy=-beta*(1.-2.*y)*mu;
      const double dmudz=-beta*(1.-2.*z)*mu;

      Tensor<1,dim> g;

      g[0]=((y*z+3.*Utilities::fixed_power<2>(x)*Utilities::fixed_power<3>(y)*z)- mu*(2.+6.*x*y))
           -dmudx*(2.+4.*x+2.*y+6.*Utilities::fixed_power<2>(x)*y)
           -dmudy*(x+Utilities::fixed_power<3>(x)+y+2.*x*Utilities::fixed_power<2>(y))
           -dmudz*(-3.*z-10.*x*y*z);

      g[1]=((x*z+3.*Utilities::fixed_power<3>(x)*Utilities::fixed_power<2>(y)*z)- mu*(2.+2.*Utilities::fixed_power<2>(x)+2.*Utilities::fixed_power<2>(y)))
           -dmudx*(x+Utilities::fixed_power<3>(x)+y+2.*x*Utilities::fixed_power<2>(y))
           -dmudy*(2.+2.*x+4.*y+4.*Utilities::fixed_power<2>(x)*y)
           -dmudz*(-3.*z-5.*Utilities::fixed_power<2>(x)*z);

      g[2]=((x*y+Utilities::fixed_power<3>(x)*Utilities::fixed_power<3>(y)) - mu*(-10.*y*z))
           -dmudx*(-3.*z-10.*x*y*z)
           -dmudy*(-3.*z-5.*Utilities::fixed_power<2>(x)*z)
           -dmudz*(-4.-6.*x-6.*y-10.*Utilities::fixed_power<2>(x)*y);

      return g;
    }

    template <int dim>
    void
    BursteddeGravity<dim>::declare_parameters (ParameterHandler &)
    {
      //nothing to declare here. This plugin will however, read parameters
      //declared by the material model in the "Burstedde benchmark" section
    }

    template <int dim>
    void
    BursteddeGravity<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Burstedde benchmark");
      {
        beta = prm.get_double ("Viscosity parameter");
      }
      prm.leave_subsection();
    }


    /**
      * A postprocessor that evaluates the accuracy of the solution.
      *
      * The implementation of error evaluators that correspond to the
      * benchmarks defined in the paper Duretz et al. reference above.
      */
    template <int dim>
    class BursteddePostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;
    };

    template <int dim>
    std::pair<std::string,std::string>
    BursteddePostprocessor<dim>::execute (TableHandler &)
    {
      std::unique_ptr<Function<dim>> ref_func;
      {
        const BursteddeMaterial<dim> &
        material_model
          = Plugins::get_plugin_as_type<const BursteddeMaterial<dim>>(this->get_material_model());

        ref_func = std::make_unique<AnalyticSolutions::FunctionBurstedde<dim>>(material_model.get_beta());
      }

      const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities+2);

      Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_ul2 (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_pl2 (this->get_triangulation().n_active_cells());

      ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                          dim+2);
      ComponentSelectFunction<dim> comp_p(dim, dim+2);

      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_u,
                                         quadrature_formula,
                                         VectorTools::L1_norm,
                                         &comp_u);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_p,
                                         quadrature_formula,
                                         VectorTools::L1_norm,
                                         &comp_p);
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

      const double u_l1 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_u, VectorTools::L1_norm);
      const double p_l1 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p, VectorTools::L1_norm);
      const double u_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_ul2, VectorTools::L2_norm);
      const double p_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_pl2, VectorTools::L2_norm);

      std::ostringstream os;

      os << std::scientific <<  u_l1
         << ", " << p_l1
         << ", " << u_l2
         << ", " << p_l2;

      return std::make_pair("Errors u_L1, p_L1, u_L2, p_L2:", os.str());
    }

  }
}



// explicit instantiations
namespace aspect
{
  namespace BursteddeBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(BursteddeMaterial,
                                   "BursteddeMaterial",
                                   "A material model that corresponds to the `Burstedde' benchmark. "
                                   "See the manual for more information.")

    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(BursteddeBoundary,
                                            "BursteddeBoundary",
                                            "Implementation of the velocity boundary conditions for the "
                                            "`Burstedde' benchmark. See the manual for more information about this "
                                            "benchmark.")

    ASPECT_REGISTER_POSTPROCESSOR(BursteddePostprocessor,
                                  "BursteddePostprocessor",
                                  "A postprocessor that compares the solution of the `Burstedde' benchmark "
                                  "with the one computed by ASPECT "
                                  "and reports the error. See the manual for more information.")
    ASPECT_REGISTER_GRAVITY_MODEL(BursteddeGravity,
                                  "BursteddeGravity",
                                  "A gravity model in corresponding to the `Burstedde' benchmark. "
                                  "See the manual for more information.")
  }
}
