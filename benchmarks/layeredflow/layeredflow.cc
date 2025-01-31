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
#include <aspect/simulator.h>
#include <aspect/material_model/simple.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/utilities.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace LayeredFlowBenchmark
  {
    namespace AnalyticSolutions
    {

      Tensor<1,2>
      LayeredFlow_velocity (const Point<2> &pos,
                            const double beta,
                            const double epsilon)
      {
        const double y0 = 1./3.;
        const double yplus=(1.+y0)/beta;
        const double yminus=(1.-y0)/beta;
        const double C1 = 2.*numbers::PI /
                          ( beta*std::log(beta*beta+Utilities::fixed_power<2>(1.+y0))-2.*(1.+y0)*std::atan(yplus)
                            - beta*std::log(beta*beta+Utilities::fixed_power<2>(1.-y0))+2.*(1.-y0)*std::atan(yminus)
                            + 2.*numbers::PI*(1.+2.*epsilon) );

        const double C2 = ( beta*std::log( beta*beta+Utilities::fixed_power<2>(1+y0) )- 2.*(1.+y0)*std::atan(yplus)
                            + numbers::PI*(1.+2.*epsilon) )*C1 ;
        const double y = pos[1]-1.0;
        const double v_x = (-beta*C1*std::log(beta*beta+Utilities::fixed_power<2>(y-y0))+2.*(y-y0)*C1*std::atan((y-y0)/beta)
                            + numbers::PI*(1.+2.*epsilon)*y*C1+C2)/(2.*numbers::PI);
        const double v_y = 0;
        return Point<2> (v_x,v_y);
      }

      double
      LayeredFlow_pressure (const Point<2> &,
                            const double &,
                            const double)
      {
        return 0;
      }


      template <int dim>
      class FunctionLayeredFlow : public Function<dim>
      {
        public:
          FunctionLayeredFlow (const double beta, const double epsilon)
            :
            Function<dim>(dim+2),
            beta_(beta),
            epsilon_(epsilon)
          {}

          void vector_value (const Point<dim>   &pos,
                             Vector<double>   &values) const override
          {
            Assert (dim == 2, ExcNotImplemented());
            Assert (values.size() >= 4, ExcInternalError());

            const Point<2> p (pos[0], pos[1]);

            const Tensor<1,2> v = AnalyticSolutions::LayeredFlow_velocity (p, beta_, epsilon_);
            values[0] = v[0];
            values[1] = v[1];
            values[2] = AnalyticSolutions::LayeredFlow_pressure (p, beta_, epsilon_);
          }

        private:
          const double beta_;
          const double epsilon_;
      };
    }


    template <int dim>
    class LayeredFlowBoundary : public BoundaryVelocity::Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        LayeredFlowBoundary();

        /**
         * Return the boundary velocity as a function of position.
         */
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id ,
                           const Point<dim> &position) const override;

      private:
        const double beta;
        const double epsilon;
    };


    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     * @ingroup MaterialModels
     */
    template <int dim>
    class LayeredFlowMaterial : public MaterialModel::Interface<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              const Point<dim> &pos = in.position[i];

              const double y = pos[1]-1.0;
              const double y0 = 1./3.;
              out.viscosities[i] = 1.0/(std::atan((y-y0)/beta)/numbers::PI + 0.5 + epsilon);
              out.densities[i] = 1;
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.thermal_conductivities[i] = 0.0;
            }
        }
        /**
         * @}
         */

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
        double get_epsilon() const;

        /**
         * viscosity value in the inclusion
         */
        double beta;
        double epsilon;
    };


    template <int dim>
    bool
    LayeredFlowMaterial<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    LayeredFlowMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      //create a global section in the parameter file for parameters
      //that describe this benchmark. note that we declare them here
      //in the material model, but other kinds of plugins (e.g., the gravity
      //model below) may also read these parameters even though they do not
      //declare them
      prm.enter_subsection("LayeredFlow benchmark");
      {
        prm.declare_entry("Viscosity parameter beta", "1",
                          Patterns::Double (),
                          "Viscosity parameter beta in the LayeredFlow benchmark.");
        prm.declare_entry("Viscosity parameter epsilon", "0.1",
                          Patterns::Double (),
                          "Viscosity parameter epsilon in the LayeredFlow benchmark.");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    LayeredFlowMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("LayeredFlow benchmark");
      {
        beta = prm.get_double ("Viscosity parameter beta");
        epsilon = prm.get_double ("Viscosity parameter epsilon");
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
    LayeredFlowMaterial<dim>::get_beta() const
    {
      return beta;
    }

    template <int dim>
    double
    LayeredFlowMaterial<dim>::get_epsilon() const
    {
      return epsilon;
    }

    template <int dim>
    LayeredFlowBoundary<dim>::LayeredFlowBoundary ()
      :
      beta (0),
      epsilon (0)
    {}



    template <>
    Tensor<1,2>
    LayeredFlowBoundary<2>::
    boundary_velocity (const types::boundary_id ,
                       const Point<2> &p) const
    {
      const LayeredFlowMaterial<2> &
      material_model
        = Plugins::get_plugin_as_type<const LayeredFlowMaterial<2>>(this->get_material_model());

      return AnalyticSolutions::LayeredFlow_velocity (p, material_model.get_beta(),
                                                      material_model.get_epsilon());
    }



    template <>
    Tensor<1,3>
    LayeredFlowBoundary<3>::
    boundary_velocity (const types::boundary_id ,
                       const Point<3> &) const
    {
      Assert (false, ExcNotImplemented());
      return Tensor<1,3>();
    }


    /**
      * A postprocessor that evaluates the accuracy of the solution.
      *
      * The implementation of error evaluators that correspond to the
      * benchmarks defined in the manual.
      */
    template <int dim>
    class LayeredFlowPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
    LayeredFlowPostprocessor<dim>::execute (TableHandler &)
    {
      std::unique_ptr<Function<dim>> ref_func;
      {
        const LayeredFlowMaterial<dim> &
        material_model
          = Plugins::get_plugin_as_type<const LayeredFlowMaterial<dim>>(this->get_material_model());

        ref_func = std::make_unique<AnalyticSolutions::FunctionLayeredFlow<dim>>(material_model.get_beta(),
                                                                                  material_model.get_epsilon());
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

      const double u_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_ul2, VectorTools::L2_norm);
      const double p_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_pl2, VectorTools::L2_norm);

      std::ostringstream os;
      os << std::scientific <<  u_l2
         << ", " << p_l2;

      return std::make_pair("Errors u_L2, p_L2:", os.str());
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace LayeredFlowBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(LayeredFlowMaterial,
                                   "LayeredFlowMaterial",
                                   "A material model that corresponds to the `LayeredFlow' benchmark. "
                                   "See the manual for more information.")

    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(LayeredFlowBoundary,
                                            "LayeredFlowBoundary",
                                            "Implementation of the velocity boundary conditions for the "
                                            "`LayeredFlow' benchmark. See the manual for more information about this "
                                            "benchmark.")

    ASPECT_REGISTER_POSTPROCESSOR(LayeredFlowPostprocessor,
                                  "LayeredFlowPostprocessor",
                                  "A postprocessor that compares the solution of the `LayeredFlow' benchmark "
                                  "with the one computed by ASPECT "
                                  "and reports the error. See the manual for more information.")
  }
}
