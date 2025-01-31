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

namespace aspect
{
  namespace Benchmark
  {
    /**
     * u = cos(y), sin(x)+xy
     * p = 2/3 eta x
     * grad p = 2/3 eta
     * R = div u = x
     * => grad p and compressibility term cancel
     */

    namespace AnalyticSolutions
    {
      Tensor<1,2>
      velocity (const Point<2> &pos,
                const double /*eta*/)
      {

        const double x = pos[0];
        const double y = pos[1];

        return Point<2> (std::cos(y), std::sin(x)+x*y);
      }

      double
      pressure (const Point<2> &pos,
                const double eta)
      {
        const double x = pos[0];
        //        const double y = pos[1];

        return 2./3.*eta*(x-1.0);
      }


      /**
       * The exact solution for the benchmark.
       */
      template <int dim>
      class Exact : public Function<dim>
      {
        public:
          Exact (const double eta)
            :
            Function<dim>(dim+2),
            eta(eta)
          {}

          virtual void vector_value (const Point<dim>   &pos,
                                     Vector<double>   &values) const
          {
            Assert (dim == 2, ExcNotImplemented());
            Assert (values.size() >= 3, ExcInternalError());

            const Point<2> p (pos[0], pos[1]);

            const Tensor<1,2> v = AnalyticSolutions::velocity (p, eta);
            values[0] = v[0];
            values[1] = v[1];

            values[2] = AnalyticSolutions::pressure (p, eta);
          }

        private:
          const double eta;
      };
    }



    template <int dim>
    class MyBoundary : public BoundaryVelocity::Interface<dim>
    {
      public:
        /**
         * Constructor.
         */
        MyBoundary();

        /**
         * Return the boundary velocity as a function of position.
         */
        virtual
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const;

      private:
        const double eta;
    };

    template <int dim>
    MyBoundary<dim>::MyBoundary ()
      :
      eta (0)
    {}



    template <>
    Tensor<1,2>
    MyBoundary<2>::
    boundary_velocity (const types::boundary_id ,
                       const Point<2> &p) const
    {
      return AnalyticSolutions::velocity (p, eta);
    }



    template <>
    Tensor<1,3>
    MyBoundary<3>::
    boundary_velocity (const types::boundary_id ,
                       const Point<3> &p) const
    {
      Assert (false, ExcNotImplemented());
      return Tensor<1,3>();
    }




    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MyMaterial : public MaterialModel::Interface<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          MaterialModel::PrescribedPlasticDilation<dim>
          *prescribed_dilation = out.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim>>();

          MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
          *force = out.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>>();

          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              const Point<dim> &p = in.position[i];

              const double x=p[0];
              const double y=p[1];

              out.viscosities[i] = eta;
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

              if (force)
                {
                  force->rhs_u[i][0] = -eta*(1-std::cos(y));
                  force->rhs_u[i][1] = -eta*(-std::sin(x));
                  force->rhs_p[i] = 0.;
                }
              if (prescribed_dilation)
                {
                  prescribed_dilation->dilation[i] = x;
                }

            }

        }

        virtual bool is_compressible () const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        double get_eta () const;

        /**
         * viscosity value in the inclusion
         */
        double eta;
    };



    template <int dim>
    double
    MyMaterial<dim>::
    get_eta () const
    {
      return eta;
    }

    template <int dim>
    bool
    MyMaterial<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    MyMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      // create a global section in the parameter file for parameters
      // that describe this benchmark. note that we declare them here
      // in the material model, but other kinds of plugins (e.g., the gravity
      // model below) may also read these parameters even though they do not
      // declare them
      prm.enter_subsection("My benchmark");

      {
        prm.declare_entry("Viscosity parameter", "10",
                          Patterns::Double (0),
                          "Viscosity the benchmark.");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    MyMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("My benchmark");
      {
        eta = prm.get_double ("Viscosity parameter");
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
      * A postprocessor that evaluates the accuracy of the solution.
      *
      * The implementation of error evaluators that correspond to the
      * benchmarks defined in the paper Duretz et al. reference above.
      */
    template <int dim>
    class MyPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
    MyPostprocessor<dim>::execute (TableHandler &statistics)
    {
      std::unique_ptr<Function<dim>> ref_func;
      {
        const MyMaterial<dim> *
        material_model
          = dynamic_cast<const MyMaterial<dim> *>(&this->get_material_model());

        ref_func = std::make_unique<AnalyticSolutions::Exact<dim>>(material_model->get_eta());
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

      const double u_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_ul2, VectorTools::L2_norm);
      const double p_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_pl2, VectorTools::L2_norm);

      std::ostringstream os;
      os << std::scientific << u_l2 << ", " << p_l2;

      return std::make_pair("Errors u_L2, p_L2:", os.str());
    }

  }
}



// explicit instantiations
namespace aspect
{
  namespace Benchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MyMaterial,
                                   "MyMaterial",
                                   "A material model that corresponds to the `My' benchmark. "
                                   "See the manual for more information.")

    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(MyBoundary,
                                            "MyBoundary",
                                            "Implementation of the velocity boundary conditions for the "
                                            "`My' benchmark. See the manual for more information about this "
                                            "benchmark.")

    ASPECT_REGISTER_POSTPROCESSOR(MyPostprocessor,
                                  "MyPostprocessor",
                                  "A postprocessor that compares the solution of the `My' benchmark "
                                  "with the one computed by ASPECT "
                                  "and reports the error. See the manual for more information.")
  }
}
