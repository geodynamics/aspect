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

// Check the Donea-Huerta benchmark using the Q1 x Q1 element

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
  namespace DoneaHuertaBenchmark
  {
    namespace AnalyticSolutions
    {

      Tensor<1,2>
      DoneaHuerta_velocity (const Point<2> &pos,
                            const double /*beta*/)
      {
        const double x = pos[0];
        const double y = pos[1];
        const double v_x =  x*x*(1.-x)*(1.-x)*(2.*y-6.*y*y+4.*y*y*y) ;
        const double v_y = -y*y*(1.-y)*(1.-y)*(2.*x-6.*x*x+4.*x*x*x) ;
        return Point<2> (v_x,v_y);
      }

      double
      DoneaHuerta_pressure (const Point<2> &pos,
                            const double /*beta*/)
      {
        const double x = pos[0];
        return x*(1.-x)-1./6.;
      }


      template <int dim>
      class FunctionDoneaHuerta : public Function<dim>
      {
        public:
          FunctionDoneaHuerta (const double beta)
            :
            Function<dim>(dim+2),
            beta_(beta)
          {}

          virtual void vector_value (const Point<dim>   &pos,
                                     Vector<double>   &values) const
          {
            Assert (dim == 2, ExcNotImplemented());
            Assert (values.size() >= 4, ExcInternalError());

            const Point<2> p (pos[0], pos[1]);

            const Tensor<1,2> v = AnalyticSolutions::DoneaHuerta_velocity (p, beta_);
            values[0] = v[0];
            values[1] = v[1];
            values[2] = AnalyticSolutions::DoneaHuerta_pressure (p, beta_);
          }

        private:
          const double beta_;
      };
    }



    template <int dim>
    class DoneaHuertaBoundary : public BoundaryVelocity::Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        DoneaHuertaBoundary();

        /**
         * Return the boundary velocity as a function of position.
         */
        virtual
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id ,
                           const Point<dim> &position) const;

      private:
        const double beta;
    };


    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     _r*
     * @ingroup MaterialModels
     */
    template <int dim>
    class DoneaHuertaMaterial : public MaterialModel::Interface<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              out.viscosities[i] = 1.;
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
        virtual bool is_compressible () const;
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
        virtual
        void
        parse_parameters (ParameterHandler &prm);



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
    DoneaHuertaMaterial<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    DoneaHuertaMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      //create a global section in the parameter file for parameters
      //that describe this benchmark. note that we declare them here
      //in the material model, but other kinds of plugins (e.g., the gravity
      //model below) may also read these parameters even though they do not
      //declare them
      prm.enter_subsection("DoneaHuerta benchmark");
      {
        prm.declare_entry("Viscosity parameter", "1",
                          Patterns::Double (),
                          "Viscosity parameter beta in the DoneaHuerta benchmark.");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    DoneaHuertaMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("DoneaHuerta benchmark");
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
    DoneaHuertaMaterial<dim>::get_beta() const
    {
      return beta;
    }

    template <int dim>
    DoneaHuertaBoundary<dim>::DoneaHuertaBoundary ()
      :
      beta (0)
    {}



    template <>
    Tensor<1,2>
    DoneaHuertaBoundary<2>::
    boundary_velocity (const types::boundary_id ,
                       const Point<2> &p) const
    {
      const DoneaHuertaMaterial<2> *
      material_model
        = dynamic_cast<const DoneaHuertaMaterial<2> *>(&this->get_material_model());
      return AnalyticSolutions::DoneaHuerta_velocity (p, material_model->get_beta());
    }



    template <>
    Tensor<1,3>
    DoneaHuertaBoundary<3>::
    boundary_velocity (const types::boundary_id ,
                       const Point<3> &) const
    {
      Assert (false, ExcNotImplemented());
      return Tensor<1,3>();
    }




    /**
     *gravity model for the DoneaHuerta benchmark
    */

    template <int dim>
    class DoneaHuertaGravity : public aspect::GravityModel::Interface<dim>
    {
      public:
        virtual Tensor<1,dim> gravity_vector (const Point<dim> &pos) const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
        */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        double beta;
    };


    template <int dim>
    Tensor<1,dim>
    DoneaHuertaGravity<dim>::
    gravity_vector(const Point<dim> &pos) const
    {
      const double x=pos[0];
      const double y=pos[1];
      Tensor<1,dim> g;
      g[0]= ( (12.-24.*y)*x*x*x*x + (-24.+48.*y)*x*x*x + (-48.*y+72.*y*y-48.*y*y*y+12.)*x*x
              + (-2.+24.*y-72.*y*y+48.*y*y*y)*x + 1.-4.*y+12.*y*y-8.*y*y*y );

      g[1]= ( (8.-48.*y+48.*y*y)*x*x*x + (-12.+72.*y-72*y*y)*x*x
              +  (4.-24.*y+48.*y*y-48.*y*y*y+24.*y*y*y*y)*x - 12.*y*y + 24.*y*y*y -12.*y*y*y*y);
      return g;
    }

    template <int dim>
    void
    DoneaHuertaGravity<dim>::declare_parameters (ParameterHandler &/*prm*/)
    {
      //nothing to declare here. This plugin will however, read parameters
      //declared by the material model in the "Burstedde benchmark" section
    }

    template <int dim>
    void
    DoneaHuertaGravity<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("DoneaHuerta benchmark");
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
    class DoneaHuertaPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

      private:
        /**
         * Calculate the L2 dynamic topography error.
         */
        //double
        //compute_dynamic_topography_error() const;
    };

    template <int dim>
    std::pair<std::string,std::string>
    DoneaHuertaPostprocessor<dim>::execute (TableHandler &)
    {
      std::unique_ptr<Function<dim>> ref_func;
      {
        const DoneaHuertaMaterial<dim> *
        material_model
          = dynamic_cast<const DoneaHuertaMaterial<dim> *>(&this->get_material_model());

        ref_func.reset (new AnalyticSolutions::FunctionDoneaHuerta<dim>(material_model->get_beta()));
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
  namespace DoneaHuertaBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(DoneaHuertaMaterial,
                                   "DoneaHuertaMaterial",
                                   "A material model that corresponds to the `DoneaHuerta' benchmark. "
                                   "See the manual for more information.")

    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(DoneaHuertaBoundary,
                                            "DoneaHuertaBoundary",
                                            "Implementation of the velocity boundary conditions for the "
                                            "`DoneaHuerta' benchmark. See the manual for more information about this "
                                            "benchmark.")

    ASPECT_REGISTER_POSTPROCESSOR(DoneaHuertaPostprocessor,
                                  "DoneaHuertaPostprocessor",
                                  "A postprocessor that compares the solution of the `DoneaHuerta' benchmark "
                                  "with the one computed by ASPECT "
                                  "and reports the error. See the manual for more information.")

    ASPECT_REGISTER_GRAVITY_MODEL(DoneaHuertaGravity,
                                  "DoneaHuertaGravity",
                                  "A gravity model in corresponding to the `DoneaHuerta' benchmark. "
                                  "See the manual for more information.")
  }
}
