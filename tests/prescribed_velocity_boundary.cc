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
  namespace InclusionBenchmark
  {
    namespace AnalyticSolutions
    {
      // based on http://geodynamics.org/hg/cs/AMR/Discontinuous_Stokes with permission
      void _Inclusion(double pos[2], double r_inclusion, double eta, double *vx, double *vy, double *p)
      {
        const double min_eta = 1.0;
        const double max_eta = eta;
        const double epsilon = 1; // strain rate
        const double A(min_eta*(max_eta-min_eta)/(max_eta+min_eta));
        std::complex<double> phi, psi, dphi;
        const double offset[2]= {1.0, 1.0};
        double r2_inclusion = r_inclusion * r_inclusion;

        double x = pos[0]-offset[0];
        double y = pos[1]-offset[1];
        double r2 = x*x+y*y;

        std::complex<double> z(x,y);
        if (r2<r2_inclusion)
          {
            // inside the inclusion
            phi=0;
            dphi=0;
            psi=-4*epsilon*(max_eta*min_eta/(min_eta+max_eta))*z;
          }
        else
          {
            // outside the inclusion
            phi=-2*epsilon*A*r2_inclusion/z;
            dphi=-phi/z;
            psi=-2*epsilon*(min_eta*z+A*r2_inclusion*r2_inclusion/(z*z*z));
          }
        double visc = (r2<r2_inclusion)? max_eta : 1.0;
        std::complex<double> v = (phi - z*conj(dphi) - conj(psi))/(2.0*visc);
        *vx=v.real();
        *vy=v.imag();
        *p=-2*epsilon*dphi.real();
      }





      /**
       * The exact solution for the Inclusion benchmark.
       */
      template <int dim>
      class FunctionInclusion : public Function<dim>
      {
        public:
          FunctionInclusion (double eta_B) : Function<dim>(dim+2), eta_B_(eta_B) {}
          virtual void vector_value (const Point<dim>   &p,
                                     Vector<double>   &values) const
          {
            double pos[2]= {p(0),p(1)};
            AnalyticSolutions::_Inclusion
            (pos,0.2,eta_B_, &values[0], &values[1], &values[2]);
          }

        private:
          double eta_B_;
      };
    }



    template <int dim>
    class InclusionBoundary : public BoundaryVelocity::Interface<dim>
    {
      public:
        /**
         * Constructor.
         */
        InclusionBoundary(): eta_B (1e3)
        {}


        /**
         * Return the boundary velocity as a function of position.
         */
        virtual
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_id,
                           const Point<dim> &position) const;

        virtual void initialize ()
        { }


      private:

        double eta_B;
    };


    template <int dim>
    Tensor<1,dim>
    InclusionBoundary<dim>::
    boundary_velocity (const types::boundary_id boundary_id,
                       const Point<dim> &p) const
    {
      Assert (dim == 2, ExcNotImplemented());

      std::cout << "Boundary_velocity called with boundary_id ="
                << static_cast<int>(boundary_id) << std::endl;

      return Tensor<1,dim>();
      double pos[2]= {p(0),p(1)};

      Tensor<1,dim> velocity;
      double pressure;
      AnalyticSolutions::_Inclusion
      (pos,0.2,eta_B, &velocity[0], &velocity[1], &pressure);

      return velocity;
    }




    /**
     * A material model that describes the "Pure shear/Inclusion" benchmark
     * of the paper cited in the documentation of the DuretzEtAl namespace.
     *
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class InclusionMaterial : public MaterialModel::Interface<dim>
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
              const Point<dim> &pos = in.position[i];
              const double r2 = ( pos[0] - 1.0 ) * ( pos[0] - 1.0 ) + ( pos[1] - 1.0 ) * ( pos[1] - 1.0 );
              out.viscosities[i] = ( r2 < 0.2 * 0.2 ) ? eta_B : 1.0;
              out.densities[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.compressibilities[i] = 0;
              out.thermal_conductivities[i] = 0.0;
              out.specific_heat[i] = 0;
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
        double get_eta_B() const;

      private:
        /**
         * viscosity value in the inclusion
         */
        double eta_B;
    };

    template <int dim>
    bool
    InclusionMaterial<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    InclusionMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Inclusion");
        {
          prm.declare_entry ("Viscosity jump", "1e3",
                             Patterns::Double (0),
                             "Viscosity in the Inclusion.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    InclusionMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Inclusion");
        {
          eta_B = prm.get_double ("Viscosity jump");
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

    template <int dim>
    double
    InclusionMaterial<dim>::get_eta_B() const
    {
      return eta_B;
    }



    /**
      * A postprocessor that evaluates the accuracy of the solution of the
      * aspect::MaterialModel::DuretzEtAl::Inclusion material models.
      *
      * The implementation of error evaluators that correspond to the
      * benchmarks defined in the following paper:
      * @code
      *  @Article{DMGT11,
      *    author =       {T. Duretz and D. A. May and T. V. Gerya and P. J. Tackley},
      *    title =        {Discretization errors and free surface stabilization in the
      *                  finite difference and marker-in-cell method for applied
      *                  geodynamics: {A} numerical study},
      *    journal =      {Geochemistry Geophysics Geosystems},
      *    year =         2011,
      *    volume =       12,
      *    pages =        {Q07004/1--26}}
      * @endcode
      *
      * @note While this paper summarizes the benchmarks used here, some of the
      * benchmarks actually originate in earlier papers. For the original
      * references, see the bibliography of the paper above.
      * @ingroup Postprocessing
      */
    template <int dim>
    class InclusionPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
    InclusionPostprocessor<dim>::execute (TableHandler &statistics)
    {
      std::unique_ptr<Function<dim>> ref_func;
      if (dynamic_cast<const InclusionMaterial<dim> *>(&this->get_material_model()) != nullptr)
        {
          const InclusionMaterial<dim> *
          material_model
            = dynamic_cast<const InclusionMaterial<dim> *>(&this->get_material_model());

          ref_func = std::make_unique<AnalyticSolutions::FunctionInclusion<dim>>(material_model->get_eta_B());
        }
      else
        {
          AssertThrow(false,
                      ExcMessage("Postprocessor DuretzEtAl only works with the material model SolCx, SolKz, and Inclusion."));
        }

      const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities +2);

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
      os << std::scientific << u_l1
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
  namespace InclusionBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(InclusionMaterial,
                                   "InclusionMaterial",
                                   "A material model that corresponds to the 'Inclusion' benchmark "
                                   "defined in Duretz et al., G-Cubed, 2011.")

    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(InclusionBoundary,
                                            "InclusionBoundary",
                                            "Implementation of the velocity boundary conditions for the "
                                            "``inclusion'' benchmark. See the manual and the Kronbichler, Heister "
                                            "and Bangerth paper on ASPECT for more information about this "
                                            "benchmark.")

    ASPECT_REGISTER_POSTPROCESSOR(InclusionPostprocessor,
                                  "InclusionPostprocessor",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "the Duretz et al., G-Cubed, 2011, paper with the one computed by ASPECT "
                                  "and reports the error. Specifically, it can compute the errors for "
                                  "the SolCx, SolKz and inclusion benchmarks. The postprocessor inquires "
                                  "which material model is currently being used and adjusts "
                                  "which exact solution to use accordingly.")
  }
}
