#ifndef ASPECT_INCLUSION_H
#define ASPECT_INCLUSION_H

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
  /**
   * This is the "Pure shear/Inclusion" benchmark defined in the following paper:
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
   * The results are published in Kronbichler, Heister and Bangerth paper.
   */
  namespace InclusionBenchmark
  {
    namespace AnalyticSolutions
    {
      // based on http://geodynamics.org/hg/cs/AMR/Discontinuous_Stokes with permission
      void _Inclusion(double pos[2], double r_inclusion, double eta, double *vx, double *vy, double *p)
      {
        const double min_eta = 1.0;
        const double max_eta = eta;
        const double epsilon = 1; //strain rate
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
            //inside the inclusion
            phi=0;
            dphi=0;
            psi=-4*epsilon*(max_eta*min_eta/(min_eta+max_eta))*z;
          }
        else
          {
            //outside the inclusion
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
          FunctionInclusion (double eta_B,
                             const unsigned int n_compositional_fields)
            : Function<dim>(dim+2+n_compositional_fields), eta_B_(eta_B) {}

          void vector_value (const Point<dim>   &p,
                             Vector<double>   &values) const override
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
        InclusionBoundary() : eta_B(1e3) {}

        /**
         * Return the boundary velocity as a function of position.
         */
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id ,
                           const Point<dim> &position) const override
        {
          Assert (dim == 2, ExcNotImplemented());

          double pos[2]= {position(0),position(1)};

          Tensor<1,dim> velocity;
          double pressure;
          AnalyticSolutions::_Inclusion
          (pos,0.2,eta_B, &velocity[0], &velocity[1], &pressure);

          return velocity;
        };

      private:
        double eta_B;
    };



    /**
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
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              const Point<dim> &pos = in.position[i];
              double r2 = (pos[0]-1.0)*(pos[0]-1.0) + (pos[1]-1.0)*(pos[1]-1.0);

              out.viscosities[i] = (r2<0.2*0.2)? eta_B : 1.0;

              out.densities[i] = 0;
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.thermal_conductivities[i] = 0.0;
              out.entropy_derivative_pressure[i] = 0.0;
              out.entropy_derivative_temperature[i] = 0.0;
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
        bool is_compressible () const override
        {
          return false;
        }
        /**
         * @}
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Material model");
          {
            prm.enter_subsection("Inclusion");
            {
              prm.declare_entry ("Viscosity jump", "1e3",
                                 Patterns::Double (0.),
                                 "Viscosity in the Inclusion.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override
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



        /**
         * Returns the viscosity value in the inclusion
         */
        double get_eta_B() const
        {
          return eta_B;
        }

      private:
        /**
         * viscosity value in the inclusion
         */
        double eta_B;
    };



    /**
     * A postprocessor that evaluates the accuracy of the solution.
     *
     * The implementation of error evaluators that correspond to the
     * benchmarks defined in the paper Duretz et al. reference above.
     */
    template <int dim>
    class InclusionPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &/*statistics*/) override
        {
          AssertThrow(Plugins::plugin_type_matches<const InclusionMaterial<dim>>(this->get_material_model()),
                      ExcMessage("Postprocessor only works with the inclusion material model."));

          const InclusionMaterial<dim> &
          material_model
            = Plugins::get_plugin_as_type<const InclusionMaterial<dim>>(this->get_material_model());

          std::unique_ptr<Function<dim>> ref_func
            = std::make_unique<AnalyticSolutions::FunctionInclusion<dim>>(
                material_model.get_eta_B(),
                this->n_compositional_fields());

          const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities+2);

          Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
          Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
          Vector<float> cellwise_errors_ul2 (this->get_triangulation().n_active_cells());
          Vector<float> cellwise_errors_pl2 (this->get_triangulation().n_active_cells());

          ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                              this->get_fe().n_components());
          ComponentSelectFunction<dim> comp_p(dim,
                                              this->get_fe().n_components());

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

          // Compute stokes unknowns, do not include temperature
          types::global_dof_index n_stokes_dofs = this->introspection().system_dofs_per_block[this->introspection().block_indices.velocities];
          if (this->introspection().block_indices.velocities != this->introspection().block_indices.pressure)
            n_stokes_dofs  += this->introspection().system_dofs_per_block[this->introspection().block_indices.pressure];

          std::ostringstream os;
          os << n_stokes_dofs << "; "
             << std::scientific << u_l1
             << ", " << p_l1
             << ", " << u_l2
             << ", " << p_l2;

          return std::make_pair("DoFs; Errors u_L1, p_L1, u_L2, p_L2:", os.str());
        }
    };
  }
}
#endif
