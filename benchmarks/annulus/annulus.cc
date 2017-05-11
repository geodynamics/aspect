#include <aspect/material_model/simple.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/gravity_model/interface.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  /**
   * This is the "annulus" benchmark defined in the following paper:
   * @code
   *  @Article{thph17,
   *    author =       {C. Thieulot and E.G.P. Puckett and H. Lokavarapu},
   *    title =        {Stokes flow in an annulus: analytical solution and numerical benchmark},
   *    journal =      {xxx},
   *    year =         2017,
   *    volume =       xxx,
   *    number =       {x},
   *    publisher =    {xxx},
   *    pages =        {xxx--xxx}}
   * @endcode
   *
   */
  namespace AnnulusBenchmark
  {
    using namespace dealii;

    namespace AnalyticSolutions
    {
      Tensor<1,2>
      Annulus_velocity (const Point<2> &pos,
                        const double k)
      {
        const double x = pos[0];
        const double y = pos[1];
        const double r = sqrt(x*x+y*y);
        const double theta = atan2(y,x);
        const double A=2.0, B=-3/log(2), C=-1;
        const double f_r = A*r + B/r;
        const double g_r = A*r/2 + B*log(r)/r + C/r;
        const double v_r = g_r*k*sin(k*theta);
        const double v_theta = f_r*cos(k*theta);
        const double v_x = cos(theta)*v_r - sin(theta)*v_theta;
        const double v_y = sin(theta)*v_r + cos(theta)*v_theta;
        return Point<2> (v_x,v_y);
      }

      double
      Annulus_pressure (const Point<2> &pos,
                        const double k)
      {
        const double x = pos[0];
        const double y = pos[1];
        const double r=sqrt(x*x+y*y);
        const double theta = atan2(y,x);
        const double A=2.0, B=-3/log(2), C=-1;
        const double f_r = 2*r + B/r;
        const double g_r = A*r/2 + B*log(r)/r + C/r;
        const double h_r=(2*g_r-f_r)/r;
        const double l_r=0;
        return k*h_r*sin(k*theta)+l_r;
      }


      /**
       * The exact solution for the Annulus benchmark.
       */
      template <int dim>
      class FunctionAnnulus : public Function<dim>
      {
        public:
          FunctionAnnulus (const double beta)
            :
            Function<dim>(dim+2),
            beta_(beta)
          {}

          virtual void vector_value (const Point< dim >   &pos,
                                     Vector< double >   &values) const
          {
            Assert (dim == 3, ExcNotImplemented());
            Assert (values.size() >= 4, ExcInternalError());

            const Point<2> p (pos[0], pos[1]);

            const Tensor<1,2> v = AnalyticSolutions::Annulus_velocity (p, beta_);
            values[0] = v[0];
            values[1] = v[1];
            values[2] = AnalyticSolutions::Annulus_pressure (p, beta_);
          }

        private:
          const double beta_;
      };
    }



    template <int dim>
    class AnnulusBoundary : public BoundaryVelocity::Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        AnnulusBoundary();

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
    class AnnulusMaterial : public MaterialModel::InterfaceCompatibility<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;

        virtual double density (const double temperature,
                                const double pressure,
                                const std::vector<double> &compositional_fields,
                                const Point<dim> &position) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position) const;

        virtual double thermal_expansion_coefficient (const double      temperature,
                                                      const double      pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const;

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const std::vector<double> &compositional_fields,
                                             const Point<dim> &position) const;
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
         * solve the contuity equation as $\nabla \cdot (\rho \mathbf u)=0$
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
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;
        /**
         * @}
         */
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
    double
    AnnulusMaterial<dim>::
    viscosity (const double,
               const double,
               const std::vector<double> &,       /*composition*/
               const SymmetricTensor<2,dim> &,
               const Point<dim> &p) const
    {
      return 1.;
    }


    template <int dim>
    double
    AnnulusMaterial<dim>::
    reference_viscosity () const
    {
      return 1.;
    }


    template <int dim>
    double
    AnnulusMaterial<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return 0;
    }


    template <int dim>
    double
    AnnulusMaterial<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return 0;
    }


    template <int dim>
    double
    AnnulusMaterial<dim>::
    density (const double,
             const double,
             const std::vector<double> &, /*composition*/
             const Point<dim> &pos) const
    {
      const double k=beta;
      const double x = pos[0];
      const double y = pos[1];
      const double r=sqrt(x*x+y*y);
      const double theta = atan2(y,x);
      const double  A=2.0, B=-3/log(2), C=-1;
      const double f = A*r + B/r;
      const double f_prime = 2 - B/std::pow(r,2.0);
      const double g = A*r/2 + B*log(r)/r + C/r;
      const double g_prime = A/2 - B*log(r)/std::pow(r,2.0) + B/std::pow(r,2.0) - C/std::pow(r,2.0);
      const double g_prime_prime = -B/std::pow(r,3)*(3-2*log(r)) - 2./std::pow(r,3);
      const double N = g_prime_prime - g_prime/r - (g*(std::pow(k,2) - 1))/std::pow(r,2.0) + f/std::pow(r,2.0) + f_prime/r;
      const double density = N*k*sin(k*theta);
      return density;
    }


    template <int dim>
    double
    AnnulusMaterial<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return 0;
    }


    template <int dim>
    double
    AnnulusMaterial<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0;
    }


    template <int dim>
    bool
    AnnulusMaterial<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    AnnulusMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      //create a global section in the parameter file for parameters
      //that describe this benchmark. note that we declare them here
      //in the material model, but other kinds of plugins (e.g., the gravity
      //model below) may also read these parameters even though they do not
      //declare them
      prm.enter_subsection("Annulus benchmark");
      {
        prm.declare_entry("Viscosity parameter", "0",
                          Patterns::Double (0),
                          "Viscosity parameter k in the Annulus benchmark.");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AnnulusMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Annulus benchmark");
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
    AnnulusMaterial<dim>::get_beta() const
    {
      return beta;
    }




    template <int dim>
    AnnulusBoundary<dim>::AnnulusBoundary ()
      :
      beta (0)
    {}



    template <>
    Tensor<1,2>
    AnnulusBoundary<2>::
    boundary_velocity (const types::boundary_id ,
                       const Point<2> &p) const
    {

      const AnnulusMaterial<2> *
      material_model
        = dynamic_cast<const AnnulusMaterial<2> *>(&this->get_material_model());

      return AnalyticSolutions::Annulus_velocity (p, material_model->get_beta());
    }



    template <>
    Tensor<1,3>
    AnnulusBoundary<3>::
    boundary_velocity (const types::boundary_id ,
                       const Point<3> &p) const
    {
      Assert (false, ExcNotImplemented());
      return Tensor<1,3>();
    }




    /**
     *gravity model for the Annulus benchmark
    */

    template <int dim>
    class AnnulusGravity : public aspect::GravityModel::Interface<dim>
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


    /**
      * A postprocessor that evaluates the accuracy of the solution.
      *
      * The implementation of error evaluators that correspond to the
      * benchmarks defined in the paper Duretz et al. reference above.
      */
    template <int dim>
    class AnnulusPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
    AnnulusPostprocessor<dim>::execute (TableHandler &statistics)
    {
      std_cxx1x::shared_ptr<Function<dim> > ref_func;
      {
        const AnnulusMaterial<dim> *
        material_model
          = dynamic_cast<const AnnulusMaterial<dim> *>(&this->get_material_model());

        ref_func.reset (new AnalyticSolutions::FunctionAnnulus<dim>(material_model->get_beta()));
      }

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

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

      const double u_l1 =  Utilities::MPI::sum(cellwise_errors_u.l1_norm(),this->get_mpi_communicator());
      const double p_l1 =  Utilities::MPI::sum(cellwise_errors_p.l1_norm(),this->get_mpi_communicator());
      const double u_l2 =  std::sqrt(Utilities::MPI::sum(cellwise_errors_ul2.norm_sqr(),this->get_mpi_communicator()));
      const double p_l2 =  std::sqrt(Utilities::MPI::sum(cellwise_errors_pl2.norm_sqr(),this->get_mpi_communicator()));

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
  namespace AnnulusBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(AnnulusMaterial,
                                   "AnnulusMaterial",
                                   "A material model that corresponds to the `Annulus' benchmark. "
                                   "See the manual for more information.")

    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(AnnulusBoundary,
                                            "AnnulusBoundary",
                                            "Implementation of the velocity boundary conditions for the "
                                            "`Annulus' benchmark. See the manual for more information about this "
                                            "benchmark.")

    ASPECT_REGISTER_POSTPROCESSOR(AnnulusPostprocessor,
                                  "AnnulusPostprocessor",
                                  "A postprocessor that compares the solution of the `Annulus' benchmark "
                                  "with the one computed by ASPECT "
                                  "and reports the error. See the manual for more information.")
  }
}

