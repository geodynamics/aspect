#include <aspect/material_model/simple.h>
#include <aspect/velocity_boundary_conditions/interface.h>
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
    using namespace dealii;

    namespace AnalyticSolutions
    {
      Tensor<1,3>
      burstedde_velocity (const Point<3> &pos,
                          const double eta)
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
                          const double eta)
      {
        const double x = pos[0];
        const double y = pos[1];
        const double z = pos[2];

        const double min_eta = 1.0;
        const double max_eta = eta;
        const double A(min_eta*(max_eta-min_eta)/(max_eta+min_eta));

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

          virtual void vector_value (const Point< dim >   &pos,
                                     Vector< double >   &values) const
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
    class BursteddeBoundary : public VelocityBoundaryConditions::Interface<dim>
    {
      public:
        /**
         * Constructor.
         */
        BursteddeBoundary();

        /**
         * Return the boundary velocity as a function of position.
         */
        virtual
        Tensor<1,dim>
        boundary_velocity (const Point<dim> &position) const;

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
    boundary_velocity (const Point<2> &p) const
    {
      Assert (false, ExcNotImplemented());
      return Tensor<1,2>();
    }



    template <>
    Tensor<1,3>
    BursteddeBoundary<3>::
    boundary_velocity (const Point<3> &p) const
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
    class BursteddeMaterial : public MaterialModel::InterfaceCompatibility<dim>
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
         * Return true if the viscosity() function returns something that
         * may depend on the variable identifies by the argument.
         */
        virtual bool
        viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the density() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        density_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the compressibility() function returns something
         * that may depend on the variable identifies by the argument.
         *
         * This function must return false for all possible arguments if the
         * is_compressible() function returns false.
         */
        virtual bool
        compressibility_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the specific_heat() function returns something
         * that may depend on the variable identifies by the argument.
         */
        virtual bool
        specific_heat_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the thermal_conductivity() function returns
         * something that may depend on the variable identifies by the
         * argument.
         */
        virtual bool
        thermal_conductivity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

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

        virtual double reference_density () const;

        virtual double reference_thermal_expansion_coefficient () const;

        double reference_thermal_diffusivity () const;

        double reference_cp () const;
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
    BursteddeMaterial<dim>::
    viscosity (const double,
               const double,
               const std::vector<double> &,       /*composition*/
               const SymmetricTensor<2,dim> &,
               const Point<dim> &p) const
    {
      const double mu = exp(1. - beta * (p(0)*(1.-p(0))+p(1)*(1.-p(1)) + p(2)*(1.-p(2))));
      return mu;
    }


    template <int dim>
    double
    BursteddeMaterial<dim>::
    reference_viscosity () const
    {
      return 1.;
    }

    template <int dim>
    double
    BursteddeMaterial<dim>::
    reference_density () const
    {
      return 1.;
    }

    template <int dim>
    double
    BursteddeMaterial<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return 0;
    }

    template <int dim>
    double
    BursteddeMaterial<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return 0;
    }

    template <int dim>
    double
    BursteddeMaterial<dim>::
    reference_cp () const
    {
      return 0;
    }

    template <int dim>
    double
    BursteddeMaterial<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return 0;
    }

    template <int dim>
    double
    BursteddeMaterial<dim>::
    reference_thermal_diffusivity () const
    {
      return 0;
    }

    template <int dim>
    double
    BursteddeMaterial<dim>::
    density (const double,
             const double,
             const std::vector<double> &, /*composition*/
             const Point<dim> &p) const
    {
      return 1;
    }


    template <int dim>
    double
    BursteddeMaterial<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return 0;
    }


    template <int dim>
    double
    BursteddeMaterial<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0;
    }



    template <int dim>
    bool
    BursteddeMaterial<dim>::
    viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    BursteddeMaterial<dim>::
    density_depends_on (const MaterialModel::NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    BursteddeMaterial<dim>::
    compressibility_depends_on (const MaterialModel::NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    BursteddeMaterial<dim>::
    specific_heat_depends_on (const MaterialModel::NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    BursteddeMaterial<dim>::
    thermal_conductivity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
    {
      return false;
    }

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
    }


    template <int dim>
    double
    BursteddeMaterial<dim>::get_beta() const
    {
      return beta;
    }

    /**
     *gravity model for the Burstedde benchmark
    */

    template <int dim>
    class BursteddeGravity : public aspect::GravityModel::Interface<dim>
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
    BursteddeGravity<dim>::
    gravity_vector(const Point<dim> &pos) const
    {

      const double x=pos[0];
      const double y=pos[1];
      const double z=pos[2];
      const double mu=exp(1. - beta * (x*(1.-x)+y*(1.-y) + z*(1.-z)));

      const double dmudx=-beta*(1.-2.*x)*mu;
      const double dmudy=-beta*(1.-2.*y)*mu;
      const double dmudz=-beta*(1.-2.*z)*mu;

      Tensor<1,dim> g;

      g[0]=((y*z+3.*std::pow(x,2)*std::pow(y,3)*z)- mu*(2.+6.*x*y))
           -dmudx*(2.+4.*x+2.*y+6.*std::pow(x,2)*y)
           -dmudy*(x+std::pow(x,3)+y+2.*x*std::pow(y,2))
           -dmudz*(-3.*z-10.*x*y*z);

      g[1]=((x*z+3.*std::pow(x,3)*std::pow(y,2)*z)- mu*(2.+2.*std::pow(x,2)+2.*std::pow(y,2)))
           -dmudx*(x+std::pow(x,3)+y+2.*x*std::pow(y,2))
           -dmudy*(2.+2.*x+4.*y+4.*std::pow(x,2)*y)
           -dmudz*(-3.*z-5.*std::pow(x,2)*z);

      g[2]=((x*y+std::pow(x,3)*std::pow(y,3)) - mu*(-10.*y*z))
           -dmudx*(-3.*z-10.*x*y*z)
           -dmudy*(-3.*z-5.*std::pow(x,2)*z)
           -dmudz*(-4.-6.*x-6.*y-10.*std::pow(x,2)*y);

      return g;
    }

    template <int dim>
    void
    BursteddeGravity<dim>::declare_parameters (ParameterHandler &prm)
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
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };

    template <int dim>
    std::pair<std::string,std::string>
    BursteddePostprocessor<dim>::execute (TableHandler &statistics)
    {
      std_cxx1x::shared_ptr<Function<dim> > ref_func;
      {
        const BursteddeMaterial<dim> *
        material_model
          = dynamic_cast<const BursteddeMaterial<dim> *>(&this->get_material_model());

        ref_func.reset (new AnalyticSolutions::FunctionBurstedde<dim>(material_model->get_beta()));
      }

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

      Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_ul2 (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_pl2 (this->get_triangulation().n_active_cells());

      double u_l1;
      double p_l1;
      double u_l2;
      double p_l2;

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

      u_l1 =  Utilities::MPI::sum(cellwise_errors_u.l1_norm(),MPI_COMM_WORLD);
      p_l1 =  Utilities::MPI::sum(cellwise_errors_p.l1_norm(),MPI_COMM_WORLD);
      u_l2 =  std::sqrt(Utilities::MPI::sum(cellwise_errors_ul2.norm_sqr(),MPI_COMM_WORLD));
      p_l2 =  std::sqrt(Utilities::MPI::sum(cellwise_errors_pl2.norm_sqr(),MPI_COMM_WORLD));

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

    ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS(BursteddeBoundary,
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

