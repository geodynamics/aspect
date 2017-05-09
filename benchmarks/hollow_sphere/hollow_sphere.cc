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
  /**
   * This is the "HollowSphere" benchmark defined in the following paper:
   * Analytical solution for viscous incompressible Stokes flow in a spherical shell
   * C. Thieulot, in prep.
   */
  namespace HollowSphereBenchmark
  {
    using namespace dealii;

    namespace AnalyticSolutions
    {
      Tensor<1,3>
      hollow_sphere_velocity (const Point<3> &pos,
                              const double mmm)
      {
        const double gammma = 1.0;

        const std_cxx11::array<double,3> spos =
          aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(pos);

        const double r=spos[0];
        const double phi=spos[1];
        const double theta=spos[2];

        const double R1 = 0.5;
        const double R2 = 1.0;

        double alpha,beta,fr,gr;

        if (mmm == -1)
          {
            alpha=-gammma*(pow(R2,3)-pow(R1,3))/(pow(R2,3)*log(R1)-pow(R1,3)*log(R2));
            beta=-3*gammma*(log(R2)-log(R1))/(pow(R1,3)*log(R2)-pow(R2,3)*log(R1)) ;
            fr=alpha/(r*r)+beta*r;
            gr=-2/(r*r)*(alpha*log(r)+beta/3*pow(r,3)+gammma);
          }
        else
          {
            alpha=gammma*(mmm+1)*(pow(R1,-3)-pow(R2,-3))/(pow(R1,-mmm-4)-pow(R2,-mmm-4));
            beta=-3*gammma*(pow(R1,mmm+1)-pow(R2,mmm+1))/(pow(R1,mmm+4)-pow(R2,mmm+4));
            fr=alpha/pow(r,mmm+3)+beta*r;
            gr=-2/(r*r)*(-alpha/(mmm+1)*pow(r,-mmm-1)+beta/3*pow(r,3)+gammma);
          }

        const double v_r    =gr*cos(theta);
        const double v_theta=fr*sin(theta);
        const double v_phi  =fr*sin(theta);
        const double v_x=sin(theta)*cos(phi)*v_r + cos(theta)*cos(phi)*v_theta-sin(phi)*v_phi;
        const double v_y=sin(theta)*sin(phi)*v_r + cos(theta)*sin(phi)*v_theta+cos(phi)*v_phi;
        const double v_z=cos(theta)*v_r - sin(theta)*v_theta;

        // create a Point<3> (because it has a constructor that takes
        // three doubles) and return it (it automatically converts to
        // the necessary Tensor<1,3>).
        return Point<3> (v_x,v_y,v_z);
      }

      double
      hollow_sphere_pressure (const Point<3> &pos,
                              const double mmm)
      {
        const double gammma = 1.0;
        const double mu0=1;

        const std_cxx11::array<double,3> spos =
          aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(pos);

        const double r=spos[0];
        const double phi=spos[1];
        const double theta=spos[2];

        const double R1 = 0.5;
        const double R2 = 1.0;

        double alpha,beta,fr,gr,hr,mur;


        if (mmm == -1)
          {
            mur=mu0;
            alpha=-gammma*(pow(R2,3)-pow(R1,3))/(pow(R2,3)*log(R1)-pow(R1,3)*log(R2));
            beta=-3*gammma*(log(R2)-log(R1))/(pow(R1,3)*log(R2)-pow(R2,3)*log(R1)) ;
            fr=alpha/(r*r)+beta*r;
            gr=-2/(r*r)*(alpha*log(r)+beta/3*pow(r,3)+gammma);
            hr=2/r*gr*mur;
          }
        else
          {
            mur=mu0*pow(r,mmm+1);
            alpha=gammma*(mmm+1)*(pow(R1,-3)-pow(R2,-3))/(pow(R1,-mmm-4)-pow(R2,-mmm-4));
            beta=-3*gammma*(pow(R1,mmm+1)-pow(R2,mmm+1))/(pow(R1,mmm+4)-pow(R2,mmm+4));
            fr=alpha/pow(r,mmm+3)+beta*r;
            gr=-2/(r*r)*(-alpha/(mmm+1)*pow(r,-mmm-1)+beta/3*pow(r,3)+gammma);
            hr=(mmm+3)/r*gr*mur;
          }

        return hr*cos(theta);
      }


      /**
       * The exact solution for the HollowSphere benchmark.
       */
      template <int dim>
      class FunctionHollowSphere : public Function<dim>
      {
        public:
          FunctionHollowSphere (const double mmm)
            :
            Function<dim>(dim+2),
            mmm_(mmm)
          {}

          virtual void vector_value (const Point< dim >   &pos,
                                     Vector< double >   &values) const
          {
            Assert (dim == 3, ExcNotImplemented());
            Assert (values.size() >= 4, ExcInternalError());

            const Point<3> p (pos[0], pos[1], pos[2]);

            const Tensor<1,3> v = AnalyticSolutions::hollow_sphere_velocity (p, mmm_);
            values[0] = v[0];
            values[1] = v[1];
            values[2] = v[2];

            values[3] = AnalyticSolutions::hollow_sphere_pressure (p, mmm_);
          }

        private:
          const double mmm_;
      };
    }



    template <int dim>
    class HollowSphereBoundary : public BoundaryVelocity::Interface<dim>
    {
      public:
        /**
         * Constructor.
         */
        HollowSphereBoundary();

        /**
         * Return the boundary velocity as a function of position.
         */
        virtual
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id ,
                           const Point<dim> &position) const;


        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
        */
        virtual
        void
        parse_parameters (ParameterHandler &prm);



      private:
        double mmm;
    };

    template <int dim>
    HollowSphereBoundary<dim>::HollowSphereBoundary ()
      :
      mmm (0)
    {}



    template <>
    Tensor<1,2>
    HollowSphereBoundary<2>::
    boundary_velocity (const types::boundary_id ,
                       const Point<2> &p) const
    {
      Assert (false, ExcNotImplemented());
      return Tensor<1,2>();
    }



    template <>
    Tensor<1,3>
    HollowSphereBoundary<3>::
    boundary_velocity (const types::boundary_id ,
                       const Point<3> &p) const
    {
      return AnalyticSolutions::hollow_sphere_velocity (p, mmm);
    }




    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class HollowSphereMaterial : public MaterialModel::InterfaceCompatibility<dim>
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
        double get_mmm() const;

        /**
         * viscosity value in the inclusion
         */
        double mmm;
    };

    template <int dim>
    double
    HollowSphereMaterial<dim>::
    viscosity (const double,
               const double,
               const std::vector<double> &,       /*composition*/
               const SymmetricTensor<2,dim> &,
               const Point<dim> &pos) const
    {
      const std_cxx11::array<double,dim> spos =
        aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(pos);

      const double r=spos[0];
      const double mu = pow(r,mmm+1);
      return mu;
    }


    template <int dim>
    double
    HollowSphereMaterial<dim>::
    reference_viscosity () const
    {
      return 1.;
    }


    template <int dim>
    double
    HollowSphereMaterial<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return 0;
    }


    template <int dim>
    double
    HollowSphereMaterial<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return 0;
    }


    template <int dim>
    double
    HollowSphereMaterial<dim>::
    density (const double,
             const double,
             const std::vector<double> &, /*composition*/
             const Point<dim> &p) const
    {
      return 1;
    }


    template <int dim>
    double
    HollowSphereMaterial<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return 0;
    }


    template <int dim>
    double
    HollowSphereMaterial<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0;
    }


    template <int dim>
    bool
    HollowSphereMaterial<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    HollowSphereMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      //create a global section in the parameter file for parameters
      //that describe this benchmark. note that we declare them here
      //in the material model, but other kinds of plugins (e.g., the gravity
      //model below) may also read these parameters even though they do not
      //declare them
      prm.enter_subsection("HollowSphere benchmark");
      {
        prm.declare_entry("Viscosity parameter", "-1",
                          Patterns::Double (),
                          "Viscosity in the HollowSphere benchmark.");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    HollowSphereMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("HollowSphere benchmark");
      {
        mmm = prm.get_double ("Viscosity parameter");
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
    void
    HollowSphereBoundary<dim>::declare_parameters (ParameterHandler &prm)
    {
      //nothing to declare here. This plugin will however, read parameters
      //declared by the material model in the "HollowSphere benchmark" section
    }


    template <int dim>
    void
    HollowSphereBoundary<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("HollowSphere benchmark");
      {
        mmm = prm.get_double ("Viscosity parameter");
      }
      prm.leave_subsection();
    }



    template <int dim>
    double
    HollowSphereMaterial<dim>::get_mmm() const
    {
      return mmm;
    }

    /**
     *gravity model for the HollowSphere benchmark
    */

    template <int dim>
    class HollowSphereGravity : public aspect::GravityModel::Interface<dim>
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

        double mmm;
    };

    template <int dim>
    Tensor<1,dim>
    HollowSphereGravity<dim>::
    gravity_vector(const Point<dim> &pos) const
    {

      const double x=pos[0];
      const double y=pos[1];
      const double z=pos[2];

      const std_cxx11::array<double,dim> spos =
        aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(pos);

      const double r=spos[0];
      const double phi=spos[1];
      const double theta=spos[2];

      Tensor<1,dim> g;

      const double gammma = 1.0;
      const double R1 = 0.5;
      const double R2 = 1.0;

      double alpha,beta,rho;

      if (mmm == -1)
        {
          alpha = -gammma*(pow(R2,3)-pow(R1,3))/(pow(R2,3)*log(R1)-pow(R1,3)*log(R2));
          beta  = -3*gammma*(log(R2)-log(R1))/(pow(R1,3)*log(R2)-pow(R2,3)*log(R1)) ;
          rho=-(alpha/pow(r,4)*(8*log(r)-6) + 8./3.*beta/r+8*gammma/pow(r,4))*cos(theta);
        }
      else
        {
          alpha=gammma*(mmm+1)*(pow(R1,-3)-pow(R2,-3))/(pow(R1,-mmm-4)-pow(R2,-mmm-4));
          beta=-3*gammma*(pow(R1,mmm+1)-pow(R2,mmm+1))/(pow(R1,mmm+4)-pow(R2,mmm+4));
          rho=-(2*alpha*pow(r,-4)*(mmm+3)/(mmm+1)*(mmm-1)-2*beta/3*(mmm-1)*(mmm+3)*pow(r,mmm)-mmm*(mmm+5)*2*gammma*pow(r,mmm-3) )*cos(theta);
        }

      const double g_x= -x/r;
      const double g_y= -y/r;
      const double g_z= -z/r;

      g[0]=rho*g_x;
      g[1]=rho*g_y;
      g[2]=rho*g_z;

      return g;
    }

    template <int dim>
    void
    HollowSphereGravity<dim>::declare_parameters (ParameterHandler &prm)
    {
      //nothing to declare here. This plugin will however, read parameters
      //declared by the material model in the "HollowSphere benchmark" section
    }

    template <int dim>
    void
    HollowSphereGravity<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("HollowSphere benchmark");
      {
        mmm = prm.get_double ("Viscosity parameter");
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
    class HollowSpherePostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
    HollowSpherePostprocessor<dim>::execute (TableHandler &statistics)
    {
      std_cxx1x::shared_ptr<Function<dim> > ref_func;
      {
        const HollowSphereMaterial<dim> *
        material_model
          = dynamic_cast<const HollowSphereMaterial<dim> *>(&this->get_material_model());

        ref_func.reset (new AnalyticSolutions::FunctionHollowSphere<dim>(material_model->get_mmm()));
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
  namespace HollowSphereBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(HollowSphereMaterial,
                                   "HollowSphereMaterial",
                                   "A material model that corresponds to the `HollowSphere' benchmark. "
                                   "See the manual for more information.")

    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(HollowSphereBoundary,
                                            "HollowSphereBoundary",
                                            "Implementation of the velocity boundary conditions for the "
                                            "`HollowSphere' benchmark. See the manual for more information about this "
                                            "benchmark.")

    ASPECT_REGISTER_POSTPROCESSOR(HollowSpherePostprocessor,
                                  "HollowSpherePostprocessor",
                                  "A postprocessor that compares the solution of the `HollowSphere' benchmark "
                                  "with the one computed by ASPECT "
                                  "and reports the error. See the manual for more information.")
    ASPECT_REGISTER_GRAVITY_MODEL(HollowSphereGravity,
                                  "HollowSphereGravity",
                                  "A gravity model in corresponding to the `HollowSphere' benchmark. "
                                  "See the manual for more information.")
  }
}

