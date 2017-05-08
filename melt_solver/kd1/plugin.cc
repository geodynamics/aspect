#include <aspect/melt.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/mesh_refinement/interface.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

double g = 100.0;

namespace aspect
{
    template <int dim>
    class MeltVelocityRefinement : public MeshRefinement::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Execute this mesh refinement criterion.
         *
         * @param[out] error_indicators A vector that for every active cell of
         * the current mesh (which may be a partition of a distributed mesh)
         * provides an error indicator. This vector will already have the
         * correct size when the function is called.
         */
        virtual
        void
        execute (Vector<float> &indicators) const
	  {
	    indicators = 0;

 //commented out so that we are using only global refinement
	    /*      KellyErrorEstimator<dim>::estimate (this->get_mapping(),
                                          this->get_dof_handler(),
                                          QGauss<dim-1>(this->get_fe().base_element(this->introspection().base_elements.velocities).degree+1),
                                          typename FunctionMap<dim>::type(),
                                          this->get_solution(),
                                          indicators,
                                          this->introspection().component_masks.fluid_velocities,
                                          0,
                                          0,
                                          this->get_triangulation().locally_owned_subdomain());
	    */
           }
	
    };
 
    template <int dim>
    class PCRefinement : public MeshRefinement::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Execute this mesh refinement criterion.
         *
         * @param[out] error_indicators A vector that for every active cell of
         * the current mesh (which may be a partition of a distributed mesh)
         * provides an error indicator. This vector will already have the
         * correct size when the function is called.
         */
        virtual
        void
        execute (Vector<float> &indicators) const
	  {
	    indicators = 0;
	    /*
      KellyErrorEstimator<dim>::estimate (this->get_mapping(),
                                          this->get_dof_handler(),
                                          QGauss<dim-1>(this->get_fe().base_element(this->introspection().base_elements.velocities).degree+1),
                                          typename FunctionMap<dim>::type(),
                                          this->get_solution(),
                                          indicators,
                                          this->introspection().component_masks.compaction_pressure,
                                          0,
                                          0,
                                          this->get_triangulation().locally_owned_subdomain());
	    */
	  }
	
    };




  
  template <int dim>
  class TestMeltMaterial:
      public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
      public:
      virtual bool
      viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        if ((dependence & MaterialModel::NonlinearDependence::compositional_fields) != MaterialModel::NonlinearDependence::none)
          return true;
        return false;
      }

      virtual bool
      density_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      compressibility_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      specific_heat_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }


      virtual bool
      thermal_conductivity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }

      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_viscosity () const
      {
        return 1.0;
      }

      virtual double reference_density () const
      {
        return 1.0;
      }

      virtual double reference_darcy_coefficient () const
      {
	//	double permeability = pow(exp(-0.50e2 * pow(z / 0.100000e6 - exp(-x / 0.100000e6) - 0.1e1, 0.2e1)), 0.3e1) * pow(cos(pi * y / 0.200000e6), 0.3e1) / 0.1250000000000e13;

        return 1;//0.01;
      }
      
      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                 typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        double c = 1.0;
        MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
         *force = out.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >();

        for (unsigned int i=0;i<in.position.size();++i)
          {
            const double porosity = in.composition[i][porosity_idx];
            const double x = in.position[i](0);
            const double z = in.position[i](1);
            out.viscosities[i] = 1.0;//exp(c*porosity);
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
            out.densities[i] = 1.0;
            // This is the RHS we need to use for the manufactured solution.
            // We calculated it by subtracting the first term of the RHS (\Nabla (K_D \rho_f g)) from the LHS
            // we computed using our analytical solution.
	    double reactionterm;
	    reactionterm = 0.0; // -0.100e3 * exp(-0.200e2 * pow(x + 0.2e1 * z, 0.2e1)) * (0.7840400000e1 * x + 0.1568080000e2 * z + 0.4e1 * cos(z) * x + 0.8e1 * cos(z) * z + 0.8e1 * sin(x) * x + 0.16e2 * sin(x) * z - 0.3184000000e1 * x * exp(-0.200e2 * pow(x + 0.2e1 * z, 0.2e1)) + 0.1600000000e0 * exp(-0.400e2 * pow(x + 0.2e1 * z, 0.2e1)) * x - 0.6368000000e1 * z * exp(-0.200e2 * pow(x + 0.2e1 * z, 0.2e1)) + 0.3200000000e0 * exp(-0.400e2 * pow(x + 0.2e1 * z, 0.2e1)) * z) * pow(-0.9950000000e1 + exp(-0.200e2 * pow(x + 0.2e1 * z, 0.2e1)), -0.2e1);

        out.reaction_terms[i][porosity_idx] = reactionterm;
//**********
// copy and paste here (add "out.")
//out.force_vector[i][0] = cos(z) - cos(x + z) + cos(x * z) * z;
//out.force_vector[i][1] = sin(x) - cos(x + z) + cos(x * z) * x;
//out.force_vector[i][2] = 0.4e1 * sin(x + z) - sin(x * z) * z * z - sin(x * z) * x * x;
//out.force_vector[i][3] = -0.10e2 * sin(x + z) / (0.1e1 + exp(-0.20e2 * x * x - 0.20e2 * z * z + 0.1e1));
	if (force)
	  {
	    
force->rhs_u[i][0] = cos(z) - cos(x + z) + z * cos(x * z);
force->rhs_u[i][1] = sin(x) - cos(x + z) + x * cos(x * z);
force->rhs_p[i] = 0.4e1 * sin(x + z) - z * z * sin(x * z) - x * x * sin(x * z);
force->rhs_melt_pc[i] = -0.1e1 / g * sin(x + z);
}
//***********
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim> >();

        if (melt_out != NULL)
          for (unsigned int i=0;i<in.position.size();++i)
            {
              double porosity = in.composition[i][porosity_idx];

              const double x = in.position[i](0);
              const double z = in.position[i](1);
              //porosity = 0.1000000000e-1 + 0.1000000000e0 * exp(-0.40e1 * pow(x + 0.20e1 * z, 0.2e1));
              //porosity = 0.1000000000e-1 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.2e1 * z, 0.2e1));
              melt_out->compaction_viscosities[i] = g; //0.1e0 + 0.1e0 * exp(-0.20e2 * x * x - 0.20e2 * z * z + 0.1e1);
              // xi
              melt_out->fluid_viscosities[i] = 1.0;
              melt_out->permeabilities[i] = 1.0; //1.e-6; //1.0; //0.1000000000e-1 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)); //porosity;// K_D
              melt_out->fluid_density_gradients[i] = 0.0;
              melt_out->fluid_densities[i] = 0.5;
            }

      }

  };
  

  template <int dim>
    class Gravity : public aspect::GravityModel::Interface<dim>
    {
      public:
        virtual Tensor<1,dim> gravity_vector (const Point<dim> &pos) const
	  {
	    const double x=pos[0];
	    const double z=pos[1];
	    Tensor<1,dim> gravity;
// zero when we have the force vector
gravity[0] = 0.0;
gravity[1] = 0.0;
	    
	    return gravity;
	  }
	

    }; 


      template <int dim>
      class RefFunction : public Function<dim>
      {
        public:
          RefFunction () : Function<dim>(dim+2) {}
          virtual void vector_value (const Point< dim >   &p,
                                     Vector< double >   &values) const
          {
            double x = p(0);
            double z = p(1);
//**********
// copy and paste here (add "out.")
values[0] = cos(z);
values[1] = sin(x);
values[2] = -0.2e1 * sin(x + z) + sin(x * z);
values[3] = sin(x + z);
values[4] = 1;
values[5] = 1;
values[6] = sin(x * z);
values[7] = 0;
values[8] = 0.1000000000e-1 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1));
//**********
          }
      };



    /**
      * A postprocessor that evaluates the accuracy of the solution
      * by using the L2 norm.
      */
    template <int dim>
    class ConvergenceMeltPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
    ConvergenceMeltPostprocessor<dim>::execute (TableHandler &statistics)
    {
      RefFunction<dim> ref_func;
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

      const unsigned int n_total_comp = this->introspection().n_components;

      Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p_f (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p_c (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_u_f (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_porosity (this->get_triangulation().n_active_cells());

      ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                                n_total_comp);
            ComponentSelectFunction<dim> comp_p_f(dim, n_total_comp);
            ComponentSelectFunction<dim> comp_p_c(dim+1, n_total_comp);
            ComponentSelectFunction<dim> comp_u_f(std::pair<unsigned int, unsigned int>(dim+2,dim+2+
      dim),
                                                n_total_comp);
            ComponentSelectFunction<dim> comp_p(dim+2+dim, n_total_comp);
            ComponentSelectFunction<dim> comp_porosity(dim+2+dim+2, n_total_comp);

      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         ref_func,
                                         cellwise_errors_u,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_u);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         ref_func,
                                         cellwise_errors_p_f,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_p_f);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         ref_func,
                                         cellwise_errors_p,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_p);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         ref_func,
                                         cellwise_errors_p_c,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_p_c);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         ref_func,
                                         cellwise_errors_porosity,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_porosity);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         ref_func,
                                         cellwise_errors_u_f,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_u_f);

      std::ostringstream os;
      os << std::scientific
         << " h = " << this->get_triangulation().begin_active()->diameter()
	 << " ndofs= " << this->get_solution().size()
	 << " u_L2= " << std::sqrt(Utilities::MPI::sum(cellwise_errors_u.norm_sqr(),MPI_COMM_WORLD))
	 << " p_L2= "  << std::sqrt(Utilities::MPI::sum(cellwise_errors_p.norm_sqr(),MPI_COMM_WORLD))
         << " p_f_L2= " << std::sqrt(Utilities::MPI::sum(cellwise_errors_p_f.norm_sqr(),MPI_COMM_WORLD))
         << " p_c_L= " << std::sqrt(Utilities::MPI::sum(cellwise_errors_p_c.norm_sqr(),MPI_COMM_WORLD))
         << " phi_L2= " << std::sqrt(Utilities::MPI::sum(cellwise_errors_porosity.norm_sqr(),MPI_COMM_WORLD))
         << " u_f_L2= " << std::sqrt(Utilities::MPI::sum(cellwise_errors_u_f.norm_sqr(),MPI_COMM_WORLD))
	;

      return std::make_pair("Errors", os.str());
    }

  
  template <int dim>
  class PressureBdry:
      
      public BoundaryFluidPressure::Interface<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (const dealii::types::boundary_id,
    const typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs,
    const typename MaterialModel::Interface<dim>::MaterialModelOutputs &material_model_outputs,
			const std::vector<dealii::Tensor<1,dim> >& normal_vectors,
	std::vector<double> & output
      ) const
	{
	  for (unsigned int q=0; q<output.size(); ++q)
	    {
	      const double x = material_model_inputs.position[q][0];
	      const double z = material_model_inputs.position[q][1];
	      Tensor<1,dim> gradient;
//**********
// copy and paste here (add "out.")
gradient[0] = -0.2e1 * cos(x + z) + cos(x * z) * z;
gradient[1] = -0.2e1 * cos(x + z) + cos(x * z) * x;
//**********
	      output[q] = gradient * normal_vectors[q];
	    }	  
	}
      

      
  };

}

// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_GRAVITY_MODEL(Gravity,
                                  "MyGravity",
                                  "")
  
    ASPECT_REGISTER_MATERIAL_MODEL(TestMeltMaterial,
                                   "test melt material",
				   "")
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(MeltVelocityRefinement,
                                              "melt velocity",
                                              "")
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(PCRefinement,
                                              "pc",
                                              "")

    ASPECT_REGISTER_POSTPROCESSOR(ConvergenceMeltPostprocessor,
                                  "melt error calculation",
                                  "A postprocessor that compares the numerical solution to the analytical "
                                  "solution derived for incompressible melt transport in a 2D box as described "
                                  "in the manuscript and reports the error.")

    ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(PressureBdry,
						       "PressureBdry",
						       "A fluid pressure boundary condition that prescribes the "
						       "gradient of the fluid pressure at the boundaries as "
						       "calculated in the analytical solution. ")
						       
}
