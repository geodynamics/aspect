#include <aspect/material_model/melt_interface.h>
#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/fluid_pressure_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

const double c = 1.0;
const double pi = 3.14159265359;

namespace aspect
{
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
      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                 typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        for (unsigned int i=0;i<in.position.size();++i)
          {
            const double porosity = in.composition[i][porosity_idx];
            const double x = in.position[i](0);
            const double z = in.position[i](1);
            out.viscosities[i] = std::exp(c * porosity);
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
            out.densities[i] = 1.0;
            // This is the RHS we need to use for the manufactured solution.
            // We calculated it by subtracting the first term of the RHS (\Nabla (K_D \rho_f g)) from the LHS
            // we computed using our analytical solution.
            for (unsigned int c=0;c<in.composition[i].size();++c)
              out.reaction_terms[i][c] = 0.1e1 - (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2) * exp((double) z / 0.10e2) / (double) (1 + (int) pow((double) (x + 2 * z), (double) 2)) / 0.50e2 - pow(0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2, 0.2e1) * exp((double) z / 0.10e2) / 0.100e3 + (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2) * (-0.1166666667e0 * c / (double) (1 + (int) pow((double) (x + 2 * z), (double) 2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) - 0.1000000000e0 * c / (double) (1 + (int) pow((double) (x + 2 * z), (double) 2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) * (cos((double) z) + 0.1e1) + 0.1000000000e1 * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) * sin((double) z)) / (double) (1 + (int) pow((double) (x + 2 * z), (double) 2)) / 0.10e2 + pow(0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2, 0.2e1) * (0.1166666667e0 * c * (double) (int) pow((double) (1 + (int) pow((double) (x + 2 * z), (double) 2)), (double) (-2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) * (double) (2 * x + 4 * z) - 0.5833333335e-2 * c * c * (double) (int) pow((double) (1 + (int) pow((double) (x + 2 * z), (double) 2)), (double) (-2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) + 0.1000000000e0 * c * (double) (int) pow((double) (1 + (int) pow((double) (x + 2 * z), (double) 2)), (double) (-2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) * (cos((double) z) + 0.1e1) * (double) (2 * x + 4 * z) - 0.5000000000e-2 * c * c * (double) (int) pow((double) (1 + (int) pow((double) (x + 2 * z), (double) 2)), (double) (-2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) * (cos((double) z) + 0.1e1) + 0.5000000000e-1 * c / (double) (1 + (int) pow((double) (x + 2 * z), (double) 2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) * sin((double) z)) + (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2) * (-0.3333333333e-1 * c / (double) (1 + (int) pow((double) (x + 2 * z), (double) 2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) - 0.5000000000e-1 * c / (double) (1 + (int) pow((double) (x + 2 * z), (double) 2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) * (cos((double) z) + 0.1e1) + 0.1000000000e0 * exp((double) z / 0.10e2)) / (double) (1 + (int) pow((double) (x + 2 * z), (double) 2)) / 0.5e1 + pow(0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2, 0.2e1) * (0.3333333333e-1 * c * (double) (int) pow((double) (1 + (int) pow((double) (x + 2 * z), (double) 2)), (double) (-2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) * (double) (4 * x + 8 * z) - 0.3333333333e-2 * c * c * (double) (int) pow((double) (1 + (int) pow((double) (x + 2 * z), (double) 2)), (double) (-2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) + 0.5000000000e-1 * c * (double) (int) pow((double) (1 + (int) pow((double) (x + 2 * z), (double) 2)), (double) (-2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) * (cos((double) z) + 0.1e1) * (double) (4 * x + 8 * z) - 0.5000000000e-2 * c * c * (double) (int) pow((double) (1 + (int) pow((double) (x + 2 * z), (double) 2)), (double) (-2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) * (cos((double) z) + 0.1e1) + 0.5000000000e-1 * c / (double) (1 + (int) pow((double) (x + 2 * z), (double) 2)) * exp(c * (0.15e0 + pi / 0.40e2 + atan((double) (x + 2 * z)) / 0.20e2)) * sin((double) z) + 0.1000000000e-1 * exp((double) z / 0.10e2));
          }
      }

      virtual void evaluate_with_melt(const typename MaterialModel::MeltInterface<dim>::MaterialModelInputs &in,
          typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs &out) const
      {
        evaluate(in, out);
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

        for (unsigned int i=0;i<in.position.size();++i)
          {
            double porosity = in.composition[i][porosity_idx];
            out.compaction_viscosities[i] = std::exp(c * porosity);
            out.fluid_viscosities[i] = 1.0;
            out.permeabilities[i] = porosity * porosity;
            out.fluid_compressibilities[i] = 0.0;
            out.fluid_densities[i] = 1.0;
          }

      }


      private:


  };
  
  


      template <int dim>
      class RefFunction : public Function<dim>
      {
        public:
          RefFunction () : Function<dim>(2*dim+5) {}
          virtual void vector_value (const Point< dim >   &p,
                                     Vector< double >   &values) const
          {
            double x = p(0);
            double z = p(1);
            double porosity = 0.15+1.0/20.0 * (pi/2.0 * std::atan(x + 2*z));

            values[0]= x + std::sin(z);       //x vel
            values[1]= x;    //z vel
            values[2]= std::exp(z/10.0);  // p_f
            values[3]= - std::exp(c * porosity);  // p_c
            values[4]= 0;       //x melt vel
            values[5]= 0;    //y melt vel
            values[6]= 0;  // p_s
            values[7]= 0; // T
            values[8]= porosity;
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
      AssertThrow(Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) == 1,
                  ExcNotImplemented());

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
      ComponentSelectFunction<dim> comp_u_f(std::pair<unsigned int, unsigned int>(dim+2,dim+2+dim),
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
      os << std::scientific << cellwise_errors_u.l2_norm()
         << ", " << cellwise_errors_p.l2_norm()
         << ", " << cellwise_errors_p_f.l2_norm()
         << ", " << cellwise_errors_p_c.l2_norm()
         << ", " << cellwise_errors_porosity.l2_norm()
         << ", " << cellwise_errors_u_f.l2_norm();

      return std::make_pair("Errors u_L2, p_L2, p_f_L2, p_c_L2, porosity_L2, u_f_L2:", os.str());
    }

  
  template <int dim>
  class PressureBdry:
      
      public FluidPressureBoundaryConditions::Interface<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
	const typename MaterialModel::MeltInterface<dim>::MaterialModelInputs &material_model_inputs,
	const typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs &material_model_outputs,
	std::vector<Tensor<1,dim> > & output
      ) const
	{
	  for (unsigned int q=0; q<output.size(); ++q)
	    {
	      const double z = material_model_inputs.position[q][1];
	      Tensor<1,dim> gravity;
	      gravity[0] = 0.0;
	      gravity[dim-1] = 1.0;
	      output[q] = 0.1 * std::exp(z/10.0) * gravity;
	    }	  
	}
      

      
  };

}

// explicit instantiations
namespace aspect
{

    ASPECT_REGISTER_MATERIAL_MODEL(TestMeltMaterial,
                                   "test melt material",
				   "")


    ASPECT_REGISTER_POSTPROCESSOR(ConvergenceMeltPostprocessor,
                                  "melt error calculation",
                                  "A postprocessor that compares the numerical solution to the analytical "
                                  "solution derived for incompressible melt transport in a 2D box as described "
                                  "in the manuscript and reports the error.")

    ASPECT_REGISTER_FLUID_PRESSURE_BOUNDARY_CONDITIONS(PressureBdry,
						       "PressureBdry",
						       "A fluid pressure boundary condition that prescribes the "
						       "gradient of the fluid pressure at the boundaries as "
						       "calculated in the analytical solution. ")
						       
}
