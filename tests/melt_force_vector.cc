#include <aspect/melt.h>
#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/gravity_model/interface.h>
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

  template <int dim>
  class TestMeltMaterial:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:

      virtual double reference_darcy_coefficient () const
      {
        return 1.0;
      }


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
        MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
        *force = out.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >();

        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            const double x = in.position[i](0);
            const double z = in.position[i](1);
            out.viscosities[i] = 1.0;
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
            out.densities[i] = 1.0;

            out.reaction_terms[i][porosity_idx] = 0;
//**********
// copy and paste here
            if (force)
              {
                force->rhs_u[i][0] = cos(z) - cos(x + z) + cos(x * z) * z;
                force->rhs_u[i][1] = sin(x) - cos(x + z) + cos(x * z) * x;
                force->rhs_p[i] = 0.4e1 * sin(x + z) - sin(x * z) * z * z - sin(x * z) * x * x;
                force->rhs_melt_pc[i] = -0.10e2 * sin(x + z) / (0.1e1 + exp(-0.20e2 * x * x - 0.20e2 * z * z + 0.1e1));
              }
//***********
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim> >();

        if (melt_out != nullptr)
          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              double porosity = in.composition[i][porosity_idx];
              // porosity = 0.1000000000e-1 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.2e1 * z, 0.2e1));
              const double x = in.position[i](0);
              const double z = in.position[i](1);
              melt_out->compaction_viscosities[i] = 0.1e0 + 0.1e0 * exp(-0.20e2 * x * x - 0.20e2 * z * z + 0.1e1); // xi
              melt_out->fluid_viscosities[i] = 1.0;
              melt_out->permeabilities[i] = 1.0; // K_D
              melt_out->fluid_density_gradients[i] = 0.0;
              melt_out->fluid_densities[i] = 0.5;
            }

      }

  };



  template <int dim>
  class RefFunction : public Function<dim>
  {
    public:
      RefFunction () : Function<dim>(2*dim+3+2) {}
      virtual void vector_value (const Point< dim >   &p,
                                 Vector< double >   &values) const
      {
        double x = p(0);
        double z = p(1);
//**********
// copy and paste here (add "out.")
        const double phi = 0.1000000000e-1 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1));

        values[0] = cos(z);                                                 //ux
        values[1] = sin(x);                                                 //uz
        values[2] = -0.2e1 * sin(x + z) + sin(x * z);                       //p_f
        values[3] = sin(x + z);                                             //p_c
        values[4] = values[0] - 1.0/phi * (-2.0 * cos(x+z) + cos(x*z)*z);   //u_f_x
        values[5] = values[1] - 1.0/phi * (-2.0 * cos(x+z) + cos(x*z)*x);   //u_f_z
        values[6] = values[2] + values[3] / (1.0 - phi);                    //p_s = p_f + p_c / (1-phi)
        values[7] = 0;
        values[8] = phi;                                                    //porosity
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
    const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities +2);

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

    const double u_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_u, VectorTools::L2_norm);
    const double p_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p, VectorTools::L2_norm);
    const double p_f_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p_f, VectorTools::L2_norm);
    const double p_c_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p_c, VectorTools::L2_norm);
    const double phi_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_porosity, VectorTools::L2_norm);
    const double u_f_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_u_f, VectorTools::L2_norm);


    std::ostringstream os;
    os << std::scientific
       << " h = " << this->get_triangulation().begin_active()->diameter()
       << " ndofs= " << this->get_solution().size()
       << " u_L2= " << u_l2
       << " p_L2= "  << p_l2
       << " p_f_L2= " << p_f_l2
       << " p_c_L2= " << p_c_l2
       << " phi_L2= " << phi_l2
       << " u_f_L2= " << u_f_l2
       ;

    return std::make_pair("Errors", os.str());
  }


  template <int dim>
  class PressureBdry:
    public BoundaryFluidPressure::Interface<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const dealii::types::boundary_id,
        const typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs,
        const typename MaterialModel::Interface<dim>::MaterialModelOutputs &material_model_outputs,
        const std::vector<Tensor<1,dim> > &normal_vectors,
        std::vector<double> &output
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
  ASPECT_REGISTER_MATERIAL_MODEL(TestMeltMaterial,
                                 "test melt material",
                                 "")

  ASPECT_REGISTER_POSTPROCESSOR(ConvergenceMeltPostprocessor,
                                "melt error calculation",
                                "A postprocessor that compares the numerical solution to the analytical solution")

  ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(PressureBdry,
                                                "PressureBdry",
                                                "A fluid pressure boundary condition that prescribes the "
                                                "gradient of the fluid pressure at the boundaries as "
                                                "calculated in the analytical solution. ")

}
