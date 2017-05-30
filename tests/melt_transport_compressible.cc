#include <aspect/material_model/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>


namespace aspect
{
  template <int dim>
  class CompressibleMeltMaterial:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      virtual bool is_compressible () const
      {
        return true;
      }

      virtual double reference_viscosity () const
      {
        return 1.0;
      }

      virtual double reference_darcy_coefficient () const
      {
        const double permeability = K_D_0 + 2.0 * B / E - rho_s_0 * B * D / E * (1.0/rho_s_0 - 1.0/rho_f_0) * std::exp(0.5);
        return permeability / 1.0;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        for (unsigned int i=0; i<in.position.size(); ++i)
          {
            double porosity = in.composition[i][porosity_idx];
            out.viscosities[i] = 0.5 * std::exp(2.0 * in.position[i][0]);
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 1.0 / (rho_s_0 * C);
            out.densities[i] = rho_s_0 * std::exp(-in.position[i][1]);
            for (unsigned int c=0; c<in.composition[i].size(); ++c)
              out.reaction_terms[i][c] = - rho_s_0 * B * D * std::exp(in.position[i][1]);
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim> >();

        if (melt_out != NULL)
          {
            for (unsigned int i=0; i<in.position.size(); ++i)
              {
                double porosity = in.composition[i][porosity_idx];
                melt_out->compaction_viscosities[i] = xi_1 * std::exp(-in.position[i][1]) + 2.0/3.0 * std::exp(2.0 * in.position[i][0]) + xi_0;
                melt_out->fluid_viscosities[i] = 1.0;
                melt_out->permeabilities[i] = K_D_0 + 2.0 * B / E - rho_s_0 * B * D / E * (1.0/rho_s_0 - 1.0/rho_f_0) * std::exp(in.position[i][1]);
                const double fluid_compressibility = 1.0 / (rho_f_0 * C);
                melt_out->fluid_densities[i] = rho_f_0 * std::exp(-in.position[i][1]);
                melt_out->fluid_density_gradients[i] = melt_out->fluid_densities[i] * melt_out->fluid_densities[i]
                                                       * fluid_compressibility
                                                       * this->get_gravity_model().gravity_vector(in.position[i]);
              }
          }
      }

      virtual void initialize ()
      {
        rho_s_0 = 1.2;
        rho_f_0 = 1.0;
        xi_0 = 1.0;
        xi_1 = 1.0;

        // A, B and C are constants from the velocity boundary conditions and gravity model
        // they have to be consistent!
        A = 0.1;
        B = -3.0/4.0 * A;
        C = 1.0;
        D = 0.3;
        E = - 3.0/4.0 * xi_0 * A + C * D *(rho_f_0 - rho_s_0);

        K_D_0 = 2.2;
      }


    private:
      double rho_s_0;
      double rho_f_0;
      double xi_0;
      double xi_1;
      double K_D_0;
      double A;
      double B;
      double C;
      double D;
      double E;


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
        double y = p(1);
        double porosity = 1.0 - 0.3 * std::exp(y);
        double K_D = 2.2 + 2.0 * 0.075/0.135 + (1.0 - 5.0/6.0) * 0.075 * 0.3 * 1.2 / 0.135 * std::exp(y);

        values[0]=0.1 * std::exp(y);       // x vel
        values[1]=-0.075 * std::exp(y);    // y vel
        values[2]=-0.135*(std::exp(y) - std::exp(1)) + 1.0 - y;  // p_f
        values[3]=0.75 * (std::exp(-y) + 2.0/3.0 * std::exp(2.0*x) + 1.0) * 0.1 * std::exp(y);  // p_c
        values[4]=0.1 * std::exp(y);       // x melt vel
        values[5]=-0.075 * std::exp(y) + 0.135 * std::exp(y) * K_D / porosity;    // y melt vel

        values[6]=values[2] + values[3] / (1.0 - porosity);  // p_s
        values[7]=0; // T
        values[8]=porosity; // porosity
      }
  };

  /**
    * A postprocessor that evaluates the accuracy of the solution
    * by using the L2 norm.
    */
  template <int dim>
  class CompressibleMeltPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      /**
       * Generate graphical output from the current solution.
       */
      virtual
      std::pair<std::string,std::string>
      execute (TableHandler &statistics);


    private:
      double rho_s_0;
      double rho_f_0;
      double xi_0;
      double xi_1;
      double K_D_0;
      double A;
      double B;
      double C;
      double D;
      double E;
  };

  template <int dim>
  std::pair<std::string,std::string>
  CompressibleMeltPostprocessor<dim>::execute (TableHandler &statistics)
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
    os << std::scientific << std::sqrt(Utilities::MPI::sum(cellwise_errors_u.norm_sqr(),MPI_COMM_WORLD))
       << ", " << std::sqrt(Utilities::MPI::sum(cellwise_errors_p.norm_sqr(),MPI_COMM_WORLD))
       << ", " << std::sqrt(Utilities::MPI::sum(cellwise_errors_p_f.norm_sqr(),MPI_COMM_WORLD))
       << ", " << std::sqrt(Utilities::MPI::sum(cellwise_errors_p_c.norm_sqr(),MPI_COMM_WORLD))
       << ", " << std::sqrt(Utilities::MPI::sum(cellwise_errors_porosity.norm_sqr(),MPI_COMM_WORLD))
       << ", " << std::sqrt(Utilities::MPI::sum(cellwise_errors_u_f.norm_sqr(),MPI_COMM_WORLD));

    return std::make_pair("Errors u_L2, p_L2, p_f_L2, p_c_L2, porosity_L2, u_f_L2:", os.str());
  }


  template <int dim>
  class PressureBdry:

    public BoundaryFluidPressure::Interface<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const types::boundary_id boundary_indicator,
        const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
        const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
        const std::vector<Tensor<1,dim> > &normal_vectors,
        std::vector<double> &output
      ) const
      {
        for (unsigned int q=0; q<output.size(); ++q)
          {
            const double rho_s_0 = 1.2;
            const double rho_f_0 = 1.0;
            const double xi_0 = 1.0;
            const double xi_1 = 1.0;
            const double A = 0.1;
            const double B = -3.0/4.0 * A;
            const double C = 1.0;
            const double D = 0.3;
            const double E = - 3.0/4.0 * xi_0 * A + C * D *(rho_f_0 - rho_s_0);
            const double y = material_model_inputs.position[q][1];
            Tensor<1,dim> gravity;
            gravity[dim-1] = 1.0;
            output[q] = (E * std::exp(y) - rho_f_0 * C) * gravity * normal_vectors[q];
          }
      }



  };

}

// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(CompressibleMeltMaterial,
                                 "compressible melt material",
                                 "")


  ASPECT_REGISTER_POSTPROCESSOR(CompressibleMeltPostprocessor,
                                "compressible melt error",
                                "A postprocessor that compares the numerical solution to the analytical "
                                "solution derived for compressible melt transport in a 2D box as described "
                                "in the manuscript and reports the error.")

  ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(PressureBdry,
                                                "PressureBdry",
                                                "A fluid pressure boundary condition that prescribes the "
                                                "gradient of the fluid pressure at the boundaries as "
                                                "calculated in the analytical solution. ")

}
