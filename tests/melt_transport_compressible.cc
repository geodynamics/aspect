#include <aspect/material_model/melt_interface.h>
#include <aspect/velocity_boundary_conditions/interface.h>
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
  class CompressibleMeltMaterial:
      public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
      public:
      virtual bool
      viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }

      virtual bool
      density_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        if ((dependence & MaterialModel::NonlinearDependence::compositional_fields) != MaterialModel::NonlinearDependence::none)
          return true;
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
        return true;
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
            double porosity = in.composition[i][porosity_idx];
            out.viscosities[i] = 0.5 * std::exp(2.0 * in.position[i][0]);
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 1.0 / (rho_s_0 * C);
            out.densities[i] = rho_s_0 * std::exp(-in.position[i][1]);
            for (unsigned int c=0;c<in.composition[i].size();++c)
              out.reaction_terms[i][c] = - rho_s_0 * B * D * std::exp(in.position[i][1]);
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
            out.compaction_viscosities[i] = xi_1 * std::exp(-in.position[i][1]) + 2.0/3.0 * std::exp(2.0 * in.position[i][0]) + xi_0;
            out.fluid_viscosities[i] = 1.0;
            out.permeabilities[i] = K_D_0 + 2.0 * B / E - rho_s_0 * B * D / E * (1.0/rho_s_0 - 1.0/rho_f_0) * std::exp(in.position[i][1]);
            out.fluid_compressibilities[i] = 1.0 / (rho_f_0 * C);
            out.fluid_densities[i] = rho_f_0 * std::exp(-in.position[i][1]);
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
          RefFunction () : Function<dim>(dim+2) {}
          virtual void vector_value (const Point< dim >   &p,
                                     Vector< double >   &values) const
          {
            double x = p(0);
            double y = p(1);

            values[0]=0.1 * std::exp(y);       //x vel
            values[1]=-0.075 * std::exp(y);    //y vel
            values[2]=0;  // p_s
            values[3]=-0.135*(std::exp(y) - std::exp(1)) + 1.0 - y;  // p_f
            values[4]=0; // T
            values[5]=1.0 - 0.3 * std::exp(y); // porosity
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
      AssertThrow(Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) == 1,
                  ExcNotImplemented());

      RefFunction<dim> ref_func;
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

      Vector<float> cellwise_errors_porosity (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());

      ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                          dim+4);
      ComponentSelectFunction<dim> comp_pf(dim+1, dim+4);
      ComponentSelectFunction<dim> comp_porosity(dim+3, dim+4);

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
                                         cellwise_errors_p,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_pf);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         ref_func,
                                         cellwise_errors_porosity,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_porosity);
      std::ostringstream os;
      os << std::scientific << cellwise_errors_u.l2_norm()
         << ", " << cellwise_errors_p.l2_norm()
         << ", " << cellwise_errors_porosity.l2_norm();

      return std::make_pair("Errors u_L2, p_fL2, porosity_L2:", os.str());
    }

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
                                  "solution derived for compressible melt transport in a 2D box described in "
                                  "the manuscript and reports the error.")
}
