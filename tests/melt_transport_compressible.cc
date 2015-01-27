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

      template <int dim>
      class RefFunctionFluid : public Function<dim>
      {
        public:
          RefFunctionFluid () : Function<dim>(dim+2) {}
          virtual void vector_value (const Point< dim >   &p,
                                     Vector< double >   &values) const
          {
            double x = p(0);
            double y = p(1);
            double porosity = 1.0 - 0.3 * std::exp(y);
            double K_D = 2.2 + 2.0 * 0.075/0.135 + (1.0 - 5.0/6.0) * 0.075 * 0.3 * 1.2 / 0.135 * std::exp(y);

            values[0]=0.1 * std::exp(y);       //x melt vel
            values[1]=-0.075 * std::exp(y) + 0.135 * std::exp(y) * K_D / porosity;    //y melt vel
            values[2]=0;  // p_s
            values[3]=0.75 * (std::exp(-y) + 2.0/3.0 * std::exp(2.0*x) + 1.0) * 0.1 * std::exp(y);  // p_c
            values[4]=0; // T
            values[5]=porosity;
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
      RefFunctionFluid<dim> ref_func_fluid;
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

      // we need the compaction pressure and the melt velocity, but we only have
      // the solid and the fluid pressure stored in the solution vector.
      // Hence, we create a new vector only with the compaction pressure.
      LinearAlgebra::BlockVector compaction_pressure(this->get_solution());
      LinearAlgebra::BlockVector melt_velocity(this->get_solution());

      const unsigned int por_idx = this->introspection().compositional_index_for_name("porosity");
      const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.velocities).get_unit_support_points());
      std::vector<double> porosity_values(quadrature.size());
      std::vector<Tensor<1,dim> > pressure_gradients(quadrature.size());
      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature.size()));

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values | update_gradients);

      typename MaterialModel::MeltInterface<dim>::MaterialModelInputs in(fe_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs out(fe_values.n_quadrature_points, this->n_compositional_fields());

      std::vector<types::global_dof_index> local_dof_indices (this->get_fe().dofs_per_cell);
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            cell->get_dof_indices (local_dof_indices);

            fe_values[this->introspection().extractors.compositional_fields[por_idx]].get_function_values (
            		this->get_solution(), porosity_values);
            fe_values[this->introspection().extractors.compaction_pressure].get_function_gradients (
            		this->get_solution(), pressure_gradients);

            // get the various components of the solution, then
            // evaluate the material properties there
            fe_values[this->introspection().extractors.temperature]
            .get_function_values (this->get_solution(), in.temperature);
            fe_values[this->introspection().extractors.pressure]
            .get_function_values (this->get_solution(), in.pressure);
            fe_values[this->introspection().extractors.velocities]
            .get_function_symmetric_gradients (this->get_solution(), in.strain_rate);

            in.position = fe_values.get_quadrature_points();

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              fe_values[this->introspection().extractors.compositional_fields[c]]
              .get_function_values(this->get_solution(),
                                   composition_values[c]);
            for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
              {
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  in.composition[i][c] = composition_values[c][i];
              }

            const typename MaterialModel::MeltInterface<dim> * melt_mat = dynamic_cast<const MaterialModel::MeltInterface<dim>*> (&this->get_material_model());
            AssertThrow(melt_mat != NULL, ExcMessage("Need MeltMaterial if include_melt_transport is on."));
            melt_mat->evaluate_with_melt(in, out);

            for (unsigned int j=0; j<this->get_fe().base_element(this->introspection().base_elements.pressure).dofs_per_cell; ++j)
              {
                unsigned int pressure_idx
                = this->get_fe().component_to_system_index(this->introspection().component_indices.pressure,
                    /*dof index within component=*/ j);

                // skip entries that are not locally owned:
                if (!this->get_dof_handler().locally_owned_dofs().is_element(pressure_idx))
                  continue;

                unsigned int p_f_idx
                = this->get_fe().component_to_system_index(this->introspection().component_indices.compaction_pressure,
                    /*dof index within component=*/ j);

                double p_s = this->get_solution()(local_dof_indices[pressure_idx]);
                double p_f = this->get_solution()(local_dof_indices[p_f_idx]);

                double phi = porosity_values[j];
                double p_c;
                p_c = (1.0-phi) * (p_s - p_f);

                compaction_pressure(local_dof_indices[p_f_idx]) = p_c;
              }

            unsigned int end = this->get_fe().base_element(this->introspection().base_elements.velocities).dofs_per_cell;
            for (unsigned int j=0; j<end; ++j)
              {
                unsigned int velocity_idx[dim];
                for (unsigned int d=0; d<dim; ++d)
                  velocity_idx[d] = this->get_fe().component_to_system_index(this->introspection().component_indices.velocities[d],
                     j);

                Tensor<1,dim> u_f;
                for (unsigned int d=0; d<dim; ++d)
                  u_f[d] = this->get_solution()(local_dof_indices[velocity_idx[d]]);

                double phi = porosity_values[j];

                if (phi > 1e-7)
                  {
                    double K_D = out.permeabilities[j] / out.fluid_viscosities[j];
                    Tensor<1,dim> grad_p_f = pressure_gradients[j];
                    const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(in.position[j]);

                    // v_f =  v_s - K_D (nabla p_f - rho_f g) / phi
                    Tensor<1,dim> correction = - K_D * (grad_p_f - out.fluid_densities[j] * gravity) / phi;
                    u_f += correction;
//                    std::cout << K_D << " " << grad_p_f << " " << out.fluid_densities[j] << " " << gravity << std::endl;
//                    std::cout << j << " "<< in.position[j](1) << " " << u_f << std::endl;
//                    std::cout << j << " " << out.permeabilities[j] << " " << out.fluid_viscosities[j] << std::endl;
                  }

                for (unsigned int d=0; d<dim; ++d)
                  melt_velocity(local_dof_indices[velocity_idx[d]]) = u_f[d];
              }
          }

      Vector<float> cellwise_errors_porosity (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p_c (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_u_f (this->get_triangulation().n_active_cells());

      ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                          dim+4);
      ComponentSelectFunction<dim> comp_p(dim+1, dim+4);
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
                                         &comp_p);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
    		                             compaction_pressure,
                                         ref_func_fluid,
                                         cellwise_errors_p_c,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_p);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         ref_func,
                                         cellwise_errors_porosity,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_porosity);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         melt_velocity,
                                         ref_func_fluid,
                                         cellwise_errors_u_f,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_u);

      std::ostringstream os;
      os << std::scientific << cellwise_errors_u.l2_norm()
         << ", " << cellwise_errors_p.l2_norm()
         << ", " << cellwise_errors_p_c.l2_norm()
         << ", " << cellwise_errors_porosity.l2_norm()
         << ", " << cellwise_errors_u_f.l2_norm();

      return std::make_pair("Errors u_L2, p_fL2, p_cL2, porosity_L2, u_fL2:", os.str());
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
	std::vector<double> & output
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
	      output[q] = (E * std::exp(y) - rho_f_0 * C) * std::exp(-y)/ (-C );
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

    ASPECT_REGISTER_FLUID_PRESSURE_BOUNDARY_CONDITIONS(PressureBdry,
						       "PressureBdry",
						       "A fluid pressure boundary condition that prescribes the "
						       "gradient of the fluid pressure at the boundaries as "
						       "calculated in the analytical solution. ")
						       
}
