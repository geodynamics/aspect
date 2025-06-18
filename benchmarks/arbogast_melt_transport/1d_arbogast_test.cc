#include <aspect/material_model/interface.h>
#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>
#include <aspect/simulator.h>
#include <aspect/postprocess/visualization.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_postprocessor.h>

const double phi_plus = 0.001, L = 2.0;

namespace aspect
{
  template <int dim>
  class Arbogast:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
      virtual bool is_compressible () const override
      {
        return false;
      }

      virtual double reference_darcy_coefficient () const override
      {
        return std::pow(0.001, 2.0) / 1.0;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const override
      {
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        for (unsigned int i=0; i<in.position.size(); ++i)
          {
            const double porosity = std::max(0.0, in.composition[i][porosity_idx]);

            for (unsigned int c=0; c<in.composition[i].size(); ++c)
              out.reaction_terms[i][c] = 0.0;

            out.viscosities[i] = 1.0 * (1.0 - porosity);
            out.densities[i] = 2.0;
            out.thermal_expansion_coefficients[i] = 1.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim>>();

        if (melt_out != NULL)
          {
            for (unsigned int i=0; i<in.position.size(); ++i)
              {
                const double porosity = std::max(0.0, in.composition[i][porosity_idx]);

                melt_out->compaction_viscosities[i] = porosity>1e-9 ? 1.0 / porosity : std::numeric_limits<double>::max();
                melt_out->inverse_compaction_viscosities[i] = porosity;
                melt_out->fluid_viscosities[i]= 1.0;
                melt_out->permeabilities[i]= std::pow(porosity,2);
                melt_out->fluid_densities[i]= 1.0;
                melt_out->fluid_density_gradients[i] = Tensor<1,dim>();
              }
          }
      }

  };


  template <int dim>
  class RefFunction : public Function<dim>
  {
    public:
      RefFunction () : Function<dim>(2*dim+5) {}
      virtual void vector_value (const Point< dim >   &p,
                                 Vector< double >   &values) const override
      {
        const double z = p(1);
        const double r1 = (3.0 + sqrt(9.0 + 4.0 / phi_plus)) / 2.0;
        // const double r2 = (3.0 - sqrt(9.0 + 4.0 / phi_plus)) / 2.0;
        const double v_s = (z>0
                            ?
                            - phi_plus * phi_plus / (1.0 - 4.0 * phi_plus) * (std::pow(L, 4.0-r1) * std::pow(z, r1) - std::pow(z, 4))
                            :
                            0.0);
        const double q_f = (z>0
                            ?
                            1.0 / (1.0 - 4.0 * phi_plus) * (z - (std::pow(L, 4.0-r1) * std::pow(z, r1-3.0)) / (r1-3.0))
                            :
                            z);
        const double q_s = z;
        const double phi = (z>0 ? phi_plus*z *z : 0);

        values[0] = 0.0;                        // x velocity
        values[1] = v_s;                        // y velocity
        values[2] = q_f + z;                    // p_f
        values[3] = values[2] + (1.0 - phi) * (q_s - q_f);  // p_t = p_f + p_c
        values[4] = 0.0;                        // x melt vel
        values[5] = (phi>0 ? v_s * (1.0 - 1.0/phi) : v_s);              // y melt vel TODO: fix

        values[6] = q_s + z;  // p_s
        values[7] = 0; // T
        values[8] = phi; // porosity
      }
  };

  /**
    * A postprocessor that evaluates the accuracy of the solution
    * by using the L2 norm.
    */
  template <int dim>
  class ConvergencePostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      /**
       * Generate graphical output from the current solution.
       */
      virtual
      std::pair<std::string,std::string>
      execute (TableHandler &statistics) override;

  };

  template <int dim>
  std::pair<std::string,std::string>
  ConvergencePostprocessor<dim>::execute (TableHandler &)
  {
    RefFunction<dim> ref_func;
    const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

    const unsigned int n_total_comp = this->introspection().n_components;

    Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_c (this->get_triangulation().n_active_cells());

    ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                        n_total_comp);
    ComponentSelectFunction<dim> comp_p(dim, n_total_comp);
    ComponentSelectFunction<dim> comp_p_c(dim+1, n_total_comp);

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
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p_c,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p_c);

    std::ostringstream os;
    os << std::scientific
       << "ndofs= " << this->get_solution().size()
       << " u_L2= " << std::sqrt(Utilities::MPI::sum(cellwise_errors_u.norm_sqr(),MPI_COMM_WORLD))
       << " p_L2= "  << std::sqrt(Utilities::MPI::sum(cellwise_errors_p.norm_sqr(),MPI_COMM_WORLD))
       << " p_c_L2= "  << std::sqrt(Utilities::MPI::sum(cellwise_errors_p_c.norm_sqr(),MPI_COMM_WORLD))
       ;

    return std::make_pair("Errors", os.str());
  }


  /**
    * A postprocessor that evaluates reference function.
    */
  template <int dim>
  class VisualizeRefFunction : public Postprocess::VisualizationPostprocessors::Interface<dim>, public ::aspect::SimulatorAccess<dim>, public DataPostprocessor<dim>
  {
    public:
      virtual
      void
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const override;

      virtual std::vector<std::string> get_names () const override;

      virtual
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      get_data_component_interpretation () const override;

      virtual UpdateFlags get_needed_update_flags () const override;
  };

  template <int dim>
  void
  VisualizeRefFunction<dim>::
  evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                        std::vector<Vector<double>> &computed_quantities) const
  {
    RefFunction<dim> ref_func;

    const unsigned int n_quadrature_points = input_data.solution_values.size();
    Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
    Assert ((computed_quantities[0].size() == 9),
            ExcInternalError());
    Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());
    Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,  ExcInternalError());

    MaterialModel::MaterialModelInputs<dim> in(input_data,
                                               this->introspection());

    for (unsigned int q=0; q<n_quadrature_points; ++q)
      {
        Vector<double> function_values (9);
        ref_func.vector_value(in.position[q],function_values);

        for (unsigned int i=0; i<9; ++i)
          computed_quantities[q](i) = function_values[i];
      }
  }


  template <int dim>
  std::vector<std::string>
  VisualizeRefFunction<dim>::get_names () const
  {
    std::vector<std::string> names;
    names.push_back ("x_velocity_solution");
    names.push_back ("y_velocity_solution");
    names.push_back ("fluid_pressure_solution");
    names.push_back ("total_pressure_solution");
    names.push_back ("x_melt_velocity_solution");
    names.push_back ("y_melt_velocity_solution");
    names.push_back ("pressure_solution");
    names.push_back ("temperature_solution");
    names.push_back ("porosity_solution");
    return names;
  }


  template <int dim>
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  VisualizeRefFunction<dim>::get_data_component_interpretation () const
  {
    return
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      (9, DataComponentInterpretation::component_is_scalar);
  }



  template <int dim>
  UpdateFlags
  VisualizeRefFunction<dim>::get_needed_update_flags () const
  {
    return update_gradients | update_values | update_quadrature_points;
  }



  template <int dim>
  class PressureBdry:

    public BoundaryFluidPressure::Interface<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const types::boundary_id /*boundary_indicator*/,
        const MaterialModel::MaterialModelInputs<dim> &/*material_model_inputs*/,
        const MaterialModel::MaterialModelOutputs<dim> &/*material_model_outputs*/,
        const std::vector<Tensor<1,dim>> &normal_vectors,
        std::vector<double> &output
      ) const override
      {
        for (unsigned int q=0; q<output.size(); ++q)
          {
            Tensor<1,dim> gradient;
            output[q] = gradient * normal_vectors[q];
          }
      }



  };
}


// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_MATERIAL_MODEL(Arbogast,
                                 "1d arbogast material",
                                 "")

  ASPECT_REGISTER_POSTPROCESSOR(ConvergencePostprocessor,
                                "error calculation",
                                "")

  ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(VisualizeRefFunction,
                                              "reference function",
                                              "A postprocessor that visualizes the reference function.")

  ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(PressureBdry,
                                                "PressureBdry",
                                                "A fluid pressure boundary condition that prescribes the "
                                                "gradient of the fluid pressure at the boundaries as "
                                                "calculated in the analytical solution. ")
}
