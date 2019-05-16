#ifndef ASPECT_RIGID_SHEAR_H
#define ASPECT_RIGID_SHEAR_H

#include <aspect/material_model/interface.h>
#include <aspect/postprocess/particles.h>
#include <aspect/particle/property/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/postprocess/interface.h>
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
  namespace RigidShearBenchmark
  {
  /**
   * This is the "Sol Kz" benchmark defined in the following paper:
   * @code
   *  @Article{DMGT11,
   *    author =       {T. Duretz and D. A. May and T. V. Gerya and P. J. Tackley},
   *    title =        {Discretization errors and free surface stabilization in the
   *                  finite difference and marker-in-cell method for applied
   *                  geodynamics: {A} numerical study},
   *    journal =      {Geochemistry Geophysics Geosystems},
   *    year =         2011,
   *    volume =       12,
   *    pages =        {Q07004/1--26}}
   * @endcode
   *
   * The results are published in Kronbichler, Heister and Bangerth paper.
   */
  namespace AnalyticSolutions
  {
    using namespace dealii;

      /**
       * The exact solution for the SolKz benchmark.
       */
      template<int dim>
      class FunctionRigidShear : public Function<dim>
      {
        public:
          FunctionRigidShear(unsigned int n_components) : Function<dim>(n_components) {}

          virtual void vector_value(const Point<dim> &p,
                                    Vector<double> &values) const
          {
            const double pi = dealii::numbers::PI;
            values[0] = std::sin(pi*p[0]) * std::cos(pi*p[1]);
            values[1] = -std::cos(pi*p[0]) * std::sin(pi*p[1]);
            values[2] = 2.0 * pi * std::cos(pi*p[0]) * std::cos(pi*p[1]);
            values[4] = std::sin(pi*p[0]) * std::sin(pi*p[1]);
            return;
          }
      };
    }



    template<int dim>
    class RigidShearMaterial : public MaterialModel::Interface<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.position.size(); ++i)
            {
              const double x = in.position[i][0];
              const double y = in.position[i][1];
              const double pi = numbers::PI;

              out.densities[i] = std::sin(pi*x) * std::sin(pi*y);
              out.viscosities[i] = 1.0;
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0.0;
              out.thermal_conductivities[i] = 0.0;
            }
        }

        virtual bool is_compressible() const
        {
          return false;
        }

        virtual double reference_viscosity() const
        {
          return 1;
        }

        void
        parse_parameters(ParameterHandler &/*prm*/)
        {
          // Declare dependencies on solution variables
          this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
          this->model_dependence.density = MaterialModel::NonlinearDependence::none;
          this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
          this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
          this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
        }
    };



    template<int dim>
    class RigidShearMaterialTimeDependent : public RigidShearMaterial<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.position.size(); ++i)
            {
              out.densities[i] = in.composition[i][1];
              out.viscosities[i] = 1.0;
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0.0;
              out.thermal_conductivities[i] = 0.0;
            }
        }

        virtual bool is_compressible() const
        {
          return false;
        }

        virtual double reference_viscosity() const
        {
          return 1;
        }

        void
        parse_parameters(ParameterHandler &/*prm*/)
        {
          // Declare dependencies on solution variables
          this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
          this->model_dependence.density = MaterialModel::NonlinearDependence::none;
          this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
          this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
          this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
        }
    };



    template <int dim>
    class RigidShearForcingTerm : public aspect::Particle::Property::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual
        void
        initialize_one_particle_property (const Point<dim> &position,
                                          std::vector<double> &particle_properties) const
        {
          const auto &property_manager = this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Particles<dim> >().
              get_particle_world().
              get_property_manager();
          const unsigned density_index = property_manager.get_data_info().get_field_index_by_name("function");

          const double g_y = -4.0 * numbers::PI * numbers::PI * std::cos(numbers::PI * position[0]) / std::sin(numbers::PI * position[0]);
          particle_properties.push_back(particle_properties[density_index]*g_y);
        }

        virtual
        void
        update_one_particle_property (const unsigned int data_position,
                                      const Point<dim> &position,
                                      const Vector<double> &/*solution*/,
                                      const std::vector<Tensor<1,dim> > &/*gradients*/,
                                      const ArrayView<double> &particle_properties) const
        {
          const auto &property_manager = this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Particles<dim> >().
              get_particle_world().
              get_property_manager();
          const unsigned density_index = property_manager.get_data_info().get_field_index_by_name("function");

          const double g_y = -4.0 * numbers::PI * numbers::PI * std::cos(numbers::PI * position[0]) / std::sin(numbers::PI * position[0]);
          particle_properties[data_position] = particle_properties[density_index]*g_y;
        }

        aspect::Particle::Property::UpdateTimeFlags
        need_update () const
        {
          return aspect::Particle::Property::UpdateTimeFlags::update_time_step;
        }

        virtual
        UpdateFlags
        get_needed_update_flags () const
        {
          return update_default;
        }

        virtual
        std::vector<std::pair<std::string, unsigned int> >
        get_property_information() const
        {
          std::vector<std::pair<std::string,unsigned int> > property_information;

          const std::string field_name = "forcing_term_y";
          property_information.emplace_back(field_name,1);
          return property_information;
        }
    };




    /**
     * A postprocessor that evaluates the accuracy of the solution.
     *
     * The implementation of error evaluators that correspond to the
     * benchmarks defined in the paper Duretz et al. reference above.
     */
    template<int dim>
    class RigidShearPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        virtual
        std::pair<std::string, std::string>
        execute(TableHandler &/*statistics*/)
        {
          AnalyticSolutions::FunctionRigidShear<dim> ref_func(this->introspection().n_components);

          if (dynamic_cast<const RigidShearMaterial<dim> *>(&this->get_material_model()) == NULL)
            {
              AssertThrow(false,
                          ExcMessage(
                            "Postprocessor RigidShearPostprocessor only works with the material model RigidShearMaterial."));
            }

          const QGauss<dim> quadrature_formula(this->introspection().polynomial_degree.velocities + 2);

          Vector<float> cellwise_errors_u(this->get_triangulation().n_active_cells());
          Vector<float> cellwise_errors_p(this->get_triangulation().n_active_cells());
          Vector<float> cellwise_errors_ul2(this->get_triangulation().n_active_cells());
          Vector<float> cellwise_errors_pl2(this->get_triangulation().n_active_cells());
          Vector<float> cellwise_errors_rhol2(this->get_triangulation().n_active_cells());

          ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0, dim),
                                              this->get_fe().n_components());
          ComponentSelectFunction<dim> comp_p(dim,
                                              this->get_fe().n_components());
          ComponentSelectFunction<dim> comp_rho(dim+2,
                                                this->get_fe().n_components());

          VectorTools::integrate_difference(this->get_mapping(), this->get_dof_handler(),
                                            this->get_solution(),
                                            ref_func,
                                            cellwise_errors_u,
                                            quadrature_formula,
                                            VectorTools::L1_norm,
                                            &comp_u);
          VectorTools::integrate_difference(this->get_mapping(), this->get_dof_handler(),
                                            this->get_solution(),
                                            ref_func,
                                            cellwise_errors_p,
                                            quadrature_formula,
                                            VectorTools::L1_norm,
                                            &comp_p);
          VectorTools::integrate_difference(this->get_mapping(), this->get_dof_handler(),
                                            this->get_solution(),
                                            ref_func,
                                            cellwise_errors_ul2,
                                            quadrature_formula,
                                            VectorTools::L2_norm,
                                            &comp_u);
          VectorTools::integrate_difference(this->get_mapping(), this->get_dof_handler(),
                                            this->get_solution(),
                                            ref_func,
                                            cellwise_errors_pl2,
                                            quadrature_formula,
                                            VectorTools::L2_norm,
                                            &comp_p);
          VectorTools::integrate_difference(this->get_mapping(), this->get_dof_handler(),
                                            this->get_solution(),
                                            ref_func,
                                            cellwise_errors_rhol2,
                                            quadrature_formula,
                                            VectorTools::L2_norm,
                                            &comp_rho);

          const double u_l1 = Utilities::MPI::sum(cellwise_errors_u.l1_norm(), this->get_mpi_communicator());
          const double p_l1 = Utilities::MPI::sum(cellwise_errors_p.l1_norm(), this->get_mpi_communicator());
          const double u_l2 = std::sqrt(
                                Utilities::MPI::sum(cellwise_errors_ul2.norm_sqr(), this->get_mpi_communicator()));
          const double p_l2 = std::sqrt(
                                Utilities::MPI::sum(cellwise_errors_pl2.norm_sqr(), this->get_mpi_communicator()));
          const double rho_l2 = std::sqrt(
                                Utilities::MPI::sum(cellwise_errors_rhol2.norm_sqr(), this->get_mpi_communicator()));

          std::ostringstream os;
          os << std::scientific << u_l1
             << ", " << p_l1
             << ", " << u_l2
             << ", " << p_l2
             << ", " << rho_l2;

          return std::make_pair("Errors u_L1, p_L1, u_L2, p_L2, rho_L2:", os.str());
        }
    };
  }
}

// explicit instantiations
namespace aspect
{
  namespace RigidShearBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(RigidShearMaterial,
                                   "rigid shear",
                                   "")

    ASPECT_REGISTER_MATERIAL_MODEL(RigidShearMaterialTimeDependent,
        "rigid shear time dependent",
        "")

    ASPECT_REGISTER_PARTICLE_PROPERTY(RigidShearForcingTerm,
        "rigid shear forcing term",
        "")

    ASPECT_REGISTER_POSTPROCESSOR(RigidShearPostprocessor,
                                  "rigid shear",
                                  "")
  }
}
#endif
