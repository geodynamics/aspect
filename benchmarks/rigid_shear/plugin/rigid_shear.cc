/*
  Copyright (C) 2019 - 2022 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef aspect_rigid_shear_h
#define aspect_rigid_shear_h

#include <aspect/material_model/interface.h>
#include <aspect/postprocess/particles.h>
#include <aspect/particle/property/interface.h>
#include <aspect/particle/world.h>
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
  using namespace dealii;

  namespace RigidShearBenchmark
  {
    /**
     * This is the "Rigid shear" benchmark based on a suggestion in the following paper:
     * @code
     *  @Article{KKSCND97,
     *  author =       {P. E. van Keken and S. D. King and H. Schmeling and U. R. Christensen and D. Neumeister and M.-P. Doin},
     *  title =        {A comparison of methods for the modeling of thermochemical convection},
     *  journal =      {J. Geoph. Res.},
     *  year =         1997,
     *  volume =       102,
     *  pages =        {22477--22495}}
     * @endcode
     *
     * The results of the modification are published in Gassmoeller,
     * Lokavarapu, Bangerth, Puckett, 2019, "Evaluating the Accuracy of Hybrid
     * Finite Element/Particle-In-Cell Methods for Modeling Incompressible
     * Stokes Flow".
     */
    namespace AnalyticSolutions
    {
      /**
       * The exact solution for the Rigid Shear benchmark.
       */
      template<int dim>
      class FunctionRigidShear : public Function<dim>
      {
        public:
          FunctionRigidShear(unsigned int n_components) : Function<dim>(n_components) {}

          virtual void vector_value(const Point<dim> &p,
                                    Vector<double> &values) const
          {
            const double pi = numbers::PI;
            values[0] = std::sin(pi*p[0]) * std::cos(pi*p[1]);
            values[1] = -std::cos(pi*p[0]) * std::sin(pi*p[1]);
            values[2] = 2.0 * pi * std::cos(pi*p[0]) * std::cos(pi*p[1]);
            values[4] = std::sin(pi*p[0]) * std::sin(pi*p[1]);
            return;
          }
      };
    }



    /**
     * A material model for the stationary form of the rigid shear benchmark. All properties
     * are defined in dependence of position.
     */
    template<int dim>
    class RigidShearMaterial : public MaterialModel::Interface<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
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

        void
        parse_parameters(ParameterHandler &/*prm*/)
        {
          this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
          this->model_dependence.density = MaterialModel::NonlinearDependence::none;
          this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
          this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
          this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
        }
    };



    /**
     * A material model for the time-dependent form of the benchmark.
     * The density depends on the composition (e.g. advected by particles).
     */
    template<int dim>
    class RigidShearMaterialTimeDependent : public RigidShearMaterial<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
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

        void
        parse_parameters(ParameterHandler &/*prm*/)
        {
          this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
          this->model_dependence.density = MaterialModel::NonlinearDependence::none;
          this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
          this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
          this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
        }
    };



    /**
     * A particle property that represents the force term in the Stokes equation (density * gravity).
     * The reason this is implemented as a separate particle property is that the gravity becomes
     * infinite at x=0 and x=1, while (density * gravity) always remains finite. Therefore,
     * interpolating (density * gravity) does not break down at the model edges, while interpolating
     * density and then multiplying with the analytical density does.
     */
    template <int dim>
    class RigidShearForcingTerm : public aspect::Particle::Property::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual
        void
        initialize_one_particle_property (const Point<dim> &position,
                                          std::vector<double> &particle_properties) const
        {
          const auto &property_manager = this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Particles<dim>>().
                                         get_particle_world().
                                         get_property_manager();
          const unsigned density_index = property_manager.get_data_info().get_field_index_by_name("function");

          const double g_y = -4.0 * numbers::PI * numbers::PI * std::cos(numbers::PI * position[0]) / std::sin(numbers::PI * position[0]);
          particle_properties.push_back(particle_properties[density_index]*g_y);
        }

        virtual
        void
        update_particle_property (const unsigned int data_position,
                                  const Vector<double> &/*solution*/,
                                  const std::vector<Tensor<1,dim>> &/*gradients*/,
                                  typename Particles::ParticleHandler<dim>::particle_iterator &particle) const
        {
          const auto &property_manager = this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Particles<dim>>().
                                         get_particle_world().
                                         get_property_manager();
          const unsigned density_index = property_manager.get_data_info().get_field_index_by_name("function");

          const Point<dim> position = particle->get_location();
          const double g_y = -4.0 * numbers::PI * numbers::PI * std::cos(numbers::PI * position[0]) / std::sin(numbers::PI * position[0]);
          particle->get_properties()[data_position] = particle->get_properties()[density_index]*g_y;
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
        std::vector<std::pair<std::string, unsigned int>>
        get_property_information() const
        {
          std::vector<std::pair<std::string,unsigned int>> property_information;

          const std::string field_name = "forcing_term_y";
          property_information.emplace_back(field_name,1);
          return property_information;
        }
    };




    /**
     * A postprocessor that evaluates the accuracy of the solution.
     *
     * The implementation of error evaluators that correspond to the
     * benchmarks defined in the paper Gassmoeller et al. referenced above.
     */
    template<int dim>
    class RigidShearPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate error output for velocity, pressure, and density.
         */
        virtual
        std::pair<std::string, std::string>
        execute(TableHandler &/*statistics*/)
        {
          AnalyticSolutions::FunctionRigidShear<dim> ref_func(this->introspection().n_components);

          AssertThrow(Plugins::plugin_type_matches<const RigidShearMaterial<dim>>(this->get_material_model()),
                      ExcMessage(
                        "Postprocessor RigidShearPostprocessor only works with the material model RigidShearMaterial."));

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
                                   "A material model for the stationary form of the rigid shear benchmark. All properties "
                                   "are defined in dependence of position.")

    ASPECT_REGISTER_MATERIAL_MODEL(RigidShearMaterialTimeDependent,
                                   "rigid shear time dependent",
                                   "A material model for the time-dependent form of the rigid shear benchmark. "
                                   "The density depends on the composition (e.g. advected by particles).")

    ASPECT_REGISTER_PARTICLE_PROPERTY(RigidShearForcingTerm,
                                      "rigid shear forcing term",
                                      "A particle property that represents the force term in the Stokes equation (density * gravity). "
                                      "The reason this is implemented as a separate particle property is that the gravity becomes "
                                      "infinite at x=0 and x=1, while (density * gravity) always remains finite. Therefore, "
                                      "interpolating (density * gravity) does not break down at the model edges, while interpolating "
                                      "density and then multiplying with the analytical gravity does.")

    ASPECT_REGISTER_POSTPROCESSOR(RigidShearPostprocessor,
                                  "rigid shear",
                                  "The implementation of error evaluators that correspond to the "
                                  "benchmarks defined in the paper Gassmoeller et al. referenced above.")
  }
}
#endif
