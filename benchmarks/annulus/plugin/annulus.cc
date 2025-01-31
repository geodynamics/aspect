/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
#include <aspect/simulator.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/utilities.h>

#include <aspect/gravity_model/interface.h>
#include <aspect/postprocess/dynamic_topography.h>
#include <aspect/postprocess/visualization.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  /**
   * This is the "annulus" benchmark defined in the following paper:
   * @code
   *  @Article{thph17,
   *    author =       {Rene Gassm\"{o}ller, Juliane Dannberg, Wolfgang Bangerth, Elbridge Gerry Puckett, Cedric Thieulot},
   *    title =        {Benchmarking the accuracy of higher order particle methods in geodynamic models of transient flow},
   *    journal =      {xxx},
   *    year =         2023,
   *    volume =       xxx,
   *    number =       {x},
   *    publisher =    {xxx},
   *    pages =        {xxx--xxx}}
   * @endcode
   *
   */
  namespace AnnulusBenchmark
  {
    using namespace aspect::Utilities::Coordinates;

    /**
     * The exact solution for the Annulus benchmark.
     */
    namespace AnalyticSolutions
    {
      // Inner and outer radius are hard-coded to guarantee
      // tangential flow at the boundaries
      const double R1 = 1.0;
      const double R2 = 2.0;
      const double A = 2.0;
      const double B = -3.0/std::log(2.0);
      const double C=-1;
      const double rho_0 = 1000.;
      const double gravity = 1.;

      double
      phase(const double t)
      {
        return std::exp(t)-1;
      }



      double
      f(const double r)
      {
        return A*r + B/r;
      }



      double
      g(const double r)
      {
        return A*r/2 + B*std::log(r)/r + C/r;
      }



      double
      h(const double r)
      {
        return (2*g(r)-f(r))/r;
      }



      double
      m(const double r, const double k)
      {
        const double f_prime = A - B/(r*r);
        const double g_prime = A/2 - B*std::log(r)/(r*r) + B/(r*r) - C/(r*r);
        const double g_prime_prime = -B/(r*r*r)*(3-2*std::log(r)) - 2./(r*r*r);
        return g_prime_prime - g_prime/r - (g(r)*(k*k - 1))/(r*r) + f(r)/(r*r) + f_prime/r;
      }



      template <int dim>
      Tensor<1,dim>
      Annulus_velocity (const Point<dim> &pos,
                        const double k,
                        const double t,
                        const bool transient)
      {
        Assert (dim == 2, ExcNotImplemented());

        const std::array<double,dim> spherical_position = cartesian_to_spherical_coordinates(pos);
        const double r = spherical_position[0];
        const double theta = spherical_position[1];

        Tensor<1,dim> spherical_velocity;
        spherical_velocity[0] = g(r) * k * std::sin(k*(theta-phase(t)));
        spherical_velocity[1] = f(r) * std::cos(k*(theta-phase(t)));

        if (transient == true)
          spherical_velocity[1] += r * std::exp(t);

        return spherical_to_cartesian_vector(spherical_velocity, pos);
      }



      template <int dim>
      double
      Annulus_pressure (const Point<dim> &pos,
                        const double k,
                        const double t)
      {
        Assert (dim == 2, ExcNotImplemented());
        const std::array<double,dim> spherical_position = cartesian_to_spherical_coordinates(pos);
        const double r = spherical_position[0];
        const double theta = spherical_position[1];

        return k * h(r) * std::sin(k*(theta-phase(t)))
               + rho_0 * gravity * (R2-r);
      }



      template <int dim>
      double
      Annulus_normal_traction (const Point<dim> &pos,
                               const double k,
                               const double t)
      {
        Assert (dim == 2, ExcNotImplemented());
        const std::array<double,dim> spherical_position = cartesian_to_spherical_coordinates(pos);
        const double r = spherical_position[0];
        const double theta = spherical_position[1];

        return k * 3.*f(r)/r * std::sin(k*(theta-phase(t))) - rho_0 * g(r) * (R2 - r);
      }



      template <int dim>
      double
      Annulus_density (const Point<dim> &pos,
                       const double k,
                       const double t,
                       const bool transient)
      {
        Assert (dim == 2, ExcNotImplemented());
        const std::array<double,dim> spherical_position = cartesian_to_spherical_coordinates(pos);
        const double r = spherical_position[0];
        const double theta = spherical_position[1];

        double density = 0.0;
        if (transient == true)
          {
            const double stream_function = -((A/2)*r*r + B * std::log(r) + C) * std::cos(k*(theta-phase(t)));
            density = stream_function + rho_0;
          }
        else
          {
            density = m(r,k) * k * std::sin(k*theta) + rho_0;
          }

        return density;
      }



      template <int dim>
      Tensor<1,dim>
      Annulus_gravity (const Point<dim> &pos,
                       const double k,
                       const double t,
                       const bool transient)
      {
        Assert (dim == 2, ExcNotImplemented());

        double gravity_magnitude = 1.0;

        // see benchmark description for a justification. we want to change
        // the density, but keep the forcing term of the Stokes equation constant
        if (transient == true)
          {
            const std::array<double,dim> spherical_position = cartesian_to_spherical_coordinates(pos);
            const double r = spherical_position[0];
            const double theta = spherical_position[1];

            const double forcing_term = m(r,k)*k*std::sin(k*(theta-phase(t))) + rho_0;
            gravity_magnitude = forcing_term / Annulus_density(pos,k,t, transient);
          }

        const Tensor<1,dim> gravity_direction = -pos / pos.norm();
        return gravity_magnitude * gravity_direction;
      }



      template <int dim>
      class FunctionAnnulus : public Function<dim>
      {
        public:
          FunctionAnnulus (const double k,
                           const unsigned int n_comp,
                           const bool transient)
            :
            Function<dim>(n_comp),
            k_(k),
            n_components(n_comp),
            transient_(transient)
          {
            Assert(n_components >= 4, ExcInternalError());
          }

          void
          vector_value (const Point<dim>   &pos,
                        Vector<double>   &values) const override
          {
            Assert (dim == 2, ExcNotImplemented());
            Assert (values.size() >= n_components, ExcInternalError());

            const double t = this->get_time();

            const Tensor<1,dim> v = AnalyticSolutions::Annulus_velocity (pos, k_, t, transient_);
            values[0] = v[0];
            values[1] = v[1];
            values[2] = AnalyticSolutions::Annulus_pressure (pos, k_, t);
            values[3] = 0.0;

            if (n_components >= 5)
              values[4] = AnalyticSolutions::Annulus_density (pos, k_, t, transient_);
          }

        private:
          const double k_;
          const unsigned int n_components;
          const bool transient_;
      };
    }



    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class AnnulusMaterial : public MaterialModel::Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                 MaterialModel::MaterialModelOutputs<dim> &out) const override
        {
          double t = 0.0;
          if (this->simulator_is_past_initialization() == true)
            t = this->get_time();

          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              out.viscosities[i] = 1;
              if (use_analytical_density == true)
                {
                  out.densities[i] = AnalyticSolutions::Annulus_density(in.position[i], k, t, use_transient_flow_solution);
                }
              else
                {
                  out.densities[i] = in.composition[i][density_index];
                }

              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.thermal_conductivities[i] = 0.0;

              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                out.reaction_terms[i][c] = 0.0;
            }
        }



        bool is_compressible () const override
        {
          return false;
        }



        static
        void
        declare_parameters (ParameterHandler &prm)
        {
          // create a global section in the parameter file for parameters
          // that describe this benchmark. note that we declare them here
          // in the material model, but other kinds of plugins
          // may also read these parameters even though they do not
          // declare them.
          prm.enter_subsection("Annulus benchmark");
          {
            prm.declare_entry("k", "0",
                              Patterns::Double (0),
                              "Parameter k in the Annulus benchmark solution. This "
                              "parameter controls the number of convection cells in "
                              "the solution (which is 2*k).");

            prm.declare_entry("Use analytical density", "true",
                              Patterns::Bool(),
                              "Whether to use the analytical density solution, or to look for a "
                              "field named <density_field> to use as density.");

            prm.declare_entry("Use transient solution", "false",
                              Patterns::Bool(),
                              "Whether to use the transient flow solution, or to use the "
                              "default instantaneous solution.");
          }
          prm.leave_subsection();
        }



        void
        parse_parameters (ParameterHandler &prm) override
        {
          prm.enter_subsection("Annulus benchmark");
          {
            k = prm.get_double ("k");
            use_analytical_density = prm.get_bool ("Use analytical density");
            use_transient_flow_solution = prm.get_bool ("Use transient solution");
          }
          prm.leave_subsection();

          density_index = numbers::invalid_unsigned_int;
          if (use_analytical_density == false)
            {
              Assert(this->introspection().compositional_name_exists("density_field"),
                     ExcInternalError());
              density_index = this->introspection().compositional_index_for_name("density_field");
            }

          // Declare dependencies on solution variables
          this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
          this->model_dependence.density = MaterialModel::NonlinearDependence::none;
          this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
          this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
          this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
        }



        double get_k() const
        {
          return k;
        }



        bool use_transient_solution() const
        {
          return use_transient_flow_solution;
        }



        bool analytical_density() const
        {
          return use_analytical_density;
        }

      private:
        /**
         * Whether to use the analytical density.
         */
        bool use_analytical_density;

        /**
         * If not using the analytical density, store the
         * field index of the density field.
        */
        unsigned int density_index;

        /**
         * Number of positive and negative density anomalies,
         * resulting in 2*k convection cells.
         */
        double k;

        /**
         * Whether to use the transient flow solution.
        */
        bool use_transient_flow_solution;
    };



    /**
     * Velocity boundary condition for the Annulus benchmark
     */
    template <int dim>
    class AnnulusBoundary : public BoundaryVelocity::Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id ,
                           const Point<dim> &position) const override
        {
          Assert (dim==2, ExcNotImplemented());

          const AnnulusMaterial<dim> &
          material_model
            = Plugins::get_plugin_as_type<const AnnulusMaterial<dim>>(this->get_material_model());

          return AnalyticSolutions::Annulus_velocity (position,
                                                      material_model.get_k(),
                                                      this->get_time(),
                                                      material_model.use_transient_solution());
        }
    };



    /**
     * Gravity model for the Annulus benchmark
     */
    template <int dim>
    class AnnulusGravity : public aspect::GravityModel::Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        Tensor<1,dim>
        gravity_vector (const Point<dim> &pos) const override
        {
          const AnnulusMaterial<dim> &material_model
            = Plugins::get_plugin_as_type<const AnnulusMaterial<dim>>(this->get_material_model());

          double t = 0.0;
          if (this->simulator_is_past_initialization() == true)
            t = this->get_time();

          return AnalyticSolutions::Annulus_gravity(pos,
                                                    material_model.get_k(),
                                                    t,
                                                    material_model.use_transient_solution());
        }
    };


    /**
      * A postprocessor that evaluates the accuracy of the solution.
      */
    template <int dim>
    class AnnulusPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        void initialize() override
        {
          const AnnulusMaterial<dim> &material_model
            = Plugins::get_plugin_as_type<const AnnulusMaterial<dim>>(this->get_material_model());
          const unsigned int n_components = this->introspection().n_components;

          compute_density_error = (material_model.analytical_density() == false);
          reference_solution = std::make_unique<AnalyticSolutions::FunctionAnnulus<dim>>(material_model.get_k(),
                                                                                          n_components,
                                                                                          material_model.use_transient_solution());
        }



        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override
        {
          const unsigned int n_components = this->introspection().n_components;
          reference_solution->set_time(this->get_time());

          const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities+2);

          Vector<double> cellwise_errors_ul2 (this->get_triangulation().n_active_cells());
          Vector<double> cellwise_errors_pl2 (this->get_triangulation().n_active_cells());

          ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                              n_components);
          ComponentSelectFunction<dim> comp_p(dim, n_components);

          VectorTools::integrate_difference (this->get_mapping(),
                                             this->get_dof_handler(),
                                             this->get_solution(),
                                             *reference_solution,
                                             cellwise_errors_ul2,
                                             quadrature_formula,
                                             VectorTools::L2_norm,
                                             &comp_u);
          VectorTools::integrate_difference (this->get_mapping(),
                                             this->get_dof_handler(),
                                             this->get_solution(),
                                             *reference_solution,
                                             cellwise_errors_pl2,
                                             quadrature_formula,
                                             VectorTools::L2_norm,
                                             &comp_p);


          const double u_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_ul2, VectorTools::L2_norm);
          const double p_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_pl2, VectorTools::L2_norm);

          // we output the density error either way, but set it to 0 if we do not compute it.
          // this simplifies postprocessing the output.
          double rho_l2 = 0.0;
          if (compute_density_error == true)
            {
              Vector<double> cellwise_errors_rhol2 (this->get_triangulation().n_active_cells());
              ComponentSelectFunction<dim> comp_rho(dim+2, n_components);
              VectorTools::integrate_difference (this->get_mapping(),
                                                 this->get_dof_handler(),
                                                 this->get_solution(),
                                                 *reference_solution,
                                                 cellwise_errors_rhol2,
                                                 quadrature_formula,
                                                 VectorTools::L2_norm,
                                                 &comp_rho);
              rho_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_rhol2, VectorTools::L2_norm);
            }

          const double topo_l2 = compute_dynamic_topography_error();

          statistics.add_value ("u_L2",
                                u_l2);
          statistics.set_precision ("u_L2", 14);
          statistics.set_scientific ("u_L2", true);

          statistics.add_value ("p_L2",
                                p_l2);
          statistics.set_precision ("p_L2", 14);
          statistics.set_scientific ("p_L2", true);

          statistics.add_value ("rho_L2",
                                rho_l2);
          statistics.set_precision ("rho_L2", 14);
          statistics.set_scientific ("rho_L2", true);

          statistics.add_value ("topo_L2",
                                topo_l2);
          statistics.set_precision ("topo_L2", 14);
          statistics.set_scientific ("topo_L2", true);

          std::ostringstream os;

          os << std::scientific <<  u_l2
             << ", " << p_l2
             << ", " << rho_l2
             << ", " << topo_l2;

          return std::make_pair("Errors u_L2, p_L2, rho_L2, topo_L2:", os.str());
        }



        std::list<std::string>
        required_other_postprocessors() const override
        {
          return std::list<std::string> (1, "dynamic topography");
        }

      private:
        /**
         * Function that defines the reference solution to compare against.
         */
        std::unique_ptr<Function<dim>> reference_solution;

        /**
         * Whether to compute the density error. This only makes sense
         * if the material model does not use the analytical density solution,
         * otherwise the error is always 0.
         */
        bool compute_density_error;

        /**
         * Calculate the L2 dynamic topography error.
         */
        double
        compute_dynamic_topography_error() const
        {
          const Postprocess::DynamicTopography<dim> &dynamic_topography =
            this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::DynamicTopography<dim>>();

          const AnnulusMaterial<dim> &material_model
            = Plugins::get_plugin_as_type<const AnnulusMaterial<dim>>(this->get_material_model());
          const double k = material_model.get_k();

          const QGauss<dim-1> quadrature_formula (this->introspection().polynomial_degree.velocities+2);
          FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                           this->get_fe(),
                                           quadrature_formula,
                                           update_values | update_gradients |
                                           update_quadrature_points | update_JxW_values);
          LinearAlgebra::BlockVector topo_vector = dynamic_topography.topography_vector();
          std::vector<double> topo_values(quadrature_formula.size());

          double l2_error = 0.;

          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();

          for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
              if (cell->at_boundary())
                for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                  if (cell->face(f)->at_boundary())
                    {
                      fe_face_values.reinit(cell, f);
                      MaterialModel::MaterialModelInputs<dim> in_face(fe_face_values, cell, this->introspection(), this->get_solution());
                      MaterialModel::MaterialModelOutputs<dim> out_face(fe_face_values.n_quadrature_points, this->n_compositional_fields());
                      fe_face_values[this->introspection().extractors.temperature].get_function_values(topo_vector, topo_values);
                      in_face.requested_properties = MaterialModel::MaterialProperties::density;
                      material_model.evaluate(in_face, out_face);

                      for (unsigned int q=0; q < quadrature_formula.size(); ++q)
                        {
                          const Point<dim> p = fe_face_values.quadrature_point(q);
                          const double analytic_normal_stress = AnalyticSolutions::Annulus_normal_traction<dim>(p, k, this->get_time());
                          const double gravity = this->get_gravity_model().gravity_vector(p).norm();
                          const double density = out_face.densities[q];
                          const double diff = - analytic_normal_stress / gravity/ density - topo_values[q];
                          l2_error += (diff*diff) * fe_face_values.JxW(q);
                        }
                    }
          const double total_l2_error =  Utilities::MPI::sum(l2_error,this->get_mpi_communicator());
          return std::sqrt(total_l2_error);
        }
    };

    /**
      * A postprocessor that visualizes the analytical solution.
      */
    template <int dim>
    class AnnulusVisualizationPostprocessor : public DataPostprocessor<dim>,
      public Postprocess::VisualizationPostprocessors::Interface<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:
        AnnulusVisualizationPostprocessor ()
          :
          DataPostprocessor<dim> ()
        {};

        std::vector<std::string>
        get_names () const override
        {
          std::vector<std::string> solution_names;

          solution_names.push_back ("analytic_pressure");
          solution_names.push_back ("analytic_density");
          return solution_names;
        };

        std::vector<DataComponentInterpretation::DataComponentInterpretation>
        get_data_component_interpretation () const override
        {
          std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
          interpretation.push_back (DataComponentInterpretation::component_is_scalar);
          interpretation.push_back (DataComponentInterpretation::component_is_scalar);

          return interpretation;
        };

        UpdateFlags
        get_needed_update_flags () const override
        {
          return update_quadrature_points;
        };

        void
        evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                              std::vector<Vector<double>> &computed_quantities) const override
        {
          const unsigned int n_quadrature_points = input_data.solution_values.size();
          Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
          Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());


          const AnnulusMaterial<dim> &material_model
            = Plugins::get_plugin_as_type<const AnnulusMaterial<dim>>(this->get_material_model());

          for (unsigned int q=0; q<n_quadrature_points; ++q)
            {
              computed_quantities[q][0] = AnalyticSolutions::Annulus_pressure(input_data.evaluation_points[q],
                                                                              material_model.get_k(),
                                                                              this->get_time());
              computed_quantities[q][1] = AnalyticSolutions::Annulus_density(input_data.evaluation_points[q],
                                                                             material_model.get_k(),
                                                                             this->get_time(),
                                                                             material_model.use_transient_solution());
            }
        };
    };



    /**
     * A particle property that represents the density on particles.
     */
    template <int dim>
    class ParticleDensity : public aspect::Particle::Property::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        void
        initialize_one_particle_property (const Point<dim> &position,
                                          std::vector<double> &particle_properties) const override
        {
          const AnnulusMaterial<dim> &material_model
            = Plugins::get_plugin_as_type<const AnnulusMaterial<dim>>(this->get_material_model());

          const double density = AnalyticSolutions::Annulus_density(position,
                                                                    material_model.get_k(),
                                                                    /*time*/ 0.0,
                                                                    material_model.use_transient_solution());

          particle_properties.push_back(density);
        }



        std::vector<std::pair<std::string, unsigned int>>
        get_property_information() const override
        {
          return {{"particle_density",1}};
        }
    };



    /**
    * A initial condition for compositional fields that represents the initial density.
    */
    template <int dim>
    class AnnulusInitialDensity : public aspect::InitialComposition::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        double initial_composition (const Point<dim> &position,
                                    const unsigned int n_comp) const override
        {
          const double density_index = this->introspection().compositional_index_for_name("density_field");

          if (n_comp == density_index)
            {
              const AnnulusMaterial<dim> &material_model
                = Plugins::get_plugin_as_type<const AnnulusMaterial<dim>>(this->get_material_model());

              return AnalyticSolutions::Annulus_density(position,
                                                        material_model.get_k(),
                                                        /*time*/ 0.0,
                                                        material_model.use_transient_solution());
            }

          return 0.0;
        }
    };
  }
}



// explicit instantiations
namespace aspect
{
  namespace AnnulusBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(AnnulusMaterial,
                                   "AnnulusMaterial",
                                   "A material model that corresponds to the `Annulus' benchmark. "
                                   "See the manual for more information.")

    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(AnnulusBoundary,
                                            "AnnulusBoundary",
                                            "Implementation of the velocity boundary conditions for the "
                                            "`Annulus' benchmark. See the manual for more information about this "
                                            "benchmark.")

    ASPECT_REGISTER_POSTPROCESSOR(AnnulusPostprocessor,
                                  "AnnulusPostprocessor",
                                  "A postprocessor that compares the solution of the `Annulus' benchmark "
                                  "with the one computed by ASPECT "
                                  "and reports the error. See the manual for more information.")

    ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(AnnulusVisualizationPostprocessor,
                                                "AnnulusVisualizationPostprocessor",
                                                "A visualization output object that visualizes "
                                                "the analytical solution of the Annulus benchmark.")

    ASPECT_REGISTER_PARTICLE_PROPERTY(ParticleDensity,
                                      "particle density",
                                      "A particle property that represents the density.")

    ASPECT_REGISTER_GRAVITY_MODEL(AnnulusGravity,
                                  "annulus gravity",
                                  "A gravity model for the `Annulus' benchmark.")

    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(AnnulusInitialDensity,
                                              "initial density",
                                              "An initial condition for composition following the initial density.")
  }
}
