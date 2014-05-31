#include <aspect/material_model/melt_interface.h>
#include <aspect/compositional_initial_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>



namespace aspect
{
  /**
   * This is the "Solitary wave" benchmark defined in the following paper:
   * @code
   *  @Article{DMGT11,
   *    author =       {T. Keller and D. A. May and B. J. P. Kaus},
   *    title =        {Numerical modelling of magma dynamics coupled
   *                    to tectonic deformation of lithosphere and crust},
   *    journal =      {Geophysical Journal International},
   *    year =         2013,
   *    volume =       195(3),
   *    pages =        {1406-1442}}
   * @endcode
   *
   */
  namespace SolitaryWaveBenchmark
  {
    using namespace dealii;

    namespace AnalyticSolutions
    {
      /**
       * @note The solitary wave solution only exists as a function x = func(phi)
       * and not phi = func(x), which is what we would like to have for describing
       * the shape of the wave. Thus, we calculate x = func(phi) for a range of phis
       * between the background porosity and the amplitude of the wave. In a next
       * step, we interpolate these values to the grid.
       *
       * @param phi The characteristic shape of the wave, with phi --> 1
       * for x --> +- infinity
       *
       * @param amplitude The amplitude of the solitary wave, which is always
       * greater than 1.
       */
       double solitary_wave_solution (const double phi, const double amplitude)
       {
         AssertThrow(phi > 1.0 && phi <= amplitude,
                     ExcMessage("The solitary wave solution can only be computed "
                         "for porosities larger than the background porosity of 1"
                         "and smaller than or equal to the amplitude of the wave."));
         AssertThrow(amplitude > 1,
                     ExcMessage("Amplitude of the solitary wave must be larger than 1!"));
         const double A_1   = std::sqrt(amplitude - 1.0);
         const double A_phi = std::sqrt(amplitude - phi);
         return std::sqrt(amplitude + 0.5)
                * (2 * A_phi - 1.0/A_1 * std::log((A_1 - A_phi)/(A_1 + A_phi)));
       }

       /**
        * This function gets the coordinate as an input parameters and gives
        * back the porosity of the solitary wave. As this function is only defined
        * implicitly, we have to interpolate from the coordinates where we have the
        * porosity to our mesh.
        *
        * @param x The coordinate vector where we want to interpolate the solitary
        * wave function to.
        * @param amplitude The amplitude of the solitary wave, which is always
        * greater than 1.
        * @param offset The offset of the center of the solitary wave from the
        * boundary of the domain.
        *
        * @param phi The porosity of the solitary wave function.
        */
       void compute_porosity (const std::vector<double> &grid_coordinate,
                              const double amplitude,
                              const double offset,
                              std::vector<double> &grid_porosity)
       {
         const unsigned int max_points = 1000;
         std::vector<double> porosity(max_points), coordinate(max_points);

         // get the coordinates where we have the solution
         for (unsigned int i=0;i<max_points;++i)
           {
             porosity[i] = 1.0 + 1e-10*amplitude
                           + double(i)/double(max_points-1) * (amplitude * (1.0 - 1e-10) - 1.0);
             coordinate[i] = solitary_wave_solution(porosity[i], amplitude);
           }

         // interpolate from the solution grid to the mesh used in the simulation
         // solitary wave is a monotonically decreasing function, so the coordinates
         // should be in descending order
         for (unsigned int i=0;i<grid_coordinate.size();++i)
           {
             // we only have the solution of the solitary wave for
             // coordinates larger than 0 (one half of the wave)
             const double x = (grid_coordinate[i] > offset
                 ?
                     grid_coordinate[i] - offset
                     :
                     offset - grid_coordinate[i]);

             AssertThrow(x < coordinate[0],
                         ExcMessage("The solitary wave solution can not be computed "
                             "for this coordinate."));

             unsigned int j= max_points-2;
             while (x < coordinate[j] && j>=0)
               j--;

             const double distance = (grid_coordinate[i] - coordinate[j+1])
                                      /(coordinate[j] - coordinate[j+1]);
             grid_porosity[i] = porosity[j+1] + distance /(porosity[j] - porosity[j+1]);
             }
         // TODO rescale
       }

       // vectors to store the porosity field and the corresponding coordinate in
       std::vector<double> porosity_vector;
       std::vector<double> position_vector;
     }


    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class SolitaryWaveMaterial : public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
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

      /**
       * Declare the parameters this class takes through input files.
       */
      static
      void
      declare_parameters (ParameterHandler &prm);

      /**
       * Read the parameters this class declares from the parameter file.
       */
      virtual
      void
      parse_parameters (ParameterHandler &prm);

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {

        for (unsigned int i=0;i<in.position.size();++i)
          {
            out.viscosities[i] = eta_0;
            out.densities[i] = reference_rho_s;
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 0.0;
            out.thermal_conductivities[i] = 0.0;
            out.compressibilities[i] = 0.0;
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

            out.compaction_viscosities[i] = xi_0;
            out.fluid_viscosities[i]= eta_f;
            out.permeabilities[i]= reference_permeability * std::pow(porosity,3);
            out.fluid_densities[i]= reference_rho_f;
            out.fluid_compressibilities[i] = 0.0;
          }

      }

      private:
        double reference_rho_s;
        double reference_rho_f;
        double eta_0;
        double xi_0;
        double eta_f;
        double reference_permeability;
    };

    template <int dim>
    void
    SolitaryWaveMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Solitary wave");
        {
          prm.declare_entry ("Reference solid density", "3000",
              Patterns::Double (0),
              "Reference density of the solid $\\rho_{s,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference melt density", "2500",
              Patterns::Double (0),
              "Reference density of the melt/fluid$\\rho_{f,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference shear viscosity", "1e20",
              Patterns::Double (0),
              "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
              "Units: $Pa s$.");
          prm.declare_entry ("Reference compaction viscosity", "1e20",
              Patterns::Double (0),
              "The value of the constant volumetric viscosity $\\xi_0$ of the solid matrix. "
              "Units: $Pa s$.");
          prm.declare_entry ("Reference melt viscosity", "100.0",
              Patterns::Double (0),
              "The value of the constant melt viscosity $\\eta_f$. Units: $Pa s$.");

          prm.declare_entry ("Reference permeability", "5e-9",
              Patterns::Double(),
              "Reference permeability of the solid host rock."
              "Units: $m^2$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SolitaryWaveMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Solitary wave");
        {
          reference_rho_s            = prm.get_double ("Reference solid density");
          reference_rho_f            = prm.get_double ("Reference melt density");
          eta_0                      = prm.get_double ("Reference shear viscosity");
          xi_0                       = prm.get_double ("Reference compaction viscosity");
          eta_f                      = prm.get_double ("Reference melt viscosity");
          reference_permeability     = prm.get_double ("Reference permeability");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    /**
     * An initial conditions model for the solitary waves benchmark.
     */
    template <int dim>
    class SolitaryWaveInitialCondition : public CompositionalInitialConditions::Interface<dim>,
                                         public ::aspect::SimulatorAccess<dim>
    {
      public:

      /**
       * Initialization function. Take references to the geometry model, the
       * object that describes the temperature boundary values, and the
       * adiabatic conditions and store them so that derived classes can
       * access them.
       */
      void
      initialize (const GeometryModel::Interface<dim>       &geometry_model);

      /**
       * Return the boundary velocity as a function of position.
       */
      virtual
      double
      initial_composition (const Point<dim> &position, const unsigned int n_comp) const;

      static
      void
      declare_parameters (ParameterHandler &prm);

      virtual
      void
      parse_parameters (ParameterHandler &prm);

      private:
        double amplitude;
        double background_porosity;
        double offset;
    };

    template <int dim>
    void
    SolitaryWaveInitialCondition<dim>::initialize (const GeometryModel::Interface<dim> &)
    {
      // initialize position and porosity vectors
      const unsigned int base_element = this->introspection().base_elements.compositional_fields;
      const FiniteElement<dim> &finite_element = this->get_fe();

      // get the composition support points
      const std::vector<Point<dim> > support_points
        = finite_element.base_element(base_element).get_unit_support_points();
      Assert (support_points.size() != 0,
              ExcInternalError());

      // create an FEValues object with just the composition element
      FEValues<dim> fe_values (this->get_mapping(), finite_element,
                               support_points,
                               update_quadrature_points);

      std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);
      const DoFHandler<dim> &dof_handler = this->get_dof_handler();

      unsigned int j=0;
      for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
           cell != dof_handler.end(); ++cell)
        {
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            // go through the composition dofs and set their global values
            // to the porosity wave solution field interpolated at these points
            cell->get_dof_indices (local_dof_indices);
            const unsigned int dofs_per_cell = finite_element.base_element(base_element).dofs_per_cell;
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              AnalyticSolutions::position_vector.push_back(fe_values.quadrature_point(i)[dim-1]);
          }
        ++j;
        }

      AnalyticSolutions::porosity_vector.resize(AnalyticSolutions::position_vector.size());

      // non-dimensionalize the amplitude and compute the porosity for every point
      const double non_dim_amplitude = amplitude / background_porosity;

      AnalyticSolutions::compute_porosity(AnalyticSolutions::position_vector,
                                          non_dim_amplitude,
                                          offset,
                                          AnalyticSolutions::porosity_vector);
      // TODO dimensionalize again
    }


    template <int dim>
    double
    SolitaryWaveInitialCondition<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      for (unsigned int i=0;i<AnalyticSolutions::position_vector.size();++i)
        if(position[dim-1] == AnalyticSolutions::position_vector[i])
          return AnalyticSolutions::porosity_vector[i];

      AssertThrow(false,
                  ExcMessage("No solitary wave solution was found for this position."));
      return std::numeric_limits<double>::quiet_NaN();
    }


    template <int dim>
    void
    SolitaryWaveInitialCondition<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Compositional initial conditions");
      {
        prm.enter_subsection("Solitary wave initial condition");
        {
          prm.declare_entry ("Amplitude", "0.01",
              Patterns::Double (0),
              "Amplitude of the solitary wave. Units: none.");
          prm.declare_entry ("Background porosity", "0.001",
              Patterns::Double (0),
              "Background porosity of the solitary wave. Units: none.");
          prm.declare_entry ("Offset", "1000",
              Patterns::Double (0),
              "Offset of the center of the solitary wave from the boundary"
              "of the domain. "
              "Units: $m$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SolitaryWaveInitialCondition<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Compositional initial conditions");
      {
        prm.enter_subsection("Solitary wave initial condition");
        {
          amplitude            = prm.get_double ("Amplitude");
          background_porosity  = prm.get_double ("Background porosity");
          offset               = prm.get_double ("Offset");

          AssertThrow(amplitude > background_porosity,
                      ExcMessage("Amplitude of the solitary wave must be larger "
                          "than the background porosity."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    /**
      * A postprocessor that evaluates the accuracy of the solution.
      *
      * The implementation of error evaluators that correspond to the
      * benchmarks defined in the paper Keller et al. reference above.
      */
    template <int dim>
    class SolitaryWavePostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
    SolitaryWavePostprocessor<dim>::execute (TableHandler &statistics)
    {
      AssertThrow(Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) == 1,
                  ExcNotImplemented());

      std_cxx1x::shared_ptr<Function<dim> > ref_func;
      if (dynamic_cast<const SolitaryWaveMaterial<dim> *>(&this->get_material_model()) != NULL)
        {
          const SolitaryWaveMaterial<dim> *
          material_model
            = dynamic_cast<const SolitaryWaveMaterial<dim> *>(&this->get_material_model());
        }
      else
        {
          AssertThrow(false,
                      ExcMessage("Postprocessor Solitary Wave only works with the material model Solitary wave."));
        }

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

      Vector<float> cellwise_errors_f (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());

      std::ostringstream os;
      os << std::scientific << 0.0
         << ", " << 0.0;

      return std::make_pair("Errors e_f, e_p:", os.str());
    }

  }
}



// explicit instantiations
namespace aspect
{
  namespace SolitaryWaveBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SolitaryWaveMaterial,
                                   "Solitary Wave",
                                   "A material model that corresponds to the 'SolitaryWave' benchmark "
                                   "defined in Keller et al., JGI, 2013.")

    ASPECT_REGISTER_POSTPROCESSOR(SolitaryWavePostprocessor,
                                  "SolitaryWavePostprocessor",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "the Keller et al., JGI, 2013, paper with the one computed by ASPECT "
                                  "and reports the error.")

   ASPECT_REGISTER_COMPOSITIONAL_INITIAL_CONDITIONS(SolitaryWaveInitialCondition,
                                                    "Solitary wave initial condition",
                                                    "Composition is set to a solitary wave function.")
  }
}
