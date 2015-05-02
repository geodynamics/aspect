#include <aspect/material_model/melt_interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/fe_field_function.h>


namespace aspect
{
  /**
   * This is the melt transport benchmark defined in the following paper:
   * @code
   *  @Article{DMGT11,
   *    author =       {H. Schmeling et al.},
   *    title =        {A community benchmark on mantle convection with melting and 
                        melt segregation},
   *    journal =      {},
   *    year =         2015,
   *    volume =       XXX(X),
   *    pages =        {XX}
   * @endcode
   */
  namespace SchmelingBenchmark
  {
    using namespace dealii;

    /**
     * @note 
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class SchmelingMeltMaterial : public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
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
          else if ((dependence & MaterialModel::NonlinearDependence::temperature) != MaterialModel::NonlinearDependence::none)
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
          return false;
        }

        virtual double reference_viscosity () const
        {
          return eta_0;
        }

        virtual double reference_density () const
        {
          return reference_rho_s;
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
          std::vector<double> maximum_melt_fractions(in.position.size());

          // we want to get the depletion field from the old solution here,
          // because it tells us how much of the material was already molten
          if(this->include_melt_transport() && in.cell != this->get_dof_handler().end())
            {
              // Prepare the field function
              Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector>
                fe_value(this->get_dof_handler(), this->get_old_solution(), this->get_mapping());

              // only get depletion field from the old the solution
              Assert(this->introspection().compositional_name_exists("depletion"),
                     ExcMessage("Material model schmeling melt only works with melt transport if there is a "
                                "compositional field called depletion."));
              const unsigned int depletion_idx = this->introspection().compositional_index_for_name("depletion");

              fe_value.set_active_cell(in.cell);
              fe_value.value_list(in.position,
                                  maximum_melt_fractions,
                                  this->introspection().component_indices.compositional_fields[depletion_idx]);
            }

          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              // get the porosity and depletion of material if they exist
              double porosity = 0.0;
              if(this->introspection().compositional_name_exists("porosity"))
                {
                  const double porosity_idx = this->introspection().compositional_index_for_name("porosity");
                  porosity = in.composition[i][porosity_idx];
                }

              double depletion = 0.0;
              if(this->introspection().compositional_name_exists("depletion"))
                {
                  const double depletion_idx = this->introspection().compositional_index_for_name("depletion");
                  depletion = in.composition[i][depletion_idx];
                }

              // calculate all the simple material parameters
              // we have to scale the density to account for the use of the Boussinesq approximation
              // (we want to have an almost constant density everywhere except in the buoyancy term)
              out.viscosities[i] = eta_0;
              out.densities[i] = reference_rho_s
                                 + density_scaling * (- thermal_expansivity * (in.temperature[i] - reference_T)
                                                      + depletion * (1.0 - porosity) * (rho_depletion - reference_rho_s));
              if (!this->include_melt_transport())
                out.densities[i] += density_scaling * porosity * (reference_rho_f - reference_rho_s);

              out.thermal_expansion_coefficients[i] = thermal_expansivity;
              out.specific_heat[i] = 1.0;
              out.thermal_conductivities[i] = 1.0;
              out.compressibilities[i] = 0.0;

              // calculate the solidus for the melt fraction and latent heat
              const double delta_T_solidus = (false //this->include_melt_transport()
                                              ?
                                              (depletion - porosity) / (Delta_T_liquidus_solidus)
                                              :
                                              0.0);

              const double T_solidus = A1 + A2 * (1.0 - in.position[i](1)) + delta_T_solidus;
              const double c_melt = ((in.temperature[i] > T_solidus) && (in.temperature[i] < T_solidus + Delta_T_liquidus_solidus))
                                    ?
                                    1.0
                                    :
                                    0.0;
              const double gravity = this->get_gravity_model().gravity_vector(in.position[i]).norm();

              // reaction terms are the melting rates of material
              double melt_fraction = 0.0;
              if ((in.temperature[i] > T_solidus) && (in.temperature[i] < T_solidus + Delta_T_liquidus_solidus))
                melt_fraction = (in.temperature[i] - T_solidus) / Delta_T_liquidus_solidus;
              else if (in.temperature[i] > T_solidus + Delta_T_liquidus_solidus)
                melt_fraction = 1.0;


              // we have to use the depletion here for both cases, and this doesn't work with iterated IMPES
              // (for that have a look at melt_simple material model)
              for (unsigned int c=0;c<in.composition[i].size();++c)
                {
                  if (this->introspection().name_for_compositional_index(c) == "porosity")
                    {
                      if(this->introspection().compositional_name_exists("depletion"))
                        {
                          const double porosity_idx = this->introspection().compositional_index_for_name("depletion");
                          out.reaction_terms[i][c] = melt_fraction - in.composition[i][porosity_idx];

                          if (this->include_melt_transport() && this->get_timestep_number() > 0)
                            out.reaction_terms[i][c] = c_melt * (melt_fraction - maximum_melt_fractions[i])
                                                       * reference_rho_s  / this->get_timestep();
                        }
                      else
                        out.reaction_terms[i][c] = melt_fraction - in.composition[i][c];
                    }
                  else if (this->introspection().name_for_compositional_index(c) == "depletion")
                    {
                      if (this->include_melt_transport() && this->get_timestep_number() > 0)
                        out.reaction_terms[i][c] = c_melt * (melt_fraction - maximum_melt_fractions[i]);
                      else
                        out.reaction_terms[i][c] = melt_fraction - in.composition[i][c];
                    }
                  else
                    out.reaction_terms[i][c] = 0.0;
                }

              // we have to scale the latent heat terms with the temperature and gravity, respectively, because the
              // terms in the benchmark are defined in terms of depth and in nondimensional form, so we have to scale
              // scale them back to the dimensions we use
              // also note that we have a different sign in the pressure derivative compared to the manuscript,
              // resulting from the opposite directions of the z axis
              out.entropy_derivative_temperature[i] = - c_melt * latent_heat / (Delta_T_liquidus_solidus * in.temperature[i]);
              out.entropy_derivative_pressure[i] = A2 * c_melt * latent_heat / (Delta_T_liquidus_solidus * in.temperature[i] * reference_rho_s * gravity);
            }
        }

        virtual void evaluate_with_melt(const typename MaterialModel::MeltInterface<dim>::MaterialModelInputs &in,
                                        typename MaterialModel::MeltInterface<dim>::MaterialModelOutputs &out) const
        {
          evaluate(in, out);
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              double porosity = in.composition[i][porosity_idx];

              out.compaction_viscosities[i] = xi_0;
              out.fluid_viscosities[i]= eta_f;
              out.permeabilities[i]= std::max(reference_permeability * std::pow(porosity,3),0.0);
              out.fluid_densities[i]= reference_rho_f - thermal_expansivity * (in.temperature[i] - reference_T);
              out.fluid_compressibilities[i] = 0.0;
            }

        }

      private:
        double reference_rho_s;
        double reference_rho_f;
        double rho_depletion;
        double eta_0;
        double xi_0;
        double eta_f;
        double reference_permeability;
        double thermal_expansivity;
        double reference_T;
        double Delta_T_liquidus_solidus;
        double A1;
        double A2;
        double density_scaling;
        double latent_heat;
    };

    template <int dim>
    void
    SchmelingMeltMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Schmeling melt");
        {
          prm.declare_entry ("Reference solid density", "1.0",
                             Patterns::Double (0),
                             "Reference density of the solid $\\rho_{s,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference melt density", "1.0",
                             Patterns::Double (),
                             "Reference density of the melt $\\rho_{f,0}$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference depleted density", "1.0",
                             Patterns::Double (),
                             "Reference density of the depleted material. Units: $kg/m^3$.");
          prm.declare_entry ("Density scaling", "1e-3",
                             Patterns::Double (0),
                             "The scaling factor the buoyancy term is divided by in order to"
                             "account for the use of Boussinesq approximation instead of variable"
                             "density everywhere. Units: $None$.");
          prm.declare_entry ("Thermal expansion coefficient", "1.0",
                             Patterns::Double (0),
                             "The value of the constant thermal expansivity $\\alpha$. "
                             "Units: $1/Ks$.");
          prm.declare_entry ("Reference temperature", "0.0",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: $K$.");
          prm.declare_entry ("Reference shear viscosity", "1e-4",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
                             "Units: $Pa s$.");
          prm.declare_entry ("Reference compaction viscosity", "1e-4",
                             Patterns::Double (0),
                             "The value of the constant volumetric viscosity $\\xi_0$ of the solid matrix. "
                             "Units: $Pa s$.");
          prm.declare_entry ("Reference melt viscosity", "1.0",
                             Patterns::Double (0),
                             "The value of the constant melt viscosity $\\eta_f$. Units: $Pa s$.");
          prm.declare_entry ("Reference permeability", "0.02",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: $m^2$.");
          prm.declare_entry ("Temperature difference liquidus solidus", "0.4",
                             Patterns::Double(),
                             "Fixed temperature difference between liquidus and solidus "
                             "(does not depend on composition / depletion)."
                             "Units: $T$.");
          prm.declare_entry ("A1", "0.5769",
                             Patterns::Double (0),
                             "The constant A1 in the solidus law. "
                             "Units: None.");
          prm.declare_entry ("A2", "0.4731",
                             Patterns::Double (0),
                             "The prefactor A2 in the solidus law. "
                             "Units: None.");
          prm.declare_entry ("Latent heat of melting", "0.256016385",
                             Patterns::Double (0),
                             "The latent heat change L across the phase transition from "
                             "solid to fluid. This is related to the Stephan number by "
                             "$St = \frac{c_p Delta T}{L}$. "
                             "Units: $J / kg$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SchmelingMeltMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Schmeling melt");
        {
          reference_rho_s            = prm.get_double ("Reference solid density");
          reference_rho_f            = prm.get_double ("Reference melt density");
          rho_depletion              = prm.get_double ("Reference depleted density");
          thermal_expansivity        = prm.get_double ("Thermal expansion coefficient");
          reference_T                = prm.get_double ("Reference temperature");
          eta_0                      = prm.get_double ("Reference shear viscosity");
          xi_0                       = prm.get_double ("Reference compaction viscosity");
          eta_f                      = prm.get_double ("Reference melt viscosity");
          reference_permeability     = prm.get_double ("Reference permeability");
          Delta_T_liquidus_solidus   = prm.get_double ("Temperature difference liquidus solidus");
          A1                         = prm.get_double ("A1");
          A2                         = prm.get_double ("A2");
          density_scaling            = prm.get_double ("Density scaling");
          latent_heat                = prm.get_double ("Latent heat of melting");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace SchmelingBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SchmelingMeltMaterial,
                                   "schmeling melt",
                                   "A material model that corresponds to the melt transport benchmark "
                                   "defined in Schmeling et al., 2015.")

  }
}
