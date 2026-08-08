/*
  Copyright (C) 2018 - 2024 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
 */

#include <aspect/geometry_model/interface.h>
#include <aspect/global.h>
#include <aspect/geometry_model/box.h>
#include <aspect/postprocess/interface.h>
#include <aspect/newton.h>

#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace aspect
{
  /**
   * This is the material model used for the viscoplastic convection
   * benchmark models that were presented in the open access Tosi et al. (2015) paper:
   * @code
   *  @Article{T15,
   *    Author = {Tosi, N. and Stein, C. and Noack, L. and H\"uttig, C. and Maierova, P. and Samuel, H. and Davies, D. R. and Wilson, C. R. and Kramer, S. C. and Thieulot, C. and Glerum, A. and Fraters, M. and Spakman, W. and Rozel, A. and Tackley, P. J.},
   *    Title = {A community benchmark for viscoplastic thermal convection in a 2-D square box},
   *    Journal = {Geochemistry, Geophysics, Geosystems},
   *    Year = {2015}
   *    Pages = {2175-2196},
   *    Volume = {16},
   *    Doi = {10.1002/2015GC005807}}
   * @endcode
   *
   * It is a material model that consists of globally constant values for all
   * material parameters except density and viscosity. Density is temperature dependent,
   * but practically constant through a small value of the thermal expansivity,
   * and viscosity depends on temperature, depth and strain rate.
   *
   * The model is considered incompressible, following the definition
   * described in Interface::is_compressible.
   *
   * @ingroup MaterialModels
   */
  namespace TosiRathmannBenchmark
  {
    /**
     * @ingroup MaterialModels
     */
    template <int dim>
    class TosiRathmannMaterial : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * also initializes cpo av model to set the assembler which can handle anisotropic vicosity 
         */
        void
        initialize() override; 

        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override
        {
          /**
           * As described in Tosi et al (2015), the viscosity \eta is computed as the
           * harmonic average of a linear and nonlinear plastic part. The rheology is ammended to use and anistrotropic non-linear anisotropic viscosity by Rathmann et. al. (2024)
           *
           * The linear part is calculated as follows
           * (see below for the meaning of the used parameter names):
           * $\eta_{lin}(T,z) = \exp(-\ln(\text{eta\_T} * T + \ln(\text{eta\_Z}) * z)$
           * while the strain rate dependent nonlinear part is computed as:
           * $\eta_{plast}(\dot\epsilon) = \text{eta\_asterisk} + \text{sigma\_yield} * \text{anisotropic\_viscosity}$
           */
          
          //set up additional output for the derivatives
          const std::shared_ptr<MaterialModel::MaterialModelDerivatives<dim>> derivatives
            = out.template get_additional_output_object<MaterialModel::MaterialModelDerivatives<dim>>();
          
          cpo_av->evaluate(in, out); 

          const std::vector<double> anisotropic_viscosity = out.viscosities;

          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              if (in.requests_property(MaterialModel::MaterialProperties::viscosity))
                {
                  out.viscosities[i] = viscosity (in.temperature[i],
                                                  in.pressure[i],
                                                  in.composition[i],
                                                  in.strain_rate[i],
                                                  in.position[i], 
                                                  anisotropic_viscosity[i]);
                }

              out.densities[i] = reference_rho * (1.0 - thermal_alpha * (in.temperature[i] - reference_T));
              out.thermal_expansion_coefficients[i] = thermal_alpha;
              out.specific_heat[i] = reference_specific_heat;
              out.thermal_conductivities[i] = thermal_k;
              out.compressibilities[i] = 0.0;
              // Pressure derivative of entropy at the given positions.
              out.entropy_derivative_pressure[i] = 0.0;
              // Temperature derivative of entropy at the given positions.
              out.entropy_derivative_temperature[i] = 0.0;
              // Change in composition due to chemical reactions at the
              // given positions. The term reaction_terms[i][c] is the
              // change in compositional field c at point i.
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                out.reaction_terms[i][c] = 0.0;

            }

        }

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.
         * Incompressibility does not necessarily imply that the density is
         * constant; rather, it may still depend on temperature or pressure.
         * In the current context, compressibility means whether we should
         * solve the continuity equation as $\nabla \cdot (\rho \mathbf u)=0$
         * (compressible Stokes) or as $\nabla \cdot \mathbf{u}=0$
         * (incompressible Stokes).
         */
        bool is_compressible () const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Wraps CPO_AV_3D::create_additional_named_outputs and initialize to generate
         * Those generate the anisotropic viscosity tensor also called stress strain director as an additional output
         * choose the correct assembler which can use this anisotropic viscosity tensor.  
         */
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;


      private:

        double viscosity (const double                  temperature,
                          const double                  pressure,
                          const std::vector<double>    &compositional_fields,
                          const SymmetricTensor<2,dim> &strain_rate,
                          const Point<dim>             &position,
                          const double                  aniso_viscosity) const;

        /*
         * Function to compute the linear viscosity
         * according to equation (7) of Tosi et al. 2015.
         */
        double viscolin(const double etaT,
                        const double etaZ,
                        const double T,
                        const double depth) const;

        /*
         * Function to compute the plastic viscosity
         * according to equation (19) of Rathmann et.al..
         */
        double viscoplast(const double eta_asterisk,
                          const double aniso_visc) const;

        /**
         * The density at reference temperature
         */
        double reference_rho;

        /*
         * The reference temperature
         */
        double reference_T;

        /*
         * The thermal expansivity
         */
        double thermal_alpha;

        /*
         * The reference specific heat
         */
        double reference_specific_heat;

        /**
         * The thermal conductivity.
         */
        double thermal_k;

        /*
         * The linear viscosity parameter pertaining to
         * the viscosity contrast due to temperature
         */
        double eta_T;

        /*
         * The linear viscosity parameter pertaining to
         * the viscosity contrast due to pressure (depth)
         */
        double eta_Z;

        /*
         * The effective viscosity at high stresses that is
         * part of the plastic viscosity
         */
        double eta_asterisk;

        /*
         * The lower viscosity cut-off value
         */
        double eta_minimum;

        /*
         * The upper viscosity cut-off value
         */
        double eta_maximum;

        /*
         * The initial guess of the viscosity
         */
        double eta_initial;

        /**
         * Pointer to cpo_av material model used as the base model
         */
        std::unique_ptr<MaterialModel::Interface<dim>> cpo_av;

    };

    /*
     * Function to calculate the viscosity according
     * to equation (6) of Tosi et.al.
     */
    template <int dim>
    double
    TosiRathmannMaterial<dim>::
    viscosity (const double temperature,
               const double,
               const std::vector<double> &,
               const SymmetricTensor<2,dim> &strain_rate,
               const Point<dim> &,
               const double aniso_viscosity) const
    {

      // In the first nonlinear iteration of the (pre-refinement steps of the) first time step,
      // strain rate is zero, so we set viscosity to eta_initial, a user-defined guess of the viscosity.
      if (strain_rate.norm() == 0)
        {
          return eta_initial;
        }

      // Otherwise we compute the linear viscosity and the plastic viscosity.
      double viscosity = 0.0;
      const double visc_linear = viscolin(eta_T,eta_Z,temperature,0);

      const double visc_plastic = viscoplast(eta_asterisk, aniso_viscosity);

      // Compute the harmonic average (equation (6) of Tosi et.al.)
      viscosity = 2.0 / ((1.0 / visc_linear) + (1.0 / visc_plastic));
  
      // Cut-off the viscosity by user-defined values to avoid possible very large viscosity ratios
      viscosity = std::max(std::min(viscosity,eta_maximum),eta_minimum);

      return viscosity;
    }

    /**
     * Function to compute the linear viscosity
     * according to equation (7) of Tosi et al. 2015.
     */
    template <int dim>
    double
    TosiRathmannMaterial<dim>::
    viscolin(const double etaT,
             const double etaZ,
             const double T,
             const double z) const
    {

      return std::exp((-1.0 * std::log(etaT) * T ) + (std::log(etaZ) * z));
    }

    /**
     * Function to compute the plastic viscosity
     * according to equation (19) of Rathmann et.al.
     */
    template <int dim>
    double
    TosiRathmannMaterial<dim>::
    viscoplast(const double etaasterisk,
               const double aniso_visc) const
    { 

      return etaasterisk + aniso_visc;
    }

    /**
     * create additional named outputs for Anisotropic tensor
     */
    template <int dim>
    void
    TosiRathmannMaterial<dim>::initialize()
    {
      cpo_av->initialize();
    }
    
    template <int dim>
    bool
    TosiRathmannMaterial<dim>::
    is_compressible () const
    {
      AssertThrow(!(cpo_av->is_compressible()), 
        ExcMessage("anisotropic material model is not allowed to be compressible"));
      return false;
    }

    template <int dim>
    void
    TosiRathmannMaterial<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      cpo_av->create_additional_named_outputs(out);
    }

    template <int dim>
    void
    TosiRathmannMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      // Default values are for Case 1 of Tosi et al. (2015).
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Tosi Rathmann benchmark");
        {
          prm.declare_entry ("Reference density", "1",
                             Patterns::Double (0),
                             "The value of the reference density $\\rho_0$.");
          prm.declare_entry ("Reference temperature", "0",
                             Patterns::Double (0),
                             "The value of the reference temperature $T_0$. The reference temperature is used "
                             "in the density calculation.");
          prm.declare_entry ("Minimum viscosity", "1e-6",
                             Patterns::Double (0),
                             "The value of the minimum cut-off viscosity $\\eta_min$.");
          prm.declare_entry ("Maximum viscosity", "1e1",
                             Patterns::Double (0),
                             "The value of the maximum cut-off viscosity $\\eta_max$.");
          prm.declare_entry ("Initial viscosity", "1e-1",
                             Patterns::Double (0),
                             "The value of the initial viscosity guess $\\eta_init$.");
          prm.declare_entry ("Thermal conductivity", "1",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$.");
          prm.declare_entry ("Reference specific heat", "1",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$.");
          prm.declare_entry ("Thermal expansion coefficient", "1e-6",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\alpha$.");
          prm.declare_entry ("Thermal viscosity parameter", "1e5",
                             Patterns::Double (0),
                             "The value of the thermal viscosity parameter $\\eta_T$, "
                             "as used in equation (7) of the paper.");
          prm.declare_entry ("Pressure viscosity parameter", "1e0",
                             Patterns::Double (0),
                             "The value of the pressure viscosity parameter $\\eta_Z$, "
                             "as used in equation (7) of the paper.");
          prm.declare_entry ("Nonlinear viscosity constant", "0",
                             Patterns::Double (0),
                             "The value of the plastic viscosity constant $\\eta_asterisk$, "
                             "as used in equation (8) of the paper.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    TosiRathmannMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Tosi Rathmann benchmark");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          thermal_k                  = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          eta_T                      = prm.get_double ("Thermal viscosity parameter");
          eta_Z                      = prm.get_double ("Pressure viscosity parameter");
          eta_asterisk               = prm.get_double ("Nonlinear viscosity constant");
          eta_minimum                = prm.get_double ("Minimum viscosity");
          eta_maximum                = prm.get_double ("Maximum viscosity");
          eta_initial                = prm.get_double ("Initial viscosity");

          
          //create cpo_av base model and initialize its SimulatorAccess base class 
          cpo_av = MaterialModel::create_material_model<dim>("CPO-induced anisotropic viscosity");
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(cpo_av.get()))
            sim->initialize_simulator (this->get_simulator());

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = MaterialModel::NonlinearDependence::strain_rate | MaterialModel::NonlinearDependence::temperature | MaterialModel::NonlinearDependence::pressure;
      this->model_dependence.density = MaterialModel::NonlinearDependence::temperature;
      this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
      this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;

      /**
      * base model syntax -> get model dependencies
      */
      cpo_av->parse_parameters(prm);
      this->model_dependence = cpo_av->get_model_dependence();

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace TosiRathmannBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(TosiRathmannMaterial,
                                   "TosiRathmannMaterial",
                                   "A material model that has constant values "
                                   "for all coefficients except the density and viscosity as described in the open access paper of Tosi et al. 2015. "
                                   "Ammended by including an anisotropic non-linear viscosity and an anisotropic tensor in the paper by Rathmann et.al. 2024"
                                   "The default parameter values are chosen according to Case 1 of the paper by Tosi et.al."
                                   "All of the values that define this model are read "
                                   "from a section ``Material model/Tosi Rathmann model'' in the input file, see "
                                   "Section~\\ref{parameters:Material_model/Tosi_Rathmann_model}."
                                   "\n\n"
                                   "This model uses the following set of equations for the two coefficients that "
                                   "are non-constant (see equation (6) - (10) of Tosi et.al. and (19) of Rathmann et.al.: "
                                   "\\begin{align}"
                                   "  \\eta(T,z,\\dot \\epsilon) &= 2\\frac{1}{\\frac{1}{\\eta_{lin}(T,z)}\\frac{1}{\\eta_{plast}(\\dot\\epsilon}}, \\\\"
                                   "  \\rho(T) &= \\left(1-\\alpha (T-T_0)\\right)\\rho_0,"
                                   "\\end{align}"
                                   "where $z$ represents depth."
                                   "\n\n"
                                   "The linear and plastic viscosity parts are defined as follows:"
                                   "\\begin{align}"
                                   "  \\eta_{lin}(T,z) &= \\exp(-\\ln(\\eta_T)T+\\ln(\\eta_z)z), \\\\"
                                   "  \\eta_{plast}(\\dot\\epsilon) &= \\eta^{*}+\\eta_0(\\dot\\epsilon) "
                                   "\\end{align} "
                                   "\n\n"
                                   "The constitutive relationship is ammended according to Rathmann et.al. equation (21), which one can reformulate to:"
                                   "\\sigma_{ij} = \\eta(T,z,\\dot \\epsilon) A_{ijkl}\\dot\\epsilon_{kl}"
                                   "with the anisotropic tensor from the CPO_AV_3D material model"
                                   "Note that this model uses the formulation that assumes an incompressible "
                                   "medium despite the fact that the density follows the law "
                                   "$\\rho(T)=\\rho_0(1-\\alpha(T-T_0))$. ")
  }
}
