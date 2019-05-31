/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#include <aspect/material_model/viscoelastic_plastic.h>
#include <aspect/utilities.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <numeric>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    ViscoelasticPlastic<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {

      // Store which components to exclude during volume fraction computation.
      ComponentMask composition_mask(this->n_compositional_fields(),true);
      // Elastic stress fields
      for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
        composition_mask.set(i,false);

      std::vector<double> average_elastic_shear_moduli (in.temperature.size());
      EquationOfStateOutputs<dim> eos_outputs (this->n_compositional_fields()+1);


      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const std::vector<double> composition = in.composition[i];
          const double pressure = in.pressure[i];
          const SymmetricTensor<2,dim> strain_rate = in.strain_rate[i];
          const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(composition, composition_mask);

          equation_of_state.evaluate(in, i, eos_outputs);

          // Arithmetic averaging of specific heat.
          // This may not be strictly the most reasonable thing, but for most Earth materials we hope
          // that they do not vary so much that it is a big problem. This statement also applies to
          // the arithmetic averaging of density and thermal conductivity below.
          out.densities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);
          out.specific_heat[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);

          out.thermal_conductivities[i] = MaterialUtilities::average_value(volume_fractions, thermal_conductivities, MaterialUtilities::arithmetic);

          // Set properties that are not relevant for this material model to 0
          out.compressibilities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);
          out.entropy_derivative_pressure[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_pressure, MaterialUtilities::arithmetic);
          out.entropy_derivative_temperature[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_temperature, MaterialUtilities::arithmetic);

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

          if (in.strain_rate.size())
            {
              // Calculate the square root of the second moment invariant for the deviatoric strain rate tensor.
              // The first time this function is called (first iteration of first time step)
              // a specified "reference" strain rate is used as the returned value would
              // otherwise be zero.
              const double edot_ii = ( (this->get_timestep_number() == 0 && strain_rate.norm() <= std::numeric_limits<double>::min() )
                                       ?
                                       reference_strain_rate
                                       :
                                       std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
                                                minimum_strain_rate) );

              const std::vector<double> viscosities_pre_yield = linear_viscosities;

              // TODO: Add strain-weakening of cohesion and friction
              const std::vector<double> coh = cohesions;
              const std::vector<double> phi = angles_internal_friction;

              // Initialize variables
              std::vector<double> stresses_viscous(volume_fractions.size());
              std::vector<double> stresses_yield(volume_fractions.size());
              std::vector<double> viscosities_viscoplastic(volume_fractions.size());
              std::vector<double> viscosities_viscoelastic(volume_fractions.size());

              // Loop through all compositions
              for (unsigned int j=0; j < volume_fractions.size(); ++j)
                {

                  // Note that edot_ii is the full computed strain rate, which includes elastic stresses
                  stresses_viscous[j] = 2. * viscosities_pre_yield[j] * edot_ii;

                  stresses_yield[j] = ( (dim==3)
                                        ?
                                        ( 6.0 * coh[j] * std::cos(phi[j]) + 6.0 * std::max(pressure,0.0) * std::sin(phi[j]) )
                                        / ( std::sqrt(3.0) * (3.0 + std::sin(phi[j]) ) )
                                        :
                                        coh[j] * std::cos(phi[j]) + std::max(pressure,0.0) * std::sin(phi[j]) );

                  // If the viscous stress is greater than the yield strength, rescale the viscosity back to yield surface.
                  // If the viscous stress is less than the yield stress, the yield viscosity is equal to the pre-yield value.
                  if ( stresses_viscous[j] >= stresses_yield[j]  )
                    {
                      viscosities_viscoplastic[j] = stresses_yield[j] / (2.0 * edot_ii);
                    }
                  else
                    {
                      viscosities_viscoplastic[j] = viscosities_pre_yield[j];
                    }
                }

              // Average viscosity
              const double average_viscosity = MaterialUtilities::calculate_average_vector<dim>(composition,
                                               viscosities_viscoplastic,
                                               viscosity_averaging);

              // Average elastic shear modulus
              std::vector<double> elastic_shear_moduli(elastic_rheology.get_elastic_shear_moduli());
              average_elastic_shear_moduli[i] = MaterialUtilities::calculate_average_vector<dim>(composition,
                                                elastic_shear_moduli,
                                                viscosity_averaging);

              // Average viscoelastic (e.g., effective) viscosity (equation 28 in Moresi et al., 2003, J. Comp. Phys.)
              out.viscosities[i] = elastic_rheology.calculate_viscoelastic_viscosity(average_viscosity,
                                                                                     average_elastic_shear_moduli[i]);

              // Fill the material properties that are part of the elastic additional outputs
              if (Rheology::ElasticAdditionalOutputs<dim> *elastic_out = out.template get_additional_output<Rheology::ElasticAdditionalOutputs<dim> >())
                {
                  elastic_out->elastic_shear_moduli[i] = average_elastic_shear_moduli[i];
                }
            }
        }

      elastic_rheology.fill_elastic_force_outputs(in, average_elastic_shear_moduli, out);
      elastic_rheology.fill_reaction_outputs(in, average_elastic_shear_moduli, out);
    }

    template <int dim>
    double
    ViscoelasticPlastic<dim>::
    reference_viscosity () const
    {
      return input_reference_viscosity;
    }

    template <int dim>
    bool
    ViscoelasticPlastic<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }

    template <int dim>
    void
    ViscoelasticPlastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Viscoelastic Plastic");
        {
          // Equation of state parameters
          EquationOfState::MulticomponentIncompressible<dim>::declare_parameters (prm);
          Rheology::Elasticity<dim>::declare_parameters (prm);

          // Reference and minimum/maximum values
          prm.declare_entry ("Minimum strain rate", "1.0e-20", Patterns::Double(0),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double(0),
                             "Reference strain rate for first time step. Units: $1 / s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0),
                             "Lower cutoff for effective viscosity. Units: $Pa \\, s$");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0),
                             "Upper cutoff for effective viscosity. Units: $Pa \\, s$");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double(0),
                             "Reference viscosity for nondimensionalization. "
                             "To understand how pressure scaling works, take a look at "
                             "\\cite{KHB12}. In particular, the value of this parameter "
                             "would not affect the solution computed by \\aspect{} if "
                             "we could do arithmetic exactly; however, computers do "
                             "arithmetic in finite precision, and consequently we need to "
                             "scale quantities in ways so that their magnitudes are "
                             "roughly the same. As explained in \\cite{KHB12}, we scale "
                             "the pressure during some computations (never visible by "
                             "users) by a factor that involves a reference viscosity. This "
                             "parameter describes this reference viscosity."
                             "\n\n"
                             "For problems with a constant viscosity, you will generally want "
                             "to choose the reference viscosity equal to the actual viscosity. "
                             "For problems with a variable viscosity, the reference viscosity "
                             "should be a value that adequately represents the order of "
                             "magnitude of the viscosities that appear, such as an average "
                             "value or the value one would use to compute a Rayleigh number."
                             "\n\n"
                             "Units: $Pa \\, s$");

          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $W/m/K$ ");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition "),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          prm.declare_entry ("Linear viscosities", "1.e21",
                             Patterns::List(Patterns::Double(0)),
                             "List of linear (fixed) viscosities for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $Pa s$");

          // Plasticity parameters
          prm.declare_entry ("Angles of internal friction", "0",
                             Patterns::List(Patterns::Double(0)),
                             "List of angles of internal friction, $\\phi$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "For a value of zero, in 2D the von Mises criterion is retrieved. "
                             "Angles higher than 30 degrees are harder to solve numerically. Units: degrees.");
          prm.declare_entry ("Cohesions", "1e20",
                             Patterns::List(Patterns::Double(0)),
                             "List of cohesions, $C$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The extremely large default cohesion value (1e20 Pa) prevents the viscous stress from "
                             "exceeding the yield stress. Units: $Pa$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    ViscoelasticPlastic<dim>::parse_parameters (ParameterHandler &prm)
    {

      // Get the number of fields for composition-dependent material properties
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Viscoelastic Plastic");
        {
          elastic_rheology.initialize_simulator (this->get_simulator());
          elastic_rheology.parse_parameters(prm);

          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm);

          AssertThrow(this->get_parameters().enable_elasticity == true,
                      ExcMessage ("Material model Viscoelastic plastic only works if 'Enable elasticity' is set to true"));

          minimum_strain_rate = prm.get_double("Minimum strain rate");
          reference_strain_rate = prm.get_double("Reference strain rate");
          minimum_viscosity = prm.get_double ("Minimum viscosity");
          maximum_viscosity = prm.get_double ("Maximum viscosity");
          input_reference_viscosity = prm.get_double ("Reference viscosity");

          viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);


          // Plasticity parameters
          angles_internal_friction = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Angles of internal friction"))),
                                                                             n_fields,
                                                                             "Angles of internal friction");
          // Convert angles from degrees to radians
          for (unsigned int i = 0; i<n_fields; ++i)
            angles_internal_friction[i] *= numbers::PI/180.0;
          cohesions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesions"))),
                                                              n_fields,
                                                              "Cohesions");

          // Parse additional material properties
          linear_viscosities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Linear viscosities"))),
                                                                       n_fields,
                                                                       "Viscosities");
          thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                                                                           n_fields,
                                                                           "Thermal conductivities");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();



      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::strain_rate | NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::compositional_fields;
      this->model_dependence.thermal_conductivity = NonlinearDependence::compositional_fields;
    }



    template <int dim>
    void
    ViscoelasticPlastic<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      elastic_rheology.create_elastic_outputs(out);
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ViscoelasticPlastic,
                                   "viscoelastic plastic",
                                   "A material model that combines non-linear plasticity with a simple "
                                   "linear viscoelastic material behavior. The model is incompressible. "
                                   "Note that this material model is based heavily on and combines "
                                   "functionality from the following material models: "
                                   "DiffusionDislocation, DruckerPrager, ViscoPlastic and Viscoelastic. "
                                   "\n\n"
                                   "Plasticity limits viscous stress through a Drucker Prager "
                                   "yield criterion, where the yield stress in 3D is  "
                                   "$\\sigma_y = \\frac{6*C*\\cos(\\phi) + 2*P*\\sin(\\phi)} "
                                   "{\\sqrt(3)*(3+\\sin(\\phi))}$ "
                                   "and "
                                   "$\\sigma_y = C\\cos(\\phi) + P\\sin(\\phi)$ "
                                   "in 2D. Above, $C$ is cohesion and $\\phi$  is the angle of "
                                   "internal friction.  Note that the 2D form is equivalent to the "
                                   "Mohr Coulomb yield surface.  If $\\phi$ is 0, the yield stress "
                                   "is fixed and equal to the cohesion (Von Mises yield criterion). "
                                   "When the viscous stress ($2v{\\varepsilon}_{ii}$) exceeds "
                                   "the yield stress, the viscosity is rescaled back to the yield "
                                   "surface: $v_{y}=\\sigma_{y}/(2{\\varepsilon}_{ii})$. "
                                   "This form of plasticity is commonly used in geodynamic models. "
                                   "See, for example, Thieulot, C. (2011), PEPI 188, pp. 47-68. "
                                   "\n\n"
                                   "The viscoelastic rheology behavior takes into account the elastic shear "
                                   "strength (e.g., shear modulus), while the tensile and volumetric "
                                   "strength (e.g., Young's and bulk modulus) are not considered. The "
                                   "model is incompressible and allows specifying an arbitrary number "
                                   "of compositional fields, where each field represents a different "
                                   "rock type or component of the viscoelastic stress tensor. The stress "
                                   "tensor in 2D and 3D, respectively, contains 3 or 6 components. The "
                                   "compositional fields representing these components must be named "
                                   "and listed in a very specific format, which is designed to minimize "
                                   "mislabeling stress tensor components as distinct 'compositional "
                                   "rock types' (or vice versa). For 2D models, the first three "
                                   "compositional fields must be labeled 'stress\\_xx', 'stress\\_yy' and 'stress\\_xy'. "
                                   "In 3D, the first six compositional fields must be labeled 'stress\\_xx', "
                                   "'stress\\_yy', 'stress\\_zz', 'stress\\_xy', 'stress\\_xz', 'stress\\_yz'. "
                                   "\n\n "
                                   "Combining this viscoelasticity implementation with non-linear viscous flow "
                                   "and plasticity produces a constitutive relationship commonly referred to "
                                   "as partial elastoviscoplastic (e.g., pEVP) in the geodynamics community. "
                                   "While extensively discussed and applied within the geodynamics "
                                   "literature, notable references include: "
                                   "Moresi et al. (2003), J. Comp. Phys., v. 184, p. 476-497. "
                                   "Gerya and Yuen (2007), Phys. Earth. Planet. Inter., v. 163, p. 83-105. "
                                   "Gerya (2010), Introduction to Numerical Geodynamic Modeling. "
                                   "Kaus (2010), Tectonophysics, v. 484, p. 36-47. "
                                   "Choi et al. (2013), J. Geophys. Res., v. 118, p. 2429-2444. "
                                   "Keller et al. (2013), Geophys. J. Int., v. 195, p. 1406-1442. "
                                   "\n\n "
                                   "The overview below directly follows Moresi et al. (2003) eqns. 23-32. "
                                   "However, an important distinction between this material model and "
                                   "the studies above is the use of compositional fields, rather than "
                                   "tracers, to track individual components of the viscoelastic stress "
                                   "tensor. The material model will be updated when an option to track "
                                   "and calculate viscoelastic stresses with tracers is implemented. "
                                   "\n\n "
                                   "Moresi et al. (2003) begins (eqn. 23) by writing the deviatoric "
                                   "rate of deformation ($\\hat{D}$) as the sum of elastic "
                                   "($\\hat{D_{e}}$) and viscous ($\\hat{D_{v}}$) components: "
                                   "$\\hat{D} = \\hat{D_{e}} + \\hat{D_{v}}$.  "
                                   "These terms further decompose into "
                                   "$\\hat{D_{v}} = \\frac{\\tau}{2\\eta}$ and "
                                   "$\\hat{D_{e}} = \\frac{\\overset{\\triangledown}{\\tau}}{2\\mu}$, where "
                                   "$\\tau$ is the viscous deviatoric stress, $\\eta$ is the shear viscosity, "
                                   "$\\mu$ is the shear modulus and $\\overset{\\triangledown}{\\tau}$ is the "
                                   "Jaumann corotational stress rate. If plasticity is included the deviatoric "
                                   "rate of deformation may be written as: "
                                   "$\\hat{D} = \\hat{D_{e}} + \\hat{D_{v}} + \\hat{D_{p}}$, where $\\hat{D_{p}}$ "
                                   "is the plastic component. As defined in the second paragraph, $\\hat{D_{p}}$ "
                                   "decomposes to $\\frac{\\tau_{y}}{2\\eta_{y}}$, where $\\tau_{y}$ is the yield "
                                   "stress and $\\eta_{y}$ is the viscosity rescaled to the yield surface. "
                                   "\n\n "
                                   "Above, the Jaimann corotational stress rate (eqn. 24) from the elastic "
                                   "component contains the time derivative of the deviatoric stress ($\\dot{\\tau}$) "
                                   "and terms that account for material spin (e.g., rotation) due to advection: "
                                   "$\\overset{\\triangledown}{\\tau} = \\dot{\\tau} + {\\tau}W -W\\tau$. "
                                   "Above, $W$ is the material spin tensor (eqn. 25): "
                                   "$W_{ij} = \\frac{1}{2} \\left (\\frac{\\partial V_{i}}{\\partial x_{j}} - "
                                   "\\frac{\\partial V_{j}}{\\partial x_{i}} \\right )$. "
                                   "\n\n "
                                   "The Jaumann stress-rate can also be approximated using terms from the time "
                                   "at the previous time step ($t$) and current time step ($t + \\Delta t^{e}$): "
                                   "$\\smash[t]{\\overset{\\triangledown}{\\tau}}^{t + \\Delta t^{e}} \\approx "
                                   "\\frac{\\tau^{t + \\Delta t^{e} - \\tau^{t}}}{\\Delta t^{e}} - "
                                   "W^{t}\\tau^{t} + \\tau^{t}W^{t}$. "
                                   "In this material model, the size of the time step above ($\\Delta t^{e}$) "
                                   "can be specified as the numerical time step size or an independent fixed time "
                                   "step. If the latter case is a selected, the user has an option to apply a "
                                   "stress averaging scheme to account for the differences between the numerical "
                                   "and fixed elastic time step (eqn. 32). If one selects to use a fixed elastic time "
                                   "step throughout the model run, this can still be achieved by using CFL and "
                                   "maximum time step values that restrict the numerical time step to a specific time. "
                                   "\n\n "
                                   "The formulation above allows rewriting the total rate of deformation (eqn. 29) as "
                                   "$\\tau^{t + \\Delta t^{e}} = \\eta_{eff} \\left ( "
                                   "2\\hat{D}^{t + \\triangle t^{e}} + \\frac{\\tau^{t}}{\\mu \\Delta t^{e}} + "
                                   "\\frac{W^{t}\\tau^{t} - \\tau^{t}W^{t}}{\\mu}  \\right )$. "
                                   "\n\n "
                                   "The effective viscosity (eqn. 28) is a function of the viscosity ($\\eta$), "
                                   "elastic time step size ($\\Delta t^{e}$) and shear relaxation time "
                                   "($ \\alpha = \\frac{\\eta}{\\mu} $): "
                                   "$\\eta_{eff} = \\eta \\frac{\\Delta t^{e}}{\\Delta t^{e} + \\alpha}$ "
                                   "The magnitude of the shear modulus thus controls how much the effective "
                                   "viscosity is reduced relative to the initial viscosity. "
                                   "\n\n "
                                   "Elastic effects are introduced into the governing Stokes equations through "
                                   "an elastic force term (eqn. 30) using stresses from the previous time step: "
                                   "$F^{e,t} = -\\frac{\\eta_{eff}}{\\mu \\Delta t^{e}} \\tau^{t}$. "
                                   "This force term is added onto the right-hand side force vector in the "
                                   "system of equations. "
                                   "\n\n "
                                   "The value of each compositional field representing distinct rock types at a "
                                   "point is interpreted to be a volume fraction of that rock type. If the sum of "
                                   "the compositional field volume fractions is less than one, then the remainder "
                                   "of the volume is assumed to be 'background material'."
                                   "\n\n "
                                   "Several model parameters (densities, elastic shear moduli, thermal expansivities, "
                                   "plasticity parameters, viscosity terms, etc) can be defined per-compositional field. "
                                   "For each material parameter the user supplies a comma delimited list of length "
                                   "N+1, where N is the number of compositional fields. The additional field corresponds "
                                   "to the value for background material. They should be ordered ''background, "
                                   "composition1, composition2...''. However, the first 3 (2D) or 6 (3D) composition "
                                   "fields correspond to components of the elastic stress tensor and their material "
                                   "values will not contribute to the volume fractions. If a single value is given, then "
                                   "all the compositional fields are given that value. Other lengths of lists are not "
                                   "allowed. For a given compositional field the material parameters are treated as "
                                   "constant, except density, which varies linearly with temperature according to the "
                                   "thermal expansivity. "
                                   "\n\n "
                                   "When more than one compositional field is present at a point, they are averaged "
                                   "arithmetically. An exception is viscosity, which may be averaged arithmetically, "
                                   "harmonically, geometrically, or by selecting the viscosity of the composition field "
                                   "with the greatest volume fraction.")
  }
}
