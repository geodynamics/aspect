/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


#include <aspect/material_model/viscoelastic.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <aspect/global.h>
#include <numeric>

//using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    const std::vector<double>
    Viscoelastic<dim>::
    compute_volume_fractions( const std::vector<double> &compositional_fields) const
    {
      std::vector<double> volume_fractions( compositional_fields.size()+1);

      //clip the compositional fields so they are between zero and one
      std::vector<double> x_comp = compositional_fields;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);

      // assign compositional fields associated with strain a value of 0
      for ( unsigned int i=0; i < SymmetricTensor<2,dim>::n_independent_components; ++i)
        x_comp[i] = 0.0;

      //sum the compositional fields for normalization purposes
      double sum_composition = 0.0;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        sum_composition += x_comp[i];

      if (sum_composition >= 1.0)
        {
          volume_fractions[0] = 0.0;  //background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1]/sum_composition;
        }
      else
        {
          volume_fractions[0] = 1.0 - sum_composition; //background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1];
        }
      return volume_fractions;
    }

    template <int dim>
    double
    Viscoelastic<dim>::
    average_value ( const std::vector<double> &volume_fractions,
                    const std::vector<double> &parameter_values,
                    const enum AveragingScheme &average_type) const
    {
      double averaged_parameter = 0.0;

      switch (average_type)
        {
          case arithmetic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*parameter_values[i];
            break;
          }
          case harmonic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]/(parameter_values[i]);
            averaged_parameter = 1.0/averaged_parameter;
            break;
          }
          case geometric:
          {
            for (unsigned int i=0; i < volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*std::log(parameter_values[i]);
            averaged_parameter = std::exp(averaged_parameter);
            break;
          }
          case maximum_composition:
          {
            const unsigned int i = (unsigned int)(std::max_element( volume_fractions.begin(),
                                                                    volume_fractions.end() )
                                                  - volume_fractions.begin());
            averaged_parameter = parameter_values[i];
            break;
          }
          default:
          {
            AssertThrow( false, ExcNotImplemented() );
            break;
          }
        }
      return averaged_parameter;
    }


    template <int dim>
    void
    Viscoelastic<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {

      // Define elastic time step
      const double dte = ( ( this->get_timestep_number() > 0 && use_fixed_elastic_time_step == false )
                           ?
                           this->get_timestep()
                           :
                           fixed_elastic_time_step * year_in_seconds );

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const double temperature = in.temperature[i];
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = compute_volume_fractions(composition);

          out.specific_heat[i] = average_value ( volume_fractions, specific_heats, arithmetic);


          // Arithmetic averaging of thermal conductivities
          // This may not be strictly the most reasonable thing, but for most Earth materials we hope
          // that they do not vary so much that it is a big problem.
          out.thermal_conductivities[i] = average_value ( volume_fractions, thermal_conductivities, arithmetic);

          double density = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              // not strictly correct if thermal expansivities are different, since we are interpreting
              // these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor= (1.0 - thermal_expansivities[j] * (temperature - reference_T));
              density += volume_fractions[j] * densities[j] * temperature_factor;
            }
          out.densities[i] = density;

          out.thermal_expansion_coefficients[i] = average_value ( volume_fractions, thermal_expansivities, arithmetic);

          // Compressibility at the given positions.
          // The compressibility is given as
          // $\frac 1\rho \frac{\partial\rho}{\partial p}$.
          // (here we use an incompressible medium)
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

          // Compute average values of viscosity and elastic shear modulus
          const double vis_avg = average_value ( volume_fractions, viscosities, viscosity_averaging);
          const double esm_avg = average_value ( volume_fractions, elastic_shear_moduli, viscosity_averaging);

          //Compute viscoelastic (e.g., effective) viscosity (equation 28 in Moresi et al., 2003, J. Comp. Phys.) 
          const double vev_avg = ( vis_avg * dte ) / ( dte + ( vis_avg / esm_avg ) );

          out.viscosities[i] = vev_avg;

        }

      // Viscoelasticity section
      if (in.cell && this->get_timestep_number() > 0)
        {

          // Get current and old (previous time step) velocity gradients
          const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+1);
          FEValues<dim> fe_values (this->get_mapping(),
                                   this->get_fe(),
                                   quadrature_formula,
                                   update_gradients);
          std::vector<Tensor<2,dim> > velocity_gradients (quadrature_formula.size(), Tensor<2,dim>());
          fe_values.reinit (*in.cell);
          fe_values[this->introspection().extractors.velocities].get_function_gradients (this->get_solution(),
                                                                                         velocity_gradients);

          std::vector<Tensor<2,dim> > old_velocity_gradients (quadrature_formula.size(), Tensor<2,dim>());
          fe_values[this->introspection().extractors.velocities].get_function_gradients (this->get_old_solution(),
                                                                                         old_velocity_gradients);


          MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
          *force = out.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >();

          for (unsigned int i=0; i < in.position.size(); ++i)
            {

              // Get old stresses from compositional fields
              SymmetricTensor<2,dim> stress_old;
              for (unsigned int j=0; j < SymmetricTensor<2,dim>::n_independent_components; ++j)
                {
                  stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)] = in.composition[i][j];
                }

              // Calculate the rotated stresses
              // Rotation (vorticity) tensor (equation 25 in Moresi et al., 2003, J. Comp. Phys.)
              const Tensor<2,dim> rotation = 0.5 * ( old_velocity_gradients[i] - transpose(old_velocity_gradients[i]) );

              // Recalculate average values of viscosity, elastic shear modulus and viscoelastic (effective) viscosity 
              const std::vector<double> composition = in.composition[i];
              const std::vector<double> volume_fractions = compute_volume_fractions(composition);
              const double vis_avg = average_value ( volume_fractions, viscosities, viscosity_averaging);
              const double esm_avg = average_value ( volume_fractions, elastic_shear_moduli, viscosity_averaging);
              const double vev_avg = ( vis_avg * dte ) / ( dte + ( vis_avg / esm_avg ) );

              // Calculate new stress terms (see equation 29 in Moresi et al., 2003, J. Comp. Phys.)
              SymmetricTensor<2,dim> stress_new = ( 2. * vev_avg * deviator(in.strain_rate[i]) ) +
                                                  ( ( vev_avg / ( esm_avg * dte ) ) * stress_old ) +
                                                  ( ( vev_avg / esm_avg ) *
                                                    ( symmetrize(rotation * Tensor<2,dim>(stress_old) ) - symmetrize(Tensor<2,dim>(stress_old) * rotation) ) );

              // Stress averaging scheme to account for difference betweed fixed elastic time step
              // and numerical time step (see equation 32 in Moresi et al., 2003, J. Comp. Phys.)
              const double dt = this->get_timestep();
              if (use_fixed_elastic_time_step == true && use_stress_averaging == true)
                {
                  stress_new = ( ( 1. - ( dt / dte ) ) * stress_old ) + ( ( dt / dte ) * stress_new ) ; 
                }

              // Fill reaction terms
              for (unsigned int j = 0; j < SymmetricTensor<2,dim>::n_independent_components ; ++j)
                out.reaction_terms[i][j] = -stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)]
                                           + stress_new[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)];

              // Fill elastic outputs (see equation 30 in Moresi et al., 2003, J. Comp. Phys.) 
              if (force)
                {
                  force->rhs_e[i] = -1. * ( ( vev_avg / ( esm_avg * dte  ) ) * stress_old );
                }

            }
        }

    }

    template <int dim>
    double
    Viscoelastic<dim>::
    reference_viscosity () const
    {
      return viscosities[0]; //background
    }

    template <int dim>
    bool
    Viscoelastic<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    Viscoelastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Viscoelastic");
        {
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Viscosities", "1.e21",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $Pa s$");
          prm.declare_entry ("Thermal expansivities", "4.e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $1/K$");
          prm.declare_entry ("Specific heats", "1250.",
                             Patterns::List(Patterns::Double(0)),
                             "List of specific heats $C_p$ for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $J /kg /K$");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $W/m/K$ ");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          prm.declare_entry ("Elastic shear moduli", "75.0e9",
                             Patterns::List(Patterns::Double(0)),
                             "List of elastic shear moduli, $G$, "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The default value of 75 GPa is representative of mantle rocks. Units: none.");
          prm.declare_entry ("Use fixed elastic time step", "false",
                             Patterns::Bool (),
                             "Select whether material time scale in the viscoelastic constitutive"
                             "relationship uses the regular numerical time step or a separate fixed"
                             "elastic time step throughout the model run. The fixed elastic time step"
                             "is always used during the initial time step.");
          prm.declare_entry ("Fixed elastic time step", "1.e3",
                             Patterns::Double (0),
                             "The fixed elastic time step $dte$. Units: $yr$.");
          prm.declare_entry ("Use stress averaging","false",
                             Patterns::Bool (),
                             "Apply a stress averaging scheme to account for differences between the"
                             "fixed elastic time step and numerical time step.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Viscoelastic<dim>::parse_parameters (ParameterHandler &prm)
    {
      //not pretty, but we need to get the number of compositional fields before
      //simulatoraccess has been initialized here...
      unsigned int n_foreground_fields;
      prm.enter_subsection ("Compositional fields");
      {
        n_foreground_fields = prm.get_integer ("Number of fields");
      }
      prm.leave_subsection();

      const unsigned int n_fields= n_foreground_fields + 1;


      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Viscoelastic");
        {
          reference_T = prm.get_double ("Reference temperature");

          if (prm.get ("Viscosity averaging scheme") == "harmonic")
            viscosity_averaging = harmonic;
          else if (prm.get ("Viscosity averaging scheme") == "arithmetic")
            viscosity_averaging = arithmetic;
          else if (prm.get ("Viscosity averaging scheme") == "geometric")
            viscosity_averaging = geometric;
          else if (prm.get ("Viscosity averaging scheme") == "maximum composition")
            viscosity_averaging = maximum_composition;
          else
            AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));

          // Parse viscoelastic properties
          densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                              n_fields,
                                                              "Densities");
          viscosities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosities"))),
                                                                n_fields,
                                                                "Viscosities");
          thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                                                                           n_fields,
                                                                           "Thermal conductivities");
          thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),
                                                                          n_fields,
                                                                          "Thermal expansivities");
          specific_heats = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Specific heats"))),
                                                                   n_fields,
                                                                   "Specific heats");
          elastic_shear_moduli = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Elastic shear moduli"))),
                                                                         n_fields,
                                                                         "Elastic shear moduli");
          
          use_fixed_elastic_time_step = prm.get_bool ("Use fixed elastic time step");

          use_stress_averaging = prm.get_bool ("Use stress averaging");
          if (use_stress_averaging)
            AssertThrow(use_fixed_elastic_time_step == true,
                        ExcMessage("A fixed elastic time step must also be used with stress averaging"));
              
          fixed_elastic_time_step = prm.get_double ("Fixed elastic time step");


        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::compositional_fields;
      this->model_dependence.thermal_conductivity = NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Viscoelastic,
                                   "viscoelastic",
                                   "An implementation of a simple linear viscoelastic rheology that "
                                   "only includes the deviatoric components of elasticity. Specifically, "
                                   "the viscoelastic rheology only takes into account the elastic shear "
                                   "strength (e.g., shear modulus), while the tensile and volumetric "
                                   "strength (e.g., Young's and bulk modulus) are not considered. The "
                                   "model is incompressible and allows specifying an arbitrary number "
                                   "of compositional fields, where each field represents a different "
                                   "rock type or component of the viscoelastic stress tensor. The stress "
                                   "tensor in 2D and 3D, respectively, contains 3 or 6 components. The "
                                   "compositional fields representing these components must be the first " 
                                   "listed compositional fields in the parameter file."
                                   "\n\n "
                                   "Expanding the model to include non-linear viscous flow (e.g., "
                                   "diffusion/dislocation creep) and plasticity would produce a "
                                   "of constitutive relationship commonly referred to as partial "
                                   "elastoviscoplastic (e.g., pEVP) in the geodynamics community. "
                                   "While extensively discussed and applied within the geodynamics "
                                   "literature, notable references include: "
                                   "Moresi et al. (2003), J. Comp. Phys., v. 184, p. 476-497. "
                                   "Gerya and Yuen (2007), Phys. Earth. Planet. Inter., v. 163, p. 83-105. "
                                   "Gerya (2010), Introduction to Numerical Geodynamic Modeling."
                                   "Kaus (2010), Tectonophysics, v. 484, p. 36-47. "
                                   "Choi et al. (2013), J. Geophys. Res., v. 118, p. 2429-2444. "
                                   "Keller et al. (2013), Geophys. J. Int., v. 195, p. 1406-1442. "
                                   "\n\n"
                                   "The overview below directly follows Moresi et al. (2003) eqns. 23-32. "
                                   "However, an important distinction between this material model and "
                                   "the studies above is the use of compositional fields, rather than "
                                   "tracers, to track individual components of the viscoelastic stress "
                                   "tensor. The material model will be udpated when an option to track "
                                   "and calculate viscoelastic stresses with tracers is implemented. "
                                   "\n\n"
                                   "Moresi et al. (2003) begins (eqn. 23) by writing the deviatoric "
                                   "rate of deformation ($\\hat{D}$) as the sum of elastic "
                                   "(($\\hat{D_{e}}$) and viscous (($\\hat{D_{v}}$)) components: "
                                   "$\\hat{D} = \\hat{D_{e}} + \\hat{D_{v}}$  "
                                   "These terms further decompose into"
                                   "$\\hat{D_{v}} = \\frac{\\tau}{2\\eta}$ and " 
                                   "$\\hat{D_{e}} = \\frac{\\overset{\\triangledown}{\\tau}}{2\\mu}$, where "
                                   "$\\tau$ is the viscous deviatoric stress, $\\eta$ is the shear viscosity, "
                                   "$\\mu$ is the shear modulus and $\\overset{\\triangledown}{\\tau}$ is the "
                                   "Jaumann corotational stress rate. This later term (eqn. 24) contains the "
                                   "time derivative of the deviatoric stress ($\\dot{\\tau}$) and terms that "
                                   "account for material spin (e.g., rotation) due to advection: "
                                   "$\\overset{\\triangledown}{\\tau} = \\dot{\\tau} + {\\tau}W -W\\tau$. "
                                   "Above, $W$ is the material spin tensor (eqn. 25): "
                                   "$W_{ij} = \\frac{1}{2} \\left (\\frac{\\partial V_{i}}{\\partial x_{j}} -$ " 
                                   "$\\frac{\\partial V_{j}}{\\partial x_{i}} \\right )$ "
                                   "\n\n"
                                   "The Jaumann stress-rate can also be approximated using terms from the time "
                                   "at the previous time step ($t$) and current time step ($t + \\triangle t_^{e}$): "
                                   "$\\smash[t]{\\overset{\\triangledown}{\\tau}}^{t + \\triangle t^{e}} \\approx $ "
                                   "$\\frac{\\tau^{t + \\triangle t^{e} - \\tau^{t}}}{\\triangle t^{e}} - $ "
                                   "$W^{t}\\tau^{t} + \\tau^{t}W^{t}$. "
                                   "In this material model, the size of the time step above ($\\triangle t^{e}$) "
                                   "can be specified as the numerical time step size or an independent fixed time "
                                   "step. If the latter case is a selected, the user has an option to apply a "
                                   "stress averaging scheme to account for the differences between the numerical "
                                   "and fixed elastic time step (eqn. 32). "
                                   "\n\n"
                                   "This formulation allows rewriting the total rate of deformation (eqn. 29) as "
                                   "$\\tau^{t + \\triangle t^{e}} = \\eta_{eff} \\left ( $ "
                                   "$2\\hat{D}^{t + \\triangle t^{e}} + \\frac{\\tau^{t}}{\\mu \\triangle t^{e}} +$ "
                                   "$\\frac{W^{t}\\tau^{t} - \\tau^{t}W^{t}}{\\mu}  \\right ) $ "
                                   "\n\n"
                                   "The effective viscosity (eqn. 28) is a function of the viscosity ($\\eta$), "
                                   "elastic time step size ($\\triangle t^{e}$) and shear relaxation time "
                                   "($ \\alpha = \\frac{\\eta}{\\mu} $): "
                                   "$\\eta_{eff} = \\eta \\frac{\\triangle t^{e}}{\\triangle t^{e} + \\alpha}$ "
                                   "The magnitude of the shear modulus thus controls how much the effective "
                                   "viscosity is reduced relative to the initial viscosity."
                                   "\n\n"
                                   "Elastic effects are introduced into the governing stokes equations through "
                                   "an elastic force term (eqn. 30) using stresses from the previous time step: "
                                   "$F^{e,t} = -\\frac{\\eta_{eff}}{\\mu \\triangle t^{e}} \\tau^{t}$. "
                                   "This force term is added onto the right-hand side force vector in the "
                                   "system of equations. ")
  }
}
