/*
  Copyright (C) 2018 - 2022 by the authors of the ASPECT code.

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
  namespace TosiBenchmark
  {
    using namespace dealii;


    /**
     * @ingroup MaterialModels
     */
    template <int dim>
    class TosiMaterial : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          /**
           * As described in Tosi et al (2015), the viscosity \eta is computed as the
           * harmonic average of a linear and nonlinear part.
           *
           * The linear part is calculated as follows
           * (see below for the meaning of the used parameter names):
           * $\eta_{lin}(T,z) = \exp(-\ln(\text{eta\_T} * T + \ln(\text{eta\_Z}) * z)$
           * while the strain rate dependent nonlinear part is computed as:
           * $\eta_{plast}(\dot\epsilon) = \text{eta\_asterisk} + \frac{\text{sigma\_yield}}{\sqrt(\dot\epsilon:\dot\epsilon)}$
           */

          //set up additional output for the derivatives
          MaterialModel::MaterialModelDerivatives<dim> *derivatives;
          derivatives = out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim>>();

          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              if (in.requests_property(MaterialModel::MaterialProperties::viscosity))
                {
                  out.viscosities[i] = viscosity (in.temperature[i],
                                                  in.pressure[i],
                                                  in.composition[i],
                                                  in.strain_rate[i],
                                                  in.position[i]);
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

              // If requested compute viscosity derivatives. This is only important
              // if using the Newton solver.
              if (derivatives != NULL && in.requests_property(MaterialModel::MaterialProperties::viscosity))
                {
                  if (use_analytical_derivative)
                    {
                      if (in.strain_rate[i].norm() == 0)
                        {
                          derivatives->viscosity_derivative_wrt_strain_rate[i] = 0;
                          derivatives->viscosity_derivative_wrt_pressure[i] = 0;
                        }
                      else
                        {
                          const double deviator_strain_rate_norm = in.strain_rate[i].norm();
                          const double part1 = (eta_asterisk * deviator_strain_rate_norm + sigma_yield);

                          derivatives->viscosity_derivative_wrt_strain_rate[i] = -(0.5*out.viscosities[i]*out.viscosities[i])
                                                                                 * (sigma_yield/(part1 * part1 * deviator_strain_rate_norm)) * in.strain_rate[i];
                          derivatives->viscosity_derivative_wrt_pressure[i] = 0;
                        }
                    }
                  else
                    {
                      // finite difference derivative
                      const double finite_difference_accuracy = 1e-7;

                      SymmetricTensor<2,dim> &deta = derivatives->viscosity_derivative_wrt_strain_rate[i];

                      // derivative in xx direction
                      SymmetricTensor<2,dim> dstrain_rate = in.strain_rate[i];
                      dstrain_rate[0][0] += std::fabs(in.strain_rate[i][0][0]) * finite_difference_accuracy;
                      const double eta_zero_zero = viscosity(in.temperature[i], in.pressure[i], in.composition[i], dstrain_rate, in.position[i]);
                      deta[0][0] = eta_zero_zero - out.viscosities[i];

                      if (dstrain_rate[0][0] != 0)
                        deta[0][0] /= std::fabs(in.strain_rate[i][0][0]) * finite_difference_accuracy;
                      else
                        deta[0][0] = 0;

                      // derivative in xy direction
                      dstrain_rate = in.strain_rate[i];
                      // dstrain_rate in yx direction is multiplied by 0.5, because the symmetric tensor
                      // is modified by 0.5 in xy and yx direction simultaneously and we compute the combined
                      // derivative
                      dstrain_rate[1][0] += 0.5 * std::fabs(in.strain_rate[i][1][0]) * finite_difference_accuracy;
                      const double eta_one_zero = viscosity(in.temperature[i], in.pressure[i], in.composition[i], dstrain_rate, in.position[i]);
                      deta[1][0] = eta_one_zero - out.viscosities[i];

                      if (dstrain_rate[1][0] != 0)
                        deta[1][0] /= std::fabs(in.strain_rate[i][1][0]) * finite_difference_accuracy;
                      else
                        deta[1][0] = 0;

                      // derivative in yy direction
                      dstrain_rate = in.strain_rate[i];
                      dstrain_rate[1][1] += std::fabs(in.strain_rate[i][1][1]) * finite_difference_accuracy;
                      const double eta_one_one = viscosity(in.temperature[i], in.pressure[i], in.composition[i], dstrain_rate, in.position[i]);
                      deta[1][1] = eta_one_one - out.viscosities[i];

                      if (dstrain_rate[1][1] != 0)
                        deta[1][1] /= std::fabs(in.strain_rate[i][1][1]) * finite_difference_accuracy;
                      else
                        deta[1][1] = 0;

                      derivatives->viscosity_derivative_wrt_pressure[i] = 0;
                    }
                }

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
        virtual bool is_compressible () const;
        /**
         * @}
         */
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

      private:

        double viscosity (const double                  temperature,
                          const double                  pressure,
                          const std::vector<double>    &compositional_fields,
                          const SymmetricTensor<2,dim> &strain_rate,
                          const Point<dim>             &position) const;

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
         * according to equation (8) of the paper.
         */
        double viscoplast(const double eta_asterisk,
                          const double stress_y,
                          const double strain_rate_norm) const;

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
         * The constant yield stress that is used in the
         * plastic viscosity formula
         */
        double sigma_yield;

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
         * Whether to use analytical or finite-difference viscosity derivatives
         * if the Newton solver is used.
         */
        bool use_analytical_derivative;

    };

    /*
     * Function to calculate the viscosity according
     * to equation (6) of the paper.
     */
    template <int dim>
    double
    TosiMaterial<dim>::
    viscosity (const double temperature,
               const double,
               const std::vector<double> &,
               const SymmetricTensor<2,dim> &strain_rate,
               const Point<dim> &) const
    {

      // In the first nonlinear iteration of the (pre-refinement steps of the) first time step,
      // strain rate is zero, so we set viscosity to eta_initial, a user-defined guess of the viscosity.
      if (strain_rate.norm() == 0)
        {
          return eta_initial;
        }

      // Otherwise we compute the linear viscosity and, if needed, the plastic viscosity.
      double viscosity = 0.0;
      const double visc_linear = viscolin(eta_T,eta_Z,temperature,0);

      if (eta_asterisk == 0.0 && sigma_yield == 0.0)
        {
          viscosity = visc_linear;
        }
      else
        {
          const double visc_plastic = viscoplast(eta_asterisk,sigma_yield,strain_rate.norm());

          // Compute the harmonic average (equation (6) of the paper)
          viscosity = 2.0 / ((1.0 / visc_linear) + (1.0 / visc_plastic));
        }

      // Cut-off the viscosity by user-defined values to avoid possible very large viscosity ratios
      viscosity = std::max(std::min(viscosity,eta_maximum),eta_minimum);

      return viscosity;
    }

    /*
     * Function to compute the linear viscosity
     * according to equation (7) of Tosi et al. 2015.
     */
    template <int dim>
    double
    TosiMaterial<dim>::
    viscolin(const double etaT,
             const double etaZ,
             const double T,
             const double z) const
    {

      return std::exp((-1.0 * std::log(etaT) * T ) + (std::log(etaZ) * z));
    }

    /*
     * Function to compute the plastic viscosity
     * according to equation (8) of the paper.
     */
    template <int dim>
    double
    TosiMaterial<dim>::
    viscoplast(const double etaasterisk,
               const double stressy,
               const double strainratenorm) const
    {
      return etaasterisk + (stressy/strainratenorm);
    }



    template <int dim>
    bool
    TosiMaterial<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    TosiMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      // Default values are for Case 1 of Tosi et al. (2015).
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Tosi benchmark");
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
          prm.declare_entry ("Yield stress", "0",
                             Patterns::Double (0),
                             "The value of the plastic viscosity yield stress $\\sigma_yield$, "
                             "as used in equation (8) of the paper.");
          prm.declare_entry ("Nonlinear viscosity constant", "0",
                             Patterns::Double (0),
                             "The value of the plastic viscosity constant $\\eta_asterisk$, "
                             "as used in equation (8) of the paper.");
          prm.declare_entry ("Use analytical derivative", "false",
                             Patterns::Bool (),
                             "Wether to use the analytical or the finite difference derivative for the Newton method.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    TosiMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Tosi benchmark");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          thermal_k                  = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          eta_T                      = prm.get_double ("Thermal viscosity parameter");
          eta_Z                      = prm.get_double ("Pressure viscosity parameter");
          sigma_yield                = prm.get_double ("Yield stress");
          eta_asterisk               = prm.get_double ("Nonlinear viscosity constant");
          eta_minimum                = prm.get_double ("Minimum viscosity");
          eta_maximum                = prm.get_double ("Maximum viscosity");
          eta_initial                = prm.get_double ("Initial viscosity");
          use_analytical_derivative  = prm.get_bool ("Use analytical derivative");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = (sigma_yield != 0.0) ? MaterialModel::NonlinearDependence::strain_rate | MaterialModel::NonlinearDependence::temperature | MaterialModel::NonlinearDependence::pressure : MaterialModel::NonlinearDependence::temperature | MaterialModel::NonlinearDependence::pressure;
      this->model_dependence.density = MaterialModel::NonlinearDependence::temperature;
      this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
      this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
    }


    /**
     * A postprocessor that computes some statistics about the
     * rate of viscous dissipation, the rate of work done against gravity
     * and the error between the two. These diagnostic quantities
     * are reported in the Tosi et al. (2015) paper for a visco-plastic
     * 2D thermal convection model.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class TosiPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for statistics on the rate of viscous dissipation,
         * rate of work and the error between the two.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };

    template <int dim>
    std::pair<std::string,std::string>
    TosiPostprocessor<dim>::execute (TableHandler &statistics)
    {
      AssertThrow(this->get_geometry_model().natural_coordinate_system() == aspect::Utilities::Coordinates::CoordinateSystem::cartesian,
                  ExcMessage("The current calculation of rate of work only makes sense in a Cartesian geometry."));

      AssertThrow(Plugins::plugin_type_matches<const TosiMaterial<dim>>(this->get_material_model()),
                  ExcMessage("The current calculation of viscous dissipation is only for incompressible models "
                             "and specifically computes the difference between work and dissipation as requested "
                             "in the paper of Tosi et al. 2015."));

      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(this->introspection().base_elements.velocities)
                                            .degree+1);

      const unsigned int n_q_points = quadrature_formula.size();
      std::vector<Tensor<1,dim>> velocities(n_q_points);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_gradients |
                               update_JxW_values);

      // the local integral values of rate of viscous dissipation and work against gravity
      double local_dissipation_integral = 0.0;
      double local_work = 0.0;

      // the values of the compositional fields are stored as blockvectors for each field
      // we have to extract them in this structure
      std::vector<std::vector<double>> prelim_composition_values (this->n_compositional_fields(),
                                                                   std::vector<double> (n_q_points));

      typename MaterialModel::Interface<dim>::MaterialModelInputs in(n_q_points,
                                                                     this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(n_q_points,
                                                                       this->n_compositional_fields());
      in.requested_properties = MaterialModel::MaterialProperties::viscosity;

      // loop over active, locally owned cells and
      // extract material model input and compute integrals
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            // retrieve velocities
            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocities);

            // retrieve the input for the material model
            in.position = fe_values.get_quadrature_points();

            fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                      in.pressure);

            fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                         in.temperature);

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values
              (this->get_solution(),prelim_composition_values[c]);

            for (unsigned int i=0; i<n_q_points; ++i)
              {
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  in.composition[i][c] = prelim_composition_values[c][i];
              }

            fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(),
                in.strain_rate);

            // get the output from the material model
            this->get_material_model().evaluate(in, out);

            // calculate the local viscous dissipation integral and local rate of work against gravity
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const SymmetricTensor<2,dim> strain_rate_dev = deviator(in.strain_rate[q]);
                local_dissipation_integral += 2.0 * out.viscosities[q] * strain_rate_dev * strain_rate_dev * fe_values.JxW(q);
                local_work += in.temperature[q] * velocities[q][dim-1] * fe_values.JxW(q);
              }
          }

      // compute the integrals on the whole domain
      const double viscous_dissipation
        = Utilities::MPI::sum (local_dissipation_integral, this->get_mpi_communicator());

      const double work
        = Utilities::MPI::sum(local_work, this->get_mpi_communicator());

      // compute the percentage error between rate of dissipation (divided by the surface rayleigh number of the Tosi et al. (2015)
      // benchmark (Ra=100)) and rate of work (see equation (21) of the paper)
      const double error = 100.0 * (std::abs(work - (viscous_dissipation / 100.0)) / std::max(work, viscous_dissipation / 100.0));

      if (this->convert_output_to_years() == true)
        {
          // fill statistics file
          // make sure that the columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.add_value ("Viscous dissipation (J/yr)", viscous_dissipation * year_in_seconds);
          statistics.set_precision ("Viscous dissipation (J/yr)", 8);
          statistics.set_scientific ("Viscous dissipation (J/yr)", true);
          statistics.add_value ("Rate of work (Km/yr)", work * year_in_seconds);
          statistics.set_precision ("Rate of work (Km/yr)", 8);
          statistics.set_scientific ("Rate of work (Km/yr)", true);
          statistics.add_value ("Error", error);
          statistics.set_precision ("Error", 8);
          statistics.set_scientific ("Error", true);
        }
      else
        {
          // fill statistics file
          // make sure that the columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.add_value ("Viscous dissipation (W)", viscous_dissipation);
          statistics.set_precision ("Viscous dissipation (W)", 8);
          statistics.set_scientific ("Viscous dissipation (W)", true);
          statistics.add_value ("Rate of work (W)", work);
          statistics.set_precision ("Rate of work (W)", 8);
          statistics.set_scientific ("Rate of work (W)", true);
          statistics.add_value ("Error", error);
          statistics.set_precision ("Error", 8);
          statistics.set_scientific ("Error", true);
        }

      std::ostringstream output;
      output.precision(3);
      if (this->convert_output_to_years() == true)
        output << viscous_dissipation *year_in_seconds
               << " J/yr "
               << work *year_in_seconds
               << " J/yr "
               << error
               << " error";
      else
        output << viscous_dissipation
               << " W "
               << work
               << " W "
               << error
               << " error";

      return std::pair<std::string, std::string> ("Viscous dissipation, rate of work, error:",
                                                  output.str());
    }


  }
}



// explicit instantiations
namespace aspect
{
  namespace TosiBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(TosiMaterial,
                                   "TosiMaterial",
                                   "A material model that has constant values "
                                   "for all coefficients except the density and viscosity as described in the "
                                   "open access paper of Tosi et al. 2015. "
                                   "The default parameter values are chosen according to Case 1 of the paper. "
                                   "All of the values that define this model are read "
                                   "from a section ``Material model/Tosi model'' in the input file, see "
                                   "Section~\\ref{parameters:Material_model/Tosi_model}."
                                   "\n\n"
                                   "This model uses the following set of equations for the two coefficients that "
                                   "are non-constant (see equation (6) - (10) of the paper): "
                                   "\\begin{align}"
                                   "  \\eta(T,z,\\dot \\epsilon) &= 2\\frac{1}{\\frac{1}{\\eta_{lin}(T,z)}\\frac{1}{\\eta_{plast}(\\dot\\epsilon}}, \\\\"
                                   "  \\rho(T) &= \\left(1-\\alpha (T-T_0)\\right)\\rho_0,"
                                   "\\end{align}"
                                   "where $z$ represents depth."
                                   "\n\n"
                                   "The linear and plastic viscosity parts are defined as follows:"
                                   "\\begin{align}"
                                   "  \\eta_{lin}(T,z) &= \\exp(-\\ln(\\eta_T)T+\\ln(\\eta_z)z), \\\\"
                                   "  \\eta_{plast}(\\dot\\epsilon) &= \\eta^{*}+\\frac{\\sigma_y}{\\sqrt(\\dot\\epsilon : \\dot\\epsilon)} "
                                   "\\end{align} "
                                   "\n\n"
                                   "Note that this model uses the formulation that assumes an incompressible "
                                   "medium despite the fact that the density follows the law "
                                   "$\\rho(T)=\\rho_0(1-\\alpha(T-T_0))$. ")
    ASPECT_REGISTER_POSTPROCESSOR(TosiPostprocessor,
                                  "TosiPostprocessor",
                                  "A postprocessor that computes the viscous dissipation"
                                  "for the whole domain as: "
                                  "$\\left<\\Phi\\right>=\\int_{V} \\tau : \\dot{\\epsilon}dV$ "
                                  "= $\\int_{V} 2\\mu\\dot{\\epsilon}:\\dot{\\epsilon} dV$. "
                                  "Besides the dissipation, for the Tosi et al. (2015) benchmark "
                                  "the rate of work against gravity is "
                                  "calculated, as well as the percentage error between the two: "
                                  "$W$ = $\\int_{V} T u_{y} dV$ and "
                                  "$\\delta$ = $ \\frac{|W - \\frac{\\left<\\Phi\\right>}{Ra}|}{\\max\\left(\\left<W\\right>,\\frac{\\left<\\Phi\\right>}{Ra}\\right)}$. "
                                  "This error W should tend to zero if in steady state the thermal energy is accurately preserved."
                                  "Note that this postprocessor only makes sense for box geometries.")
  }
}
