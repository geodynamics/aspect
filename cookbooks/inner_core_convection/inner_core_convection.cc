#include <aspect/material_model/simple.h>
#include <aspect/heating_model/interface.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class InnerCore : public MaterialModel::Simple<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the boundary values will
         * next be evaluated. For the current class, the function passes to
         * the parsed function what the current time is.
         */
        virtual
        void
        update ();

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * A function object representing resistance to phase change at the
         * inner core boundary as a function the position (and, optionally,
         * the model time).
         */
        Functions::ParsedFunction<dim> resistance_to_phase_change;

      private:
        /**
         * Parameters related to the phase transition.
         */
        double transition_depth;
        double transition_width;
        double transition_temperature;
        double transition_clapeyron_slope;
        double transition_density_change;

        /**
         * Percentage of material that has already undergone the phase
         * transition to the higher-pressure material.
         */
        virtual
        double
        phase_function (const Point<dim> &position,
                        const double temperature) const;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    double
    InnerCore<dim>::
    phase_function (const Point<dim> &position,
                    const double temperature) const
    {
      // We need to convert the depth to pressure,
      // keeping in mind that for this model, the pressure is only the dynamic pressure,
      // so we have to compute the hydrostatic pressure explicitly as rho * g.
      const double depth = this->get_geometry_model().depth(position);
      const double reference_rho = 1.0; // formulation is nondimensionalized
      const double pressure_gradient = this->get_gravity_model().gravity_vector(position).norm() * reference_rho;
      double depth_deviation = depth - transition_depth;

      if (std::abs(pressure_gradient) > 0)
        depth_deviation -= transition_clapeyron_slope * (temperature - transition_temperature) / pressure_gradient;

      double phase_func;
      // use delta function for width = 0
      if (transition_width == 0)
        (depth_deviation > 0) ? phase_func = 1 : phase_func = 0;
      else
        phase_func = 0.5*(1.0 + std::tanh(depth_deviation / transition_width));
      return phase_func;
    }


    template <int dim>
    void
    InnerCore<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // First, we use the material descriptions of the 'simple' material model to fill all of the material
      // model outputs. Below, we will then overwrite selected properties (the specific heat) to make the
      // product of density and specific heat a constant.
      Simple<dim>::evaluate(in, out);

      // We want the right-hand side of the momentum equation to be (- Ra T gravity) and
      // density * cp to be 1
      for (unsigned int q=0; q < in.position.size(); ++q)
        {
          out.densities[q] = - out.thermal_expansion_coefficients[q] * in.temperature[q]
                             + phase_function (in.position[q], in.temperature[q]) * transition_density_change;
          if (std::abs(out.densities[q]) > 0.0)
            out.specific_heat[q] /= out.densities[q];
        }
    }


    template <int dim>
    void
    InnerCore<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        resistance_to_phase_change.set_time (this->get_time() / year_in_seconds);
      else
        resistance_to_phase_change.set_time (this->get_time());
    }


    template <int dim>
    void
    InnerCore<dim>::declare_parameters (ParameterHandler &prm)
    {
      Simple<dim>::declare_parameters (prm);

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Inner core");
        {
          prm.enter_subsection("Phase change resistance function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
          }
          prm.leave_subsection();

          prm.declare_entry ("Phase transition depth", "0.0",
                             Patterns::Double (0),
                             "The depth where the phase transition occurs. "
                             "Units: m.");
          prm.declare_entry ("Phase transition width", "0.0",
                             Patterns::Double (0),
                             "The width of the phase transition. The argument of the phase function "
                             "is scaled with this value, leading to a jump between phases "
                             "for a value of zero and a gradual transition for larger values. "
                             "Units: m.");
          prm.declare_entry ("Phase transition temperature", "0.0",
                             Patterns::Double (0),
                             "The temperature at which the phase transition occurs in the depth "
                             "given by the 'Phase transition depth' parameter.  Higher or lower "
                             "temperatures lead to phase transition occurring in shallower or greater "
                             "depths, depending on the Clapeyron slope given in 'Phase transition "
                             "Clapeyron slope'. "
                             "Units: K.");
          prm.declare_entry ("Phase transition Clapeyron slope", "0.0",
                             Patterns::Double (),
                             "The Clapeyron slope of the phase transition. A positive "
                             "Clapeyron slope indicates that the phase transition will occur in "
                             "a greater depth if the temperature is higher than the one given in "
                             "Phase transition temperatures (and in a smaller depth if the "
                             "temperature is smaller than the one given in Phase transition temperatures). "
                             "For negative Clapeyron slopes, the effect is in the opposite direction. "
                             "Units: Pa/K.");
          prm.declare_entry ("Phase transition density change", "0.0",
                             Patterns::Double (),
                             "The density change that occurs across the phase transition. "
                             "A positive value means that the density increases with depth. "
                             "Units: kg/m$^3$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    InnerCore<dim>::parse_parameters (ParameterHandler &prm)
    {
      Simple<dim>::parse_parameters (prm);

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Inner core");
        {
          prm.enter_subsection("Phase change resistance function");
          try
            {
              resistance_to_phase_change.parse_parameters (prm);
            }
          catch (...)
            {
              std::cerr << "ERROR: FunctionParser failed to parse\n"
                        << "\t'Phase boundary model.Function'\n"
                        << "with expression\n"
                        << "\t'" << prm.get("Function expression") << "'"
                        << "More information about the cause of the parse error \n"
                        << "is shown below.\n";
              throw;
            }
          prm.leave_subsection();
        }

        transition_depth           = prm.get_double ("Phase transition depth");
        transition_width           = prm.get_double ("Phase transition width");
        transition_temperature     = prm.get_double ("Phase transition temperature");
        transition_clapeyron_slope = prm.get_double ("Phase transition Clapeyron slope");
        transition_density_change  = prm.get_double ("Phase transition density change");

        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}



namespace aspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements a constant radiogenic heating rate.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class ConstantCoreHeating : public Interface<dim>
    {
      public:
        /**
         * Return the heating terms. For the current class, this
         * function obviously simply returns a constant value.
         */
        virtual
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
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
        /**
         * @}
         */

      private:
        double radiogenic_heating_rate;
    };
  }
}

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    ConstantCoreHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &/*material_model_inputs*/,
              const MaterialModel::MaterialModelOutputs<dim> &/*material_model_outputs*/,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          // return a constant value
          heating_model_outputs.heating_source_terms[q] = radiogenic_heating_rate;
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    void
    ConstantCoreHeating<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Constant core heating");
        {
          prm.declare_entry ("Radiogenic heating rate", "0e0",
                             Patterns::Double (0),
                             "The specific rate of heating due to radioactive decay (or other bulk sources "
                             "you may want to describe). This parameter corresponds to the variable "
                             "$H$ in the temperature equation stated in the manual, and the heating "
                             "term is $\rho H$. Units: W/kg.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ConstantCoreHeating<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Constant core heating");
        {
          radiogenic_heating_rate    = prm.get_double ("Radiogenic heating rate");
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
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(InnerCore,
                                   "inner core material",
                                   "A simple material model that is like the "
                                   "'Simple' model, but has a constant $\rho c_p$, "
                                   "and implements a function that characterizes the "
                                   "resistance to melting/freezing at the inner core "
                                   "boundary.")
  }

  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ConstantCoreHeating,
                                  "constant core heating",
                                  "Implementation of a model in which the heating "
                                  "rate is constant.")
  }
}
