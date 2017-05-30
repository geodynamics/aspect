

#include <aspect/postprocess/interface.h>
#include <aspect/material_model/simple.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>



namespace aspect
{
  using namespace dealii;



  /**
   * This benchmark is from the article
   * @code
   * @endcode
   *
   * @ingroup Postprocessing
   */


  namespace MaterialModel
  {

    template <int dim>
    class Material : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:



        /**
        * Evaluate material properties.
        */
        virtual void evaluate(const MaterialModelInputs<dim> &in,
                              MaterialModelOutputs<dim> &out) const
        {

          for (unsigned int i=0; i < in.temperature.size(); ++i)
            {
              const Point<dim> position = in.position[i];
              const double temperature = in.temperature[i];
              const double pressure = in.pressure[i];

              const double depth = 1.0-position(dim-1);

              out.viscosities[i] = ((Di==0.0)?1.0:Di)/Ra*exp(-b*(temperature- this->get_adiabatic_conditions().temperature(position))+c*depth);

              out.specific_heat[i] = reference_specific_heat;
              out.thermal_conductivities[i] = 1.0;
              out.thermal_expansion_coefficients[i] = (Di==0.0)?1.0:Di;

              double rho = reference_rho * exp(depth * Di/gamma);
              rho *= 1.0 - out.thermal_expansion_coefficients[i] * (temperature - this->get_adiabatic_conditions().temperature(position))
                     + (tala?0.0:1.0)*Di*gamma
                     *  (pressure - this->get_adiabatic_conditions().pressure(position));

              out.densities[i] = rho;

              out.compressibilities[i] = 0.0;
              out.entropy_derivative_pressure[i] = 0.0;
              out.entropy_derivative_temperature[i] = 0.0;
              // Change in composition due to chemical reactions at the
              // given positions. The term reaction_terms[i][c] is the
              // change in compositional field c at point i.
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                out.reaction_terms[i][c] = 0.0;
            }
        }

        virtual bool is_compressible () const
        {
          return true;
        }



        virtual double reference_viscosity () const
        {
          return (Di==0.0?1.0:Di)/Ra;
        }

        /**
         * @}
         */

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

        bool tala;

        /**
         * blankenbach parameters
         */
        double b, c;


        /**
         * The reference surface temperature
         */
        double reference_rho;

        /**
         * The constant viscosity
         */
        double Di, Ra, gamma;

        /**
         * The constant specific heat
         */
        double reference_specific_heat;

        /**
         * The constant compressibility.
         */
        double reference_compressibility;


    };



    template <int dim>
    void
    Material<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Tan Gurnis model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Ra", "1e4",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("Di", "0.0",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("gamma", "1.0",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("b", "6.907755279",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("c", "0",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("Use TALA", "false",
                             Patterns::Bool (),
                             "Whether to use the TALA instead of the ALA "
                             "approximation.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Material<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Tan Gurnis model");
        {
          reference_rho   = prm.get_double ("Reference density");
          //          reference_T = prm.get_double ("Reference temperature");
          // eta                   = prm.get_double ("Viscosity");
          b               = prm.get_double ("b");
          c               = prm.get_double ("c");
          Di              = prm.get_double ("Di");
          Ra              = prm.get_double ("Ra");
          gamma           = prm.get_double ("gamma");

          tala            = prm.get_bool ("Use TALA");

          reference_specific_heat = prm.get_double ("Reference specific heat");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature;
      this->model_dependence.density = NonlinearDependence::none;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
    }
  }

}



// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Material,
                                   "Material",
                                   "A simple compressible material model based on a benchmark"
                                   " from the paper of Tan/Gurnis (2007). This does not use the"
                                   " temperature equation, but has a hardcoded temperature.")
  }

}

