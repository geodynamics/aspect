/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_cpp_h
#define _aspect_material_model_cpp_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <regex>

// This plugin requires that we can load shared libraries.
#if ASPECT_USE_SHARED_LIBS==1
#include <dlfcn.h>
#endif

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A structure for organizing user code snippets.
     */
    struct UserCode
    {
      std::vector<std::string> includes;
      std::string variable_defs;
      std::string update_function;
      std::string viscosity_function;
      std::string density_function;
      std::string thermal_conductivity_function;
      std::string thermal_expansivity_function;
      std::string specific_heat_function;
      std::string compressibility_function;
      std::string entropy_derivative_p_function;
      std::string entropy_derivative_t_function;
      std::string reaction_function;
    };

    /**
     * Converts a string describing a nonlinear dependency
     * to its equivalent `Dependence'.
     */
    inline
    NonlinearDependence::Dependence
    str2dep (const std::string depstr)
    {
      if (depstr == "none")
        {
          return NonlinearDependence::none;
        }
      else if (depstr == "compositional fields")
        {
          return NonlinearDependence::compositional_fields;
        }
      else if (depstr == "pressure")
        {
          return NonlinearDependence::pressure;
        }
      else if (depstr == "strain rate")
        {
          return NonlinearDependence::strain_rate;
        }
      else if (depstr == "temperature")
        {
          return NonlinearDependence::temperature;
        }
      else
        {
          AssertThrow(false, ExcMessage("Nonlinear dependencies must be one of "
                                        "non|compositional fields|pressure|strain rate|temperature."));
        }
      return NonlinearDependence::uninitialized;
    }

    /**
     * Execute a command using the system() command processor. Return the
     * exit code.
     */
    int
    execute (const std::string &cmd)
    {
      AssertThrow(system((char *)0) != 0,
                  ExcMessage("The 'CPP' material model requires a command-processor, "
                             "which appears to be unavailable on this system."));

      return system(cmd.c_str());
    }

    /**
     * A material model that compiles user-defined C++ code at runtime.
     * It exposes plugin functionality at the input parameter level,
     * allowing quick deployment of new material models without requiring
     * knowledge of Aspect's plugin architecture.
     *
     * It aims to make models more easily auditable, as more of the relevant
     * model code can be contained in a single parameter file.
     *
     * The generated code can be assumed to have the form:
     * ```
        #include <USER>
        #include <INCLUDES>
        #include <aspect/material_model/interface.h>
        #include <aspect/simulator.h>
        #include <aspect/simulator_access.h>
        namespace aspect {
          namespace MaterialModel {
            using namespace dealii;
            class Local_
            {
              public:
                const unsigned int dim = dim;

                USER GLOBAL VARIABLES

                void update (position, temperature, ...)
                  { USER_UPDATE_CODE }
                double viscosity (position, temperature, pressure, ...)
                  { USER VISCOSITY CODE }
                double density (position, temperature, pressure, ...)
                  { USER DENSITY CODE }
                ...
            } local;

            eval (in, out, simulator)
            {
              for (unsigned int _i=0; _i<in.position.size(); ++_i)
                {
                  local.update(in.position[_i], in.temperature[_i], ...);
                  out.viscosities[_i] = local.viscosity (in.position[_i], ...);
                  out.densities[_i]   = local.density (in.position[_i], ...);
                  ...
                }
            }
          }
        }
     * ```
     * where USER_* are defined by input parameters.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class CPP : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        virtual bool is_compressible () const;

        virtual double reference_viscosity () const;

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;


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
         * Reads parameters from user input, then coordinates generating,
         * compiling, and linking functions to be used when evaluating
         * material properties.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

      private:
        double reference_viscosity_param;
        bool is_compressible_param;

        /**
         * Define function pointer for the new evaluate function.
         */
        typedef void (*eval_t)(const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out,
                               const ::aspect::SimulatorAccess<dim> *simulator);

        /**
         * Pointer to the user's evaluation function.
         */
        eval_t user_eval;
    };

  }
}

#endif
