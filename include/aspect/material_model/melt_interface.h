/*
  Copyright (C) 2011, 2012, 2013 by the authors of the ASPECT code.

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


#ifndef __aspect__material_model_melt_interface_h
#define __aspect__material_model_melt_interface_h

#include <aspect/plugins.h>
#include <aspect/material_model/interface.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with modeling
   * convecting material, including descriptions of material parameters such
   * as viscosities, densities, etc.
   *
   * @ingroup MaterialModels
   */
  namespace MaterialModel
  {
    using namespace dealii;
    /**
     * A base class for parameterizations of material models. Classes derived
     * from this class will need to implement functions that provide material
     * parameters such as the viscosity, density, etc, typically as a function
     * of position, temperature and pressure at that location.
     *
     * There is two ways to implement a material model and they can not be
     * mixed: Option one is to override all the virtual functions like
     * viscosity(), density(), etc. but not change evaluate().
     *
     * Option two only requires you to override evaluate() and fill the output
     * argument struct instead of implementing the functions viscosity(),
     * density(), etc.. In this case, all other functions are being ignored.
     *
     * The second option is more efficient in general, but it is okay to use
     * option one for simple material models.
     *
     * In all cases, *_depends_on(), is_compressible(), reference_viscosity(),
     * reference_density() need to be implemented.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MeltInterface: public Interface<dim>
    {
      public:
        /**
         * Destructor. Made virtual to enforce that derived classes also have
         * virtual destructors.
         */
        virtual ~MeltInterface();

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in. If MaterialModelInputs.strain_rate has the length
         * 0, then the viscosity does not need to be computed.
         */
        virtual void evaluate(const MaterialModelInputs &in, MaterialModelOutputs &out) const = 0;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
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
         * @}
         */
    };
  }
}


#endif
