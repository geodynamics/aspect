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

        struct MaterialModelInputs: public Interface<dim>::MaterialModelInputs
         {
            MaterialModelInputs (const unsigned int n_points,
                                              const unsigned int n_comp);
         };
        struct MaterialModelOutputs: public Interface<dim>::MaterialModelOutputs
         {
            MaterialModelOutputs (const unsigned int n_points,
                                  const unsigned int n_comp);

            std::vector<double> compaction_viscosities;
            std::vector<double> fluid_viscosities;
            std::vector<double> permeabilities;
            std::vector<double> fluid_densities;
            std::vector<double> fluid_compressibilities;
         };


        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in. If MaterialModelInputs.strain_rate has the length
         * 0, then the viscosity does not need to be computed.
         */
        virtual void evaluate_with_melt(const MaterialModelInputs &in, MaterialModelOutputs &out) const = 0;

        virtual void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
            typename Interface<dim>::MaterialModelOutputs &out) const = 0;

    };
  }
}


#endif
