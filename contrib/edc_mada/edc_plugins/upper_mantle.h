/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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


#ifndef __aspect__upper_mantle_h
#define __aspect__upper_mantle_h

#include <aspect/geometry_model/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/drucker_prager.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace internal
    {
      namespace Boundary
      {
        class BoundaryLookup;
      }
    }

    template <int dim>
    class UpperMantle : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        virtual double reference_thermal_expansion_coefficient () const;

        double reference_thermal_diffusivity () const;

        double reference_cp () const;
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

      private:

        double reference_rho;
        double reference_T;
        double eta;
        double composition_viscosity_prefactor;
        double thermal_viscosity_exponent;
        double thermal_alpha;
        double reference_specific_heat;

        /*
         * Rheology parameters
         */
        double grain_size;
        double activation_energie_diffusion;
        double activation_volume_diffusion;
        double stress_exponent_diffusion;
        double grain_size_exponent_diffusion;
        double activation_energie_dislocation;
        double activation_volume_dislocation;
        double stress_exponent_dislocation;
        double prefactor_diffusion;
        double prefactor_dislocation;
        double min_strain_rate;
        double ref_strain_rate;
        double min_visc;
        double max_visc;
        double C_OH;
        double reference_compressibility;

        enum averaging_scheme
        {
          harmonic,
          arithmetic,
          geometric,
          maximum_composition
        } viscosity_averaging, density_averaging;

        /**
         * The thermal conductivity.
         */
        double k_value;

        double compositional_delta_rho;

        std::vector<double> compute_volume_fractions( const std::vector<double> &compositional_fields) const;

        double average_value (const std::vector<double> &composition,
                              const std::vector<double> &parameter_values,
                              const averaging_scheme &average_type) const;

        double diffusion_creep (const double &pressure,
                                const double &temperature) const;

        double dislocation_creep (const double &pressure,
                                  const double &temperature,
                                  const SymmetricTensor<2,dim> &strain_rate) const;
        /**
         * Pointer to the material model used as the base model
         */
        std::shared_ptr<MaterialModel::Interface<dim> > base_model;

    };

  }
}

#endif
