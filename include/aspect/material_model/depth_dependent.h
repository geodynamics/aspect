/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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



#ifndef __aspect__model_depth_dependent_h
#define __aspect__model_depth_dependent_h

#include <deal.II/base/function_lib.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that applies a depth-dependent viscosity to a ''base model''
     * chosen from any of the other available material models. This depth-dependent
     * material model allows the user to specify a depth-dependent reference viscosity
     * either through a parsed function, lists of depth and viscosity values, or a file.
     * The current implementation only allows for depth-dependence of viscosity - all
     * other properties are derived from the base model.
     * @ingroup MaterialModels
     */
    template <int dim>
    class DepthDependent : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in.
         */
        virtual
        void
        evaluate (const typename Interface<dim>::MaterialModelInputs &in,
                  typename Interface<dim>::MaterialModelOutputs &out) const;
        /**
         * Method to declare parameters related to depth-dependent model
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * Method to parse parameters related to depth-dependent model
         */
        virtual void
        parse_parameters (ParameterHandler &prm);

        /**
         * Method that indicates whether material is compressible. Depth dependent model is compressible
         * if and only if base model is compressible.
         */
        virtual bool is_compressible () const;

        /**
         * Method to calculate reference viscosity for the depth-dependent model. The reference
         * viscosity is determined by evaluating the depth-dependent part of the viscosity at
         * the mean depth of the model.
         */
        virtual double reference_viscosity () const;

        /**
         * Method to calculate reference density for the depth-dependent model. Because the depth-
         * dependent model deos not modify density, the reference density is equivalent to the
         * base model's reference density.
         */
        virtual double reference_density () const;

      private:

        /**
         * An enum to describe where the depth dependency of the viscosity is coming from.
         */
        enum ViscositySource
        {
          Function,
          File,
          List,
          None
        };

        /**
         * Currently chosen source for the viscosity.
         **/
        ViscositySource viscosity_source;

        /**
         * Function to read depth-dependent lookup table and set up interpolating function
         * for File depth dependence method
         */
        void
        read_viscosity_file(const std::string &filename,
                            const MPI_Comm &comm);

        /**
         * Data structures to store depth and viscosity lookup tables as well as interpolating
         * function to calculate viscosity for File Depth dependence method
         */
        std_cxx11::shared_ptr< Functions::InterpolatedTensorProductGridData<1> > viscosity_file_function;

        /**
         * Function to calculate viscosity at depth using values provided as List input
         */
        double
        viscosity_from_list(const double &depth) const;

        /**
         * function to calculate viscosity at depth using values provided from File input
         */
        double
        viscosity_from_file(const double &depth) const;

        /**
         * Function to calculate depth-dependent multiplicative prefactor to be applied
         * to base model viscosity.
         */
        double
        calculate_depth_dependent_prefactor(const double &depth) const;

        /**
         * Values of depth specified by the user if using List depth dependence method
         */
        std::vector<double> depth_values;
        std::vector<double> viscosity_values;
        /**
         * Parsed function that specifies viscosity depth-dependence when using the Function
         * method.
         */
        Functions::ParsedFunction<1> viscosity_function;

        /**
         * Pointer to the material model used as the base model
         */
        std_cxx11::shared_ptr<MaterialModel::Interface<dim> > base_model;
    };
  }
}

#endif
