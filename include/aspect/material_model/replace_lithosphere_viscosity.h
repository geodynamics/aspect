/*
  Copyright (C) 2014 - 2018 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_replace_lithosphere_viscosity_h
#define _aspect_material_model_replace_lithosphere_viscosity_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

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
    class ReplaceLithosphereViscosity : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
    	ReplaceLithosphereViscosity ();

        /**
         * Initialize the base model at the beginning of the run.
         */
        virtual
        void initialize();

        /**
         * Directory in which the LAB depth file is present.
         */
        std::string data_directory;

        /**
         * File name of the LAB depth file.
         */
        std::string LAB_file_name;

        /**
         * Return LAB depth as a function of position (latitude and longitude). For the
         * current class, this function returns value from the text files.
         */
        double
        ascii_lab (const Point<2> &position) const;

        /**
         * This parameter gives the viscosity set within the lithosphere.
         */
        double lithosphere_viscosity;

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

      private:


        Utilities::AsciiDataLookup<2> lab_depths;

        /**
         * Pointer to the material model used as the base model
         */
        std::shared_ptr<MaterialModel::Interface<dim> > base_model;
    };
  }
}

#endif
