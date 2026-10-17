// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_initial_temperature_spherical_shell_h
#define _aspect_initial_temperature_spherical_shell_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace InitialTemperature
  {
    /**
     * A class that describes a perturbed initial temperature field for the
     * spherical shell.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class SphericalHexagonalPerturbation : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * The angular mode is the number of perturbations to apply to the
         * spherical shell. Historically, this was permanently set to 6 (hence
         * the class name SphericalHexagonalPerturbation) The default is 6 in
         * order to provide backwards compatibility.
         *
         * This variable is read from the parameter file through a parameter called 'Angular mode'.
         */
        unsigned int angular_mode;

        /**
         * The rotation offset describes the number of degrees to rotate the
         * perturbation counterclockwise. Setting the rotation offset to 0
         * will cause one of the perturbations to point north/up. Rotation
         * offset is set to -45 degrees by default in order to provide
         * backwards compatibility.
         * This variable is read from the parameter file through a parameter called 'Rotation offset'.
         */
        double rotation_offset;

        /**
         * Outer radius.
         */
        double R1;

    };


    /**
     * A class that describes a perturbed initial temperature field for the
     * spherical shell.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class SphericalGaussianPerturbation : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        /**
         * Constructor.
         */
        SphericalGaussianPerturbation();

        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * This variable is read from the parameter file through a parameter called 'Angle'.
         */
        double angle;
        /**
         * This variable is read from the parameter file through a parameter called 'Non-dimensional depth'.
         */
        double depth;
        /**
         * This variable is read from the parameter file through a parameter called 'Amplitude'.
         */
        double amplitude;
        /**
         * This variable is read from the parameter file through a parameter called 'Sigma'.
         */
        double sigma;
        /**
         * This variable is read from the parameter file through a parameter called 'Sign'.
         */
        double sign;
        unsigned int npoint;

        std::vector<double> radial_position;
        std::vector<double> geotherm;

        /**
         * Inner radius.
         */
        double R0;

        /**
         * Outer radius.
         */
        double R1;
    };
  }
}

#endif
