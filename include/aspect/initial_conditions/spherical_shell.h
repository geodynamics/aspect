//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__initial_conditions_spherical_shell_h
#define __aspect__initial_conditions_spherical_shell_h

#include <aspect/initial_conditions/interface.h>
#include <deal.II/base/table.h>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * A class that describes a perturbed initial temperature field for
     * the spherical shell.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class SphericalHexagonalPerturbation : public Interface<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;
    };


    /**
     * A class that describes a perturbed initial temperature field for
     * the spherical shell.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class SphericalGaussianPerturbation : public Interface<dim>
    {
      public:

        SphericalGaussianPerturbation<dim>();
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        double angle;
        double depth;
        double amplitude;
        double sigma;
        double sign;
        unsigned int npoint;
        dealii::Table<1,double> radial_position;
        dealii::Table<1,double> geotherm;

    };
  }
}

#endif
