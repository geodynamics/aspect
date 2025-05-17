/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_composition_rift_box_initial_plastic_strain_h
#define _aspect_initial_composition_rift_box_initial_plastic_strain_h

#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/function_lib.h>

namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements initial conditions for the compositional fields
     * based on a functional description provided in the input file.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class RiftBoxInitialPlasticStrain : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        RiftBoxInitialPlasticStrain ();

        /**
         * Initialization function.
         */
        void
        initialize () override;

        /**
         * Return the initial composition as a function of position and number
         * of compositional field. The composition varies as a Gaussian distribution
         * around a user-defined set of line-segments.
         */
        virtual
        double initial_composition (const Point<dim> &position, const unsigned int n_comp) const override;

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
        parse_parameters(ParameterHandler &prm) override;

      private:

        /**
         * The number of the compositional field representing the strain.
         */
        unsigned int strain_composition_number;

        /**
         * The value of the seed for the random number generator
         */
        unsigned int seed;

        /**
         * The maximum amplitude of the Gaussian amplitude of the noise.
         */
        double A;

        /**
         * The standard deviation of the Gaussian amplitude of the noise.
         */
        double sigma;

        /**
         * The height around and halfwidth with which the noise is smoothed out.
         */
        double strain_height;
        double strain_halfwidth;

        /**
         * The list of line segments consisting of two 2d coordinates per segment (even in 2d).
         * The segments represent the rift axis.
         */
        std::vector<std::array<Point<2>,2 >> point_list;

        /**
         * A table with random noise for the
         * second invariant of the strain.
         */
        std::array<unsigned int,dim> grid_intervals;
        Functions::InterpolatedUniformGridData<dim> *interpolate_noise;

        /**
         * A pointer to the initial topography model.
         */
        InitialTopographyModel::Interface<dim> *topo_model;

        /**
         * The min and max of the strain grid in each direction
         */
        std::array<std::pair<double,double>,dim> grid_extents;

        /**
         * Origin of the box in x, y, and z (in 3d) coordinates.
         */
        Point<dim> origin;

        /**
         * Extent of the domain in x-, y-, and z-direction (in 3d).
         */
        Point<dim> extents;

        /**
         * The number of cells in each coordinate direction (x, y, z).
         */
        std::array<unsigned int, dim> repetitions;

        /**
         * The height we add to vertical domain extent
         * to deal with any initial topography.
         */
        const double maximum_initial_topography = 10000.;
    };
  }
}


#endif
