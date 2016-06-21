/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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


#ifndef __aspect__postprocess_visualization_melt_h
#define __aspect__postprocess_visualization_melt_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessor that takes an output vector
       * and computes a variable that represents the melt density at every point.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class MeltDensity
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          MeltDensity ();

          virtual
          void
          compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                             const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                             const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                             const std::vector<Point<dim> >                  &normals,
                                             const std::vector<Point<dim> >                  &evaluation_points,
                                             std::vector<Vector<double> >                    &computed_quantities) const;
      };

      /**
       * A class derived from DataPostprocessor that takes an output vector
       * and computes a variable that represents the compaction viscosity at every point.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class CompactionViscosity
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          CompactionViscosity ();

          virtual
          void
          compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                             const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                             const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                             const std::vector<Point<dim> >                  &normals,
                                             const std::vector<Point<dim> >                  &evaluation_points,
                                             std::vector<Vector<double> >                    &computed_quantities) const;
      };

      /**
       * A class derived from DataPostprocessor that takes an output vector
       * and computes a variable that represents the permeability at every point.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class Permeability
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          Permeability ();

          virtual
          void
          compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                             const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                             const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                             const std::vector<Point<dim> >                  &normals,
                                             const std::vector<Point<dim> >                  &evaluation_points,
                                             std::vector<Vector<double> >                    &computed_quantities) const;
      };

      /**
       * A class derived from DataPostprocessor that takes an output vector
       * and computes a variable that represents the melt viscosity at every point.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class MeltViscosity
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          MeltViscosity ();

          virtual
          void
          compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                             const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                             const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                             const std::vector<Point<dim> >                  &normals,
                                             const std::vector<Point<dim> >                  &evaluation_points,
                                             std::vector<Vector<double> >                    &computed_quantities) const;
      };
    }
  }
}

#endif
