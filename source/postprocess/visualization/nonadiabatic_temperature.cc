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
/*  $Id$  */


#include <aspect/postprocess/visualization.h>
#include <aspect/simulator.h>

#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessor that takes an output vector and
       * computes a variable that represents the non-adiabatic part of the
       * temperature at every point.
       *
       * The member functions are all implementations of those declared in the base
       * class. See there for their meaning.
       */
      template <int dim>
      class NonadiabaticTemperature
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          NonadiabaticTemperature ();

          virtual
          void
          compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                             const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                             const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                             const std::vector<Point<dim> >                  &normals,
                                             const std::vector<Point<dim> >                  &evaluation_points,
                                             std::vector<Vector<double> >                    &computed_quantities) const;
      };


      template <int dim>
      NonadiabaticTemperature<dim>::
      NonadiabaticTemperature ()
        :
        DataPostprocessorScalar<dim> ("nonadiabatic_temperature",
                                      update_values | update_q_points)
      {}



      template <int dim>
      void
      NonadiabaticTemperature<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                         const std::vector<Point<dim> >                  &normals,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points,  ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                 ExcInternalError());
        Assert (uh[0].size() == dim+2,                              ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const double temperature = uh[q][dim+1];

            computed_quantities[q](0) = temperature - this->get_adiabatic_conditions().temperature(evaluation_points[q]);
          }
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(NonadiabaticTemperature,
                                                  "nonadiabatic temperature",
                                                  "A visualization output object that generates output "
                                                  "for the non-adiabatic component of the pressure.")
    }
  }
}
