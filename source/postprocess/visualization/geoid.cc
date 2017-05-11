/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/geoid.h>
#include <aspect/postprocess/geoid.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Geoid<dim>::
      Geoid ()
        :
        DataPostprocessorScalar<dim> ("geoid",
                                      update_q_points)
      {}



      template <int dim>
      void
      Geoid<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> > &,
                                         const std::vector<std::vector<Tensor<1,dim> > > &,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const Postprocess::Geoid<dim> *geoid_postprocessor = this->template find_postprocessor<Postprocess::Geoid<dim> >();
        const internal::HarmonicCoefficients geoid_coefficients = geoid_postprocessor->get_geoid_coefficients();

        for (unsigned int q=0; q<evaluation_points.size(); ++q)
          {
            const std_cxx11::array<double,dim> spherical_position =
              Utilities::spherical_coordinates(evaluation_points[q]);

            unsigned int i = 0;
            unsigned int j = 0;
            unsigned int k = 0;

            computed_quantities[q](0) = 0.0;

            while (k<geoid_coefficients.sine_coefficients.size())
              {
                computed_quantities[q](0) += geoid_coefficients.cosine_coefficients[k]
                                             * boost::math::spherical_harmonic_r(i,j,spherical_position[2],spherical_position[1])
                                             + geoid_coefficients.sine_coefficients[k]
                                             * boost::math::spherical_harmonic_i(i,j,spherical_position[2],spherical_position[1]);
                if (j<i)
                  ++j;
                else
                  {
                    ++i;
                    j = 0;
                  }
                ++k;
              }
          }
      }

      template <int dim>
      std::list<std::string>
      Geoid<dim>::required_other_postprocessors () const
      {
        std::list<std::string> names;
        names.push_back("geoid");
        return names;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Geoid,
                                                  "geoid",
                                                  "A visualization output object that generates output "
                                                  "for the geoid.")
    }
  }
}
