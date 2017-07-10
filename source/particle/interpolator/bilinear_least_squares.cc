/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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

#include <aspect/particle/interpolator/bilinear_least_squares.h>
#include <aspect/postprocess/particles.h>
#include <aspect/simulator.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/full_matrix.templates.h>

#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      template <int dim>
      std::vector<std::vector<double> >
      BilinearLeastSquares<dim>::properties_at_points(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                                      const std::vector<Point<dim> > &positions,
                                                      const ComponentMask &selected_properties,
                                                      const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        const unsigned int n_particle_properties = particles.begin()->second.get_properties().size();

        const unsigned int property_index = selected_properties.first_selected_component(selected_properties.size());

        AssertThrow(property_index != numbers::invalid_unsigned_int,
                    ExcMessage("Internal error: the particle property interpolator was "
                               "called without a specified component to interpolate."));

        AssertThrow(dim == 2,
                    ExcMessage("Currently, the particle interpolator `bilinear' is only supported for 2D models."));

        AssertThrow(selected_properties.n_selected_components(n_particle_properties) == 1,
                    ExcNotImplemented("Interpolation of multiple components is not supported."));


        const Point<dim> approximated_cell_midpoint = std::accumulate (positions.begin(), positions.end(), Point<dim>())
                                                      / static_cast<double> (positions.size());

        typename parallel::distributed::Triangulation<dim>::active_cell_iterator found_cell;

        if (cell == typename parallel::distributed::Triangulation<dim>::active_cell_iterator())
          {
            // We can not simply use one of the points as input for find_active_cell_around_point
            // because for vertices of mesh cells we might end up getting ghost_cells as return value
            // instead of the local active cell. So make sure we are well in the inside of a cell.
            Assert(positions.size() > 0,
                   ExcMessage("The particle property interpolator was not given any "
                              "positions to evaluate the particle cell_properties at."));


            found_cell =
              (GridTools::find_active_cell_around_point<> (this->get_mapping(),
                                                           this->get_triangulation(),
                                                           approximated_cell_midpoint)).first;
          }
        else
          found_cell = cell;

        const types::LevelInd cell_index = std::make_pair<unsigned int, unsigned int> (found_cell->level(),found_cell->index());
        const std::pair<typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator,
              typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator> particle_range = particles.equal_range(cell_index);


        std::vector<std::vector<double> > cell_properties(positions.size(),
                                                          std::vector<double>(n_particle_properties,
                                                                              numbers::signaling_nan<double>()));

        const unsigned int n_particles = std::distance(particle_range.first,particle_range.second);

        AssertThrow(n_particles != 0,
                    ExcMessage("At least one cell contained no particles. The `bilinear'"
                               "interpolation scheme does not support this case. "));


        // Noticed that the size of matrix A is n_particles x matrix_dimension
        // which usually is not a square matrix. Therefore, we solve Ax=r by
        // solving A^TAx= A^Tr.
        const unsigned int matrix_dimension = 4;
        dealii::LAPACKFullMatrix<double> A(n_particles, matrix_dimension);
        Vector<double> r(n_particles);
        r = 0;

        unsigned int index = 0;
        const double cell_diameter = found_cell->diameter();
        for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator particle = particle_range.first;
             particle != particle_range.second; ++particle, ++index)
          {
            const double particle_property_value = particle->second.get_properties()[property_index];
            r[index] = particle_property_value;

            const Point<dim> position = particle->second.get_location();
            A(index,0) = 1;
            A(index,1) = (position[0] - approximated_cell_midpoint[0])/cell_diameter;
            A(index,2) = (position[1] - approximated_cell_midpoint[1])/cell_diameter;
            A(index,3) = (position[0] - approximated_cell_midpoint[0]) * (position[1] - approximated_cell_midpoint[1])/std::pow(cell_diameter,2);
          }

        dealii::LAPACKFullMatrix<double> B(matrix_dimension, matrix_dimension);

        Vector<double> c_ATr(matrix_dimension);
        Vector<double> c(matrix_dimension);

        const double threshold = 1e-15;
        unsigned int index_positions = 0;

        // Matrix A can be rank deficient if it does not have full rank, therefore singular.
        // To circumvent this issue, we solve A^TAx=A^Tr by using singular value
        // decomposition (SVD).
        A.Tmmult(B, A, false);
        A.Tvmult(c_ATr,r);

        dealii::LAPACKFullMatrix<double> B_inverse(B);
        B_inverse.compute_inverse_svd(threshold);
        B_inverse.vmult(c, c_ATr);

        for (typename std::vector<Point<dim> >::const_iterator itr = positions.begin(); itr != positions.end(); ++itr, ++index_positions)
          {
            const Point<dim> support_point = *itr;
            double interpolated_value = c[0] +
                                        c[1]*(support_point[0] - approximated_cell_midpoint[0])/cell_diameter +
                                        c[2]*(support_point[1] - approximated_cell_midpoint[1])/cell_diameter +
                                        c[3]*(support_point[0] - approximated_cell_midpoint[0])*(support_point[1] - approximated_cell_midpoint[1])/std::pow(cell_diameter,2);

            // Overshoot and undershoot correction of interpolated particle property.
            if (use_global_valued_limiter)
              {
                interpolated_value = std::min(interpolated_value, global_maximum_particle_properties[property_index]);
                interpolated_value = std::max(interpolated_value, global_minimum_particle_properties[property_index]);
              }

            cell_properties[index_positions][property_index] = interpolated_value;
          }
        return cell_properties;
      }

      template <int dim>
      void
      BilinearLeastSquares<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Interpolator");
            {
              prm.enter_subsection("Bilinear least squares");
              {
                prm.declare_entry ("Global particle property maximum",
                                   boost::lexical_cast<std::string>(std::numeric_limits<double>::max()),
                                   Patterns::List(Patterns::Double ()),
                                   "The maximum global particle property values that will be used as a "
                                   "limiter for the bilinear least squares interpolation. The number of the input "
                                   "'Global particle property maximum' values separated by ',' has to be "
                                   "the same as the number of particle properties.");
                prm.declare_entry ("Global particle property minimum",
                                   boost::lexical_cast<std::string>(-std::numeric_limits<double>::max()),
                                   Patterns::List(Patterns::Double ()),
                                   "The minimum global particle property that will be used as a "
                                   "limiter for the bilinear least squares interpolation. The number of the input "
                                   "'Global particle property minimum' values separated by ',' has to be "
                                   "the same as the number of particle properties.");
                prm.declare_entry("Use limiter", "false",
                                  Patterns::Bool (),
                                  "Whether to apply a global particle property limiting scheme to the interpolated "
                                  "particle properties.");

              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      BilinearLeastSquares<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Interpolator");
            {
              prm.enter_subsection("Bilinear least squares");
              {
                use_global_valued_limiter = prm.get_bool("Use limiter");
                if (use_global_valued_limiter)
                  {
                    global_maximum_particle_properties = Utilities::string_to_double(Utilities::split_string_list(prm.get("Global particle property maximum")));
                    global_minimum_particle_properties = Utilities::string_to_double(Utilities::split_string_list(prm.get("Global particle property minimum")));

                    const Postprocess::Particles<dim> *particle_postprocessor = this->template find_postprocessor<Postprocess::Particles<dim> >();
                    const unsigned int n_property_components = particle_postprocessor->get_particle_world().get_property_manager().get_n_property_components();

                    AssertThrow(global_minimum_particle_properties.size() == n_property_components,
                                ExcMessage("Make sure that the size of list 'Global minimum particle property' "
                                           "is equivalent to the number of particle properties."));

                    AssertThrow(global_maximum_particle_properties.size() == n_property_components,
                                ExcMessage("Make sure that the size of list 'Global maximum particle property' "
                                           "is equivalent to the number of particle properties."));
                  }
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(BilinearLeastSquares,
                                            "bilinear least squares",
                                            "Interpolates particle properties onto a vector of points using a "
                                            "bilinear least squares method. Currently only 2D models are "
                                            "supported. Note that deal.II must be configured with BLAS/LAPACK.")
    }
  }
}
