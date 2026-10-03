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


#include "rift_box_initial_plastic_strain.h"
#include <aspect/postprocess/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    RiftBoxInitialPlasticStrain<dim>::RiftBoxInitialPlasticStrain ()
    {}

    template <int dim>
    void
    RiftBoxInitialPlasticStrain<dim>::initialize ()
    {
      AssertThrow(Plugins::plugin_type_matches<MaterialModel::ViscoPlastic<dim>>(this->get_material_model()),
                  ExcMessage("This initial condition only makes sense in combination with the visco_plastic material model."));


      // From shear_bands.cc
      Point<dim> extents_max;
      TableIndices<dim> size_idx;
      for (unsigned int d=0; d<dim; ++d)
        size_idx[d] = grid_intervals[d]+1;

      Table<dim,double> white_noise;
      white_noise.TableBase<dim,double>::reinit(size_idx);

      // Max of each direction (m)
      extents_max = extents + origin;

      // Add initial topography to vertical direction
      extents_max[dim-1] += maximum_initial_topography;

      for (unsigned int d=0; d<dim; ++d)
        {
          grid_extents[d].first=origin[d];
          grid_extents[d].second=extents_max[d];
        }

      // use a fixed number as seed for random generator
      // this is important if we run the code on more than 1 processor
      std::srand(seed);

      TableIndices<dim> idx;

      for (unsigned int i=0; i<white_noise.size()[0]; ++i)
        {
          idx[0] = i;
          for (unsigned int j=0; j<white_noise.size()[1]; ++j)
            {
              idx[1] = j;
              if (dim == 3)
                {
                  for (unsigned int k=0; k<white_noise.size()[dim-1]; ++k)
                    {
                      idx[dim-1] = k;
                      // std::rand will give a value between zero and RAND_MAX (usually INT_MAX).
                      // The modulus of this value and 10000, gives a value between 0 and 10000-1.
                      // Subsequently dividing by 5000.0 will give value between 0 and 2 (excluding 2).
                      // Subtracting 1 will give a range [-1,1)
                      // Because we want values [0,1), we change our white noise computation to:
                      white_noise(idx) = ((std::rand() % 10000) / 10000.0);
                    }
                }
              else
                white_noise(idx) = ((std::rand() % 10000) / 10000.0);
            }
        }

      interpolate_noise = new Functions::InterpolatedUniformGridData<dim> (grid_extents,
                                                                           grid_intervals,
                                                                           white_noise);
    }

    template <int dim>
    double
    RiftBoxInitialPlasticStrain<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      // If n_comp does not represent the strain field,
      // return 0 right away.
      if (n_comp != strain_composition_number)
        return 0.0;

      // Initiate distance with large value
      double distance_to_rift_axis = 1e23;
      double temp_distance = 0;

      const Point<dim> natural_coords = position;

      // Loop over all line segments
      for (unsigned int i_segments = 0; i_segments < point_list.size(); ++i_segments)
        {
          if (dim == 2)
            temp_distance = std::abs(natural_coords[0]-point_list[i_segments][0][0]);
          else
            {
              // Get the surface coordinates by dropping the last coordinate
              const Point<2> surface_position = Point<2>(natural_coords[0],natural_coords[1]);
              temp_distance = std::abs(Utilities::distance_to_line(point_list[i_segments], surface_position));
            }

          // Get the minimum distance
          distance_to_rift_axis = std::min(distance_to_rift_axis, temp_distance);
        }

      // Smoothing of noise with height (depth with respect to the unperturbed surface)
      const double depth_smoothing = 0.5 * (1.0 + std::tanh((position[dim-1] - strain_height) / strain_halfwidth));
      // Smoothing of noise with lateral distance to the rift axis
      const double noise_amplitude = A * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma,2)))) * depth_smoothing;

      // Add randomness
      return noise_amplitude * interpolate_noise->value(natural_coords);
    }

    template <int dim>
    void
    RiftBoxInitialPlasticStrain<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Rift box initial plastic strain");
        {
          prm.declare_entry ("Random number generator seed", "0",
                             Patterns::Integer (0),
                             "The value of the seed used for the random number generator. "
                             "Units: none.");
          prm.declare_entry ("Standard deviation of Gaussian noise amplitude distribution", "20000",
                             Patterns::Double (0),
                             "The standard deviation of the Gaussian distribution of the amplitude of the strain noise. "
                             "Note that this parameter is taken to be the same for all rift segments. "
                             "Units: $m$ or degrees.");
          prm.declare_entry ("Maximum amplitude of Gaussian noise amplitude distribution", "0.2",
                             Patterns::Double (0),
                             "The amplitude of the Gaussian distribution of the amplitude of the strain noise. "
                             "Note that this parameter is taken to be the same for all rift segments. "
                             "Units: none.");
          prm.declare_entry ("Depth around which Gaussian noise is smoothed out", "40000",
                             Patterns::Double (0),
                             "The depth around which smoothing out of the strain noise starts with a hyperbolic tangent. "
                             "Note that this parameter is taken to be the same for all rift segments. Also, depth "
                             "is with respect to the unperturbed surface of the box domain. "
                             "Units: $m$.");
          prm.declare_entry ("Halfwidth with which Gaussian noise is smoothed out in depth", "40000",
                             Patterns::Double (0),
                             "The halfwidth with which smoothing out of the strain noise is done with a hyperbolic tangent. "
                             "Note that this parameter is taken to be the same for all rift segments. "
                             "Units: $m$.");
          prm.declare_entry ("Rift axis line segments",
                             "",
                             Patterns::Anything(),
                             "Set the line segments that represent the rift axis. In 3d each segment is made up of "
                             "two points that represent horizontal coordinates (x,y). "
                             "The exact format for the point list describing the segments is "
                             "\"x1,y1>x2,y2;x2,y2>x3,y3;x4,y4>x5,y5\". In 2d, a segment is made up by 1 horizontal "
                             "x coordinate: \"x1;x2;x3\". Note that the segments can be connected "
                             "or isolated. Units: $m$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    RiftBoxInitialPlasticStrain<dim>::parse_parameters (ParameterHandler &prm)
    {
      // Check that there is a compositional field called strain and retrieve its index
      AssertThrow(this->introspection().compositional_name_exists("plastic_strain"),
                  ExcMessage("This plugin requires a compositional field named plastic_strain. "));
      strain_composition_number = this->introspection().compositional_index_for_name("plastic_strain");

      AssertThrow(Plugins::plugin_type_matches<GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("This initial condition only makes sense in combination with the box geometry model."));

      // Get the number of mesh cells in each direction
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box");
        {
          origin[0] = prm.get_double ("Box origin X coordinate");
          extents[0] = prm.get_double ("X extent");
          repetitions[0] = prm.get_integer ("X repetitions");
          origin[1] = prm.get_double ("Box origin Y coordinate");
          extents[1] = prm.get_double ("Y extent");
          repetitions[1] = prm.get_integer ("Y repetitions");
          if (dim >= 3)
            {
              origin[dim-1] = prm.get_double ("Box origin Z coordinate");
              extents[dim-1] = prm.get_double ("Z extent");
              repetitions[dim-1] = prm.get_integer ("Z repetitions");
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      unsigned int refinement;
      prm.enter_subsection("Mesh refinement");
      {
        refinement = prm.get_integer("Initial adaptive refinement") + prm.get_integer("Initial global refinement");
      }
      prm.leave_subsection();

      // The grid intervals equal the number of mesh cells.
      // refinement, repetitions and grid_intervals are integers,
      // but pow returns a double, so convert to int.
      for (unsigned int d = 0; d < dim; ++d)
        grid_intervals[d] = repetitions[d] * int(std::pow(2, refinement));

      // Both initial mesh deformation and initial topography models
      // can distort the initial mesh such that the vertical coordinate
      // of some locations is larger than the initial, unperturbed,
      // vertical extent of the model domain. Since not all initial topography
      // models (correctly) return a maximum topography, and
      // initial mesh deformation models do not at all, we here assume
      // a maximal initial topography of 10 km. This is added to the
      // vertical extent.
      //
      // In the vertical direction, we add enough grid cells
      // to cover 10 km at the same resolution. Round up the number of cells.
      grid_intervals[dim-1] += static_cast<int>(std::ceil(maximum_initial_topography / (extents[dim-1] / grid_intervals[dim-1])));

      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Rift box initial plastic strain");
        sigma                = prm.get_double ("Standard deviation of Gaussian noise amplitude distribution");
        A                    = prm.get_double ("Maximum amplitude of Gaussian noise amplitude distribution");
        seed                 = prm.get_integer ("Random number generator seed");
        strain_height        = origin[dim-1] + extents[dim-1] - prm.get_double ("Depth around which Gaussian noise is smoothed out");
        strain_halfwidth     = prm.get_double ("Halfwidth with which Gaussian noise is smoothed out in depth");

        // Read in the string of rift segments
        const std::string temp_all_segments = prm.get("Rift axis line segments");
        // Split the string into segment strings
        const std::vector<std::string> temp_segments = Utilities::split_string_list(temp_all_segments,';');
        // The number of segments, each consisting of a begin and an end point in 3d and one point in 3d
        const unsigned int n_temp_segments = temp_segments.size();
        point_list.resize(n_temp_segments);

        // Loop over the segments to extract the points
        for (unsigned int i_segment = 0; i_segment < n_temp_segments; i_segment++)
          {
            // In 3d a line segment consists of 2 points,
            // in 2d of only 1 (ridge axis orthogonal to x and y).
            // Also, a 3d point has 2 coordinates (x and y),
            // a 2d point only 1 (x).
            const std::vector<std::string> temp_segment = Utilities::split_string_list(temp_segments[i_segment],'>');

            if (dim == 3)
              {
                // const std::vector<std::string> temp_segment = Utilities::split_string_list(temp_segments[i_segment],'>');
                AssertThrow(temp_segment.size() == 2,ExcMessage ("The given coordinate '" + temp_segment[i_segment] + "' is not correct. "
                                                                 "It should only contain 2 parts: "
                                                                 "the two points of the segment, separated by a '>'."));

              }
            else
              {
                // Add the point to the list of points for this segment
                // As we're in 2d all segments correspond to 1 point consisting of 1 coordinate
                AssertThrow(temp_segment.size() == 1,ExcMessage ("The given coordinate '" + temp_segment[i_segment] + "' is not correct. "
                                                                 "It should only contain 1 part: "
                                                                 "the point representing the rift axis."));
              }

            // Loop over the 2 points of each segment
            for (unsigned int i_points = 0; i_points < dim-1; i_points++)
              {
                std::vector<double> temp_point = Utilities::string_to_double(Utilities::split_string_list(temp_segment[i_points],','));
                if (dim == 3)
                  {
                    AssertThrow(temp_point.size() == 2,ExcMessage ("The given coordinates of segment '" + temp_segment[i_points] + "' are not correct. "
                                                                   "It should only contain 2 parts: "
                                                                   "the x and y coordinates of the segment begin/end point, separated by a ','."));
                  }
                else
                  {
                    AssertThrow(temp_point.size() == 1,ExcMessage ("The given coordinates of segment '" + temp_segment[i_points] + "' are not correct. "
                                                                   "It should only contain 1 part: "
                                                                   "the one coordinate of the segment end point."));
                  }

                // Add the point to the list of points for this segment
                point_list[i_segment][i_points][0] = temp_point[0];
                point_list[i_segment][i_points][1] = temp_point[dim-2];
              }
          }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(RiftBoxInitialPlasticStrain,
                                              "rift box initial plastic strain",
                                              "Specify the plastic strain initial compositional field value based on "
                                              "the distance to a list of line segments "
                                              "and the user-defined Gaussian distribution around these segments, "
                                              "combined with random noise. "
                                              "A maximum initial topography of 10 km can be accommodated. Initial topography that "
                                              "is higher than 10 km will lead to plastic strain values set to the strain value at 10 km."
                                              "This plugin only works in combination with the box geometry model and the visco_plastic "
                                              "material model.")
  }
}
