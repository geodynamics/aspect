/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/initial_composition/ascii_data_slice.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    AsciiDataSlice<dim>::AsciiDataSlice ()
    {}


    template <int dim>
    void
    AsciiDataSlice<dim>::initialize ()
    {
      AssertThrow(dim==2,
                  ExcMessage("The ascii data slice plugin can only be used in 2d models"));

      Utilities::AsciiDataInitial<dim, 3>::initialize(this->n_compositional_fields());

      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical)
        {
          rotation_matrix[0][0] = 1.0;
          rotation_matrix[1][1] = 1.0;
          rotation_matrix[2][2] = 1.0;

          AssertThrow(slice_normal_vector.norm() > std::numeric_limits<double>::min(),
                      ExcMessage("The normal vector of the slice can not have length zero."));

          // Set up the normal vector of an unrotated 2D spherical shell
          // that by default lies in the x-y plane.
          const Tensor<1,3> unrotated_normal_vector({0.0,0.0,1.0});

          Tensor<1,3> rotated_normal_vector = slice_normal_vector / slice_normal_vector.norm();

          if ((rotated_normal_vector - unrotated_normal_vector).norm() > 1.e-3)
            {
              const Tensor<1,3> point_on_slice = cross_product_3d(rotated_normal_vector, unrotated_normal_vector);

              // Calculate the crossing line of the two normals,
              // which will be the rotation axis to transform the one
              // normal into the other
              Tensor<1,3> rotation_axis = cross_product_3d(unrotated_normal_vector,rotated_normal_vector);
              rotation_axis /= rotation_axis.norm();

              // Calculate the rotation angle from the inner product rule
              const double rotation_angle = std::acos(rotated_normal_vector*unrotated_normal_vector);
              rotation_matrix = rotation_matrix_from_axis(rotation_axis,rotation_angle);

              // Now apply the rotation that will project point_one onto the known point
              // (0,1,0).
              const Tensor<1,3> rotated_point_one = transpose(rotation_matrix) * point_on_slice;
              const Tensor<1,3> final_point_one ({0.0,1.0,0.0});

              const double second_rotation_angle = std::acos(rotated_point_one*final_point_one);
              Tensor<1,3> second_rotation_axis = cross_product_3d(final_point_one,rotated_point_one);

              if (second_rotation_axis.norm() > std::numeric_limits<double>::min())
                {
                  second_rotation_axis /= second_rotation_axis.norm();
                  const Tensor<2,3> second_rotation_matrix = rotation_matrix_from_axis(second_rotation_axis,second_rotation_angle);

                  // The final rotation used for the model will be the combined
                  // rotation of the two operation above. This is achieved by a
                  // matrix multiplication of the rotation matrices.
                  // This concatenation of rotations is the reason for using a
                  // rotation matrix instead of a combined rotation_axis + angle
                  rotation_matrix = rotation_matrix * second_rotation_matrix;
                }
            }
        }
    }


    template <int dim>
    double
    AsciiDataSlice<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int n_comp) const
    {
      // make a 3D point based on the 2D position
      Point<3> position_3d(position[0], position[1], 0.0);

      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::cartesian)
        position_3d[2] = slice_normal_vector[2];
      else if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical)
        {
          Tensor<1,3> position_tensor({position_3d[0], position_3d[1], position_3d[2]});
          position_3d = Point<3> (rotation_matrix * position_tensor);
        }
      else
        AssertThrow(false,
                    ExcMessage("The Ascii data slice model can only be used in spherical or "
                               "cartesian geometry models."));

      return Utilities::AsciiDataInitial<dim, 3>::get_data_component(position_3d, n_comp);
    }



    template <int dim>
    Tensor<2,3>
    AsciiDataSlice<dim>::rotation_matrix_from_axis (const Tensor<1,3> &rotation_axis,
                                                    const double rotation_angle) const
    {
      Tensor<2,3> rotation_matrix;
      rotation_matrix[0][0] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[0] + std::cos(rotation_angle);
      rotation_matrix[0][1] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[1] - rotation_axis[2] * std::sin(rotation_angle);
      rotation_matrix[0][2] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[2] + rotation_axis[1] * std::sin(rotation_angle);
      rotation_matrix[1][0] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[0] + rotation_axis[2] * std::sin(rotation_angle);
      rotation_matrix[1][1] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[1] + std::cos(rotation_angle);
      rotation_matrix[1][2] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[2] - rotation_axis[0] * std::sin(rotation_angle);
      rotation_matrix[2][0] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[0] - rotation_axis[1] * std::sin(rotation_angle);
      rotation_matrix[2][1] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[1] + rotation_axis[0] * std::sin(rotation_angle);
      rotation_matrix[2][2] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[2] + std::cos(rotation_angle);
      return rotation_matrix;
    }



    template <int dim>
    void
    AsciiDataSlice<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-composition/ascii-data/test/",
                                                          "box_2d.txt");
        prm.enter_subsection("Ascii data slice");
        {
          prm.declare_entry ("Slice normal vector", "0.0,0.0,1.0",
                             Patterns::Anything (),
                             "Point that determines the plane in which the 2D slice lies in. "
                             "In a spherical geometry, the slice will go though the center of the "
                             "model domain. In Cartesian geometry, the slice will be in the third "
                             "coordinate direction. The point has to be given in Cartesian "
                             "coordinates.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiDataSlice<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm);

        prm.enter_subsection("Ascii data slice");
        {
          std::vector<double> point = Utilities::string_to_double(Utilities::split_string_list(prm.get("Slice normal vector")));

          AssertThrow(point.size() == 3,
                      ExcMessage("The Point on slice in the Ascii data slice model needs "
                                 "to be given in three dimensions; in other words, as three "
                                 "numbers, separated by commas."));

          for (unsigned int d=0; d<3; d++)
            slice_normal_vector[d] = point[d];
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
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(AsciiDataSlice,
                                              "ascii data slice",
                                              "Implementation of a model in which the initial "
                                              "composition is derived from files containing data "
                                              "in ascii format. For more details on the data format, "
                                              "see the ascii data plugin. The speacial feature of "
                                              "this model is that it reads in a 3d ascii data file, "
                                              "but then only uses a slice of it in a 2d model (so it "
                                              "can only be used in 2d.")
  }
}
