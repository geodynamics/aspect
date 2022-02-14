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

      rotation_matrix = Utilities::compute_rotation_matrix_for_slice(first_point_on_slice, second_point_on_slice);
    }


    template <int dim>
    double
    AsciiDataSlice<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int n_comp) const
    {
      // Compute the coordinates of a 3d point based on the 2D position.
      Tensor<1,3> position_tensor({position[0], position[1], 0.0});
      Point<3> position_3d (rotation_matrix * position_tensor);

      return Utilities::AsciiDataInitial<dim, 3>::get_data_component(position_3d, n_comp);
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
          prm.declare_entry ("First point on slice", "0.0,1.0,0.0",
                             Patterns::Anything (),
                             "Point that determines the plane in which the 2D slice lies in. "
                             "The slice will go through this point, the point defined by the "
                             "parameter 'Second point on slice', and the center of the model "
                             "domain. After the rotation, this first point will lie along the "
                             "(0,1,0) axis of the coordinate system. The coordinates of the "
                             "point have to be given in Cartesian coordinates.");
          prm.declare_entry ("Second point on slice", "1.0,0.0,0.0",
                             Patterns::Anything (),
                             "Second point that determines the plane in which the 2D slice lies in. "
                             "The slice will go through this point, the point defined by the "
                             "parameter 'First point on slice', and the center of the model "
                             "domain. The coordinates of the point have to be given in Cartesian "
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
          std::vector<double> point_one = Utilities::string_to_double(Utilities::split_string_list(prm.get("First point on slice")));
          std::vector<double> point_two = Utilities::string_to_double(Utilities::split_string_list(prm.get("Second point on slice")));

          AssertThrow(point_one.size() == 3 && point_two.size() == 3,
                      ExcMessage("The points on the slice in the Ascii data slice model need "
                                 "to be given in three dimensions; in other words, as three "
                                 "numbers, separated by commas."));

          for (unsigned int d=0; d<3; d++)
            first_point_on_slice[d] = point_one[d];

          for (unsigned int d=0; d<3; d++)
            second_point_on_slice[d] = point_two[d];
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
                                              "see the ascii data plugin. The special feature of "
                                              "this model is that it reads in a 3d ascii data file, "
                                              "but then only uses a slice of it in a 2d model (so it "
                                              "can only be used in 2d.")
  }
}
