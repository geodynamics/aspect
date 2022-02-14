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
#include <aspect/initial_composition/ascii_data.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
    {}


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      if (slice_data == true)
        {
          rotation_matrix = Utilities::compute_rotation_matrix_for_slice(first_point_on_slice, second_point_on_slice);
          ascii_data_slice->initialize(this->n_compositional_fields());
        }
      else
        {
          ascii_data_initial->initialize(this->n_compositional_fields());
        }
    }


    template <int dim>
    double
    AsciiData<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int n_comp) const
    {
      if (slice_data == true)
        {
          // Compute the coordinates of a 3d point based on the 2D position.
          Tensor<1,3> position_tensor({position[0], position[1], 0.0});
          Point<3> position_3d (rotation_matrix * position_tensor);

          return ascii_data_slice->get_data_component(position_3d, n_comp);
        }

      return ascii_data_initial->get_data_component(position,n_comp);
    }


    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-composition/ascii-data/test/",
                                                          "box_2d.txt");

        prm.enter_subsection("Ascii data model");
        {
          prm.declare_entry("Slice dataset in 2D plane", "false",
                            Patterns::Bool (),
                            "Whether to use a 2D data slice from a 3D data file "
                            "or the entire data file. Slicing a 3D dataset is "
                            "only supported for 2D models.");
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
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Ascii data model");
        {
          slice_data = prm.get_bool ("Slice dataset in 2D plane");
          std::vector<double> point_one = Utilities::string_to_double(Utilities::split_string_list(prm.get("First point on slice")));
          std::vector<double> point_two = Utilities::string_to_double(Utilities::split_string_list(prm.get("Second point on slice")));

          AssertThrow(point_one.size() == 3 && point_two.size() == 3,
                      ExcMessage("The points on the slice in the Ascii data model need "
                                 "to be given in three dimensions; in other words, as three "
                                 "numbers, separated by commas."));

          for (unsigned int d=0; d<3; d++)
            first_point_on_slice[d] = point_one[d];

          for (unsigned int d=0; d<3; d++)
            second_point_on_slice[d] = point_two[d];
        }
        prm.leave_subsection();

        if (slice_data == true)
          {
            AssertThrow(dim==2,
                        ExcMessage("The ascii data plugin can only slice data in 2d models."));

            ascii_data_slice = std::make_unique<Utilities::AsciiDataInitial<dim,3>>();
            ascii_data_slice->initialize_simulator(this->get_simulator());
            ascii_data_slice->parse_parameters(prm);
          }
        else
          {
            ascii_data_initial = std::make_unique<Utilities::AsciiDataInitial<dim>>();
            ascii_data_initial->initialize_simulator(this->get_simulator());
            ascii_data_initial->parse_parameters(prm);
          }

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
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(AsciiData,
                                              "ascii data",
                                              "Implementation of a model in which the initial "
                                              "composition is derived from files containing data "
                                              "in ascii format. Note the required format of the "
                                              "input data: The first lines may contain any number of comments "
                                              "if they begin with `#', but one of these lines needs to "
                                              "contain the number of grid points in each dimension as "
                                              "for example `# POINTS: 3 3'. "
                                              "The order of the data columns "
                                              "has to be `x', `y', `composition1', `composition2', "
                                              "etc. in a 2d model and `x', `y', `z', `composition1', "
                                              "`composition2', etc. in a 3d model, according "
                                              "to the number of compositional fields, which means that "
                                              "there has to be a single column "
                                              "for every composition in the model."
                                              "Note that the data in the input "
                                              "files need to be sorted in a specific order: "
                                              "the first coordinate needs to ascend first, "
                                              "followed by the second and the third at last in order to "
                                              "assign the correct data to the prescribed coordinates. "
                                              "If you use a spherical model, "
                                              "then the assumed grid changes. `x' will be replaced by "
                                              "the radial distance of the point to the bottom of the model, "
                                              "`y' by the azimuth angle and `z' by the polar angle measured "
                                              "positive from the north pole. The grid will be assumed to be "
                                              "a latitude-longitude grid. Note that the order "
                                              "of spherical coordinates is `r', `phi', `theta' "
                                              "and not `r', `theta', `phi', since this allows "
                                              "for dimension independent expressions.")
  }
}
