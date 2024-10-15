/*
  Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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
#include <aspect/geometry_model/initial_topography_model/ascii_data.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/two_merged_boxes.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>

#include <deal.II/base/parameter_handler.h>
#include <array>



namespace aspect
{
  namespace InitialTopographyModel
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
      :
      surface_boundary_id(1)
    {}


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      Utilities::AsciiDataBoundary<dim>::initialize({surface_boundary_id},
                                                    1);
    }


    template <int dim>
    double
    AsciiData<dim>::value (const Point<dim-1> &surface_point) const
    {
      // In a first step, create a global 'dim'-dimensional point that we can pass to the
      // function expression as input -- because the function is a dim-dimensional
      // function.
      //
      // Different geometry models pass the surface point in in different ways.
      // In the following, we will first normalize the input to a dim-dimensional
      // point with a dummy vertical/radial coordinate that, one would hope,
      // the AsciiDataBoundary class will then simply ignore.
      Point<dim> global_point;
      if (Plugins::plugin_type_matches<const GeometryModel::Box<dim>> (this->get_geometry_model()) ||
          Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model()))
        {
          for (unsigned int d=0; d<dim-1; ++d)
            global_point[d] = surface_point[d];

          // Now for the vertical component:
          global_point[dim-1] = 0;
        }
      else if (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>> (this->get_geometry_model()) ||
               Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model()) ||
               Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>> (this->get_geometry_model()))
        {
          // AsciiDataBoundary always expects to get the input
          // parameters for its functions in Cartesian
          // coordinates. The get_data_component function then changes
          // the coordinate system, or more precisely it asks the
          // geometry model to convert the point into its natural
          // coordinate system, before doing the table lookup.
          //
          // This is of course all not very efficient (TODO: Think
          // about a better scheme), but the first step then needs to
          // be to convert what we have into Cartesian coordinates...
          std::array<double, dim> point;
          point[0] = 6371000.0;
          for (unsigned int d=0; d<dim-1; ++d)
            point[d+1] = surface_point[d];

          global_point = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(point);
        }
      else
        AssertThrow(false, ExcNotImplemented());

      const double topo = this->Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id,
                                                                                      global_point,
                                                                                      0);

      return topo;
    }



    template <int dim>
    Tensor<1,dim-1>
    AsciiData<dim>::vector_gradient(const Point<dim> &point) const
    {
      return Utilities::AsciiDataBoundary<dim>::vector_gradient(surface_boundary_id, point,0);
    }



    template <int dim>
    double
    AsciiData<dim>::max_topography () const
    {
      return Utilities::AsciiDataBoundary<dim>::get_maximum_component_value(surface_boundary_id,0);
    }



    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Initial topography model");
        {
          Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                            "$ASPECT_SOURCE_DIR/data/geometry-model/initial-topography-model/ascii-data/test/",
                                                            "box_2d_%s.0.txt");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Initial topography model");
        {
          Utilities::AsciiDataBase<dim>::parse_parameters(prm);
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
  namespace InitialTopographyModel
  {
    ASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(AsciiData,
                                             "ascii data",
                                             "Implementation of a model in which the surface "
                                             "topography is derived from a file containing data "
                                             "in ascii format. The following geometry models "
                                             "are currently supported: box, chunk. "
                                             "Note the required format of the "
                                             "input data: The first lines may contain any number of comments "
                                             "if they begin with `#', but one of these lines needs to "
                                             "contain the number of grid points in each dimension as "
                                             "for example `# POINTS: 3 3'. "
                                             "The order of the data columns "
                                             "has to be `x', `Topography [m]' in a 2d model and "
                                             " `x', `y', `Topography [m]' in a 3d model, which means that "
                                             "there has to be a single column "
                                             "containing the topography. "
                                             "Note that the data in the input "
                                             "file needs to be sorted in a specific order: "
                                             "the first coordinate needs to ascend first, "
                                             "followed by the second in order to "
                                             "assign the correct data to the prescribed coordinates. "
                                             "If you use a spherical model, "
                                             "then the assumed grid changes. "
                                             "`x' will be replaced by the azimuth angle (longitude) "
                                             "in radians (between zero and $2\\pi$, not between $-\\pi$ "
                                             "corresponding to 180 degrees west, and $+\\pi$ corresponding "
                                             "to 180 degrees east), "
                                             "and `y' by the polar angle in radians (between $0$ and $\\pi$) measured "
                                             "positive from the north pole. The grid will be assumed to be "
                                             "a longitude-colatitude grid. Note that the order "
                                             "of spherical coordinates is `phi', `theta' "
                                             "and not `theta', `phi', since this allows "
                                             "for dimension independent expressions.")
  }
}
