// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/global.h>
#include <aspect/initial_composition/ascii_data_layered.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    AsciiDataLayered<dim>::AsciiDataLayered ()
      = default;


    template <int dim>
    void
    AsciiDataLayered<dim>::initialize ()
    {
      // First additional field is the layer depth, so we need one additional non-compositional field
      Utilities::AsciiDataLayered<dim>::initialize(this->n_compositional_fields() + 1);
    }


    template <int dim>
    double
    AsciiDataLayered<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int n_comp) const
    {

      // First additional field is the layer depth, so we need the data component from the n_comp+1th column
      return Utilities::AsciiDataLayered<dim>::get_data_component(position,n_comp + 1);
    }


    template <int dim>
    void
    AsciiDataLayered<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        Utilities::AsciiDataLayered<dim>::declare_parameters(prm,
                                                             "$ASPECT_SOURCE_DIR/data/initial-composition/ascii-data/test/",
                                                             "initial_composition_top_mantle_box_3d.txt",
                                                             "Ascii data layered");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiDataLayered<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        Utilities::AsciiDataLayered<dim>::parse_parameters(prm,
                                                           "Ascii data layered");
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
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(AsciiDataLayered,
                                              "ascii data layered",
                                              "Implementation of a model in which the initial "
                                              "composition is derived from files containing data "
                                              "in ascii format. Each file defines a surface on which "
                                              "compositional fields are defined. "
                                              "Between the surfaces, the fields can be chosen to be "
                                              "constant (with a value defined by the nearest shallower "
                                              "surface), or linearly interpolated between surfaces. "
                                              "Note the required format of the input ascii data file: "
                                              "The first lines may contain any number of comments "
                                              "if they begin with `#', but one of these lines needs to "
                                              "contain the number of grid points in each dimension as "
                                              "for example `# POINTS: 3 3'. "
                                              "The order of the data columns "
                                              "has to be `x', `y', `composition1', `composition2' etc. "
                                              "in a 2d model and `x', `y', `z', `composition1', "
                                              "`composition2' etc. in a 3d model; i.e. "
                                              "the columns before the compositional field always contains "
                                              "the position of the surface along the vertical direction. "
                                              "The first column needs to ascend first, "
                                              "followed by the second in order to assign the correct data "
                                              "to the prescribed coordinates. If you use a spherical model, "
                                              "then the assumed grid changes. `x' will be replaced by the "
                                              "azimuth angle and `y' (if 3d) by the polar angle measured "
                                              "positive from the north pole. The last column will be the "
                                              "distance of the point from the origin "
                                              "(i.e. radial position). The grid in this case will be "
                                              "a latitude-longitude grid. Note that the order "
                                              "of spherical coordinates in 3d is `phi', `theta', `r', `T'"
                                              "and not `theta', `phi', `r', `T' as this is "
                                              "more consistent with other ASPECT plugins. Outside of the "
                                              "region defined by the grid, the plugin will use the value "
                                              "at the edge of the region.")
  }
}
