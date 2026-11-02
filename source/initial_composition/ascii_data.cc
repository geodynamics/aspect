// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2023 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/global.h>
#include <aspect/initial_composition/ascii_data.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
      = default;


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      ascii_data_initial->initialize(this->n_compositional_fields());
    }


    template <int dim>
    double
    AsciiData<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int n_comp) const
    {
      return ascii_data_initial->get_data_component(position,n_comp);
    }


    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        Utilities::AsciiDataInitial<dim>::declare_parameters(prm,
                                                             "$ASPECT_SOURCE_DIR/data/initial-composition/ascii-data/test/",
                                                             "box_2d.txt");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        ascii_data_initial = std::make_unique<Utilities::AsciiDataInitial<dim>>();
        ascii_data_initial->initialize_simulator(this->get_simulator());
        ascii_data_initial->parse_parameters(prm);
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
