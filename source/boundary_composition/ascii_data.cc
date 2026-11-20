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
#include <aspect/boundary_composition/ascii_data.h>

#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace BoundaryComposition
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
      = default;

    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      const unsigned int n_my_fixed_fields = this->get_boundary_composition_manager().get_fixed_compositional_fields_for_plugin("ascii data").size();
      const std::set<types::boundary_id> my_fixed_boundary_indicators = this->get_boundary_composition_manager().get_fixed_composition_boundaries_for_plugin("ascii data");
      Utilities::AsciiDataBoundary<dim>::initialize(my_fixed_boundary_indicators,
                                                    n_my_fixed_fields);
    }


    template <int dim>
    void
    AsciiData<dim>::update ()
    {
      Interface<dim>::update ();

      Utilities::AsciiDataBoundary<dim>::update();
    }


    template <int dim>
    double
    AsciiData<dim>::
    boundary_composition (const types::boundary_id boundary_indicator,
                          const Point<dim> &position,
                          const unsigned int compositional_field) const
    {
      // In case not all fields are fixed by the plugin,
      // figure out the right index of the given field.
      unsigned int field_id = compositional_field;
      if (this->get_boundary_composition_manager().boundaries_with_fixed_subset_of_fields_exist())
        {
          const std::set<unsigned int> my_fixed_fields = this->get_boundary_composition_manager().get_fixed_compositional_fields_for_plugin("ascii data");
          Assert (my_fixed_fields.find(compositional_field) != my_fixed_fields.end(),
                  ExcMessage ("Boundary composition was requested for field " +
                              Utilities::int_to_string(compositional_field) +
                              " on boundary " +
                              Utilities::int_to_string(boundary_indicator) +
                              " but this field is not prescribed by the `ascii data` plugin. "));
          field_id = std::distance(my_fixed_fields.begin(), my_fixed_fields.find(compositional_field));
        }

      return Utilities::AsciiDataBoundary<dim>::get_data_component(boundary_indicator,
                                                                   position,
                                                                   field_id);
    }

    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                              "$ASPECT_SOURCE_DIR/data/boundary-composition/ascii-data/test/",
                                                              "box_2d_%s.%d.txt");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        Utilities::AsciiDataBoundary<dim>::parse_parameters(prm);
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryComposition
  {
    ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(AsciiData,
                                               "ascii data",
                                               "Implementation of a model in which the boundary "
                                               "composition is derived from files containing data "
                                               "in ascii format. Note the required format of the "
                                               "input data: The first lines may contain any number of comments "
                                               "if they begin with `#', but one of these lines needs to "
                                               "contain the number of grid points in each dimension as "
                                               "for example `# POINTS: 3 3'. "
                                               "The order of the data columns "
                                               "has to be `x', `composition1', `composition2', "
                                               "etc. in a 2d model and "
                                               "`x', `y', `composition1', `composition2', "
                                               "etc., in a 3d model, with a column for each of "
                                               "the compositional fields that are "
                                               "prescribed by this plugin (on all boundaries). "
                                               "This means that there has to be a single column "
                                               "for every composition that is prescribed with `ascii data' "
                                               "in the model. Note that the data in the input "
                                               "files need to be sorted in a specific order: "
                                               "the first coordinate needs to ascend first, "
                                               "followed by the second in order to "
                                               "assign the correct data to the prescribed coordinates."
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
