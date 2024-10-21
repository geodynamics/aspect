/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/boundary_traction/ascii_data.h>
#include <aspect/global.h>

namespace aspect
{
  namespace BoundaryTraction
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
      = default;


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      unsigned int i=0;
      for (const auto &plugin : this->get_boundary_traction_manager().get_active_plugins())
        {
          if (plugin.get() == this)
            boundary_ids.insert(this->get_boundary_traction_manager().get_active_plugin_boundary_indicators()[i]);

          ++i;
        }

      AssertThrow(boundary_ids.empty() == false,
                  ExcMessage("Did not find the boundary indicator for the traction ascii data plugin."));

      Utilities::AsciiDataBoundary<dim>::initialize(boundary_ids,
                                                    1);
    }


    template <int dim>
    Tensor<1,dim>
    AsciiData<dim>::
    boundary_traction (const types::boundary_id boundary_indicator,
                       const Point<dim> &position,
                       const Tensor<1,dim> &normal_vector) const
    {
      const double pressure = Utilities::AsciiDataBoundary<dim>::get_data_component(boundary_indicator,
                                                                                    position,
                                                                                    0);
      return -pressure * normal_vector;
    }


    template <int dim>
    void
    AsciiData<dim>::update()
    {
      Interface<dim>::update ();
      Utilities::AsciiDataBoundary<dim>::update();
    }


    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary traction model");
      {
        Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                              "$ASPECT_SOURCE_DIR/data/boundary-traction/ascii-data/test/",
                                                              "box_2d_%s.%d.txt");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary traction model");
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
  namespace BoundaryTraction
  {
    ASPECT_REGISTER_BOUNDARY_TRACTION_MODEL(AsciiData,
                                            "ascii data",
                                            "Implementation of a model in which the boundary "
                                            "traction is derived from files containing pressure data "
                                            "in ascii format. The pressure given in the data file is "
                                            "applied as traction normal to the surface of a given boundary, "
                                            "pointing inwards. Note the required format of the "
                                            "input data: The first lines may contain any number of comments "
                                            "if they begin with `#', but one of these lines needs to "
                                            "contain the number of grid points in each dimension as "
                                            "for example `# POINTS: 3 3'. "
                                            "The order of the data columns "
                                            "has to be `x', `Pressure [Pa]' in a 2d model and "
                                            " `x', `y', `Pressure [Pa]' in a 3d model, which means that "
                                            "there has to be a single column "
                                            "containing the pressure. "
                                            "Note that the data in the input "
                                            "files need to be sorted in a specific order: "
                                            "the first coordinate needs to ascend first, "
                                            "followed by the second in order to "
                                            "assign the correct data to the prescribed coordinates. "
                                            "If you use a spherical model, "
                                            "then the data will still be handled as Cartesian, "
                                            "however the assumed grid changes. `x' will be replaced by "
                                            "the radial distance of the point to the bottom of the model, "
                                            "`y' by the azimuth angle and `z' by the polar angle measured "
                                            "positive from the north pole. The grid will be assumed to be "
                                            "a latitude-longitude grid. Note that the order "
                                            "of spherical coordinates is `r', `phi', `theta' "
                                            "and not `r', `theta', `phi', since this allows "
                                            "for dimension independent expressions.")
  }
}
