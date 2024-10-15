/*
  Copyright (C) 2023 - 2024 by the authors of the ASPECT code.

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
#include <aspect/initial_composition/slab_model.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    void
    SlabModel<dim>::initialize ()
    {
      AssertThrow(this->introspection().compositional_name_exists("slabs"),
                  ExcMessage("The initial composition plugin `slab model' did not find a "
                             "compositional field called `slabs' to initialize. Please add a "
                             "compositional field with this name."));

      slab_index = this->introspection().compositional_index_for_name("slabs");

      // The input slabs are defined from the surface of the model
      surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      // The two columns correspond to slabs depth and thickness
      slab_boundary.initialize({surface_boundary_id}, 2);
    }


    template <int dim>
    double
    SlabModel<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int compositional_index) const
    {
      if (compositional_index != slab_index)
        return 0.0;

      // The first data column corresponds to the slab depth and the second column to the slab thickness.
      // 'Slab depth' stands for the depth of the upper surface of the slab, 'Slab thickness'
      // for the vertical distance between upper and lower surface.
      const double slab_depth     = slab_boundary.get_data_component(surface_boundary_id, position, 0);
      const double slab_thickness = slab_boundary.get_data_component(surface_boundary_id, position, 1);

      // Return 0.0 if there is no slab in this location in the data. No slab is encoded in
      // the data file as a slab thickness of 0.0 and/or a depth larger than the depth of
      // the domain.
      if (slab_thickness == 0.0 ||
          slab_depth > this->get_geometry_model().maximal_depth())
        return 0.0;

      // Return 1.0 if we are inside the depth range of the slab.
      // We use a semi-closed interval so that a slab thickness of 0.0
      // means no slab exists.
      const double depth = this->get_geometry_model().depth(position);
      if ((depth >= slab_depth) && (depth < slab_depth + slab_thickness))
        return 1.0;

      return 0.0;
    }


    template <int dim>
    const Utilities::AsciiDataBoundary<dim> &
    SlabModel<dim>::get_slab_boundary () const
    {
      return slab_boundary;
    }


    template <int dim>
    void
    SlabModel<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                              "$ASPECT_SOURCE_DIR/data/initial-composition/slab-model/",
                                                              "shell_3d.txt",
                                                              "Slab model",
                                                              /*time dependent = */ false);
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    SlabModel<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        slab_boundary.initialize_simulator (this->get_simulator());
        slab_boundary.parse_parameters(prm, "Slab model", /*time dependent = */ false);
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
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(SlabModel,
                                              "slab model",
                                              "An initial composition model that implements subducted "
                                              "slab geometries as a compositional field determined from "
                                              "an input file. The file defines the depth to the top of "
                                              "the slab and the slab thickness. "
                                              "The computed compositional value is 1 within the slabs "
                                              "and zero elsewhere. "
                                              "An example model that is included is Slab2 described in "
                                              "Hayes, G. P., Moore, G. L., Portner, D. E., Hearne, M., "
                                              "Flamme, H., Furtney, M., \\& Smoczyk, G. M. (2018). Slab2, "
                                              "a comprehensive subduction zone geometry model. Science, "
                                              "362(6410), 58-61. The script to convert the Slab2 model "
                                              "into an aspect input data file is available in the directory "
                                              "data/initial-composition/slab-model/. Please note that "
                                              "Slab2 and the example data file assume spherical geometry "
                                              "(latitude, longitude coordinates), however, that is not "
                                              "necessary for this plugin, data files in "
                                              "Cartesian coordinates will work with box geometries.")
  }
}
