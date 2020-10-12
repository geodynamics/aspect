/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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
#include <aspect/gravity_model/interface.h>
#include <aspect/mesh_deformation/ascii_data.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/parameter_handler.h>
#include <array>


namespace aspect
{
  namespace MeshDeformation
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

      std::set<types::boundary_id> surface_boundary_set;
      surface_boundary_set.insert(surface_boundary_id);

      Utilities::AsciiDataBoundary<dim>::initialize(surface_boundary_set,
                                                    1);
    }



    template <int dim>
    Tensor<1,dim>
    AsciiData<dim>::compute_initial_deformation_on_boundary(const types::boundary_id boundary_indicator,
                                                            const Point<dim> &position) const
    {
      const double topo = Utilities::AsciiDataBoundary<dim>::get_data_component(boundary_indicator,
                                                                                position,
                                                                                0);

      const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(position);

      Tensor<1,dim> topography_direction;
      if (gravity.norm() > 0.0)
        topography_direction = -gravity / gravity.norm();

      return topo * topography_direction;
    }



    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh deformation");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/geometry-model/initial-topography-model/ascii-data/test/",
                                                          "box_3d_%s.0.txt");
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh deformation");
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm);
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(AsciiData,
                                           "ascii data",
                                           "Implementation of a model in which the initial mesh "
                                           "deformation (initial topography) is "
                                           "derived from a file containing data "
                                           "in ascii format. The following geometry models "
                                           "are currently supported: box, chunk, spherical. "
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
                                           "`x' will be replaced by the azimuth angle in radians "
                                           " and `y' by the polar angle in radians measured "
                                           "positive from the north pole. The grid will be assumed to be "
                                           "a longitude-colatitude grid. Note that the order "
                                           "of spherical coordinates is `phi', `theta' "
                                           "and not `theta', `phi', since this allows "
                                           "for dimension independent expressions.")
  }
}
