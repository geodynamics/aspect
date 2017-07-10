/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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

#include <aspect/gravity_model/ascii_data.h>

#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>

namespace aspect
{
  namespace GravityModel
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
      :
      gravity_index(numbers::invalid_unsigned_int)
    {}


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      this->initialize(this->get_mpi_communicator());
      gravity_index = this->get_column_index_from_name("gravity");
    }


    template <int dim>
    Tensor<1,dim>
    AsciiData<dim>::
    gravity_vector (const Point<dim> &position) const
    {
      const double depth = this->get_geometry_model().depth(position);
      const double magnitude = this->get_data_component(Point<1>(depth),gravity_index);

      // in dependence of what the geometry model is, gravity points in a different direction
      if (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model())
          || dynamic_cast<const GeometryModel::Chunk<dim>*> (&this->get_geometry_model())
          || dynamic_cast<const GeometryModel::EllipsoidalChunk<dim>*> (&this->get_geometry_model())
          || dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()))
        return - magnitude * position/position.norm();
      else if (dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model()))
        {
          Tensor<1,dim> g;
          g[dim-1] = -magnitude;
          return g;
        }
      else
        AssertThrow (false,
                     ExcMessage ("Not a valid geometry model for the gravity model"
                                 "ascii data."));
      return Tensor<1,dim>();
    }

    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      // we use the same file that is used for the adiabatic conditions model,
      // as it also contains gravity
      prm.enter_subsection("Gravity model");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/adiabatic-conditions/ascii-data/test/",
                                                          "");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
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
  namespace GravityModel
  {
    ASPECT_REGISTER_GRAVITY_MODEL(AsciiData,
                                  "ascii data",
                                  "Gravity is read from a file that describes the reference "
                                  "state. Note the required format of the "
                                  "input data: The first lines may contain any number of comments "
                                  "if they begin with '#', but one of these lines needs to "
                                  "contain the number of points in the reference state as "
                                  "for example '# POINTS: 3'. "
                                  "Following the comment lines there has to be a single line "
                                  "containing the names of all data columns, separated by arbitrarily "
                                  "many spaces. Column names are not allowed to contain spaces. "
                                  "The file can contain unnecessary columns, but for this plugin it "
                                  "needs to at least provide a column named `gravity'. "
                                  "Note that the data lines in the file need to be sorted in order "
                                  "of increasing depth from 0 to the maximal depth in the model "
                                  "domain. Points in the model that are outside of the provided "
                                  "depth range will be assigned the maximum or minimum depth values, "
                                  "respectively. Points do not need to be equidistant, "
                                  "but the computation of properties is optimized in speed "
                                  "if they are.")
  }
}
