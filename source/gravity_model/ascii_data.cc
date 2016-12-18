/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
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
    {}


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      this->initialize(6,this->get_mpi_communicator());
    }


    template <int dim>
    Tensor<1,dim>
    AsciiData<dim>::
    gravity_vector (const Point<dim> &position) const
    {
      const double depth = this->get_geometry_model().depth(position);
      const double magnitude = this->get_data_component(Point<1>(depth),3);
      Point<dim> vertical_direction;

      // in dependence of what the geometry model is, gravity points in a different direction
      if (const GeometryModel::SphericalShell<dim> *geometry = dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()))
        return - magnitude * position/position.norm();
      else if (const GeometryModel::Chunk<dim> *gm = dynamic_cast<const GeometryModel::Chunk<dim>*> (&this->get_geometry_model()))
        return - magnitude * position/position.norm();
      else if (const GeometryModel::EllipsoidalChunk<dim> *gm = dynamic_cast<const GeometryModel::EllipsoidalChunk<dim>*> (&this->get_geometry_model()))
        return - magnitude * position/position.norm();
      else if (const GeometryModel::Sphere<dim> *gm = dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()))
        return - magnitude * position/position.norm();
      else if (const GeometryModel::Box<dim> *gm = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model()))
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
      prm.enter_subsection("Adiabatic conditions model");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/adiabatic-conditions/ascii-data/test/",
                                                          "box_2d.txt");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Adiabatic conditions model");
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
                                  "Gravity is read from an ascii data file "
                                  "that contains a vertical or radial profile. "
                                  "File structure : ... TBD")
  }
}
