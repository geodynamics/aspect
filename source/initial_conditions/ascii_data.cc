/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/initial_conditions/ascii_data.h>

#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>

namespace aspect
{
  namespace InitialConditions
  {
    namespace internal
    {
      template <int dim>
      AsciiDataInitial<dim>::AsciiDataInitial ()
      {}


      template <int dim>
      void
      AsciiDataInitial<dim>::initialize (const unsigned int components)
      {
        lookup.reset(new Utilities::AsciiDataLookup<dim> (components,
                                                          Utilities::AsciiDataBase<dim>::scale_factor));

        const std::string filename = Utilities::AsciiDataBase<dim>::data_directory
                                     + Utilities::AsciiDataBase<dim>::data_file_name;

        this->get_pcout() << std::endl << "   Loading Ascii data initial file "
                          << filename << "." << std::endl << std::endl;
        lookup->load_file(filename);
      }

      template <int dim>
      double
      AsciiDataInitial<dim>::
      get_data_component (const Point<dim>                    &position,
                          const unsigned int                   component) const
      {
        Point<dim> internal_position = position;

        if (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != 0)
          {
            const std_cxx11::array<double,dim> spherical_position =
              ::aspect::Utilities::spherical_coordinates(position);

            for (unsigned int i = 0; i < dim; i++)
              internal_position[i] = spherical_position[i];
          }
        return lookup->get_data(internal_position,component);
      }
    }

    template <int dim>
    AsciiData<dim>::AsciiData ()
    {}


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      internal::AsciiDataInitial<dim>::initialize(1);
    }


    template <int dim>
    double
    AsciiData<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      return internal::AsciiDataInitial<dim>::get_data_component(position,0);
    }


    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-conditions/ascii-data/test/",
                                                          "box_2d.txt");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
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
  namespace InitialConditions
  {
    namespace internal
    {
      template class AsciiDataInitial<2>;
      template class AsciiDataInitial<3>;
    }

    ASPECT_REGISTER_INITIAL_CONDITIONS(AsciiData,
                                       "ascii data",
                                       "Implementation of a model in which the initial "
                                       "temperature is derived from files containing data "
                                       "in ascii format. Note the required format of the "
                                       "input data: The first lines may contain any number of comments"
                                       "if they begin with '#', but one of these lines needs to"
                                       "contain the number of grid points in each dimension as"
                                       "for example '# POINTS: 3 3'."
                                       "The order of the data columns "
                                       "has to be 'x', 'y', 'Temperature [K]' in a 2d model and "
                                       " 'x', 'y', 'z', 'Temperature [K]' in a 3d model, which means that "
                                       "there has to be a single column "
                                       "containing the temperature. "
                                       "Note that the data in the input "
                                       "files need to be sorted in a specific order: "
                                       "the first coordinate needs to ascend first, "
                                       "followed by the second and the third at last in order to "
                                       "assign the correct data to the prescribed coordinates."
                                       "If you use a spherical model, "
                                       "then the data will still be handled as cartesian,"
                                       "however the assumed grid changes. 'x' will be replaced by "
                                       "the radial distance of the point to the bottom of the model, "
                                       "'y' by the azimuth angle and 'z' by the polar angle measured "
                                       "positive from the north pole. The grid will be assumed to be "
                                       "a latitude-longitude grid. Note that the order "
                                       "of spherical coordinates is 'r', 'phi', 'theta' "
                                       "and not 'r', 'theta', 'phi', since this allows "
                                       "for dimension independent expressions. ")
  }
}
