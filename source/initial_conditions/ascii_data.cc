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


namespace aspect
{
  namespace InitialConditions
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
      :
    lookup()
    {}


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      lookup.reset(new Utilities::AsciiDataLookup<dim,dim> (this->get_geometry_model(),
                                                           1,
                                                           Utilities::AsciiDataBase<dim>::scale_factor));

      lookup->screen_output(this->get_pcout());

      const std::string filename = Utilities::AsciiDataBase<dim>::data_directory
                                 + Utilities::AsciiDataBase<dim>::data_file_name;

      this->get_pcout() << std::endl << "   Loading Ascii data initial file "
          << filename << "." << std::endl << std::endl;

      // We load the file twice, this is because AsciiDataLookup also performs time
      // interpolation for the boundary conditions, which is not necessary here.
      lookup->load_file(filename);
      lookup->load_file(filename);
    }


    template <int dim>
    double
    AsciiData<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      return lookup->get_data(position,0,0);
    }


    template <int dim>
     void
     AsciiData<dim>::declare_parameters (ParameterHandler &prm)
     {
       prm.enter_subsection("Initial conditions");
       {
         Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                           "$ASPECT_SOURCE_DIR/data/initial-conditions/ascii-data/test/",
                                                           "box_2d.csv");
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
