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
#include <aspect/compositional_initial_conditions/ascii_data.h>

#include <aspect/utilities.h>

namespace aspect
{
  namespace CompositionalInitialConditions
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
      :
    scale_factor(1.0),
    lookup()
    {}


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      std_cxx11::array<unsigned int,dim> boundary_dimensions;
      for (unsigned int i = 0; i < dim; i++)
        boundary_dimensions[i] = i;

      lookup.reset(new VelocityBoundaryConditions::internal::AsciiDataLookup<dim,dim>    (data_points,
                                                                                          this->get_geometry_model(),
                                                                                          this->n_compositional_fields(),
                                                                                          scale_factor));

      lookup->screen_output(this->get_pcout());

      this->get_pcout() << std::endl << "   Loading Ascii data initial file "
          << data_directory+data_file_name << "." << std::endl << std::endl;

      // We load the file twice, this is because AsciiDataLookup also performs time
      // interpolation for the boundary conditions, which is not necessary here.
      lookup->load_file(data_directory+data_file_name);
      lookup->load_file(data_directory+data_file_name);
    }


    template <int dim>
    double
    AsciiData<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int n_comp) const
    {
      return lookup->get_data(position,n_comp,0);
    }

    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Compositional initial conditions");
      {
        prm.enter_subsection ("Ascii data model");
        {
          prm.declare_entry ("Data directory",
                             "$ASPECT_SOURCE_DIR/data/compositional-initial-conditions/ascii-data/test/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the model data. This path "
                             "may either be absolute (if starting with a '/') or relative to "
                             "the current directory. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Data file name",
                             "box_2d_%s.%d.csv",
                             Patterns::Anything (),
                             "The file name of the material data. Provide file in format: "
                             "(Data file name).\\%s%d where \\\\%s is a string specifying "
                             "the boundary of the model according to the names of the boundary "
                             "indicators (of a box or a spherical shell).%d is any sprintf integer "
                             "qualifier, specifying the format of the current file number. ");
          prm.declare_entry ("Number of x grid points", "0",
                             Patterns::Integer (0),
                             "Number of grid points in x direction.");
          prm.declare_entry ("Number of y grid points", "0",
                             Patterns::Integer (0),
                             "Number of grid points in y direction.");
          prm.declare_entry ("Number of z grid points", "0",
                             Patterns::Integer (0),
                             "Number of grid points in z direction.");
          prm.declare_entry ("Scale factor", "1",
                             Patterns::Double (0),
                             "Scalar factor, which is applied to the boundary data. "
                             "You might want to use this to scale the data to a "
                             "reference model.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Compositional initial conditions");
      {
        prm.enter_subsection("Ascii data model");
        {
          // Get the path to the data files. If it contains a reference
          // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
          // as a #define
          data_directory    = prm.get ("Data directory");
          {
            const std::string      subst_text = "$ASPECT_SOURCE_DIR";
            std::string::size_type position;
            while (position = data_directory.find (subst_text),  position!=std::string::npos)
              data_directory.replace (data_directory.begin()+position,
                                      data_directory.begin()+position+subst_text.size(),
                                      ASPECT_SOURCE_DIR);
          }

          data_file_name    = prm.get ("Data file name");

          scale_factor      = prm.get_double ("Scale factor");
          data_points[0]    = prm.get_integer ("Number of x grid points");
          data_points[1]    = prm.get_integer ("Number of y grid points");
          data_points[2]    = prm.get_integer ("Number of z grid points");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace CompositionalInitialConditions
  {
    ASPECT_REGISTER_COMPOSITIONAL_INITIAL_CONDITIONS(AsciiData,
                                                     "ascii data",
                                                     "Implementation of a model in which the initial "
                                                     "composition is derived from files containing data "
                                                     "in ascii format. Note the required format of the "
                                                     "input data: The order of the columns "
                                                     "has to be 'x', 'y', 'composition_1', 'composition_2', "
                                                     "etc. in a 2d model and 'x', 'y', 'z', 'composition_1', "
                                                     "'composition_2', etc. in a 3d model, according "
                                                     "to the number of compositional fields, which means that "
                                                     "there has to be a single column "
                                                     "for every composition in the model."
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
