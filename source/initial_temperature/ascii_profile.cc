// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2022 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/global.h>
#include <aspect/initial_temperature/ascii_profile.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>


namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    AsciiProfile<dim>::AsciiProfile ()
      = default;


    template <int dim>
    void
    AsciiProfile<dim>::initialize ()
    {
      this->initialize(this->get_mpi_communicator());
      temperature_index = this->get_column_index_from_name("temperature");
    }



    template <int dim>
    double AsciiProfile<dim>::initial_temperature (const Point<dim> &p) const
    {
      const double depth = this->get_geometry_model().depth(p);
      return this->get_data_component(Point<1>(depth),temperature_index);
    }


    template <int dim>
    void
    AsciiProfile<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-profile/tests/",
                                                          "simple_test.txt",
                                                          "Ascii profile");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiProfile<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm,
                                                        "Ascii profile");
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(AsciiProfile,
                                              "ascii profile",
                                              "Implementation of a model in which the initial temperature is "
                                              "read from a file that provides these values as a function of "
                                              "depth. Note the required format of the input data: The first "
                                              "lines may contain any number of comments if they begin with "
                                              "`#', but one of these lines needs to contain the number of "
                                              "points in the temperature profile, for example `# POINTS: 10'. "
                                              "Following the comment lines, there has to be a single line "
                                              "containing the names of all data columns, separated by arbitrarily "
                                              "many spaces. Column names are not allowed to contain spaces. "
                                              "The file can contain unnecessary columns, but for this plugin it "
                                              "needs to at least provide columns named `depth' and"
                                              "`temperature'."
                                              "Note that the data lines in the file need to be sorted in order "
                                              "of increasing depth from 0 to the maximal depth in the model "
                                              "domain. Points in the model that are outside of the provided "
                                              "depth range will be assigned the maximum or minimum depth values, "
                                              "respectively. Points do not need to be equidistant, "
                                              "but the computation of properties is optimized in speed "
                                              "if they are.")
  }
}
