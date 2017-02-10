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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/global.h>
#include <aspect/adiabatic_conditions/ascii_data.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>


namespace aspect
{
  namespace AdiabaticConditions
  {
    template <int dim>
    AsciiData<dim>::AsciiData()
      :
      initialized(false)
    {}



    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      if (initialized)
        return;

      this->initialize(7,this->get_mpi_communicator());

      initialized = true;
    }



    template <int dim>
    bool
    AsciiData<dim>::is_initialized() const
    {
      return initialized;
    }



    template <int dim>
    double AsciiData<dim>::pressure (const Point<dim> &p) const
    {
      const double depth = this->get_geometry_model().depth(p);
      return this->get_data_component(Point<1>(depth),0);
    }



    template <int dim>
    double AsciiData<dim>::temperature (const Point<dim> &p) const
    {
      const double depth = this->get_geometry_model().depth(p);
      return this->get_data_component(Point<1>(depth),1);
    }



    template <int dim>
    double AsciiData<dim>::density (const Point<dim> &p) const
    {
      const double depth = this->get_geometry_model().depth(p);
      return this->get_data_component(Point<1>(depth),2);
    }



    template <int dim>
    double AsciiData<dim>::density_derivative (const Point<dim> &p) const
    {
      const double depth = this->get_geometry_model().depth(p);
      const double eps = std::sqrt(std::numeric_limits<double>::epsilon()) * this->get_geometry_model().maximal_depth();
      return (this->get_data_component(Point<1>(depth+eps),2)
              -
              this->get_data_component(Point<1>(depth),2))
             /
             eps;
    }



    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Adiabatic conditions model");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/tests/adiabatic-conditions/ascii-data/test/",
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
  namespace AdiabaticConditions
  {
    ASPECT_REGISTER_ADIABATIC_CONDITIONS_MODEL(AsciiData,
                                               "ascii data",
                                               "A model in which the adiabatic profile is "
                                               "read from a file that describes the reference "
                                               "state. Note the required format of the "
                                               "input data: The first lines may contain any number of comments "
                                               "if they begin with '#', but one of these lines needs to "
                                               "contain the number of points in the reference state as "
                                               "for example '# POINTS: 3'. "
                                               "The order of the data columns has to be 'depth (m)', "
                                               "'pressure (Pa)', 'temperature (K)', 'density (kg/m^3)', "
                                               "'gravity (m/s^2)', 'thermal expansivity (1/K)', "
                                               "'specific heat (J/K/kg)', and 'compressibility (1/Pa)'. "
                                               "For incompressible models the 'compressibility' column will "
                                               "not be used, but needs to be present in the file."
                                               "Note that the data in the file need to be sorted in order "
                                               "of increasing depth from 0 to the maximal depth in the model "
                                               "domain. Points in the model that are outside of the provided "
                                               "depth range will be assigned the maximum or minimum depth values, "
                                               "respectively. Points to do not need to be equidistant, "
                                               "but the computation of properties is optimized in speed, "
                                               "if they are.")
  }
}
