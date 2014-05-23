/*
  Copyright (C) 2012 by the authors of the ASPECT code.

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
/*  $Id: two_plates.cc Nov 19, 2013 Katrina M Arredondo Modified from boussinesq_incomp.cc $  */


#include <aspect/initial_conditions/two_plates.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <cmath>

namespace aspect
{
  namespace InitialConditions
  {
    template <int dim>
    double
    Two_plates<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // convert input ages to seconds
      // const double ovplate_age_ma = 40;
      const double trench_plate_age_sec = plate_age_ma * 1e6 * 365 * 24 * 60 * 60;
	const double ovplate_age_sec = ovplate_age_ma * 1e6 * 365 * 24 * 60 * 60;
      // const double T_surface = this->boundary_temperature->minimal_temperature(this->get_fixed_temperature_boundary_indicators()); // Not giving correct surface temp.
        const double depth = this->geometry_model->depth(position);
        const double Re =  6371e3; //meters
	//const double depth = Re - position.norm();      
	const double R1 = dynamic_cast<const GeometryModel::SphericalShell<dim>&> (*this->geometry_model).outer_radius(); 
	const double y = Re * atan(position(0)/position(1));
	//const double trench_distance = R1 / 2;
	const double plate_age_sec = (trench_plate_age_sec / trench_distance) * y;
	const double wedge_edge = trench_distance + depth/tan(wedge_angle * numbers::PI / 180);
	const double wedge_z = (y - trench_distance) * tan(wedge_angle * numbers::PI / 180);
	const double wedge_y_max = trench_distance + wedge_depth/tan(wedge_angle * numbers::PI / 180);
	double temp;

	if (y <= trench_distance)   
      		temp = plate_age_sec > 0.0 ?
                	          erf(depth/(2e0 * std::pow(K_0 * plate_age_sec,0.5))) * (T1 - T_surface) + T_surface
                        	  : T1;
	else if ((y > trench_distance) && (depth >= wedge_z) && (y <= wedge_edge) && (y <= wedge_y_max)){
		temp = plate_age_sec > 0.0 ?
                                  erf(depth/(2e0 * std::pow(K_0 * plate_age_sec,0.5))) * (T1 - T_surface) + T_surface
                                  : T1;
		// std::cout << "Depth: " << depth << " wedge_depth: " << wedge_depth << " wedge_edge: " << wedge_edge <<  std::endl; 
	}
	else
	 	temp = ovplate_age_sec > 0.0 ?
                          	erf(depth/(2e0 * std::pow(K_0 * ovplate_age_sec,0.5))) * (T1 - T_surface) + T_surface
                                  : T1;
        return temp;
    }


    template <int dim>
    void
    Two_plates<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Two_plates");
        {
          prm.declare_entry ("Overriding Plate Age Ma", "80",
                             Patterns::Double (),
                             "Age of the single, constant-age overriding plate. Units: millions of years.");
          prm.declare_entry ("Subducting Plate Age at the Trench Ma", "80",
                             Patterns::Double (),
                             "Age of the subducting plate at the trench, with a ridge at the grid edge. Units: millions of years.");
          prm.declare_entry ("Trench Location", "3000",
                             Patterns::Double (),
                             "Location of the trench, measured in meters from the left? edge of the grid. Units: meters.");
          prm.declare_entry ("Surface Temperature", "275.15",
                             Patterns::Double (0),
                             "Mantle temperature at the surface. Units: $K$.");	  
          prm.declare_entry ("Mantle Temperature Below Plate", "1673.15",
                             Patterns::Double (0),
                             "Mantle temperature directly below the plate. Units: $K$.");
          prm.declare_entry ("Reference Thermal diffusivity", "1e-6",
                             Patterns::Double (0),
                             "The value of the thermal diffusivity $K$. "
                             "Units: $m^2/s$.");
	  prm.declare_entry ("Wedge Angle", "30",
                             Patterns::Double (0),
                             "Angle of wedge of subducting plate under potential weakzone. "
                             "Units: Degrees.");
	  prm.declare_entry ("Wedge Depth", "80e3",
                             Patterns::Double (0),
                             "Depth of wedge of subducting plate under potential weakzone. "
                             "Units: meters.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Two_plates<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Two_plates");
        {
	  ovplate_age_ma    = prm.get_double ("Overriding Plate Age Ma"); 
          plate_age_ma      = prm.get_double ("Subducting Plate Age at the Trench Ma");
	  trench_distance   = prm.get_double ("Trench Location");
	  T_surface	    = prm.get_double ("Surface Temperature");
          T1                = prm.get_double ("Mantle Temperature Below Plate");
          K_0               = prm.get_double ("Reference Thermal diffusivity");
	  wedge_angle       = prm.get_double ("Wedge Angle");
	  wedge_depth	    = prm.get_double ("Wedge Depth");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(Two_plates,
                                       "two_plates",
                                       "Half Space Cooling model with a mantle of constant  "
                                       "temperature. ")
  }
}
