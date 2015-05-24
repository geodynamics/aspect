/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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

namespace aspect
{
  // The following are a set of global constants which may be used by ASPECT:

  namespace constants
  {
    // Number of seconds in a year [s]
    const double year_in_seconds  = 60*60*24*365.2425;

    // Zero degrees Celsius to Kelvin [K]
    const double celsius_to_kelvin = 273.15;

    // Gas constant (also known as R, from NIST) [J K^-1 mol^-1]
    const double gas_constant = 8.3144621;
    // Avogadro's constant (NIST) [mol^-1]
    const double avogadro = 6.02214129e23;
    // Gravitational constant (NIST) [m^3 kg^-1 s^-2]
    const double big_g = 6.67384e-11;

    namespace earth
    {
      // From Yoder (1995), masses are taken
      namespace masses
      {
        // Planet mass [kg]
        const double planet = 5.9736e24;
        // Mass of the whole core [kg]
        const double core = 1.932e24;
        // Mass of the mantle [kg]
        const double mantle = 4.043e24;
      }


      // Values taken from the IASP91 model
      namespace iasp91_radii
      {
        // Inner core radius [m], equivalent of 5150 km depth
        const double inner_core = 1.2171e6;
        // Inner core radius [m], equivalent of 2889 km depth
        const double core = 3.482e6;
        // Lower mantle radius [m], equivalent of 660 km depth
        const double lower_mantle = 5.711e6;
        // Radius [m], equivalent of 5150 km depth
        const double planet = 6.371e6;
      }

      // Gravity values taken from the PREM (Dziewonski and Anderson, 1981)
      namespace prem_gravity
      {
        // Inner core boundary gravity [ms^-2]
        const double icb = 4.4002;
        // Core-mantle boundary gravity [ms^-2]
        const double cmb = 10.6823;
        // Upper-lower mantle boundary gravity [ms^-2]
        const double ulmb = 10.0143;
        // Surface gravity [ms^-2]
        const double surface = 9.8156;
      }

      // NIST "Standard gravity" (average gravitational acceleration at surface [ms^-2]
      const double surface_gravity = 9.80665;
    }

    // Constants for Mars
    namespace mars
    {
      namespace radii
      {
        // Radius from Seidermann et al., 2007 [m]
        const double planet = 3.3895e6;
        // Core radius from Rivoldini et al., 2011 [m]
        const double core = 1.794e6;
      }

      // Surface gravity from Lodders et al., 1998 [ms^-2]
      const double surface_gravity = 3.711;
    }
  }

  // Output parallel statistics [bool]
  const bool output_parallel_statistics = false;
}
