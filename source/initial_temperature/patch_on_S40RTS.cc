/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/initial_temperature/patch_on_S40RTS.h>
#include <aspect/utilities.h>

#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    PatchOnS40RTS<dim>::PatchOnS40RTS ()
      = default;


    template <int dim>
    void
    PatchOnS40RTS<dim>::initialize ()
    {
      this->Utilities::AsciiDataInitial<dim>::initialize(1);
    }


    template <int dim>
    double
    PatchOnS40RTS<dim>::
    ascii_grid_vs (const Point<dim> &position) const
    {
      const double vs_perturbation = Utilities::AsciiDataInitial<dim>::get_data_component(position,0);
      return vs_perturbation;
    }


    template <int dim>
    double
    PatchOnS40RTS<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      const double depth = this->get_geometry_model().depth(position);

      double vs_perturbation;
      if (depth <= max_grid_depth - smoothing_length_scale)
        {
          vs_perturbation = ascii_grid_vs(position);
        }
      //add smoothing between the two models
      else if (depth > max_grid_depth - smoothing_length_scale && depth < max_grid_depth)
        {
          const double scale_factor = (depth-(max_grid_depth-smoothing_length_scale))/smoothing_length_scale;
          vs_perturbation = s40rts.get_Vs(position)*(scale_factor) + ascii_grid_vs(position)*(1.0-scale_factor);
        }
      else
        {
          vs_perturbation = s40rts.get_Vs(position);
        }

      // use either the user-input reference temperature as background temperature
      // (incompressible model) or the adiabatic temperature profile (compressible model)
      const double background_temperature = this->get_material_model().is_compressible() ?
                                            this->get_adiabatic_conditions().temperature(position) :
                                            s40rts.reference_temperature;

      // get the Vs to density conversion
      double vs_to_density = 0.0;
      if (s40rts.vs_to_density_method == s40rts.file)
        vs_to_density = s40rts.profile.get_data_component(Point<1>(depth), s40rts.vs_to_density_index);
      else if (s40rts.vs_to_density_method == s40rts.constant)
        vs_to_density = s40rts.vs_to_density_constant;
      else
        // we shouldn't get here but instead should already have been
        // kicked out when declaring the parameter, as the pattern won't match ("file|constant")
        Assert(false, ExcMessage("Unknown method for vs to density scaling."));

      // scale the perturbation in seismic velocity into a density perturbation
      // vs_to_density is read in from input file
      const double density_perturbation = vs_to_density * vs_perturbation;
      double temperature_perturbation;
      if (depth > no_perturbation_depth_patch)
        // scale the density perturbation into a temperature perturbation
        temperature_perturbation =  -1./s40rts.thermal_alpha * density_perturbation;
      else
        // set heterogeneity to zero down to a specified depth
        temperature_perturbation = 0.0;

      return background_temperature + temperature_perturbation;
    }

    template <int dim>
    void
    PatchOnS40RTS<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection ("Patch on S40RTS");
        {
          prm.declare_entry ("Maximum grid depth", "700000.0",
                             Patterns::Double (0.),
                             "The maximum depth of the Vs ascii grid. The model will read in  "
                             "Vs from S40RTS below this depth.");
          prm.declare_entry ("Smoothing length scale", "200000.0",
                             Patterns::Double (0.),
                             "The depth range (above maximum grid depth) over which to smooth. "
                             "The boundary is smoothed using a depth weighted combination of Vs "
                             "values from the ascii grid and S40RTS at each point in the region of smoothing.");
          prm.declare_entry ("Remove temperature heterogeneity down to specified depth",
                             boost::lexical_cast<std::string>(std::numeric_limits<double>::lowest()),
                             Patterns::Double (),
                             "This will set the heterogeneity prescribed by the Vs ascii grid and S40RTS to zero "
                             "down to the specified depth (in meters). Note that your resolution has "
                             "to be adequate to capture this cutoff. For example if you specify a depth "
                             "of 660 km, but your closest spherical depth layers are only at 500 km and "
                             "750 km (due to a coarse resolution) it will only zero out heterogeneities "
                             "down to 500 km. Similar caution has to be taken when using adaptive meshing.");

          Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                            "$ASPECT_SOURCE_DIR/data/initial-temperature/patch-on-S40RTS/test/",
                                                            "upper_shell_3d.txt");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    PatchOnS40RTS<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial temperature model");
      {
        prm.enter_subsection ("Patch on S40RTS");
        {
          max_grid_depth           = prm.get_double ("Maximum grid depth");
          smoothing_length_scale   = prm.get_double ("Smoothing length scale");
          no_perturbation_depth_patch   = prm.get_double ("Remove temperature heterogeneity down to specified depth");

          Utilities::AsciiDataBase<dim>::parse_parameters(prm);

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      s40rts.initialize_simulator (this->get_simulator());
      s40rts.parse_parameters(prm);
      s40rts.initialize();

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(PatchOnS40RTS,
                                              "patch on S40RTS",
                                              "Implementation of a model in which the initial "
                                              "temperature is derived from a file containing shear wave velocity perturbations "
                                              "in ascii format (e.g. a high resolution upper mantle tomography) "
                                              "combined with S40RTS. Note the required format of the "
                                              "input ascii input data: The first lines may contain any number of comments "
                                              "if they begin with '#', but one of these lines needs to "
                                              "contain the number of grid points in each dimension as "
                                              "for example '# POINTS: 3 3 3'. "
                                              "The order of the data columns has to be "
                                              " `x', `y', `z', 'Vs Perturbation' in a 3d model, which means that "
                                              "there has to be a single column "
                                              "containing the temperature. "
                                              "Note that the data in the input "
                                              "files need to be sorted in a specific order: "
                                              "the first coordinate needs to ascend first, "
                                              "followed by the second and the third at last in order to "
                                              "assign the correct data to the prescribed coordinates. "
                                              "In the spherical model data will be handled as Cartesian, "
                                              "however, `x' will be replaced by "
                                              "the radial distance of the point to the bottom of the model, "
                                              "`y' by the azimuth angle and `z' by the polar angle measured "
                                              "positive from the north pole. The grid will be assumed to be "
                                              "a latitude-longitude grid. Note that the order "
                                              "of spherical coordinates is `r', `phi', `theta' "
                                              "and not `r', `theta', `phi', since this allows "
                                              "for dimension independent expressions. "
                                              "See S40RTS documentation for details on input parameters in the "
                                              "S40RTS perturbation subsection. "
                                              "The boundary between the two tomography models is smoothed using a depth weighted "
                                              "combination of Vs values within the region of smoothing. ")
  }
}
