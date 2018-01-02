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


#include <aspect/material_model/steinberger_rheology.h>
#include <aspect/adiabatic_conditions/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <iostream>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    SteinbergerRheology<dim>::initialize()
    {
      DamageRheology<dim>::initialize();

      lateral_viscosity_lookup.reset(new internal::LateralViscosityLookup(datadirectory+lateral_viscosity_file_name, this->get_mpi_communicator()));
      radial_viscosity_lookup.reset(new internal::RadialViscosityLookup(datadirectory+radial_viscosity_file_name, this->get_mpi_communicator()));
    }

    template <int dim>
    double
    SteinbergerRheology<dim>::
    viscosity (const double temperature,
               const double,
               const std::vector<double> &,
               const SymmetricTensor<2,dim> &,
               const Point<dim> &position) const
    {
      const double depth = this->get_geometry_model().depth(position);
      const double adiabatic_temperature = this->get_adiabatic_conditions().temperature(position);

      const double delta_temperature = temperature-adiabatic_temperature;

      // For an explanation on this formula see the Steinberger & Calderwood 2006 paper
      const double vis_lateral_exp = -1.0*lateral_viscosity_lookup->lateral_viscosity(depth)*delta_temperature/(temperature*adiabatic_temperature);
      // Limit the lateral viscosity variation to a reasonable interval
      const double vis_lateral = std::max(std::min(std::exp(vis_lateral_exp),this->max_temperature_dependence_of_eta),1/this->max_temperature_dependence_of_eta);

      const double vis_radial = radial_viscosity_lookup->radial_viscosity(depth);

      return std::max(std::min(vis_lateral * vis_radial,this->max_eta),this->min_eta);
    }

    template <int dim>
    void
    SteinbergerRheology<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      DamageRheology<dim>::evaluate(in,out);

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          //Use the adiabatic pressure instead of the real one, because of oscillations
          const double pressure = (this->get_adiabatic_conditions().is_initialized())
                                  ?
                                  this->get_adiabatic_conditions().pressure(in.position[i])
                                  :
                                  in.pressure[i];

          /* We are only asked to give viscosities if strain_rate.size() > 0
           * and we can only calculate it if adiabatic_conditions are available.
           * Note that the used viscosity formulation needs the not
           * corrected temperatures in case we compare it to the lateral
           * temperature average.
           */
          if (this->get_adiabatic_conditions().is_initialized() && in.strain_rate.size())
            {
              out.viscosities[i] = viscosity (in.temperature[i], pressure, in.composition[i], in.strain_rate[i], in.position[i]);
            }
        }
    }


    template <int dim>
    void
    SteinbergerRheology<dim>::declare_parameters (ParameterHandler &prm)
    {
      DamageRheology<dim>::declare_parameters(prm);

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Damage rheology model");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/steinberger/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Radial viscosity file name", "radial-visc.txt",
                             Patterns::Anything (),
                             "The file name of the radial viscosity data. ");
          prm.declare_entry ("Lateral viscosity file name", "temp-viscosity-prefactor.txt",
                             Patterns::Anything (),
                             "The file name of the lateral viscosity data. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SteinbergerRheology<dim>::parse_parameters (ParameterHandler &prm)
    {
      DamageRheology<dim>::parse_parameters(prm);

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Damage rheology model");
        {
          // parameters for reading in tables with material properties
          datadirectory        = prm.get ("Data directory");
          {
            const std::string      subst_text = "$ASPECT_SOURCE_DIR";
            std::string::size_type position;
            while (position = datadirectory.find (subst_text),  position!=std::string::npos)
              datadirectory.replace (datadirectory.begin()+position,
                                     datadirectory.begin()+position+subst_text.size(),
                                     ASPECT_SOURCE_DIR);
          }

          radial_viscosity_file_name   = prm.get ("Radial viscosity file name");
          lateral_viscosity_file_name  = prm.get ("Lateral viscosity file name");
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
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SteinbergerRheology,
                                   "steinberger rheology",
                                   "A material model that behaves in the same way as "
                                   "the simple material model, but includes compositional "
                                   "fields that stand for average grain sizes of a mineral "
                                   "phase and source terms for them that determine the grain "
                                   "size evolution in dependence of the strain rate, "
                                   "temperature, phase transitions, and the creep regime. "
                                   "In the diffusion creep regime, the viscosity depends "
                                   "on this grain size."
                                   "We use the grain size evolution laws described in Behn "
                                   "et al., 2009. Implications of grain size evolution on the "
                                   "seismic structure of the oceanic upper mantle, "
                                   "Earth Planet. Sci. Letters, 282, 178–189.")
  }
}
