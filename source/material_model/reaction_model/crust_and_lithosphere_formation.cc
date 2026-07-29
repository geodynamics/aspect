/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/reaction_model/crust_and_lithosphere_formation.h>
#include <aspect/utilities.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      template <int dim>
      void
      CrustLithosphereFormation<dim>::calculate_reaction_terms (const typename Interface<dim>::MaterialModelInputs  &in,
                                                                typename Interface<dim>::MaterialModelOutputs       &out) const
      {
        for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
          {
            const double depth = this->get_geometry_model().depth(in.position[i]);

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              {
                // Crust and lithosphere are only generated when material approaches the surface,
                // i.e. above the reaction depth and when the velocity is predominantly upwards.
                const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[i]).norm();
                Tensor<1,dim> vertically_up = -this->get_gravity_model().gravity_vector(in.position[i]);
                if (gravity_norm > 0.0)
                  vertically_up /= gravity_norm;

                // We only convert material to crust and lithosphere if its upwelling angle is above a certain value.
                // The upwelling angle is defined as the angle between the horizontal and the velocity vector.
                const double crustal_angle = std::sin(minimum_upwelling_angle_for_crust);
                const double lithospheric_angle = std::sin(minimum_upwelling_angle_for_lithosphere);

                const bool upward_flow_crust = in.velocity[i] * vertically_up > in.velocity[i].norm() * crustal_angle;
                const bool upward_flow_lithosphere = in.velocity[i] * vertically_up > in.velocity[i].norm() * lithospheric_angle;

                const bool within_crust = depth < crust_thickness && upward_flow_crust;
                const bool within_lithosphere = depth < crust_thickness + lithosphere_thickness && upward_flow_lithosphere;

                // In the crust, we convert every material to basalt.
                if (within_crust)
                  {
                    if (c == basalt_index)
                      out.reaction_terms[i][c] = 1. - in.composition[i][c];
                    else if (c == harzburgite_index)
                      out.reaction_terms[i][c] = - in.composition[i][c];
                  }
                // In the lithosphere, we convert pyrolite to harzburgite.
                // If basalt (recycled crust) reaches lithospheric depth, it is not affected.
                // TODO: Think about what should happen to basalt when it reaches depth of melting again.
                else if (within_lithosphere)
                  {
                    if (c == basalt_index)
                      out.reaction_terms[i][c] = 0.0;
                    else if (c == harzburgite_index)
                      {
                        // Only background mantle is converted, whereas basalt is not.

                        double harzburgite_change = std::numeric_limits<double>::signaling_NaN();
                        if (harzburgite_profile == linear)
                          harzburgite_change = (crust_thickness + lithosphere_thickness - depth) / lithosphere_thickness
                                               * (1.0 - in.composition[i][basalt_index]);

                        else if (harzburgite_profile == constant)
                          harzburgite_change = 1.0 - in.composition[i][basalt_index];

                        // If we already have more harzburgite than the change, we do not change it.
                        if (in.composition[i][c] >= harzburgite_change)
                          out.reaction_terms[i][c] = 0.0;
                        else
                          out.reaction_terms[i][c] = harzburgite_change - in.composition[i][c];
                      }
                  }
                // Do not change other reaction terms.
                // There might be other reactions computed in the material model
                // (including for the basalt and harzburgite fields).
              }
          }
      }



      template <int dim>
      void
      CrustLithosphereFormation<dim>::declare_parameters (ParameterHandler &prm)
      {
        // TODO: In the future, we could make these depths depend on the solidus
        // and the temperature. However, note that technically this would affect
        // the composition of the generated melt and residual as well.
        prm.declare_entry ("Crustal thickness", "7000",
                           Patterns::Double (),
                           "Thickness of the crustal layer generated "
                           "at the surface. "
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Lithosphere thickness", "63000",
                           Patterns::Double (),
                           "Thickness of the lithosphere layer generated "
                           "below the crust. "
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Minimum upwelling angle for crust formation", "30",
                           Patterns::Double (0,90),
                           "The minimum upwelling angle required for a particle to be converted into crustal material. "
                           "This angle is measured between the horizontal direction and the particle velocity vector and ranges from 0 to 90 degrees. "
                           "A value of 0 means that all particles are converted into crust (basalt) as soon as they reach the required depth. "
                           "A value of 90 means that only particles moving vertically upward are converted into crust (basalt). "
                           "Units: \\si{\\degree}.");
        prm.declare_entry ("Minimum upwelling angle for lithosphere formation", "30",
                           Patterns::Double (0,90),
                           "The minimum upwelling angle required for a particle to be converted into lithospheric material. "
                           "This angle is measured between the horizontal direction and the particle velocity vector and ranges from 0 to 90 degrees. "
                           "A value of 0 means that all particles (except for particles with a basaltic composition) are converted into lithosphere (harzburgite) as soon as they reach the required depth. "
                           "A value of 90 means that only particles (except for particles with a basaltic composition) moving vertically upward are converted into lithosphere (harzburgite). "
                           "Units: \\si{\\degree}.");
        prm.declare_entry ("Harzburgite profile in lithosphere", "linear",
                           Patterns::Selection ("constant|linear"),
                           "Choose the profile used for the harzburgite fraction within the lithosphere. "
                           "If set to constant, the harzburgite fraction is 100\\% everywhere except where basalt is present. "
                           "If set to linear, the harzburgite fraction is 100\\% at the top of the lithosphere "
                           "and decreases linearly to 0\\% at the bottom of the lithosphere.");
      }



      template <int dim>
      void
      CrustLithosphereFormation<dim>::parse_parameters (ParameterHandler &prm)
      {
        crust_thickness     = prm.get_double ("Crustal thickness");
        lithosphere_thickness = prm.get_double ("Lithosphere thickness");
        minimum_upwelling_angle_for_crust = prm.get_double ("Minimum upwelling angle for crust formation")* constants::degree_to_radians;
        minimum_upwelling_angle_for_lithosphere = prm.get_double ("Minimum upwelling angle for lithosphere formation")* constants::degree_to_radians;

        const std::string harzburgite_profile_selection = prm.get ("Harzburgite profile in lithosphere");
        if (harzburgite_profile_selection == "constant")
          harzburgite_profile = constant;
        else if (harzburgite_profile_selection == "linear")
          harzburgite_profile = linear;
        else
          AssertThrow(false, ExcNotImplemented("The selected harzburgite profile calculation method does not exist."));

        AssertThrow(this->introspection().compositional_name_exists("basalt") &&
                    this->introspection().compositional_name_exists("harzburgite"),
                    ExcMessage("The reaction model <crust and lithosphere formation> "
                               "can only be used if there are compositional fields named "
                               "'basalt' and 'harzburgite'."));

        basalt_index = this->introspection().compositional_index_for_name("basalt");
        harzburgite_index = this->introspection().compositional_index_for_name("harzburgite");
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
#define INSTANTIATE(dim) \
  template class CrustLithosphereFormation<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)
#undef INSTANTIATE
    }
  }
}
