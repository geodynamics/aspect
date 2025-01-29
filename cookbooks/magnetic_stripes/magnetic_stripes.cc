/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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

#include <aspect/material_model/composition_reaction.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class MagneticStripes : public MaterialModel::CompositionReaction<dim>
    {
      public:

        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        static void declare_parameters (ParameterHandler &prm);
        void parse_parameters (ParameterHandler &prm) override;

      private:
        std::vector<double> reversal_times;
    };


    template <int dim>
    void
    MagneticStripes<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      CompositionReaction<dim>::evaluate(in, out);

      int magnetic_orientation = 1;
      if (this->simulator_is_past_initialization())
        for (unsigned int i=0; i < reversal_times.size(); ++i)
          if (this->get_time() < reversal_times[i])
            {
              // magnetic orientation is either 1 (normal) or -1 (reverse)
              magnetic_orientation = -2 * (i % 2) + 1;
              break;
            }

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const double depth = this->get_geometry_model().depth(in.position[i]);
          const double reaction_depth = 7000.0;

          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            {
              // Magnetic lineations are only generated when material approaches the surface,
              // i.e. above the reaction depth and when the velocity is predominantly upwards.
              if (depth < reaction_depth && in.velocity[i][1] > std::abs(in.velocity[i][0]))
                out.reaction_terms[i][c] = -in.composition[i][0] + magnetic_orientation;
              else
                out.reaction_terms[i][c] = 0.0;
            }
        }
    }

    template <int dim>
    void
    MagneticStripes<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Magnetic stripes model");
        {
          CompositionReaction<dim>::declare_parameters (prm);

          prm.declare_entry ("Reversal times", "5.0, 5.0, 2.0, 2.0, 2.092, 2.419, 2.419",
                             Patterns::List(Patterns::Double(0)),
                             "Reversal times of the magnetic field."
                             "Units: yr or s, depending on the ``Use years "
                             "in output instead of seconds'' parameter.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MagneticStripes<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Magnetic stripes model");
        {
          CompositionReaction<dim>::parse_parameters (prm);

          reversal_times = Utilities::string_to_double (Utilities::split_string_list(prm.get ("Reversal times")));
          std::sort (reversal_times.begin(), reversal_times.end());

          if (this->get_parameters().convert_to_years == true)
            for (unsigned int i=0; i<reversal_times.size(); ++i)
              reversal_times[i] *= constants::year_in_seconds;
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
    ASPECT_REGISTER_MATERIAL_MODEL(MagneticStripes,
                                   "magnetic stripes",
                                   "A simple material model that is like the 'composition "
                                   "reaction' model, but with a different reaction term. "
                                   "In areas above a depth of 7 km and with a velocity that "
                                   "does not diverge from vertically upwards by more than 45 "
                                   "degrees, the first compositional field used in the model "
                                   "is set to either one or minus one. The sign depends on a "
                                   "list of reversal times that is passed to the material "
                                   "model as an input. Before the 1st, 3rd, 5th, ... entry, "
                                   "the sign is positive (magnetic normal orientation), and "
                                   "before the 2nd, 4th, 6th, ... entry, the sign is negative "
                                   "(magnetic reverse orientation). List entries are sorted "
                                   "automatically in ascending order.")
  }
}
