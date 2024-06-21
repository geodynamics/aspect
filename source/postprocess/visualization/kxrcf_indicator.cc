/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/kxrcf_indicator.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      KXRCFIndicator<dim>::KXRCFIndicator()
        :
        CellDataVectorCreator<dim>("")
      {}


      template <int dim>
      std::pair<std::string, std::unique_ptr<Vector<float>>>
      KXRCFIndicator<dim>::execute() const
      {
        std::pair<std::string, std::unique_ptr<Vector<float>>>
        return_value ("KXRCF_indicator",
                      std::make_unique<Vector<float>>(this->get_triangulation().n_active_cells()));
        this->compute_KXRCF_indicators(*return_value.second, field_index);
        return return_value;
      }


      template <int dim>
      void
      KXRCFIndicator<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("KXRCF indicator");
            {
              prm.declare_entry ("Name of advection field", "",
                                 Patterns::Anything(),
                                 "The name of the advection field whose output "
                                 "should be visualized. ");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      KXRCFIndicator<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("KXRCF indicator");
            {
              const std::string field_name = prm.get("Name of advection field");

              if (field_name == "temperature")
                field_index = 0;
              else
                {
                  AssertThrow(this->introspection().compositional_name_exists(field_name),
                              ExcMessage("No compositional field with name <" +
                                         field_name +
                                         "> exists for which you want to visualize the KXRCF indicator."));

                  field_index = this->introspection().compositional_index_for_name(field_name) + 1;
                }
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(KXRCFIndicator,
                                                  "kxrcf indicator",
                                                  "A visualization output object that generates output "
                                                  "showing the value of the KXRCF indicator on each "
                                                  "cell."
                                                  "\n\n"
                                                  "Physical units: none.")
    }
  }
}
