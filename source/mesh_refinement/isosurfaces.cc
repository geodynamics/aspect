/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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



#include <algorithm>
#include <aspect/geometry_model/interface.h>
#include <aspect/mesh_refinement/isosurfaces.h>
#include <aspect/utilities.h>
#include <cmath>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MeshRefinement
  {
    namespace internal
    {
      bool
      Isosurface::are_all_values_in_range(const std::vector<double> &values) const
      {
        // This assumes that all the vectors in isosurfaces already have the
        // same length.
        Assert(values.size() == min_values.size(),
               ExcMessage("internal error: Vector of values passed to the isosurface class "
                          "function are_all_values_in_range, does not have the correct size."));
        for (unsigned int index = 0; index < values.size(); ++index)
          {
            if (values[index] < min_values[index] || values[index] > max_values[index])
              {
                // outside this isosurface, no need to search further
                return false;
              }
          }
        // If we made it this far, then we are inside the conditions, so return true.
        return true;
      }

      Property::Property(const std::string &property_name,
                         const std::vector<std::string> &available_compositions)
      {
        bool found = false;
        if (property_name == "Temperature")
          {
            type = PropertyType::Temperature;
            index = 0;
            found = true;
          }
        else
          {
            auto p = std::find(available_compositions.begin(), available_compositions.end(), property_name);
            if (p != available_compositions.end())
              {
                type = PropertyType::Composition;
                index = std::distance(available_compositions.begin(), p);
                found = true;
              }
          }
        if (found == false)
          {
            // The property name has not been found. This could be because the compositional field is not present.
            // Abort and warn the user.
            std::string key_list = "Temperature";
            for (auto &composition : available_compositions)
              key_list += ", " + composition;

            AssertThrow(false,
                        ExcMessage("The key given for the isosurface could not be converted. The provided key was: " + property_name + ". "
                                   "The following keys are allowed for this model: " + key_list + "."));
          }
      }

      /**
       * Convert a string describing the refinement level which potentially contains the words "min"
       * or "max" to an integer. The format is "min", "max", "min+X", "max-X", or "X", where X is a positive integer.
       *
       * For example if minimum_refinement_level  is 1, and the string is 'min+1', it will return 2.
       * The returned value will be capped  by the provided @p minimum_refinement_level and @p maximum_refinement_level.
       * The @p minimum_refinement_level and @p maximum_refinement_level are set in this input file subsection for this
       * plugin, but in the Mesh refinement subsection.
       *
       * This function will throw an assertion when the string contains `max+` or `min-` since
       * that is always incorrect.
       */
      unsigned int min_max_string_to_int(const std::string &string_value,
                                         const unsigned int minimum_refinement_level,
                                         const unsigned int maximum_refinement_level)
      {
        // start with removing the spaces so that a ' min + 1 ' would become 'min+1'.
        std::string string = string_value;
        std::string::iterator end_pos = std::remove(string.begin(), string.end(), ' ');
        string.erase(end_pos, string.end());

        // check whether the field starts with a 'min'
        if (string.compare(0,3,"min") == 0)
          {
            if (string.compare(0,4,"min+") == 0)
              {
                std::vector<std::string> tmpNumber = Utilities::split_string_list(string,'+');
                AssertThrow(tmpNumber.size() == 2,
                            ExcMessage("Could not convert value '" + string + "' to an int because it contains more than one '+' sign."));
                return std::min(std::max(minimum_refinement_level+Utilities::string_to_int(tmpNumber[1]),minimum_refinement_level),maximum_refinement_level);
              }
            else
              {
                AssertThrow(string == "min",
                            ExcMessage("A value of " + string_value + " was provided, but you can't provide a smaller value than the minimum."));
                return minimum_refinement_level;
              }
          }
        else if (string.compare(0,3,"max") == 0)
          {
            if (string.compare(0,4,"max-") == 0)
              {
                std::vector<std::string> tmpNumber = Utilities::split_string_list(string,'-');
                AssertThrow(tmpNumber.size() == 2,
                            ExcMessage("Could not convert value '" + string + "' to an int because it contains more than one '-' sign."));
                return std::min(std::max(maximum_refinement_level-Utilities::string_to_int(tmpNumber[1]),minimum_refinement_level),maximum_refinement_level);
              }
            else
              {
                AssertThrow(string=="max",
                            ExcMessage("A value of " + string_value + " was provided, but you can't provide a larger value than the maximum."));
                return maximum_refinement_level;
              }
          }
        else
          {
            return std::min(std::max(static_cast<unsigned int>(Utilities::string_to_int(string)),minimum_refinement_level),maximum_refinement_level);
          }

      }
    }

    template <int dim>
    void
    Isosurfaces<dim>::update ()
    { }

    template <int dim>
    void
    Isosurfaces<dim>::tag_additional_cells () const
    {
      if (this->get_dof_handler().n_locally_owned_dofs() == 0)
        return;

      const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.temperature).get_unit_support_points());
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values | update_gradients);

      MaterialModel::MaterialModelInputs<dim> in(quadrature.size(),
                                                 this->n_compositional_fields());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              bool coarsen = false;
              bool clear_refine = false;
              bool refine = false;
              bool clear_coarsen = false;

              fe_values.reinit(cell);
              in.reinit(fe_values, cell, this->introspection(), this->get_solution());

              for (unsigned int i_quad=0; i_quad<quadrature.size(); ++i_quad)
                {
                  for (auto &isosurface : isosurfaces)
                    {
                      // setup the vector containing the isosurfaces
                      std::vector<double> values(isosurface.min_values.size());
                      for (unsigned int index = 0; index < isosurface.properties.size(); ++index)
                        {
                          if (isosurface.properties[index].type == internal::PropertyType::Temperature)
                            {
                              values[index] = in.temperature[i_quad];
                            }
                          else if (isosurface.properties[index].type == internal::PropertyType::Composition)
                            {
                              values[index] = in.composition[i_quad][isosurface.properties[index].index];
                            }
                        }

                      if ( isosurface.are_all_values_in_range(values))
                        {
                          // If the current refinement level is smaller than or equal to the minimum
                          // refinement level, any coarsening flag should be cleared.
                          if (cell->level() <= isosurface.min_refinement)
                            {
                              clear_coarsen = true;
                            }

                          // If the current refinement level is smaller than the minimum
                          // refinement level, a refinement flag should be placed.
                          if (cell->level() <  isosurface.min_refinement)
                            {
                              refine = true;
                              break;
                            }

                          // If the current refinement level is larger than or equal to the maximum refinement
                          // level, any refinement flag should be cleared.
                          if (cell->level() >= isosurface.max_refinement)
                            {
                              clear_refine = true;
                            }

                          // If the current refinement level is larger than the maximum refinement level,
                          // a coarsening flag should be placed.
                          if (cell->level() >  isosurface.max_refinement)
                            {
                              coarsen = true;
                            }
                        }
                    }
                }

              // if both coarsen and refine are true, give preference to refinement
              if (coarsen && refine)
                coarsen = false;
              if (refine)
                clear_refine = false;
              if (clear_coarsen)
                coarsen = false;

              // Perform the actual placement of the coarsening and refinement flags.
              if (clear_coarsen)
                cell->clear_coarsen_flag ();
              if (clear_refine)
                cell->clear_refine_flag ();
              if (coarsen)
                cell->set_coarsen_flag ();
              if (refine)
                cell->set_refine_flag ();
            }
        }
    }

    template <int dim>
    void
    Isosurfaces<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Isosurfaces");
        {
          prm.declare_entry ("Isosurfaces", "",
                             Patterns::Anything(),
                             "A list of isosurfaces separated by semi-colons (;). Each isosurface entry consists of "
                             "multiple entries separated by a comma. The first two entries indicate the minimum and maximum "
                             "refinement levels respectively. The entries after the first two describe the fields the isosurface "
                             "applies to, followed by a colon (:), which is again followed by the minimum and maximum property values separated by "
                             "bar (|). An example for an isosurface is '0, 2, Temperature: 300 | 600; 2, 2, C\\_1: 0.5 | 1'. "
                             "In this example the mesh refinement is kept between level 0 and level 2 if the temperature is between "
                             "300 and 600 and at level 2 when the compositional field C\\_1 is between 0.5 and 1. If both happen at "
                             "the same location and the current refinement level is 1, it means that the first isoline will not set any flag and the "
                             "second isoline will set a refinement flag. This means the cell will be refined. If both the coarsening and refinement flags "
                             "are set, preference is given to refinement. "
                             "\n\n"
                             "The first two entries for each isosurface, describing the minimum and maximum grid levels, can be "
                             "two numbers or contain one of the key values 'min' and 'max'. This indicates the key will be replaced with "
                             "the global minimum and maximum refinement levels. The 'min' and 'max' keys also accept adding "
                             "values to be added or subtracted from them respectively. This is done by adding a '+' or '-' and a "
                             "number behind them (e.g. min+2 or max-1). Note that you can't subtract a value from a minimum value or "
                             "add a value to the maximum value. If, for example, `max-4` drops below the minimum or `min+4` goes above the "
                             "maximum, it will simply use the global minimum and maximum values respectively. The same holds for any "
                             "mesh refinement level below the global minimum or above the global maximum.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Isosurfaces<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        // lookup the minimum and maximum refinement
        const unsigned int minimum_refinement_level = this->get_parameters().min_grid_level;
        const unsigned int maximum_refinement_level = this->get_parameters().initial_global_refinement +
                                                      this->get_parameters().initial_adaptive_refinement;

        const std::vector<std::string> &compositions = this->introspection().get_composition_names();
        prm.enter_subsection("Isosurfaces");
        {
          // Split the list by comma delimited components.
          const std::vector<std::string> isosurface_entries = dealii::Utilities::split_string_list(prm.get("Isosurfaces"), ';');
          unsigned int isosurface_entry_number = 0;
          for (auto &isosurface_entry : isosurface_entries)
            {
              ++isosurface_entry_number;
              aspect::MeshRefinement::internal::Isosurface isosurface;
              std::vector<aspect::MeshRefinement::internal::Property> properties;
              std::vector<double> min_value_inputs;
              std::vector<double> max_value_inputs;
              const std::vector<std::string> field_entries = dealii::Utilities::split_string_list(isosurface_entry, ',');

              AssertThrow(field_entries.size() >= 3,
                          ExcMessage("An isosurface needs to contain at least 3 entries, but isosurface " + std::to_string(isosurface_entry_number)
                                     + " contains  only " +  std::to_string(field_entries.size()) + " entries: " + isosurface_entry + "."));

              // convert a potential min, min+1, min + 1, min+10, max, max-1, etc. to actual integers.
              isosurface.min_refinement = internal::min_max_string_to_int(field_entries[0], minimum_refinement_level, maximum_refinement_level);
              isosurface.max_refinement = internal::min_max_string_to_int(field_entries[1], minimum_refinement_level, maximum_refinement_level);
              AssertThrow(isosurface.min_refinement <= isosurface.max_refinement,
                          ExcMessage("The provided maximum refinement level has to be larger than the minimum refinement level."));
              for (auto field_entry = field_entries.begin()+2; field_entry != field_entries.end(); ++field_entry)
                {
                  AssertThrow(Patterns::Map(Patterns::Anything(),
                                            Patterns::List(Patterns::Double(), 0, std::numeric_limits<unsigned int>::max(), "|")
                                           ).match(*field_entry),
                              ExcMessage("The isosurface is not formatted correctly."));
                  std::vector<std::string> key_and_value = Utilities::split_string_list (*field_entry, ':');
                  AssertThrow(key_and_value.size() == 2,
                              ExcMessage("The isosurface property must have a key (e.g. Temperature) and two values separated by a | (e.g. (300 | 600)."));
                  properties.emplace_back(key_and_value[0], compositions); // convert key to property name
                  const std::vector<std::string> values = dealii::Utilities::split_string_list(key_and_value[1], '|');
                  AssertThrow(values.size() == 2,
                              ExcMessage("Both a maximum and a minimum value are required for each isosurface."));
                  min_value_inputs.push_back(Utilities::string_to_double(values[0]));  // get min and max values of the range
                  max_value_inputs.push_back(Utilities::string_to_double(values[1]));
                  AssertThrow(min_value_inputs.back() < max_value_inputs.back(),
                              ExcMessage("The provided maximum property value has to be larger than the provided minimum property value."));
                }
              isosurface.min_values = min_value_inputs;
              isosurface.max_values = max_value_inputs;
              isosurface.properties = properties;
              isosurfaces.push_back(isosurface);
            }
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
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Isosurfaces,
                                              "isosurfaces",
                                              "A mesh refinement criterion that computes "
                                              "coarsening and refinement indicators between two isosurfaces of "
                                              "specific field entries (e.g. temperature, composition)."
                                              "\n\n"
                                              "The way these indicators are derived between pairs of isosurfaces is by "
                                              "checking whether the solutions of specific "
                                              "fields are within the ranges of the isosurface values given. If these conditions "
                                              "hold, then coarsening and refinement indicators are set such that "
                                              "the mesh refinement levels lies within the range of levels "
                                              "given. Usage of this plugin allows the user to put a conditional "
                                              "minimum and maximum refinement function onto fields that they "
                                              "are interested in."
                                              "\n\n"
                                              "For now, only temperature and compositional fields are allowed as "
                                              "field entries. The key words could be 'Temperature' or one of the names "
                                              "of the compositional fields which are either specified by user or set up "
                                              "as C\\_0, C\\_1, etc."
                                              "\n\n"
                                              "Usage: A list of isosurfaces separated by semi-colons (;). Each "
                                              "isosurface entry consists of multiple entries separated by a comma. "
                                              "The first two entries indicate the minimum and maximum refinement "
                                              "levels respectively. The entries after the first two describe the "
                                              "fields the isosurface applies to, followed by a colon (:), which again "
                                              "is followed by the minimum and maximum field values separated by a bar (|). An "
                                              "example for two isosurface entries is '0, 2, Temperature: 300 | 600; 2, 2, C\\_1: 0.5 | 1'. "
                                              "If both isoterm entries are triggered at the same location and the current refinement level is 1, "
                                              "it means that the first isoline will not set any flag and the second isoline will set a refinement flag. "
                                              "This means the cell will be refined. If both the coarsening and refinement flags "
                                              "are set, preference is given to refinement. "
                                              "\n\n"
                                              "The minimum and maximum refinement levels per isosurface can be provided in absolute values relative to "
                                              "the global minimum and maximum refinement. This is done with the 'min' and 'max' key words. For example: "
                                              "'set Isosurfaces = max-2,  max,    Temperature: 0 | 600 ; min + 1,min+2, Temperature: 1600 | 3000,   C\\_2 : 0.0 | 0.5'."
                                             )
  }
}
