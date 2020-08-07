/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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



#include <aspect/mesh_refinement/isolines.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <math.h>
#include <algorithm>

namespace aspect
{
  namespace MeshRefinement
  {
    namespace Internal
    {
      bool
      Isoline::values_are_in_range(const std::vector<double> values) const
      {
        // This assumes that all the vectors in isolines already have the
        // same length.
        Assert(values.size() == min_values.size(),
               ExcMessage("Internal error: Vector of values passed to the isoline class function values_are_in_range, "
                          "does not have the the correct size."));
        for (unsigned int index = 0; index < values.size(); ++index)
          {
            if (values[index] < min_values[index] || values[index] > max_values[index])
              {
                // outside this isoline, no need to search further
                return false;
              }
          }
        // If we made it this far, then we are inside the conditions, so return true.
        return true;
      }

      Property::Property(const std::string &property_name,const std::vector<std::string> &available_compositions)
      {
        bool found = false;
        if (property_name == "Temperature")
          {
            name = PropertyName::Temperature;
            index = 0;
            found = true;
          }
        else if (property_name == "background")
          {
            name = PropertyName::Background;
            index = 0;
            found = true;
          }
        else
          {
            auto p = std::find(available_compositions.begin(), available_compositions.end(), property_name);
            if (p != available_compositions.end())
              {
                name = PropertyName::Composition;
                index = std::distance(available_compositions.begin(), p);
                found = true;
              }
          }
        if (found == false)
          {
            // The property name has not been found. This could come from that the compositional field was not present.
            // Abort and warn the user.
            std::string key_list = "Temperature";
            for (auto &composition : available_compositions)
              key_list += ", " + composition;

            AssertThrow(false,
                        ExcMessage("The key given for the isoline could not be converted. The provided key was: " + property_name + ". "
                                   "The following keys are allowed for this model: " + key_list + "."));
          }
      }

      unsigned int min_max_string_to_int(const std::string &string_value, const unsigned int minimum_refinement_level, const unsigned int  maximum_refinement_level)
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
                return minimum_refinement_level+Utilities::string_to_int(tmpNumber[1]);
              }
            else
              {
                AssertThrow(string.compare(0,4,"min-") != 0,
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
                return maximum_refinement_level-Utilities::string_to_int(tmpNumber[1]);
              }
            else
              {
                AssertThrow(string.compare(0,4,"maxi") != 0,
                            ExcMessage("A value of " + string_value + " was provided, but you can't provide a larger value than the maximum."));
                return maximum_refinement_level;
              }
          }
        else
          {
            return Utilities::string_to_int(string);
          }

      }
    }

    template <int dim>
    void
    Isolines<dim>::update ()
    { }

    template <int dim>
    void
    Isolines<dim>::tag_additional_cells () const
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

      unsigned i = 0;
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              bool coarsen = false;
              bool clear_refine = false;
              bool refine = false;
              bool clear_coarsen = false;

              fe_values.reinit(cell);
              in.reinit(fe_values, cell, this->introspection(), this->get_solution(), true);

              for (unsigned int i_quad=0; i_quad<quadrature.size(); ++i_quad)
                {
                  for (auto &isoline : isolines)
                    {
                      // setup the vector to check
                      std::vector<double> values(isoline.min_values.size());
                      for (unsigned int index = 0; index < isoline.properties.size(); ++index)
                        {
                          if (isoline.properties[index].name == Internal::PropertyName::Temperature)
                            {
                              values[index] = in.temperature[i_quad];
                            }
                          else if (isoline.properties[index].name == Internal::PropertyName::Composition)
                            {
                              values[index] = std::min(std::max(in.composition[i_quad][isoline.properties[index].index], 0.0), 1.0);
                            }
                          else if (isoline.properties[index].name == Internal::PropertyName::Background)
                            {
                              double sum_composition = 0.0;  // Compute background material fraction
                              for (unsigned i=0 ; i < in.composition[i_quad].size() ; ++i)
                                {
                                  sum_composition += std::min(std::max(in.composition[i_quad][i], 0.0), 1.0);
                                }
                              if (sum_composition >= 1.0)
                                values[index] = 0.0;
                              else
                                values[index] = 1.0 - sum_composition;
                            }
                        }

                      if ( isoline.values_are_in_range(values))
                        {
                          // If the current refinement level is smaller or equal to the minimum
                          // refinement level, any coarsening flags should be cleared.
                          if (cell->level() <= isoline.min_refinement)
                            {
                              clear_coarsen = true;
                            }

                          // If the current refinement level is smaller then the minimum
                          // refinement level, a refinment flag should be placed.
                          if (cell->level() <  isoline.min_refinement)
                            {
                              refine = true;
                              break;
                            }

                          // If the current refinement level is larger or equal to the maximum refinement
                          // level, any refinement flag should be cleared.
                          if (cell->level() >= isoline.max_refinement)
                            {
                              clear_refine = true;
                            }

                          // If the current refinement level is larger then the maximum refinemment level,
                          // a coarsening flag should be placed.
                          if (cell->level() >  isoline.max_refinement)
                            {
                              coarsen = true;
                            }
                        }
                    }
                }

              // if both coarsen and refine are true, give preference to refinement
              if (coarsen == true && refine == true)
                coarsen = false;
              if (refine == true)
                clear_refine = false;
              if (coarsen == true)
                clear_coarsen = false;

              // Perform the actual placement of the coarsening and refinement flags
              // We want to make sure that the refinement never goes below the minimum
              // or above the maximum, so we first check/set the coarsen/refine flag,
              // and then check/set the clear coarsen/refine flag.
              if (coarsen == true)
                {
                  cell->set_coarsen_flag ();
                }
              if (clear_coarsen == true)
                {
                  cell->clear_coarsen_flag ();
                }
              if (refine == true)
                {
                  cell->set_refine_flag ();
                }
              if (clear_refine == true)
                {
                  cell->clear_refine_flag ();
                }
            }
          ++i;
        }
    }

    template <int dim>
    void
    Isolines<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Isolines");
        {
          prm.declare_entry ("Isolines", "depth",
                             Patterns::Anything(),
                             "A list of isoline separated by semi-colins (;). Each isoline entry consists of "
                             "multiple entries separted by a comma. The first two entries indicate the minimum and maximum "
                             "refinement levels respectively. The entries after the first two describe the field the isoline "
                             "applies to, followed by a colin (:) followed by the minimum and maximum grid levels seperated by "
                             "bar (|). An example for an isoline is '0, 2, Temperature: 300 | 600'; 2, 2, C_1: 0.5 | 1");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Isolines<dim>::parse_parameters (ParameterHandler &prm)
    {
      // lookup the minimum and maximum refinment
      prm.enter_subsection("Mesh refinement");
      unsigned int minimum_refinement_level = Utilities::string_to_int(prm.get ("Minimum refinement level"));
      unsigned int maximum_refinement_level = Utilities::string_to_int(prm.get("Initial global refinement")) + Utilities::string_to_int(prm.get("Initial adaptive refinement"));
      prm.leave_subsection();

      const std::vector<std::string> compositions = this->introspection().get_composition_names();

      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Isolines");
        {
          // Split the list by comma delimited components.
          const std::vector<std::string> isoline_entries = dealii::Utilities::split_string_list(prm.get("Isolines"), ';');
          unsigned int isoline_entry_number = 0;
          for (auto &isoline_entry : isoline_entries)
            {
              isoline_entry_number++;
              aspect::MeshRefinement::Internal::Isoline isoline;  // a new object of isoline
              std::vector<aspect::MeshRefinement::Internal::Property> properties;  // a vector of Property
              std::vector<double> min_value_inputs;
              std::vector<double> max_value_inputs;
              const std::vector<std::string> field_entries = dealii::Utilities::split_string_list(isoline_entry, ',');

              AssertThrow(field_entries.size() >= 3,
                          ExcMessage("An isoline needs to contain at least 3 entries, but isoline " + std::to_string(isoline_entry_number)
                                     + " contains  only " +  std::to_string(field_entries.size()) + " entries: " + isoline_entry + "."));

              // convert a potential min, min+1, min + 1, min+10, max, max-1, etc. to actual integers.
              isoline.min_refinement = Internal::min_max_string_to_int(field_entries[0], minimum_refinement_level, maximum_refinement_level);
              isoline.max_refinement = Internal::min_max_string_to_int(field_entries[1], minimum_refinement_level, maximum_refinement_level);
              AssertThrow(isoline.min_refinement <= isoline.max_refinement,
                          ExcMessage("The provided maximum refinement level has to be larger the then the minimum refinement level."));
              for (auto field_entry = field_entries.begin()+2; field_entry < field_entries.end(); ++field_entry)
                {
                  AssertThrow(Patterns::Map(Patterns::Anything(),
                                            Patterns::List(Patterns::Double(), 0, std::numeric_limits<unsigned int>::max(), "|")
                                           ).match(*field_entry),
                              ExcMessage("The isoline is not formatted correctly."));
                  std::vector<std::string> key_and_value = Utilities::split_string_list (*field_entry, ':');
                  AssertThrow(key_and_value.size() == 2,
                              ExcMessage("The isoline property must have a key (e.g. Temperature) and two values separated by a | (e.g. (300 | 600)."));
                  properties.push_back(Internal::Property(key_and_value[0], compositions)); // convert key to property name
                  const std::vector<std::string> values = dealii::Utilities::split_string_list(key_and_value[1], '|');
                  AssertThrow(values.size() == 2,
                              ExcMessage("Both a maximum and minimum value is required for each isoline."));
                  min_value_inputs.push_back(Utilities::string_to_double(values[0]));  // get min and max values of the range
                  max_value_inputs.push_back(Utilities::string_to_double(values[1]));
                  AssertThrow(min_value_inputs.back() < max_value_inputs.back(),
                              ExcMessage("The provided maximum refinement level has to be larger than the minimum refinement level."));
                }
              isoline.min_values = min_value_inputs;
              isoline.max_values = max_value_inputs;
              isoline.properties = properties;
              isolines.push_back(isoline);
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
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Isolines,
                                              "isolines",
                                              "A mesh refinement criterion that computes "
                                              "refinement indicators between two iso-surfaces of"
                                              "specific field entries(e.g. temperature, compsitions)."
                                              "\n\n"
                                              "The way these indicators are derived on each isoline is by "
                                              "checking the conditions whether solutions of specific "
                                              "fields are within the ranges of values given. If these conditions"
                                              "hold, then an indicator is either put on or taken off"
                                              "in order to secure mesh refinement within the range of levels"
                                              "given. Usage of this plugin allows user to put an conditional"
                                              "minimum and maximum refinement function onto fields that they"
                                              "are interested in."
                                              "\n\n"
                                              "For now, only temperature and names of compositions are allowed as field"
                                              "entries"
                                             )
  }
}
