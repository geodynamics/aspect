/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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

#include <aspect/mesh_refinement/isotherms.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/derivative_approximation.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    Isotherms<dim>::tag_additional_cells () const
    {
      if (this->get_dof_handler().n_locally_owned_dofs() == 0)
        return;

      const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.compositional_fields).get_unit_support_points());
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values | update_gradients);

      MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());

      for (typename DoFHandler<dim>::active_cell_iterator
           cell = this->get_dof_handler().begin_active();
           cell != this->get_dof_handler().end(); ++cell)
        {
          if (cell->is_locally_owned())
            {
              bool coarsen = false;
              bool clear_refine = false;
              bool refine = false;
              bool clear_coarsen = false;

              fe_values.reinit(cell);
              in.reinit(fe_values, cell, this->introspection(), this->get_solution(), true);
              this->get_material_model().evaluate(in, out);

              for (unsigned int i=0; i<quadrature.size(); ++i)
                {
                  for ( unsigned int isotherm_line = 0; isotherm_line != isotherms.size(); ++isotherm_line )
                    {

                      /// determine if cell should be refined or coarsened
                      if (exclude_composition > 0 && exclude_composition <= static_cast <signed int> (this->n_compositional_fields()))  // static cast to prevent compiler warnings. Will only go wrong when there are more compositional fields then a the positive part of an int can handle ;)
                        {
                          // there is a exclude compostion (>0) and  exclude composition is smaller or equal to the current composition
                          if (in.temperature[i] <= isotherms[isotherm_line].second && in.temperature[i] >= isotherms[isotherm_line].first && in.composition[isotherm_line][exclude_composition]<0.5)
                            {
                              // the temperature is between the isotherms and the exclude composition is smaller then 0.5 at this location

                              // If the current refinement level is smaller or equal to the refinement minimum level, any coarsening flags should be cleared.
                              if (cell->level() <= isotherms_levels[isotherm_line].first)
                                {
                                  clear_coarsen = true;
                                }
                              // If the current refinement level is smaller then the minimum level, a refinment flag should be placed.
                              if (cell->level() <  isotherms_levels[isotherm_line].first)
                                {
                                  refine = true;
                                  break;
                                }

                              // If the current refinement level is larger or equal to the maximum refinement level, any refinement flag should be cleared.
                              if (cell->level() >= isotherms_levels[isotherm_line].second)
                                {
                                  clear_refine = true;
                                }
                              // If the current refinement level is larger then the maximum level, a coarsening flag should be placed.
                              if (cell->level() >  isotherms_levels[isotherm_line].second)
                                {
                                  coarsen = true;
                                }

                            }
                        }
                      else
                        {
                          // there is not a exclude compostion (>0) and/or the exclude composition is larger or equal to the current composition. Now we don't have to check the value of the composition anymore.
                          if (in.temperature[i] <= isotherms[isotherm_line].second && in.temperature[i] >= isotherms[isotherm_line].first)
                            {
                              // the temperature is between the isotherms
                              // If the current refinement level is smaller or equal to the refinement minimum level, any coarsening flags should be cleared.
                              if (cell->level() <= isotherms_levels[isotherm_line].first)
                                {
                                  clear_coarsen = true;
                                }
                              // If the current refinement level is smaller then the minimum level, a refinment flag should be placed, and we don't have to look any further.
                              if (cell->level() <  isotherms_levels[isotherm_line].first)
                                {
                                  refine = true;
                                }

                              // If the current refinement level is larger or equal to the maximum refinement level, any refinement flag should be cleared.
                              if (cell->level() >= isotherms_levels[isotherm_line].second)
                                {
                                  clear_refine = true;
                                }
                              // If the current refinement level is larger then the maximum level, a coarsening flag should be placed, and we don't have to look any further.
                              if (cell->level() >  isotherms_levels[isotherm_line].second)
                                {
                                  coarsen = true;
                                }

                            }
                        }
                    }
                }

              // if both coarsen and refine are true, give preference to refinement
              if (coarsen == true && refine == true)
                {
                  coarsen = false;
                  clear_refine = false;
                }

              // Perform the actual placement of the coarsening and refinement flags
              // We want to make sure that the refiment never goes below the minimum
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
        }
    }

    template <int dim>
    void
    Isotherms<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Isotherms");
        {
          prm.declare_entry ("Exclude composition", "-1",
                             Patterns::Double (),
                             "This option allows to exclude one composition from being condidered by this plugin. "
                             "This can for example be used to exclude the composition representing a sticky air, which "
                             "is as cold as the surface, but you do not want to be refined."
                            );
          prm.declare_entry ("List of isotherms", "",
                             Patterns::List(Patterns::List(Patterns::Anything(),0,100000000,","),0,100000000,";"),
                             "A list of isotherms separated by semi-colins (;). Each isotherm entry consists of "
                             "four values separted by a comma. The first two values indicate the minimum and maximum "
                             "refinement levels respecitvely. The can be an integer, but also the words min or max "
                             "optionally followed by a plus or minus a number. For example 'min+2' or 'max-1'. The "
                             "last two values are the mimum and maximum temperature respectively between which the "
                             "mimumun and maximum refinement are enforced."
                            );

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Isotherms<dim>::parse_parameters (ParameterHandler &prm)
    {
      std::vector<std::pair<std::string,std::string> > isotherms_raw;
      prm.enter_subsection("Mesh refinement");
      {
        minimum_refinement_level = Utilities::string_to_int(prm.get ("Minimum refinement level"));
        maximum_refinement_level = Utilities::string_to_int(prm.get("Initial global refinement")) + Utilities::string_to_int(prm.get("Initial adaptive refinement"));
        prm.enter_subsection("Isotherms");
        {
          exclude_composition = Utilities::string_to_int(prm.get("Exclude composition"));
          std::vector<std::string> isotherms_outer_loop = Utilities::split_string_list(prm.get ("List of isotherms"),';');
          // process the List of isotherms to get all data out of the structure
          isotherms_raw.resize(isotherms_outer_loop.size());
          isotherms.resize(isotherms_outer_loop.size());
          // loop through all given isotherms
          for ( unsigned int isotherm_line = 0; isotherm_line != isotherms_raw.size(); ++isotherm_line )
            {
              std::vector<std::string> isotherms_inner_loop = Utilities::split_string_list(isotherms_outer_loop[isotherm_line]);

              // The minimum and maximum refinemnt levels respecively
              isotherms_raw[isotherm_line].first  = isotherms_inner_loop[0];
              isotherms_raw[isotherm_line].second = isotherms_inner_loop[1];
              // The temperatures for the min and max isotherms respectively
              isotherms[isotherm_line].first  = Utilities::string_to_double(isotherms_inner_loop[2]);
              isotherms[isotherm_line].second = Utilities::string_to_double(isotherms_inner_loop[3]);
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      isotherms_levels.resize(isotherms_raw.size());
      // convert max,max-1,max-2 etc and min, min-1,min-2 to numbers in integer
      for ( unsigned int isotherm_line = 0; isotherm_line != isotherms_levels.size(); ++isotherm_line )
        {
          if (isotherms_raw[isotherm_line].first.compare(0,3,"max") == 0)
            {
              if (isotherms_raw[isotherm_line].first.compare(0,4,"max-") == 0)
                {
                  std::vector<std::string> tmpNumber = Utilities::split_string_list(isotherms_raw[isotherm_line].first,'-');
                  isotherms_levels[isotherm_line].first = maximum_refinement_level-Utilities::string_to_int(tmpNumber[1]);
                }
              else
                {
                  isotherms_levels[isotherm_line].first = maximum_refinement_level;
                }
            }
          else if (isotherms_raw[isotherm_line].first.compare(0,3,"min") == 0)
            {
              if (isotherms_raw[isotherm_line].first.compare(0,4,"min+") == 0)
                {
                  std::vector<std::string> tmpNumber = Utilities::split_string_list(isotherms_raw[isotherm_line].first,'+');
                  isotherms_levels[isotherm_line].first = minimum_refinement_level+Utilities::string_to_int(tmpNumber[1]);
                }
              else
                {
                  isotherms_levels[isotherm_line].first = minimum_refinement_level;
                }
            }
          else if (!isotherms_raw[isotherm_line].first.empty() && isotherms_raw[isotherm_line].first.find_first_not_of("0123456789") == std::string::npos)
            {
              isotherms_levels[isotherm_line].first = Utilities::string_to_int(isotherms_raw[isotherm_line].first);
            }
          else
            {

              AssertThrow (true,
                           ExcMessage ("Not able to read the first input at inputline "
                                       +
                                       Utilities::int_to_string(isotherm_line)
                                       +
                                       ": "
                                       +
                                       isotherms_raw[isotherm_line].first
                                       + ". Expecting a number or 'min', 'max', 'min-1', 'max-1', 'min+1', 'max+1', 'min-2', etc."
                                      )
                          );
            }

          if (isotherms_raw[isotherm_line].second.compare(0,3,"max") == 0)
            {
              if (isotherms_raw[isotherm_line].second.compare(0,4,"max-") == 0)
                {
                  std::vector<std::string> tmpNumber = Utilities::split_string_list(isotherms_raw[isotherm_line].second,'-');
                  isotherms_levels[isotherm_line].second = maximum_refinement_level-Utilities::string_to_int(tmpNumber[1]);
                }
              else
                {
                  isotherms_levels[isotherm_line].second = maximum_refinement_level;
                }
            }
          else if (isotherms_raw[isotherm_line].second.compare(0,3,"min") == 0)
            {
              if (isotherms_raw[isotherm_line].second.compare(0,4,"min+") == 0)
                {
                  std::vector<std::string> tmpNumber = Utilities::split_string_list(isotherms_raw[isotherm_line].second,'+');
                  isotherms_levels[isotherm_line].second = minimum_refinement_level-Utilities::string_to_int(tmpNumber[1]);
                }
              else
                {
                  isotherms_levels[isotherm_line].second = minimum_refinement_level;
                }
            }
          else if (!isotherms_raw[isotherm_line].second.empty() && isotherms_raw[isotherm_line].second.find_first_not_of("0123456789") == std::string::npos)
            {
              isotherms_levels[isotherm_line].second = Utilities::string_to_int(isotherms_raw[isotherm_line].second);
            }
          else
            {
              AssertThrow (true,
                           ExcMessage ("Not able to read the second input at inputline "
                                       +
                                       Utilities::int_to_string(isotherm_line)
                                       +
                                       ": "
                                       +
                                       isotherms_raw[isotherm_line].second
                                       + ". Expecting a number or 'min', 'max', 'min-1', 'max-1', 'min+1', 'max+1', 'min-2', etc."
                                      )
                          );
            }
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Isotherms,
                                              "isotherms",
                                              "A mesh refinement criterion that ensures a "
                                              "maximum and minimum refinement level between "
                                              "two temperatures (isotherms), with the possibility "
                                              "to exclude a composition from this criterion. "
                                              "To accomplish this there are two parameters "
                                              "available: 'Exclude composition' and 'List of "
                                              "isotherms'. The first parameter takes the "
                                              "number of the compositional field and excludes "
                                              "it from this plugin. "
                                              "The second parameters takes a list of isotherm "
                                              "parameters for one isotherm separated by a ';'. "
                                              "Each line of isotherm parameters contains four "
                                              "subparameters. The first subparameter sets "
                                              "the minimum refinement level for this isotherm, "
                                              "the second subparameter sets the maximum "
                                              "refinement level for this isotherm, the third "
                                              "subparameter sets the minimum temperature for "
                                              "this isotherm and the fourth subparameter sets "
                                              "the maximum temperature for the isotherm. "
                                              "The minimum and maximum level can be indicated "
                                              "by the absolute number of the refinement level "
                                              "or by using the key words 'min' or 'max', "
                                              "corresponding respectively to the set minimum "
                                              "and maximum refinement level. The keywords "
                                              "'min' and 'max' may also be used in combination "
                                              "with addition or substraction, e.g. 'max-1' "
                                              "or 'min+1'. For formatting it is useful to "
                                              "make use of the built-in option of ASPECT "
                                              "to split parameters over several lines by "
                                              "putting a backslash at the end of the line. "
                                              "\n\n"
                                              "Example input format of List of isotherms:\n"
                                              "\\begin{lstlisting}"
                                              "set List of isotherms = "
                                              "max,\t max,\t 0,\t 1000; \\ \n\t\t 5,\t max-1,\t"
                                              " 1000,\t 1500; \\ \n\t\t min,\t max-1,\t 1500,\t"
                                              " 2000"
                                              "\\end{lstlisting}"
                                              "The criterion is checked for every quaderature point "
                                              "in the cell. If there are both refine and coarsening "
                                              "flags set in a cell, the preference is given to "
                                              "refinement. ")
  }
}
