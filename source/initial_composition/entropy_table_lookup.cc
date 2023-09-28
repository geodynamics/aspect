/*
  Copyright (C) 2017 - 2023 by the authors of the ASPECT code.

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


#include <aspect/initial_composition/entropy_table_lookup.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    void
    EntropyTableLookUp<dim>::initialize()
    {
      AssertThrow (this->introspection().compositional_name_exists("entropy"),
                   ExcMessage("The 'entropy table lookup' initial composition requires the existence of a compositional field "
                              "named 'entropy'. This field does not exist."));

      // Make sure we keep track of the initial temperature manager and
      // that it continues to live beyond the time when the simulator
      // class releases its pointer to it.
      initial_temperature_manager = this->get_initial_temperature_manager_pointer();

      // Make sure we keep track of the initial composition manager and
      // that it continues to live beyond the time when the simulator
      // class releases its pointer to it.
      initial_composition_manager = this->get_initial_composition_manager_pointer();

      entropy_index = this->introspection().compositional_index_for_name("entropy");

      material_lookup = std::make_unique<Utilities::StructuredDataLookup<2>>(7,1.0);
      material_lookup->load_file(data_directory+material_file_name,
                                 this->get_mpi_communicator());
    }


    template <int dim>
    double
    EntropyTableLookUp<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int compositional_index) const
    {
      if (compositional_index == entropy_index)
        {
          const double temperature = initial_temperature_manager->initial_temperature(position);
          const double pressure = this->get_adiabatic_conditions().pressure(position);

          // Convert pressure from Pa to bar, bar is used in the table.
          Point<2> temperature_pressure(temperature, pressure / 1.e5);

          const double entropy = material_lookup->get_data(temperature_pressure, 0);

          return entropy;
        }
      return 0.0;
    }

    template <int dim>
    void
    EntropyTableLookUp<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Entropy table lookup");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/entropy-table/pyrtable/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT.");
          prm.declare_entry ("Material file name", "material_table_temperature_pressure.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file name of the material data.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    EntropyTableLookUp<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Entropy table lookup");
        {
          data_directory              = Utilities::expand_ASPECT_SOURCE_DIR(prm.get("Data directory"));
          material_file_name          = prm.get("Material file name");
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
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(EntropyTableLookUp,
                                              "entropy table lookup",
                                              "A class that implements initial conditions for the entropy field "
                                              "by converting the initial temperature field through a look up table. "
                                              "Note that this plugin only works if there is a compositional field "
                                              "called `entropy', and an additional look up table that can convert "
                                              "pressure and temperature to entropy. "
                                              "For all compositional fields except entropy this plugin returns 0.0, "
                                              "and they are therefore not changed as long as the default `add' "
                                              "operator is selected for this plugin.")
  }
}
