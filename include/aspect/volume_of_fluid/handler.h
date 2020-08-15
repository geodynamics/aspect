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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_volume_of_fluid_handler_h
#define _aspect_volume_of_fluid_handler_h

#include <aspect/simulator.h>
#include <aspect/simulator_access.h>
#include <aspect/volume_of_fluid/field.h>
#include <aspect/volume_of_fluid/assembly.h>

#include <boost/serialization/map.hpp>

using namespace dealii;

namespace aspect
{
  /**
   * A member class that isolates the functions and variables that deal
   * with the Volume of Fluid implementation. If Volume of Fluid interface
   * tracking is not active, there is no instantiation of this class at
   * all.
   */
  template <int dim>
  class VolumeOfFluidHandler : public SimulatorAccess<dim>
  {
    public:
      /**
       * Standard initial constructor
       */
      VolumeOfFluidHandler(Simulator<dim> &simulator, ParameterHandler &prm);

      /**
       * Add the Volume of Fluid field declaration to the list to be included
       * in the solution vector.
       */
      void edit_finite_element_variables (std::vector<VariableDeclaration<dim> > &vars);

      /**
       * Declare the parameters this class takes through input files.
       */
      static
      void declare_parameters (ParameterHandler &prm);

      /**
       * Read the parameters this class declares from the parameter file.
       */
      void parse_parameters (ParameterHandler &prm);

      /**
       * Get the number of volume of fluid fields in current model
       */
      unsigned int get_n_fields() const;

      /**
       * Get the name of volume of fluid field with index i
       */
      const std::string name_for_field_index(unsigned int i) const;

      /**
       * Get the structure containing the variable locations for the volume of
       * fluid field with index i.
       */
      const VolumeOfFluidField<dim> &field_struct_for_field_index(unsigned int i) const;

      /**
       * Get threshold for volume fraction
       */
      double get_volume_fraction_threshold() const;

      /**
       * Get the local index (within vof fields) for the named composition/volume of fluid field
       */
      unsigned int field_index_for_name(const std::string &fieldname) const;

      /**
       * Do necessary internal initialization that is dependent on having the
       * simulator and Finite Element initialized.
       */
      void initialize (ParameterHandler &prm);

      /**
       * Do initialization routine for all volume of fluid fields
       */
      void set_initial_volume_fractions ();

      /**
       * Initialize specified field based on a composition field initial condition
       */
      void initialize_from_composition_field (const VolumeOfFluidField<dim> &field);

      /**
       * Initialize specified field based on a level set initial condition
       */
      void initialize_from_level_set (const VolumeOfFluidField<dim> &field);

      /**
       * Do interface reconstruction for specified field and cache result in solution vector
       */
      void update_volume_of_fluid_normals (const VolumeOfFluidField<dim> &field,
                                           LinearAlgebra::BlockVector &solution);

      /**
       * Use current interface reconstruction to produce a composition field
       * approximation that is bilinear on the unit cell and write that field
       * to the specified AdvectionField
       */
      void update_volume_of_fluid_composition (const typename Simulator<dim>::AdvectionField &composition_field,
                                               const VolumeOfFluidField<dim> &volume_of_fluid_field,
                                               LinearAlgebra::BlockVector &solution);

      /**
       * Do single timestep update, includes logic for doing Strang split update
       */
      void do_volume_of_fluid_update (const typename Simulator<dim>::AdvectionField &advection_field);

      /**
       * Assemble matrix and RHS for the specified field and dimension
       * (calculation_dim).  If update_from_old_solution is true, the initial
       * values for this update step are in old_solution, otherwise the values
       * in solution are used. This allows a clean restart of the split update
       * from the last timestep if necessary without requiring the overhead of
       * copying the data.
       */
      void assemble_volume_of_fluid_system (const VolumeOfFluidField<dim> &field,
                                            const unsigned int calculation_dim,
                                            const bool update_from_old_solution);

      /**
       * Solve the diagonal matrix assembled in assemble_volume_of_fluid_system for the
       * specified field.
       */
      void solve_volume_of_fluid_system (const VolumeOfFluidField<dim> &field);


    private:
      /**
       * Parent simulator
       */
      Simulator<dim> &sim;

      /**
       * Function to copy assembled data to final system. Requires access to
       * the full matrix, so must be in this class.
       */
      void copy_local_to_global_volume_of_fluid_system (const internal::Assembly::CopyData::VolumeOfFluidSystem<dim> &data);

      /**
       * Assembler object used for doing the matrix and RHS assembly
       */
      Assemblers::VolumeOfFluidAssembler<dim> assembler;

      /**
       * Number of volume of fluid fields to calculate for
       */
      unsigned int n_volume_of_fluid_fields;

      /**
       * Structures containing the locations of the associated state data for
       * each volume of fluid field.
       */
      std::vector<VolumeOfFluidField<dim>> data;

      /**
       * Volume fraction threshold for the reconstruction and advection
       * algorithms indicating minimum relevant volume fraction.
       */
      double volume_fraction_threshold;

      /**
       * Tolerance to use for the Newton iteration in the reconstruction step
       */
      static constexpr double volume_of_fluid_reconstruct_epsilon = 1e-13;

      /**
       * Tolerance to use for the matrix solve in the timestep update
       */
      double volume_of_fluid_solver_tolerance;

      /**
       * Number of samples in each dimension to use during the Volume of Fluid
       * initialization, for a total of $n_init_samples^dim$ points sampled
       */
      unsigned int n_init_samples;

      /**
       * Vector of human readable names for the volume of fluid fields,
       * obtained from the associated composition field name
       */
      std::vector<std::string> volume_of_fluid_field_names;

      /**
       * Map relating the index of a Volume of Fluid based composition field to the index of the corresponding Volume of Fluid field.
       */
      std::map<unsigned int, unsigned int> volume_of_fluid_composition_map_index;

      /**
       * Methods to use when initializing the volume of fluid fields.  Must be
       * held here as all access to the actual methods for composition
       * initialization is handled by the manager.
       */
      std::vector<VolumeOfFluid::VolumeOfFluidInputType::Kind> initialization_data_type;

      friend class Simulator<dim>;
  };

}

#endif
