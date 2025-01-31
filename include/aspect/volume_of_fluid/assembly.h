/*
 Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_volume_of_fluid_assembly_h
#define _aspect_volume_of_fluid_assembly_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>

#include <deal.II/fe/fe_values.h>

namespace aspect
{
  template <int dim>
  class Simulator;

  template <int dim>
  struct VolumeOfFluidField;

  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        /**
         * Standard scratch data structure for matrix assembly
         */
        template <int dim>
        struct VolumeOfFluidSystem
        {
          VolumeOfFluidSystem (const FiniteElement<dim> &finite_element,
                               const FiniteElement<dim> &volume_of_fluid_element,
                               const Mapping<dim>       &mapping,
                               const Quadrature<dim>    &quadrature,
                               const Quadrature<dim-1>  &face_quadrature);

          VolumeOfFluidSystem (const VolumeOfFluidSystem &scratch);

          /**
           * Add a defaulted assignment operator because relying on it
           * implicitly is deprecated.
           */
          VolumeOfFluidSystem &operator=(const VolumeOfFluidSystem &) = default;

          FEValues<dim>          finite_element_values;
          FEValues<dim>          neighbor_finite_element_values;
          FEFaceValues<dim>      face_finite_element_values;
          FEFaceValues<dim>      neighbor_face_finite_element_values;
          FESubfaceValues<dim>   subface_finite_element_values;

          std::vector<types::global_dof_index>   local_dof_indices;

          // Field for exact cell volume (needed in some of the face flux calculations)
          double volume;

          /**
           * Variables describing the values of the shape functions at the
           * quadrature points, as they are used in the advection assembly
           * function. note that the sizes of these arrays are equal to the
           * number of shape functions corresponding to a single VolumeOfFluid field (and
           * not of all VolumeOfFluid fields!), and that they are also correspondingly
           * indexed.
           */
          std::vector<double>         phi_field;
          std::vector<double>         face_phi_field;

          std::vector<double>         old_field_values;
          /* Vector for interface normal in the unit cell */
          std::vector<Tensor<1,dim>> cell_i_n_values;
          /* "Distance" from cell center to interface as d value for the interface in the form $\vec{n}\cdot\vec{x}=d$ */
          std::vector<double>         cell_i_d_values;

          std::vector<Tensor<1,dim>> face_current_velocity_values;
          std::vector<Tensor<1,dim>> face_old_velocity_values;
          std::vector<Tensor<1,dim>> face_old_old_velocity_values;

          std::vector<double>         neighbor_old_values;
          /* Vector for interface normal in the unit cell */
          std::vector<Tensor<1,dim>> neighbor_i_n_values;
          /* "Distance" from cell center to interface as d value for the interface in the form $\vec{n}\cdot\vec{x}=d$ */
          std::vector<double>         neighbor_i_d_values;
        };
      }

      namespace CopyData
      {
        /**
         * Standard copy data structure for matrix assembly
         */
        template <int dim>
        struct VolumeOfFluidSystem
        {
          /**
           * Constructor.
           * @param finite_element The element that describes the field for which we
           *    are trying to assemble a linear system. <b>Not</b> the global finite
           *    element.
           */
          VolumeOfFluidSystem(const FiniteElement<dim> &finite_element);
          VolumeOfFluidSystem(const VolumeOfFluidSystem &data);

          /**
           * Add a defaulted assignment operator because relying on it
           * implicitly is deprecated.
           */
          VolumeOfFluidSystem &operator=(const VolumeOfFluidSystem &) = default;

          /**
           * Local contributions to the global matrix and right hand side
           * that correspond only to the variables listed in local_dof_indices
           */
          FullMatrix<double>          local_matrix;
          Vector<double>              local_rhs;

          /**
           * Local contributions to the global rhs from the face terms in the
           * discontinuous Galerkin interpretation of the VolumeOfFluid method.
           *
           * The array has a length sufficient to hold one element for each
           * possible face and sub-face of a cell.
           */
          std::vector<Vector<double>> local_face_rhs;
          std::vector<FullMatrix<double>> local_face_matrices_ext_ext;

          /**
           * Denotes which face's rhs have actually been assembled in the DG
           * field assembly. Entries not used (for example, those corresponding
           * to non-existent subfaces; or faces being assembled by the
           * neighboring cell) are set to false.
           *
           * The array has a length sufficient to hold one element for each
           * possible face and sub-face of a cell.
           */
          std::vector<bool> face_contributions_mask;

          /**
           * Indices of those degrees of freedom that actually correspond to
           * the volume_of_fluid field. Since this structure is used to represent just
           * contributions to the volume_of_fluid systems, there will be no contributions
           * to other parts of the system and consequently, we do not need to
           * list here indices that correspond to velocity or pressure degrees
           * (or, in fact any other variable outside the block we are currently
           * considering)
           */
          std::vector<types::global_dof_index>   local_dof_indices;

          /**
           * Indices of the degrees of freedom corresponding to the volume_of_fluid field
           * on all possible neighboring cells. This is used in the
           * discontinuous Galerkin interpretation of the VolumeOfFluid method.
           *
           * The array has a length sufficient to hold one element for each
           * possible face and sub-face of a cell.
           */
          std::vector<std::vector<types::global_dof_index>> neighbor_dof_indices;
        };
      }
    }
  }

  namespace Assemblers
  {

    /**
     * Class to hold VolumeOfFluid assembly logic, as analogous to that used in the main simulator.
     */
    template <int dim>
    class VolumeOfFluidAssembler : public SimulatorAccess<dim>
    {
      public:
        /**
         * Do setup and assembly on internal quadrature points and dispatch to
         * other functions for face assembly
         */
        void local_assemble_volume_of_fluid_system (const VolumeOfFluidField<dim> &field,
                                                    const unsigned int calc_dir,
                                                    const bool update_from_old,
                                                    const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                    internal::Assembly::Scratch::VolumeOfFluidSystem<dim> &scratch,
                                                    internal::Assembly::CopyData::VolumeOfFluidSystem<dim> &data) const;

        /**
         * Do assembly for cell faces on the boundary
         */
        void local_assemble_boundary_face_volume_of_fluid_system (const VolumeOfFluidField<dim> &field,
                                                                  const bool update_from_old,
                                                                  const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                                  const unsigned int face_no,
                                                                  internal::Assembly::Scratch::VolumeOfFluidSystem<dim> &scratch,
                                                                  internal::Assembly::CopyData::VolumeOfFluidSystem<dim> &data) const;

        /**
         * Function for assembling face fluxes for VolumeOfFluid system.
         */
        void local_assemble_internal_face_volume_of_fluid_system (const VolumeOfFluidField<dim> &field,
                                                                  const bool update_from_old,
                                                                  const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                                  const unsigned int face_no,
                                                                  internal::Assembly::Scratch::VolumeOfFluidSystem<dim> &scratch,
                                                                  internal::Assembly::CopyData::VolumeOfFluidSystem<dim> &data) const;

        /**
         * Set volume fraction threshold for use in assembly
         */
        void set_volume_fraction_threshold(const double value);
    };
  }
}


#endif
