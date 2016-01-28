/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#ifndef __aspect__assembly_h
#define __aspect__assembly_h

#include <deal.II/fe/fe_system.h>

#include <aspect/global.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  using namespace dealii;

  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim>      struct StokesPreconditioner;
        template <int dim>      struct StokesSystem;
        template <int dim>      struct AdvectionSystem;
      }

      namespace CopyData
      {
        template <int dim>      struct StokesPreconditioner;
        template <int dim>      struct StokesSystem;
        template <int dim>      struct AdvectionSystem;
      }

      /**
       * A namespace for the definition of sub-namespaces in which we
       * implement assemblers for the various terms in the linear systems
       * ASPECT solves.
       */
      namespace Assemblers
      {
        /**
         * A base class for objects that implement assembly
         * operations. This base class does not provide a whole lot
         * of functionality other than the fact that its destructor
         * is virtual.
         *
         * The point of this class is primarily so that we can store
         * pointers to such objects in a list. The objects are created
         * in Simulator::set_assemblers() (which is the only place where
         * we know their actual type) and destroyed in the destructor of
         * class Simulation.
         */
        template <int dim>
        class AssemblerBase
        {
          public:
            virtual ~AssemblerBase () {}
        };
      }

      /**
       * A structure that consists of member variables representing
       * each one or more functions that need to be called when
       * assembling right hand side vectors, matrices, or complete linear
       * systems. We use this approach in order to support the following
       * cases:
       * - Assembling different formulations: When assembling either the
       *   full equations or only the Boussinesq approximation (to give just
       *   two examples), one needs different terms. This could be achieved
       *   using a large number of <code>switch</code> or <code>if</code>
       *   statements in the code, or one could encapsulate each equation
       *   or approximation into a collection of functions for this particular
       *   purpose. The approach chosen here in essence allows the
       *   implementation of each set of equations in its own scope, and we
       *   then just need to store a pointer to the function that
       *   assembles the Stokes system (for example) for the selected
       *   approximation. The pointer to this function is stored in the
       *   appropriate member variable of this class.
       * - Sometimes, we want to assemble a number of terms that build on
       *   each other. An example is the addition of free boundary terms
       *   to the Stokes matrix. Rather than having to "know" in one
       *   place about all of the terms that need to be assembled,
       *   we simply add the function that computes these terms as
       *   another "slot" to the appropriate "signal" declared in this
       *   class.
       */
      template <int dim>
      struct AssemblerLists
      {
        /**
         * A signal that is called from Simulator::local_assemble_stokes_preconditioner()
         * and whose slots are supposed to assemble terms that together form the
         * Stokes preconditioner matrix.
         *
         * The arguments to the slots are as follows:
         * - The Simulator::pressure_scaling value used to scale velocity
         *   and pressure components against each other.
         * - The scratch object in which temporary data is stored that
         *   assemblers may need.
         * - The copy object into which assemblers add up their contributions.
         */
        boost::signals2::signal<void (const double,
                                      internal::Assembly::Scratch::StokesPreconditioner<dim>  &,
                                      internal::Assembly::CopyData::StokesPreconditioner<dim> &)> local_assemble_stokes_preconditioner;

        /**
         * A signal that is called from Simulator::local_assemble_stokes_system()
         * and whose slots are supposed to assemble terms that together form the
         * Stokes system matrix and right hand side.
         *
         * The arguments to the slots are as follows:
         * - The cell on which we currently assemble.
         * - The Simulator::pressure_scaling value used to scale velocity
         *   and pressure components against each other.
         * - Whether or not to actually assemble the matrix. If @p false,
         *   then only assemble the right hand side of the Stokes system.
         * - The scratch object in which temporary data is stored that
         *   assemblers may need.
         * - The copy object into which assemblers add up their contributions.
         */
        boost::signals2::signal<void (const typename DoFHandler<dim>::active_cell_iterator &,
                                      const double,
                                      const bool,
                                      internal::Assembly::Scratch::StokesSystem<dim>       &,
                                      internal::Assembly::CopyData::StokesSystem<dim>      &)> local_assemble_stokes_system;

        /**
         * A signal that is called from Simulator::local_assemble_stokes_system()
         * and whose slots are supposed to assemble terms that together form the
         * Stokes system matrix and right hand side. This signal is called
         * once for each boundary face.
         *
         * The arguments to the slots are as follows:
         * - The cell on which we currently assemble.
         * - The number of the face on which we intend to assemble. This
         *   face (of the current cell) will be at the boundary of the
         *   domain.
         * - The Simulator::pressure_scaling value used to scale velocity
         *   and pressure components against each other.
         * - Whether or not to actually assemble the matrix. If @p false,
         *   then only assemble the right hand side of the Stokes system.
         * - The scratch object in which temporary data is stored that
         *   assemblers may need.
         * - The copy object into which assemblers add up their contributions.
         */
        boost::signals2::signal<void (const typename DoFHandler<dim>::active_cell_iterator &,
                                      const unsigned int,
                                      const double,
                                      const bool,
                                      internal::Assembly::Scratch::StokesSystem<dim>       &,
                                      internal::Assembly::CopyData::StokesSystem<dim>      &)> local_assemble_stokes_system_on_boundary_face;

        /**
         * A structure that describes what information an assembler function
         * (listed as one of the signals/slots above) may need to operate.
         *
         * There are a number of pieces of information that are always
         * assumed to be needed. For example, the Stokes and advection
         * assemblers will always need to have access to the material
         * model outputs. But the Stokes assembler may or may not need
         * access to material model output for quadrature points on faces.
         *
         * These properties are all preset in a conservative way
         * (i.e., disabled) in the constructor of this class, but can
         * be enabled in Simulator::set_assemblers() when adding
         * individual assemblers. Functions such as
         * Simulator::local_assemble_stokes_preconditioner(),
         * Simulator::local_assemble_stokes_system() will then query
         * these flags to determine whether something has to be
         * initialized for at least one of the assemblers they call.
         */
        struct Properties
        {
          /**
           * Constructor. Disable all properties as described in the
           * class documentation.
           */
          Properties ();

          /**
           * Whether or not at least one of the the assembler slots in
           * a signal require the initialization and re-computation of
           * a MaterialModelOutputs object for each face. This
           * property is only relevant to assemblers that operate on
           * boundary faces.
           */
          bool need_face_material_model_data;

          /**
           * A list of FEValues UpdateFlags that are necessary for
           * a given operation. Assembler objects may add to this list
           * as necessary; it will be initialized with a set of
           * "default" flags that will always be set.
           */
          UpdateFlags needed_update_flags;
        };

        /**
         * A list of properties of the various types of assemblers.
         * These property lists are set in Simulator::set_assemblers()
         * where we add individual functions to the signals above.
         */
        Properties stokes_preconditioner_assembler_properties;
        Properties stokes_system_assembler_properties;
        Properties stokes_system_assembler_on_boundary_face_properties;
      };

    }
  }
}

#endif
