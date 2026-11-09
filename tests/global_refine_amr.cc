// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------



#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>
#include <math.h>


namespace aspect
{
  template <int dim>
  class AMRLeft : public MeshRefinement::Interface<dim>,
    public SimulatorAccess<dim>
  {
    public:
      /**
       * After cells have been marked for coarsening/refinement, apply
       * additional criteria independent of the error estimate.
       *
       */
      virtual
      void
      tag_additional_cells () const;

      /**
       * Declare the parameters this class takes through input files.
       */
      static
      void
      declare_parameters (ParameterHandler &prm);

      /**
       * Read the parameters this class declares from the parameter
       * file.
       */
      virtual
      void
      parse_parameters (ParameterHandler &prm);

    private:
  };


  template <int dim>
  void
  AMRLeft<dim>::tag_additional_cells () const
  {
    std::cout << "plugin is called!" << std::endl;

    for (typename Triangulation<dim>::active_cell_iterator
         cell = this->get_triangulation().begin_active();
         cell != this->get_triangulation().end(); ++cell)
      {
        if (cell->is_locally_owned() && cell->center(true)[0]<0.5)
          cell->clear_refine_flag ();
      }
  }

  template <int dim>
  void
  AMRLeft<dim>::
  declare_parameters (ParameterHandler &prm)
  {
  }

  template <int dim>
  void
  AMRLeft<dim>::parse_parameters (ParameterHandler &prm)
  {
  }

}




// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(AMRLeft,
                                            "AMRLeft",
                                            "TODO")
}
