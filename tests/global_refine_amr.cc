/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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



#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>
#include <math.h>


using namespace dealii;


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
          if (cell->is_locally_owned() && cell->center()[0]<0.5)
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
