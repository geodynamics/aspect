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
#include <aspect/simulator.h>
#include "strain_residual.h"

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      StrainResidual<dim>::
	  StrainResidual()
        :
	  DataPostprocessorScalar<dim> ("strain_residual",
        			update_values | update_gradients)
      {}

      template <int dim>
      void
	  StrainResidual<dim>::initialize ()
	  {
    	std::set<types::boundary_id> surface_boundary_set;
    	surface_boundary_set.insert(this->get_geometry_model().translate_symbolic_boundary_name_to_id("top"));
     // The input ascii table contains one data column (strain rates) in addition to the coordinate columns.
    	Utilities::AsciiDataBoundary<dim>::initialize(surface_boundary_set, 1);
	  }


      template <int dim>
      void
	  StrainResidual<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {

    	Assert ((computed_quantities[0].size() == dim), ExcInternalError());
        auto cell = input_data.template get_cell<DoFHandler<dim> >();

        // We only want to output dynamic topography at the top and bottom
        // boundary, so only compute it if the current cell has
        // a face at the top or bottom boundary.
        bool cell_at_top_boundary = false;
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->at_boundary(f) &&
              (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top"))
            cell_at_top_boundary = true;

        if (cell_at_top_boundary)
        {
		  for (unsigned int q=0; q<computed_quantities.size(); ++q)
		  {
		    for (unsigned int d = 0; d < dim; ++d)
			{
			  Tensor<2,dim> grad_u;

			  grad_u[d] = input_data.solution_gradients[q][d];
			  const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);
			  const SymmetricTensor<2,dim> compressible_strain_rate
			   = (this->get_material_model().is_compressible()
					?
					strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
					:
					strain_rate);

			  computed_quantities[q](0) = std::sqrt(compressible_strain_rate *
																  compressible_strain_rate);
		    }
		  }
        }
      }


      template <int dim>
      void
	  StrainResidual<dim>::declare_parameters (ParameterHandler &prm)
	  {
    	prm.enter_subsection("Postprocess");
    	{
          prm.enter_subsection("Visualization");
    	  {
    	     Utilities::AsciiDataBase<dim>::declare_parameters(prm,
    					  "$ASPECT_SOURCE_DIR/data/initial-temperature/adiabatic-boundary/",
						  "adiabatic_boundary.txt", "Surface properties");
    	  }
    	  prm.leave_subsection();
    	}
    	prm.leave_subsection();
	  }


      template <int dim>
      void
	  StrainResidual<dim>::parse_parameters (ParameterHandler &prm)
	  {
    	prm.enter_subsection("Postprocess");
        prm.enter_subsection("Visualization");
        Utilities::AsciiDataBase<dim>::parse_parameters(prm, "Surface properties");
        prm.leave_subsection();
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(StrainResidual,
                                                  "strain residual",
                                                  "A visualization for caculating the surface strain to compare"
												  "with the GPS strains.")
    }
  }
}
