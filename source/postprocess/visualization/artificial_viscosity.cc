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


#include <aspect/postprocess/visualization/artificial_viscosity.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      ArtificialViscosity<dim>::
      ArtificialViscosity ()
        :
        CellDataVectorCreator<dim>("W/m/K")
      {}



      template <int dim>
      std::pair<std::string, std::unique_ptr<Vector<float>>>
      ArtificialViscosity<dim>::execute() const
      {
        std::pair<std::string, std::unique_ptr<Vector<float>>>
        return_value ("artificial_viscosity",
                      std::make_unique<Vector<float>>(this->get_triangulation().n_active_cells()));
        this->get_artificial_viscosity(*return_value.second);

        // The function we call above sets the artificial viscosity to
        // signaling_nan on all artificial cells and, possibly, ghost cells.
        // This runs into trouble in DataOut that wants to copy this vector
        // from Vector<float> to Vector<double>, and the conversion trips
        // up over the NaNs, causing a floating point exception.
        //
        // To avoid this, strip out the NaNs and instead set these values
        // to zero -- we won't be outputting these values anyway.
        for (const auto &cell : this->get_triangulation().active_cell_iterators())
          if (cell->is_locally_owned() == false)
            (*return_value.second)[cell->active_cell_index()] = 0;

        return return_value;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(ArtificialViscosity,
                                                  "artificial viscosity",
                                                  "A visualization output object that generates output "
                                                  "showing the value of the artificial viscosity on each "
                                                  "cell."
                                                  "\n\n"
                                                  "Physical units: \\si{\\watt\\per\\meter\\per\\kelvin}.")
    }
  }
}
