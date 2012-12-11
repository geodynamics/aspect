/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/postprocess/visualization/seismic_vp_anomaly.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      std::pair<std::string, Vector<float> *>
      SeismicVpAnomaly<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("Vp_anomaly",
                      new Vector<float>(this->get_triangulation().n_active_cells()));
        this->get_Vp_anomaly(*return_value.second);

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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SeismicVpAnomaly,
                                                  "Vp anomaly",
                                                  "A visualization output object that generates output "
                                                  "showing the anomaly in the seismic compression wave "
                                                  "speed $V_p$ as a spatially variable function with one "
                                                  "value per cell. This anomaly is shown as a percentage "
                                                  "change relative to the average value of $V_p$ at "
                                                  "the depth of this cell.")
    }
  }
}
