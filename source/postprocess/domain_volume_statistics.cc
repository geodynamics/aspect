/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/postprocess/domain_volume_statistics.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    DomainVolume<dim>::execute (TableHandler &statistics)
    {
      const double global_volume = GridTools::volume (this->get_triangulation(), this->get_mapping());

      // add the volume to the statistics object
      // and create a single string that can be output to the screen
      std::ostringstream screen_text;
      const std::string name = "Model domain volume (m3) ";
      statistics.add_value (name, global_volume);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      statistics.set_precision (name, 8);
      statistics.set_scientific (name, true);

      // print to the screen
      screen_text.precision(4);
      screen_text << global_volume << " m3";

      return std::pair<std::string, std::string> ("Model domain volume:",
                                                  screen_text.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(DomainVolume,
                                  "domain volume statistics",
                                  "A postprocessor that computes the total volume "
                                  "of the computational domain. ")
  }
}
