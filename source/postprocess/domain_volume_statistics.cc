// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2023 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/domain_volume_statistics.h>

#include <deal.II/grid/grid_tools.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    DomainVolume<dim>::execute (TableHandler &statistics)
    {
      // Retrieve the current domain volume
      const double global_volume = this->get_volume();

      // add the volume to the statistics object
      const std::string unit = (dim == 2) ? "m^2" : "m^3";
      const std::string name = "Model domain volume (" + unit + ")";
      statistics.add_value (name, global_volume);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      statistics.set_precision (name, 8);
      statistics.set_scientific (name, true);

      // create a single string to output to the screen
      std::ostringstream screen_text;
      screen_text.precision(4);
      screen_text << global_volume << " " << unit;

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
                                  "A postprocessor that computes the total area (in 2d) "
                                  "or volume (in 3d) of the computational domain. ")
  }
}
