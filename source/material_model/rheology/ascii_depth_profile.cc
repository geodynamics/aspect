/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/material_model/rheology/ascii_depth_profile.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {

      template <int dim>
      AsciiDepthProfile<dim>::AsciiDepthProfile ()
      {}



      template <int dim>
      void
      AsciiDepthProfile<dim>::initialize ()
      {
        this->initialize(this->get_mpi_communicator());
        viscosity_index = this->get_column_index_from_name("viscosity");
      }



      template <int dim>
      double
      AsciiDepthProfile<dim>::compute_viscosity (const double depth) const
      {
        return this->get_data_component (Point<1>(depth), viscosity_index);
      }



      template <int dim>
      void
      AsciiDepthProfile<dim>::declare_parameters (ParameterHandler &prm,
                                                  const std::string &subsection_name)
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/material-model/rheology/",
                                                          "ascii_depth_profile.txt",
                                                          subsection_name);
      }



      template <int dim>
      void
      AsciiDepthProfile<dim>::parse_parameters (ParameterHandler &prm,
                                                const std::string &subsection_name)
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm,
                                                        subsection_name);
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class AsciiDepthProfile<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
