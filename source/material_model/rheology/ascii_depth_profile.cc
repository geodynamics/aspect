// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2020-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
        = default;



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
    namespace Rheology
    {
#define INSTANTIATE(dim) \
  template class AsciiDepthProfile<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
