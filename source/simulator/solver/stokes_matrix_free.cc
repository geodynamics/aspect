// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/simulator.h>
#include <aspect/simulator/solver/stokes_matrix_free.h>
#include <aspect/simulator/solver/matrix_free_operators.h>
#include <aspect/simulator/solver/stokes_matrix_free_local_smoothing.h>
#include <aspect/simulator/solver/stokes_matrix_free_global_coarsening.h>

namespace aspect
{
  template <int dim>
  void
  StokesMatrixFreeHandler<dim>::setup_multigrid_hierarchy()
  {}

  template <int dim>
  std::unique_ptr<StokesMatrixFreeHandler<dim>> create_matrix_free_solver(Simulator<dim> &simulator, const Parameters<dim> &parameters)
  {
    if (parameters.stokes_gmg_type == Parameters<dim>::StokesGMGType::local_smoothing)
      {
        switch (parameters.stokes_velocity_degree)
          {
            case 2:
              return std::make_unique<StokesMatrixFreeHandlerLocalSmoothingImplementation<dim,2>>(simulator, parameters);
              break;
            case 3:
              return std::make_unique<StokesMatrixFreeHandlerLocalSmoothingImplementation<dim,3>>(simulator, parameters);
              break;
            default:
              AssertThrow(false, ExcMessage("The finite element degree for the Stokes system you selected is not supported yet."));
          }
      }
    else if (parameters.stokes_gmg_type == Parameters<dim>::StokesGMGType::global_coarsening)
      {
        switch (parameters.stokes_velocity_degree)
          {
            case 2:
              return std::make_unique<StokesMatrixFreeHandlerGlobalCoarseningImplementation<dim,2>>(simulator, parameters);
              break;
            case 3:
              return std::make_unique<StokesMatrixFreeHandlerGlobalCoarseningImplementation<dim,3>>(simulator, parameters);
              break;
            default:
              AssertThrow(false, ExcMessage("The finite element degree for the Stokes system you selected is not supported yet."));
          }
      }
    else
      AssertThrow(false, ExcNotImplemented());
  }

  template <int dim>
  void StokesMatrixFreeHandler<dim>::declare_parameters(ParameterHandler &prm)
  {
    StokesMatrixFreeHandlerLocalSmoothingImplementation<dim,2>::declare_parameters(prm);
    StokesMatrixFreeHandlerGlobalCoarseningImplementation<dim,2>::declare_parameters(prm);
  }

}

// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template std::unique_ptr<StokesMatrixFreeHandler<dim>> create_matrix_free_solver(Simulator<dim> &, const Parameters<dim> &); \
  template class StokesMatrixFreeHandler<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
