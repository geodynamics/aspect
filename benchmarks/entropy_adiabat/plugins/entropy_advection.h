/*
  Copyright (C) 2016 - 2021 by the authors of the ASPECT code.

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

#ifndef _aspect_simulator_assemblers_entropy_advection_h
#define _aspect_simulator_assemblers_entropy_advection_h


#include <aspect/simulator/assemblers/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Assemblers
  {
    /**
     * This class assembles the terms for the matrix and right-hand-side of the entropy
     * advection equation for the current cell.
     */
    template <int dim>
    class EntropyAdvectionSystem : public Assemblers::Interface<dim>, public Assemblers::AdvectionStabilizationInterface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;

        virtual
        std::vector<double>
        compute_residual(internal::Assembly::Scratch::ScratchBase<dim>  &scratch) const;

        std::vector<double>
        advection_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const override;

        std::vector<double>
        diffusion_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const override;
    };
  }
}


#endif
