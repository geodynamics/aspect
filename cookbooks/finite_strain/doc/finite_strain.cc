/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

for (unsigned int q=0; q < in.position.size(); ++q)
  {
    // Convert the compositional fields into the tensor quantity they represent.
    Tensor<2,dim> strain;
    for (unsigned int i = 0; i < Tensor<2,dim>::n_independent_components ; ++i)
      strain[Tensor<2,dim>::unrolled_to_component_indices(i)] = in.composition[q][i];

    // Compute the strain accumulated in this timestep.
    const Tensor<2,dim> strain_increment = this->get_timestep() * (velocity_gradients[q] * strain);

    // Output the strain increment component-wise to its respective compositional field's reaction terms.
    for (unsigned int i = 0; i < Tensor<2,dim>::n_independent_components ; ++i)
      out.reaction_terms[q][i] = strain_increment[Tensor<2,dim>::unrolled_to_component_indices(i)];
  }
