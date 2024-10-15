/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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

for (unsigned int q=0; q < in.n_evaluation_points(); ++q)
  {
    // Convert the compositional fields into the tensor quantity they represent.
    const Tensor<2,dim> strain(make_array_view(&in.composition[q][0],
                                               &in.composition[q][0] + Tensor<2,dim>::n_independent_components));

    // Compute the strain accumulated in this timestep.
    const Tensor<2,dim> strain_increment = this->get_timestep() * (velocity_gradients[q] * strain);

    // Output the strain increment component-wise to its respective compositional field's reaction terms.
    strain_increment.unroll(&out.reaction_terms[q][0],
                            &out.reaction_terms[q][0] + Tensor<2,dim>::n_independent_components);
  }
