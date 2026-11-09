// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
