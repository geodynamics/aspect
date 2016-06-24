for (unsigned int q=0; q < in.position.size(); ++q)
  {
    // rotation tensor =
    //     asymmetric part of the displacement in this time step
    //                (= velocity gradient tensor * time step)
    //   + unit tensor
    const Tensor<2,dim> rotation = (velocity_gradients[q] - symmetrize(velocity_gradients[q]))
                                   * this->get_timestep()
                                   + unit_symmetric_tensor<dim>();

    SymmetricTensor<2,dim> accumulated_strain;
    for (unsigned int i=0; i<SymmetricTensor<2,dim>::n_independent_components; ++i)
      accumulated_strain[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)] = in.composition[q][i];

    // the new strain is the rotated old strain plus the strain of the
    // current time step
    const SymmetricTensor<2,dim> rotated_strain = symmetrize(rotation * Tensor<2,dim>(accumulated_strain) * transpose(rotation))
                                                  + in.strain_rate[q] * this->get_timestep();

    for (unsigned int c=0; c<SymmetricTensor<2,dim>::n_independent_components; ++c)
      {
        out.reaction_terms[q][c] = - in.composition[q][c]
                                   + rotated_strain[SymmetricTensor<2,dim>::unrolled_to_component_indices(c)];
      }
  }
