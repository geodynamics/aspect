for (unsigned int q=0; q < in.position.size(); ++q)
  {
    // rotation tensor =
    // asymmetric part of the displacement in this time step (= velocity gradients tensor * time step)
    // + unit tensor
    const Tensor<2,dim> rotation = (velocity_gradients[q] - symmetrize(velocity_gradients[q]))
                                   * this->get_timestep()
                                   + unit_symmetric_tensor<dim>();

    Tensor<2,dim> accumulated_strain;
    for (unsigned int i=0; i<Tensor<2,dim>::n_independent_components; ++i)
      accumulated_strain[Tensor<2,dim>::unrolled_to_component_indices(i)] = in.composition[q][i];

    // the new strain is the rotated old strain plus the strain of the current time step
    const Tensor<2,dim> rotated_strain = rotation * accumulated_strain * transpose(rotation)
                                         + in.strain_rate[q] * this->get_timestep();

    for (unsigned int c=0; c<Tensor<2,dim>::n_independent_components; ++c)
      {
        out.reaction_terms[q][c] = - in.composition[q][c]
                                   + rotated_strain[Tensor<2,dim>::unrolled_to_component_indices(c)];
      }
  }
