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
