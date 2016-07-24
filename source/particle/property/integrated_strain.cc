/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/particle/property/integrated_strain.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      IntegratedStrain<dim>::initialize_one_particle_property(const Point<dim> &,
                                                              std::vector<double> &data) const
      {
        for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
          data.push_back(0.0);
      }

      template <int dim>
      void
      IntegratedStrain<dim>::update_one_particle_property(const unsigned int data_position,
                                                          const Point<dim> &,
                                                          const Vector<double> &,
                                                          const std::vector<Tensor<1,dim> > &gradients,
                                                          std::vector<double> &data) const
      {
        SymmetricTensor<2,dim> old_strain;
        for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
          old_strain[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)] = data[data_position + i];

        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];

        // Note: We use the old_timestep, because this function is called as
        // postprocessor right now, and get_timestep() is already updated for
        // the next timestep.
        const double dt = this->get_old_timestep();

        // Here we update the integrated strain by rotating the already existing
        // strain with the rotational (asymmetric) part of the velocity gradient
        // tensor, and afterwards add the additionally experienced deformation.
        const SymmetricTensor<2,dim> deformation_rate (symmetrize(grad_u));

        // rotation tensor =
        //     asymmetric part of the displacement in this time step
        //                (= time step * (velocity gradient tensor - symmetric_part))
        //   + unit tensor
        const Tensor<2,dim> rotation (dt * (grad_u - deformation_rate) + unit_symmetric_tensor<dim>());

        // the new strain is the rotated old strain plus the
        // strain of the current time step
        const SymmetricTensor<2,dim> new_strain = symmetrize(rotation * Tensor<2,dim>(old_strain) * transpose(rotation)) + dt * deformation_rate;

        for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
          data[data_position + i] = new_strain[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];
      }

      template <int dim>
      UpdateTimeFlags
      IntegratedStrain<dim>::need_update() const
      {
        return update_time_step;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      IntegratedStrain<dim>::get_property_information() const
      {
        const unsigned int n_components = SymmetricTensor<2,dim>::n_independent_components;
        const std::vector<std::pair<std::string,unsigned int> > property_information (1,std::make_pair("integrated strain",n_components));
        return property_information;
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(IntegratedStrain,
                                        "integrated strain",
                                        "A plugin in which the tracer property tensor is "
                                        "defined as the integrated strain this "
                                        "particle has experienced.")
    }
  }
}

