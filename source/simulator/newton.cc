/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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


#include <aspect/newton.h>
#include <aspect/simulator.h>



namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    MaterialModelDerivatives<dim>::
    MaterialModelDerivatives (const unsigned int n_points)
    {
      viscosity_derivative_wrt_pressure.resize(n_points, numbers::signaling_nan<double>());
      viscosity_derivative_wrt_strain_rate.resize(n_points, numbers::signaling_nan<SymmetricTensor<2,dim> >());
    }

  }



  template <int dim>
  void
  NewtonHandler<dim>::
  create_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output)
  {
    if (output.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >() != NULL)
      return;

    const unsigned int n_points = output.viscosities.size();
    output.additional_outputs.push_back(
      std_cxx11::shared_ptr<MaterialModel::AdditionalMaterialOutputs<dim> >
      (new MaterialModel::MaterialModelDerivatives<dim> (n_points)));
  }
}




// explicit instantiation of the functions we implement in this file
namespace aspect
{

#define INSTANTIATE(dim) \
  \
  template \
  class \
  NewtonHandler<dim>; \
  \
  namespace MaterialModel \
  { \
    template \
    class \
    MaterialModelDerivatives<dim>; \
  }

  ASPECT_INSTANTIATE(INSTANTIATE)

}
