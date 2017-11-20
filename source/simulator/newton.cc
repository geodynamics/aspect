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
#include <aspect/simulator/assemblers/advection.h>
#include <aspect/simulator/assemblers/stokes.h>

#include <aspect/simulator.h>

#if DEAL_II_VERSION_GTE(9,0,0)
#include <deal.II/base/std_cxx14/memory.h>
#endif


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
  set_assemblers (Assemblers::Manager<dim> &assemblers) const
  {
    assemblers.stokes_preconditioner.push_back(std_cxx14::make_unique<aspect::Assemblers::NewtonStokesPreconditioner<dim> >());
    assemblers.stokes_system.push_back(std_cxx14::make_unique<aspect::Assemblers::NewtonStokesIncompressibleTerms<dim> >());

    if (this->get_material_model().is_compressible())
      {
        assemblers.stokes_preconditioner.push_back(
          std_cxx14::make_unique<aspect::Assemblers::StokesCompressiblePreconditioner<dim> >());

        assemblers.stokes_system.push_back(
          std_cxx14::make_unique<aspect::Assemblers::NewtonStokesCompressibleStrainRateViscosityTerm<dim> >());
      }

    if (this->get_parameters().formulation_mass_conservation ==
        Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile)
      {
        assemblers.stokes_system.push_back(
          std_cxx14::make_unique<aspect::Assemblers::NewtonStokesImplicitReferenceDensityCompressibilityTerm<dim> >());
      }
    else if (this->get_parameters().formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::reference_density_profile)
      {
        assemblers.stokes_system.push_back(
          std_cxx14::make_unique<aspect::Assemblers::NewtonStokesReferenceDensityCompressibilityTerm<dim> >());
      }
    else if (this->get_parameters().formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::incompressible)
      {
        // do nothing, because we assembled div u =0 above already
      }
    else if (this->get_parameters().formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::isothermal_compression)
      {
        assemblers.stokes_system.push_back(
          std_cxx14::make_unique<aspect::Assemblers::NewtonStokesIsothermalCompressionTerm<dim> >());
      }
    else
      AssertThrow(false,
                  ExcMessage("Unknown mass conservation equation approximation. There is no assembler"
                             " defined that handles this formulation."));

    assemblers.stokes_system_on_boundary_face.push_back(
      std_cxx14::make_unique<aspect::Assemblers::StokesBoundaryTraction<dim> >());

    // add the terms necessary to normalize the pressure
    if (this->pressure_rhs_needs_compatibility_modification())
      assemblers.stokes_system.push_back(
        std_cxx14::make_unique<aspect::Assemblers::StokesPressureRHSCompatibilityModification<dim> >());

    assemblers.advection_system.push_back(
      std_cxx14::make_unique<aspect::Assemblers::AdvectionSystem<dim> >());

    if (this->get_parameters().use_discontinuous_temperature_discretization ||
        this->get_parameters().use_discontinuous_composition_discretization)
      {
        assemblers.advection_system_on_boundary_face.push_back(
          std_cxx14::make_unique<aspect::Assemblers::AdvectionSystemBoundaryFace<dim> >());

        assemblers.advection_system_on_interior_face.push_back(
          std_cxx14::make_unique<aspect::Assemblers::AdvectionSystemInteriorFace<dim> >());
      }

    if (this->get_parameters().use_discontinuous_temperature_discretization)
      {
        assemblers.advection_system_assembler_on_face_properties[0].need_face_material_model_data = true;
        assemblers.advection_system_assembler_on_face_properties[0].need_face_finite_element_evaluation = true;
      }

    if (this->get_parameters().use_discontinuous_composition_discretization)
      {
        for (unsigned int i = 1; i<=this->introspection().n_compositional_fields; ++i)
          {
            assemblers.advection_system_assembler_on_face_properties[i].need_face_material_model_data = true;
            assemblers.advection_system_assembler_on_face_properties[i].need_face_finite_element_evaluation = true;
          }
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


  template <int dim>
  double
  NewtonHandler<dim>::
  get_newton_derivative_scaling_factor() const
  {
    return newton_derivative_scaling_factor;
  }


  template <int dim>
  void
  NewtonHandler<dim>::
  set_newton_derivative_scaling_factor(const double set_newton_derivative_scaling_factor)
  {
    newton_derivative_scaling_factor = set_newton_derivative_scaling_factor;
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
