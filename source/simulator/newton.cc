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

    // add the terms for traction boundary conditions
    if (!this->get_boundary_traction().empty())
      {
        assemblers.stokes_system_on_boundary_face.push_back(
          std_cxx14::make_unique<aspect::Assemblers::StokesBoundaryTraction<dim> >());
      }

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

  template <int dim>
  std::string
  NewtonHandler<dim>::
  get_preconditioner_stabilization() const
  {
    return preconditioner_stabilization;
  }


  template <int dim>
  void
  NewtonHandler<dim>::
  set_preconditioner_stabilization(const std::string preconditioner_stabilization_)
  {
    preconditioner_stabilization = preconditioner_stabilization_;
  }

  template <int dim>
  std::string
  NewtonHandler<dim>::
  get_velocity_block_stabilization() const
  {
    return velocity_block_stabilization;
  }


  template <int dim>
  void
  NewtonHandler<dim>::
  set_velocity_block_stabilization(const std::string velocity_block_stabilization_)
  {
    velocity_block_stabilization = velocity_block_stabilization_;
  }

  template <int dim>
  bool
  NewtonHandler<dim>::
  get_use_Newton_failsafe()
  {
    return use_Newton_failsafe;
  }

  template <int dim>
  double
  NewtonHandler<dim>::
  get_nonlinear_switch_tolerance()
  {
    return nonlinear_switch_tolerance;
  }

  template <int dim>
  unsigned int
  NewtonHandler<dim>::
  get_max_pre_newton_nonlinear_iterations()
  {
    return max_pre_newton_nonlinear_iterations;
  }

  template <int dim>
  unsigned int
  NewtonHandler<dim>::
  get_max_newton_line_search_iterations()
  {
    return max_newton_line_search_iterations;
  }

  template <int dim>
  bool
  NewtonHandler<dim>::
  get_use_newton_residual_scaling_method()
  {
    return use_newton_residual_scaling_method;
  }

  template <int dim>
  double
  NewtonHandler<dim>::
  get_maximum_linear_stokes_solver_tolerance()
  {
    return maximum_linear_stokes_solver_tolerance;
  }

  template <int dim>
  void
  NewtonHandler<dim>::
  declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("Solver parameters");
    {
      prm.enter_subsection ("Newton solver parameters");
      {
        prm.declare_entry ("Nonlinear Newton solver switch tolerance", "1e-5",
                           Patterns::Double(0,1),
                           "A relative tolerance with respect to the residual of the first "
                           "iteration, up to which the nonlinear Picard solver will iterate, "
                           "before changing to the newton solver.");

        prm.declare_entry ("Max pre-Newton nonlinear iterations", "10",
                           Patterns::Integer (0),
                           "The maximum number of Picard nonlinear iterations to be performed "
                           "before switching to Newton iterations.");

        prm.declare_entry ("Max Newton line search iterations", "5",
                           Patterns::Integer (0),
                           "The maximum number of line search iterations allowed. If the "
                           "criterion is not reached after this iteration, we apply the scaled "
                           "increment to the solution and continue.");

        prm.declare_entry ("Use Newton residual scaling method", "false",
                           Patterns::Bool (),
                           "This method allows to slowly introduce the derivatives based on the improvement "
                           "of the residual. If set to false, the scaling factor for the Newton derivatives "
                           "is set to one immediately when switching on the Newton solver.");

        prm.declare_entry ("Maximum linear Stokes solver tolerance", "0.9",
                           Patterns::Double (0,1),
                           "When the linear Stokes solver tolerance is dynamically chosen, this defines "
                           "the most loose tolerance allowed.");

        prm.declare_entry ("Stabilization preconditioner", "SPD",
                           Patterns::Selection ("SPD|PD|symmetric|none"),
                           "TODO");
        prm.declare_entry ("Stabilization velocity block", "SPD",
                           Patterns::Selection ("SPD|PD|symmetric|none"),
                           "TODO");
        prm.declare_entry ("Use Newton failsafe", "false",
                           Patterns::Bool (),
                           "Switches on SPD stabilization when solver fails.");
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();
  }

  template <int dim>
  void
  NewtonHandler<dim>::
  parse_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("Solver parameters");
    {
      prm.enter_subsection ("Newton solver parameters");
      {
        nonlinear_switch_tolerance = prm.get_double("Nonlinear Newton solver switch tolerance");
        max_pre_newton_nonlinear_iterations = prm.get_integer ("Max pre-Newton nonlinear iterations");
        max_newton_line_search_iterations = prm.get_integer ("Max Newton line search iterations");
        use_newton_residual_scaling_method = prm.get_bool("Use Newton residual scaling method");
        maximum_linear_stokes_solver_tolerance = prm.get_double("Maximum linear Stokes solver tolerance");
        preconditioner_stabilization = prm.get("Stabilization preconditioner");
        velocity_block_stabilization = prm.get("Stabilization velocity block");
        use_Newton_failsafe = prm.get_bool("Use Newton failsafe");

        AssertThrow((!DEAL_II_VERSION_GTE(9,0,0) && !use_Newton_failsafe) || DEAL_II_VERSION_GTE(9,0,0),
                    ExcMessage("The failsafe option can't be used with a deal.ii less then 9.0.0."));
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();

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
