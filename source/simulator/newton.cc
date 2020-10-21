/*
  Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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
#include <aspect/citation_info.h>

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
             Parameters<dim>::Formulation::MassConservation::isentropic_compression)
      {
        assemblers.stokes_system.push_back(
          std_cxx14::make_unique<aspect::Assemblers::NewtonStokesIsentropicCompressionTerm<dim> >());
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
  }



  template <int dim>
  void
  NewtonHandler<dim>::
  create_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output)
  {
    if (output.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >() != nullptr)
      return;

    const unsigned int n_points = output.viscosities.size();
    output.additional_outputs.push_back(
      std_cxx14::make_unique<MaterialModel::MaterialModelDerivatives<dim>>(n_points));
  }



  namespace Newton
  {
    std::string
    to_string(const Newton::Parameters::Stabilization preconditioner_stabilization)
    {
      switch (preconditioner_stabilization)
        {
          case Newton::Parameters::Stabilization::SPD:
            return "SPD";
          case Newton::Parameters::Stabilization::PD:
            return "PD";
          case Newton::Parameters::Stabilization::symmetric:
            return "symmetric";
          case Newton::Parameters::Stabilization::none:
            return "none";
          default:
            Assert(false,ExcNotImplemented());
            return "";
        }
    }



    void
    Parameters::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Solver parameters");
      {
        prm.enter_subsection ("Newton solver parameters");
        {
          prm.declare_entry ("Nonlinear Newton solver switch tolerance", "1e-5",
                             Patterns::Double(0., 1.),
                             "A relative tolerance with respect to the residual of the first "
                             "iteration, up to which the nonlinear Picard solver will iterate, "
                             "before changing to the Newton solver.");

          prm.declare_entry ("Max pre-Newton nonlinear iterations", "10",
                             Patterns::Integer (0),
                             "If the 'Nonlinear Newton solver switch tolerance' is reached before the "
                             "maximal number of Picard iterations, then the solver switches to Newton "
                             "solves anyway.");

          prm.declare_entry ("Max Newton line search iterations", "5",
                             Patterns::Integer (0),
                             "The maximum number of line search iterations allowed. If the "
                             "criterion is not reached after this number of iterations, we apply "
                             "the scaled increment even though it does not satisfy the necessary "
                             "criteria and simply continue with the next Newton iteration.");

          prm.declare_entry ("Use Newton residual scaling method", "false",
                             Patterns::Bool (),
                             "This method allows to slowly introduce the derivatives based on the improvement "
                             "of the residual. If set to false, the scaling factor for the Newton derivatives "
                             "is set to one immediately when switching on the Newton solver. When this is set to "
                             "true, the derivatives are slowly introduced by the following equation: $\\max(0.0, "
                             "(1.0-(residual/switch\\_initial\\_residual)))$, where switch\\_initial\\_residual is the "
                             "residual at the time when the Newton solver is switched on.");

          prm.declare_entry ("Maximum linear Stokes solver tolerance", "0.9",
                             Patterns::Double (0., 1.),
                             "The linear Stokes solver tolerance is dynamically chosen for the Newton solver, based "
                             "on the Eisenstat walker 1994 paper (https://doi.org/10.1137/0917003), equation 2.2. "
                             "Because this value can become larger then one, we limit this value by this parameter.");

          prm.declare_entry ("Stabilization preconditioner", "SPD",
                             Patterns::Selection ("SPD|PD|symmetric|none"),
                             "This parameters allows for the stabilization of the preconditioner. If one derives the Newton "
                             "method without any modifications, the matrix created for the preconditioning is not necessarily "
                             "Symmetric Positive Definite. This is problematic (see \\cite{FBTGS19}). When `none' is chosen, "
                             "the preconditioner is not stabilized. The `symmetric' parameters symmetrizes the matrix, and `PD' makes "
                             "the matrix Positive Definite. `SPD' is the full stabilization, where the matrix is guaranteed Symmetric "
                             "Positive Definite.");

          prm.declare_entry ("Stabilization velocity block", "SPD",
                             Patterns::Selection ("SPD|PD|symmetric|none"),
                             "This parameters allows for the stabilization of the velocity block. If one derives the Newton "
                             "method without any modifications, the matrix created for the velocity block is not necessarily "
                             "Symmetric Positive Definite. This is problematic (see \\cite{FBTGS19}). When `none' is chosen, "
                             "the velocity block is not stabilized. The `symmetric' parameters symmetrizes the matrix, and `PD' makes "
                             "the matrix Positive Definite. `SPD' is the full stabilization, where the matrix is guaranteed Symmetric "
                             "Positive Definite.");

          prm.declare_entry ("Use Newton failsafe", "false",
                             Patterns::Bool (),
                             "When this parameter is true and the linear solver fails, we try again, but now with SPD stabilization "
                             "for both the preconditioner and the velocity block. The SPD stabilization will remain active until "
                             "the next timestep, when the default values are restored.");


          prm.declare_entry ("SPD safety factor", "0.9",
                             Patterns::Double (0., 1.),
                             "When stabilizing the Newton matrix, we can encounter situations where the coefficient inside the elliptic (top-left) "
                             "block becomes negative or zero. This coefficient has the form $1+x$ where $x$ can sometimes be smaller than $-1$. In "
                             "this case, the top-left block of the matrix is no longer positive definite, and both preconditioners and iterative "
                             "solvers may fail. To prevent this, the stabilization computes an $\\alpha$ so that $1+\\alpha x$ is never negative. "
                             "This $\\alpha$ is chosen as $1$ if $x\\ge -1$, and $\\alpha=-\\frac 1x$ otherwise. (Note that this always leads to "
                             "$0\\le \\alpha \\le 1$.)  On the other hand, we also want to stay away from $1+\\alpha x=0$, and so modify the choice of "
                             "$\\alpha$ to be $1$ if $x\\ge -c$, and $\\alpha=-\\frac cx$ with a $c$ between zero and one. This way, if $c<1$, we are "
                             "assured that $1-\\alpha x>c$, i.e., bounded away from zero.");

          prm.declare_entry ("Use Eisenstat Walker method for Picard iterations", "false",
                             Patterns::Bool(),
                             "If set to true, the Picard iteration uses the Eisenstat Walker method to determine how accurately linear systems need "
                             "to be solved. The Picard iteration is used, for example, in the first few iterations of the Newton method before the "
                             "matrix is built including derivatives of the model, since the Picard iteration generally converges even from points "
                             "where Newton's method does not. "
                             "\n\n"
                             "Once derivatives are used in a Newton method, \\aspect{} always uses the Eisenstat Walker method.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    void
    Parameters::
    parse_parameters (ParameterHandler &prm)
    {
      CitationInfo::add("NewtonSolver");
      prm.enter_subsection ("Solver parameters");
      {
        prm.enter_subsection ("Newton solver parameters");
        {
          nonlinear_switch_tolerance = prm.get_double("Nonlinear Newton solver switch tolerance");
          max_pre_newton_nonlinear_iterations = prm.get_integer ("Max pre-Newton nonlinear iterations");
          max_newton_line_search_iterations = prm.get_integer ("Max Newton line search iterations");
          use_newton_residual_scaling_method = prm.get_bool("Use Newton residual scaling method");
          maximum_linear_stokes_solver_tolerance = prm.get_double("Maximum linear Stokes solver tolerance");
          std::string preconditioner_stabilization_string = prm.get("Stabilization preconditioner");
          if (preconditioner_stabilization_string == "SPD")
            preconditioner_stabilization = Stabilization::SPD;
          else if (preconditioner_stabilization_string == "PD")
            preconditioner_stabilization = Stabilization::PD;
          else if (preconditioner_stabilization_string == "symmetric")
            preconditioner_stabilization = Stabilization::symmetric;
          else if (preconditioner_stabilization_string == "none")
            preconditioner_stabilization = Stabilization::none;

          std::string velocity_block_stabilization_string = prm.get("Stabilization velocity block");
          if (velocity_block_stabilization_string == "SPD")
            velocity_block_stabilization = Stabilization::SPD;
          else if (velocity_block_stabilization_string == "PD")
            velocity_block_stabilization = Stabilization::PD;
          else if (velocity_block_stabilization_string == "symmetric")
            velocity_block_stabilization = Stabilization::symmetric;
          else if (velocity_block_stabilization_string == "none")
            velocity_block_stabilization = Stabilization::none;

          use_Newton_failsafe = prm.get_bool("Use Newton failsafe");
          SPD_safety_factor = prm.get_double("SPD safety factor");
          use_Eisenstat_Walker_method_for_Picard_iterations = prm.get_bool("Use Eisenstat Walker method for Picard iterations");

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

    }
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

#undef INSTANTIATE

}
