/*
  Copyright (C) 2022 - 2023 by the authors of the ASPECT code.

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

#include <aspect/simulator_access.h>

#include <aspect/material_model/simple.h>
#include <aspect/material_model/additional_outputs/anisotropic_viscosity.h>
#include <aspect/heating_model/shear_heating.h>
#include <aspect/heating_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/simulator/assemblers/stokes_anisotropic_viscosity.h>
#include <aspect/simulator_signals.h>

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class Anisotropic : public MaterialModel::Simple<dim>
    {
      public:
        void initialize() override;

        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;
        static void declare_parameters(ParameterHandler &prm);
        void parse_parameters(ParameterHandler &prm) override;

        /**
          * Return true if the compressibility() function returns something that
          * is not zero.
          */
        bool
        is_compressible () const override;

      private:
        void set_assemblers(const SimulatorAccess<dim> &,
                            Assemblers::Manager<dim> &assemblers) const;

        // Constitutive tensor
        SymmetricTensor<4,dim> C;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    Anisotropic<dim>::set_assemblers(const SimulatorAccess<dim> &,
                                     Assemblers::Manager<dim> &assemblers) const
    {
      for (unsigned int i=0; i<assemblers.stokes_preconditioner.size(); ++i)
        {
          if (dynamic_cast<Assemblers::StokesPreconditioner<dim> *>(assemblers.stokes_preconditioner[i].get()) != nullptr)
            assemblers.stokes_preconditioner[i] = std::make_unique<Assemblers::StokesPreconditionerAnisotropicViscosity<dim>> ();
          if (dynamic_cast<Assemblers::StokesCompressiblePreconditioner<dim> *>(assemblers.stokes_preconditioner[i].get()) != nullptr)
            assemblers.stokes_preconditioner[i] = std::make_unique<Assemblers::StokesCompressiblePreconditionerAnisotropicViscosity<dim>> ();
        }

      for (unsigned int i=0; i<assemblers.stokes_system.size(); ++i)
        {
          if (dynamic_cast<Assemblers::StokesIncompressibleTerms<dim> *>(assemblers.stokes_system[i].get()) != nullptr)
            assemblers.stokes_system[i] = std::make_unique<Assemblers::StokesIncompressibleTermsAnisotropicViscosity<dim>> ();
          if (dynamic_cast<Assemblers::StokesCompressibleStrainRateViscosityTerm<dim> *>(assemblers.stokes_system[i].get()) != nullptr)
            assemblers.stokes_system[i] = std::make_unique<Assemblers::StokesCompressibleStrainRateViscosityTermAnisotropicViscosity<dim>> ();
        }
    }

    template <int dim>
    void
    Anisotropic<dim>::
    initialize()
    {
      this->get_signals().set_assemblers.connect (std::bind(&Anisotropic<dim>::set_assemblers,
                                                            std::cref(*this),
                                                            std::placeholders::_1,
                                                            std::placeholders::_2));
    }

    template <int dim>
    void
    Anisotropic<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      const std::shared_ptr<MaterialModel::AnisotropicViscosity<dim>> anisotropic_viscosity =
        out.template get_additional_output_object<MaterialModel::AnisotropicViscosity<dim>>();

      Simple<dim>::evaluate(in, out);
      Point<dim> center;
      center[0] = 0.5;
      center[1] = 0.5;
      if (dim == 3)
        center[2] = 0.5;
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const double pressure = in.pressure[i];
          out.densities[i] = 1.0 + pressure;
          out.compressibilities[i] = 1.0 / (1. + pressure);
          if ((in.position[i]-center).norm() < 0.25)
            {
              if (anisotropic_viscosity != nullptr)
                anisotropic_viscosity->stress_strain_directors[i] = C;
            }
        }
    }

    template <int dim>
    bool
    Anisotropic<dim>::
    is_compressible () const
    {
      return true;
    }

    template <int dim>
    void
    Anisotropic<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      Simple<dim>::declare_parameters (prm);
      prm.enter_subsection ("Material model");
      {
        prm.enter_subsection ("Anisotropic");
        {
          if (dim == 2)
            prm.declare_entry ("Viscosity tensor",
                               "1, 0, 0,"
                               "0, 1, 0,"
                               "0, 0,.5",
                               Patterns::List(Patterns::Double()),
                               "Viscosity-scaling tensor in Voigt notation.");
          else
            prm.declare_entry ("Viscosity tensor",
                               "1, 0, 0, 0, 0, 0,"
                               "0, 1, 0, 0, 0, 0,"
                               "0, 0, 1, 0, 0, 0,"
                               "0, 0, 0,.5, 0, 0,"
                               "0, 0, 0, 0,.5, 0,"
                               "0, 0, 0, 0, 0,.5",
                               Patterns::List(Patterns::Double()),
                               "Viscosity-scaling tensor in Voigt notation.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Anisotropic<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      Simple<dim>::parse_parameters (prm);
      prm.enter_subsection ("Material model");
      {
        prm.enter_subsection ("Anisotropic");
        {
          const int size_voigt = (dim == 3 ? 6 : 3);
          const std::vector<double> tmp_tensor =
            Utilities::string_to_double(Utilities::split_string_list(prm.get ("Viscosity tensor")));
          Assert(tmp_tensor.size() == size_voigt*size_voigt,
                 ExcMessage("Constitutive voigt matrix must have 9 components in 2D, or 36 components in 3d"));

          std::vector<std::vector<double>> voigt_visc_tensor (size_voigt);
          for (unsigned int i=0; i<size_voigt; ++i)
            {
              voigt_visc_tensor[i].resize(size_voigt);
              for (unsigned int j=0; j<size_voigt; ++j)
                voigt_visc_tensor[i][j] = tmp_tensor[i*size_voigt+j];
            }

          // Voigt indices (For mapping back to real tensor)
          const unsigned int vi3d0[] = {0, 1, 2, 1, 0, 0};
          const unsigned int vi3d1[] = {0, 1, 2, 2, 2, 1};
          const unsigned int vi2d0[] = {0, 1, 0};
          const unsigned int vi2d1[] = {0, 1, 1};

          // Fill the constitutive tensor with values from the Voigt tensor
          for (unsigned int i=0; i<size_voigt; ++i)
            for (unsigned int j=0; j<size_voigt; ++j)
              if (dim == 2)
                C[vi2d0[i]][vi2d1[i]][vi2d0[j]][vi2d1[j]] = voigt_visc_tensor[i][j];
              else
                C[vi3d0[i]][vi3d1[i]][vi3d0[j]][vi3d1[j]] = voigt_visc_tensor[i][j];
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Anisotropic,
                                   "anisotropic",
                                   "A simple material model that is like the "
                                   "'Simple' model, but has a non-zero compressibility "
                                   "and always has a blob of anisotropic material in "
                                   "the center.")
  }
}
