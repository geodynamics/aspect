/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.
*/

#include <aspect/material_model/interface.h>
#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

#include <cstdlib>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class AdjointCellwisePerturbation : public Simple<dim>
    {
      public:
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                 MaterialModel::MaterialModelOutputs<dim> &out) const override;

        void
        parse_parameters(ParameterHandler &prm) override;

      private:
        int perturbed_active_cell_index;
        std::string perturbed_property;
        double perturbation_amplitude;
    };



    template <int dim>
    void
    AdjointCellwisePerturbation<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                                               MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      Simple<dim>::evaluate(in, out);

      if (perturbed_active_cell_index < 0
          || perturbation_amplitude == 0.0
          || in.current_cell.state() != IteratorState::valid
          || static_cast<int>(in.current_cell->active_cell_index()) != perturbed_active_cell_index)
        return;

      if (perturbed_property == "density")
        for (double &density : out.densities)
          density += perturbation_amplitude;
      else if (perturbed_property == "viscosity")
        for (double &viscosity : out.viscosities)
          viscosity += perturbation_amplitude;
    }


    template <int dim>
    void
    AdjointCellwisePerturbation<dim>::parse_parameters(ParameterHandler &prm)
    {
      Simple<dim>::parse_parameters(prm);

      const char *cell = std::getenv("ASPECT_ADJOINT_FD_CELL");
      const char *property = std::getenv("ASPECT_ADJOINT_FD_PROPERTY");
      const char *amplitude = std::getenv("ASPECT_ADJOINT_FD_AMPLITUDE");

      perturbed_active_cell_index = (cell != nullptr ? std::atoi(cell) : -1);
      perturbed_property = (property != nullptr ? property : "none");
      perturbation_amplitude = (amplitude != nullptr ? std::atof(amplitude) : 0.0);
    }
  }
}

namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(AdjointCellwisePerturbation,
                                   "adjoint cellwise perturbation",
                                   "A benchmark material model that extends the simple model "
                                   "with an optional additive density or viscosity perturbation "
                                   "on one active cell. It is used only by the dynamic-topography "
                                   "adjoint finite-difference validation script.")
  }
}
