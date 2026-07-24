/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#include "common.h"

#include <aspect/material_model/reaction_model/solid_state/cahn1956/site_saturated_n1/kinetics/eutectoid_decomposition.h>
#include <aspect/material_model/reaction_model/solid_state/cahn1956/site_saturated_n1/kinetics/interface_controlled_growth.h>
#include <aspect/material_model/reaction_model/solid_state/reaction_chain.h>

#include <deal.II/base/parameter_handler.h>
#include <cmath>

TEST_CASE("Eutectoid Decomposition Kinetics")
{
  using namespace aspect::MaterialModel::ReactionModel::SolidState::Cahn1956::SiteSaturatedN1;

  EutectoidDecomposition<2> kinetics;
  dealii::ParameterHandler prm;

  kinetics.declare_parameters(prm);

  prm.enter_subsection("Eutectoid decomposition");
  {
    prm.set("Kinetic factor", "2.0e-16");
    prm.set("Activation energy", "300000.0");
  }
  prm.leave_subsection();

  kinetics.parse_parameters(prm);

  const double temperature = 1500.0; // K
  const double pressure = 1.0e9;     // Pa
  const double R = 8.314;            // J/mol/K

  SECTION("Arrhenius factor calculation")
  {
    const double expected_arrhenius = std::exp(-300000.0 / (R * temperature));
    const double computed_arrhenius = kinetics.arrhenius_factor(temperature, pressure);

    CHECK(computed_arrhenius == Approx(expected_arrhenius));
  }

  SECTION("Thermodynamic factor calculation")
  {
    // Case 1: Forward reaction (dG < 0, B is favored)
    const double dG_forward = -5000.0;
    const double expected_tf_forward = -dG_forward * std::abs(dG_forward); // 2.5e7
    CHECK(kinetics.thermodynamic_factor(temperature, dG_forward) == Approx(expected_tf_forward));

    // Case 2: Reverse reaction (dG > 0, A is favored)
    const double dG_reverse = 5000.0;
    const double expected_tf_reverse = -dG_reverse * std::abs(dG_reverse); // -2.5e7
    CHECK(kinetics.thermodynamic_factor(temperature, dG_reverse) == Approx(expected_tf_reverse));

    // Case 3: Equilibrium (dG = 0)
    CHECK(kinetics.thermodynamic_factor(temperature, 0.0) == Approx(0.0));
  }

  SECTION("Reaction rate calculation")
  {
    const double dG = -2000.0;
    const double reactant_fraction = 0.8;

    const double arrhenius = std::exp(-300000.0 / (R * temperature));
    const double therm_factor = -dG * std::abs(dG);
    const double expected_rate = 2.0e-16 * therm_factor * arrhenius * reactant_fraction;

    const double computed_rate = kinetics.reaction_rate(temperature, pressure, dG, reactant_fraction);
    CHECK(computed_rate == Approx(expected_rate));
  }
}

TEST_CASE("Interface Controlled Growth Kinetics")
{
  using namespace aspect::MaterialModel::ReactionModel::SolidState::Cahn1956::SiteSaturatedN1;

  InterfaceControlledGrowth<2> kinetics;
  dealii::ParameterHandler prm;

  kinetics.declare_parameters(prm);

  prm.enter_subsection("Interface controlled growth");
  {
    prm.set("Kinetic factor", "1.0e8");
    prm.set("Activation enthalpy", "200000.0");
    prm.set("Activation volume", "2.0e-6");
  }
  prm.leave_subsection();

  kinetics.parse_parameters(prm);

  const double temperature = 1600.0; // K
  const double pressure = 2.0e9;     // Pa
  const double R = 8.314;            // J/mol/K

  SECTION("Arrhenius factor calculation")
  {
    const double exponent = -(200000.0 + pressure * 2.0e-6) / (R * temperature);
    const double expected_arrhenius = std::exp(exponent);
    const double computed_arrhenius = kinetics.arrhenius_factor(temperature, pressure);

    CHECK(computed_arrhenius == Approx(expected_arrhenius));
  }

  SECTION("Thermodynamic factor calculation")
  {
    // Case 1: Forward reaction (dG < 0)
    const double dG_forward = -1000.0;
    const double expected_tf_forward = 1.0 - std::exp(-std::abs(dG_forward) / (R * temperature));
    CHECK(kinetics.thermodynamic_factor(temperature, dG_forward) == Approx(expected_tf_forward));

    // Case 2: Reverse reaction (dG > 0)
    const double dG_reverse = 1000.0;
    const double expected_tf_reverse = -(1.0 - std::exp(-std::abs(dG_reverse) / (R * temperature)));
    CHECK(kinetics.thermodynamic_factor(temperature, dG_reverse) == Approx(expected_tf_reverse));

    // Case 3: Thermodynamic equilibrium (dG = 0)
    CHECK(kinetics.thermodynamic_factor(temperature, 0.0) == Approx(0.0));
  }

  SECTION("Reaction rate calculation")
  {
    const double dG = -1500.0;
    const double reactant_fraction = 0.5;

    const double arrhenius = std::exp(-(200000.0 + pressure * 2.0e-6) / (R * temperature));
    const double therm_factor = 1.0 - std::exp(-std::abs(dG) / (R * temperature));
    const double expected_rate = 1.0e8 * temperature * arrhenius * therm_factor * reactant_fraction;

    const double computed_rate = kinetics.reaction_rate(temperature, pressure, dG, reactant_fraction);
    CHECK(computed_rate == Approx(expected_rate));
  }
}

TEST_CASE("Reaction Chain Processing")
{
  using namespace aspect::MaterialModel::ReactionModel::SolidState;

  ReactionChain<2> chain;
  chain.phase_names = {"olivine", "wadsleyite", "ringwoodite", "postspinel"};
  chain.reactions.resize(chain.phase_names.size() - 1);

  SECTION("Checking active reactions status")
  {
    chain.reactions[0].enabled = false;
    chain.reactions[1].enabled = false;
    chain.reactions[2].enabled = false;
    CHECK_FALSE(chain.any_enabled());

    chain.reactions[1].enabled = true;
    CHECK(chain.any_enabled());
  }

  SECTION("Clamping cumulative reaction progress")
  {
    // Input vector with out-of-bounds and un-ordered values
    const std::vector<double> xi_raw = {1.2, 0.9, -0.2};
    const std::vector<double> xi_clamped = chain.clamp_cumulative_progress(xi_raw);

    // Expected order: 1 >= xi_0 >= xi_1 >= xi_2 >= 0
    const std::vector<double> expected = {1.0, 0.9, 0.0};
    compare_vectors_approx(xi_clamped, expected);
  }

  SECTION("Phase mass fraction extraction from cumulative progress")
  {
    SECTION("Complete conversion along chain")
    {
      // xi = {1.0, 1.0, 0.4}
      // Phase 0 (olivine): 1 - 1 = 0.0
      // Phase 1 (wadsleyite): 1 - 1 = 0.0
      // Phase 2 (ringwoodite): 1 - 0.4 = 0.6
      // Phase 3 (postspinel): 0.4
      const std::vector<double> xi = {1.0, 1.0, 0.4};
      const std::vector<double> expected_X = {0.0, 0.0, 0.6, 0.4};

      const std::vector<double> computed_X = chain.compute_phase_mass_fractions(xi);
      compare_vectors_approx(computed_X, expected_X);
    }

    SECTION("Partial reaction progress across all phases")
    {
      const std::vector<double> xi = {0.8, 0.5, 0.1};
      const std::vector<double> expected_X = {0.2, 0.3, 0.4, 0.1};

      const std::vector<double> computed_X = chain.compute_phase_mass_fractions(xi);
      compare_vectors_approx(computed_X, expected_X);
    }

    SECTION("Single phase boundary limit")
    {
      ReactionChain<2> single_reaction_chain;
      single_reaction_chain.phase_names = {"phase_A", "phase_B"};
      single_reaction_chain.reactions.resize(1);

      const std::vector<double> xi = {0.3};
      const std::vector<double> expected_X = {0.7, 0.3};

      const std::vector<double> computed_X = single_reaction_chain.compute_phase_mass_fractions(xi);
      compare_vectors_approx(computed_X, expected_X);
    }
  }
}
