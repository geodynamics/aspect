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

#include <aspect/material_model/reaction_model/kinetics/cahn1956_eutectoid_decomposition.h>
#include <aspect/material_model/reaction_model/kinetics/cahn1956_interface_controlled_growth.h>
#include <aspect/material_model/reaction_model/reaction_chain.h>

#include <deal.II/base/parameter_handler.h>
#include <cmath>

TEST_CASE("Eutectoid Decomposition Kinetics")
{
  using namespace aspect::MaterialModel::ReactionModel;

  EutectoidDecomposition<2> kinetics;
  dealii::ParameterHandler prm;

  kinetics.declare_parameters(prm);

  prm.enter_subsection("Eutectoid decomposition");
  {
    prm.set("Kinetic prefactors", "2.0e-16");
    prm.set("Activation energies", "300000.0");
  }
  prm.leave_subsection();

  kinetics.parse_parameters(prm, 1);

  const double temperature = 1500.0; // K
  const double pressure = 1.0e9;     // Pa
  const double R = aspect::constants::gas_constant;

  SECTION("Net forward reaction rate calculation")
  {
    const double arrhenius = std::exp(-300000.0 / (R * temperature));

    SECTION("Forward reaction (dG < 0)")
    {
      const double dG_forward = -2000.0;
      const double forward_reaction_progress = 0.2;
      const double reaction_progress = 1.0 - forward_reaction_progress;
      const double therm_factor = -dG_forward * std::abs(dG_forward); // 4.0e6
      const double expected_rate = 2.0e-16 * therm_factor * arrhenius * reaction_progress;

      const double computed_rate = kinetics.net_forward_reaction_rate(temperature, pressure, dG_forward, forward_reaction_progress, 0);
      CHECK(computed_rate == Approx(expected_rate));
    }

    SECTION("Forward reaction (dG < 0) from pure reactant")
    {
      const double dG_forward = -2000.0;
      const double forward_reaction_progress = 0.0;
      const double reaction_progress = 1.0 - forward_reaction_progress;
      const double therm_factor = -dG_forward * std::abs(dG_forward); // 4.0e6
      const double expected_rate = 2.0e-16 * therm_factor * arrhenius * reaction_progress;

      const double computed_rate = kinetics.net_forward_reaction_rate(temperature, pressure, dG_forward, forward_reaction_progress, 0);
      CHECK(computed_rate == Approx(expected_rate));
    }

    SECTION("Reverse reaction (dG > 0)")
    {
      const double dG_reverse = 2000.0;
      const double forward_reaction_progress = 0.2;
      const double reaction_progress = forward_reaction_progress;
      const double therm_factor = -dG_reverse * std::abs(dG_reverse); // -4.0e6
      const double expected_rate = 2.0e-16 * therm_factor * arrhenius * reaction_progress;

      const double computed_rate = kinetics.net_forward_reaction_rate(temperature, pressure, dG_reverse, forward_reaction_progress, 0);
      CHECK(computed_rate == Approx(expected_rate));
    }

    SECTION("Equilibrium (dG = 0)")
    {
      const double computed_rate = kinetics.net_forward_reaction_rate(temperature, pressure, 0.0, 0.5, 0);
      CHECK(computed_rate == Approx(0.0));
    }
  }
}

TEST_CASE("Interface Controlled Growth Kinetics")
{
  using namespace aspect::MaterialModel::ReactionModel;

  InterfaceControlledGrowth<2> kinetics;
  dealii::ParameterHandler prm;

  kinetics.declare_parameters(prm);

  prm.enter_subsection("Interface controlled growth");
  {
    prm.set("Kinetic prefactors", "1.0e8");
    prm.set("Activation enthalpies", "200000.0");
    prm.set("Activation volumes", "2.0e-6");
  }
  prm.leave_subsection();

  kinetics.parse_parameters(prm, 1);

  const double temperature = 1600.0; // K
  const double pressure = 2.0e9;     // Pa
  const double R = aspect::constants::gas_constant;

  SECTION("Net forward reaction rate calculation")
  {
    const double arrhenius = std::exp(-(200000.0 + pressure * 2.0e-6) / (R * temperature));

    SECTION("Forward reaction (dG < 0)")
    {
      const double dG_forward = -1500.0;
      const double forward_reaction_progress = 0.5;
      const double reaction_progress = 1.0 - forward_reaction_progress;
      const double therm_factor = 1.0 - std::exp(-std::abs(dG_forward) / (R * temperature));
      const double expected_rate = 1.0e8 * temperature * arrhenius * therm_factor * reaction_progress;

      const double computed_rate = kinetics.net_forward_reaction_rate(temperature, pressure, dG_forward, forward_reaction_progress, 0);
      CHECK(computed_rate == Approx(expected_rate));
    }

    SECTION("Forward reaction (dG < 0) from pure reactant")
    {
      const double dG_forward = -1500.0;
      const double forward_reaction_progress = 0.0;
      const double reaction_progress = 1.0 - forward_reaction_progress;
      const double therm_factor = 1.0 - std::exp(-std::abs(dG_forward) / (R * temperature));
      const double expected_rate = 1.0e8 * temperature * arrhenius * therm_factor * reaction_progress;

      const double computed_rate = kinetics.net_forward_reaction_rate(temperature, pressure, dG_forward, forward_reaction_progress, 0);
      CHECK(computed_rate == Approx(expected_rate));
    }

    SECTION("Reverse reaction (dG > 0)")
    {
      const double dG_reverse = 1500.0;
      const double forward_reaction_progress = 0.5;
      const double reaction_progress = forward_reaction_progress;
      const double therm_factor = -(1.0 - std::exp(-std::abs(dG_reverse) / (R * temperature)));
      const double expected_rate = 1.0e8 * temperature * arrhenius * therm_factor * reaction_progress;

      const double computed_rate = kinetics.net_forward_reaction_rate(temperature, pressure, dG_reverse, forward_reaction_progress, 0);
      CHECK(computed_rate == Approx(expected_rate));
    }

    SECTION("Equilibrium (dG = 0)")
    {
      const double computed_rate = kinetics.net_forward_reaction_rate(temperature, pressure, 0.0, 0.5, 0);
      CHECK(computed_rate == Approx(0.0));
    }
  }
}

TEST_CASE("Multiple Reactions Sharing One Kinetics Instance")
{
  using namespace aspect::MaterialModel::ReactionModel;

  InterfaceControlledGrowth<2> kinetics;
  dealii::ParameterHandler prm;

  kinetics.declare_parameters(prm);

  prm.enter_subsection("Interface controlled growth");
  {
    prm.set("Kinetic prefactors", "1.0e8|2.0e8");
    prm.set("Activation enthalpies", "200000.0|250000.0");
    prm.set("Activation volumes", "2.0e-6|3.0e-6");
  }
  prm.leave_subsection();

  kinetics.parse_parameters(prm, 2);

  const double temperature = 1600.0; // K
  const double pressure = 2.0e9;     // Pa
  const double dG = -1500.0;
  const double forward_reaction_progress = 0.5;
  const double R = aspect::constants::gas_constant;

  SECTION("Local index 0 uses the first entry of each list")
  {
    const double arrhenius = std::exp(-(200000.0 + pressure * 2.0e-6) / (R * temperature));
    const double therm_factor = 1.0 - std::exp(-std::abs(dG) / (R * temperature));
    const double expected_rate = 1.0e8 * temperature * arrhenius * therm_factor * forward_reaction_progress;

    const double computed_rate = kinetics.net_forward_reaction_rate(temperature, pressure, dG, forward_reaction_progress, 0);
    CHECK(computed_rate == Approx(expected_rate));
  }

  SECTION("Local index 1 uses the second entry of each list")
  {
    const double arrhenius = std::exp(-(250000.0 + pressure * 3.0e-6) / (R * temperature));
    const double therm_factor = 1.0 - std::exp(-std::abs(dG) / (R * temperature));
    const double expected_rate = 2.0e8 * temperature * arrhenius * therm_factor * forward_reaction_progress;

    const double computed_rate = kinetics.net_forward_reaction_rate(temperature, pressure, dG, forward_reaction_progress, 1);
    CHECK(computed_rate == Approx(expected_rate));
  }
}

TEST_CASE("Reaction Chain Processing")
{
  using namespace aspect::MaterialModel::ReactionModel;

  dealii::ParameterHandler prm;
  ReactionChain<2> reaction_chain;
  ReactionChain<2>::declare_parameters(prm);

  prm.enter_subsection("Reaction chain");
  {
    prm.set("Kinetic models", "Interface controlled growth|Eutectoid decomposition|Interface controlled growth");

    prm.enter_subsection("Interface controlled growth");
    {
      prm.set("Kinetic prefactors", "7.0e7|7.0e7");
      prm.set("Activation enthalpies", "274e3|274e3");
      prm.set("Activation volumes", "3.3e-6|3.3e-6");
    }
    prm.leave_subsection();

    prm.enter_subsection("Eutectoid decomposition");
    {
      prm.set("Kinetic prefactors", "2.7e-16");
      prm.set("Activation energies", "355e3");
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();

  reaction_chain.parse_parameters(prm);

  SECTION("Clamping cumulative reaction progress to physical range")
  {
    compare_vectors_approx(reaction_chain.clamp_cumulative_progress({1.04, 0.85, -0.04}), {1.0, 0.85, 0.0});
  }

  SECTION("Without strict tolerance, out-of-range and non-monotonic values are clamped")
  {
    // Values below 0 (reverse reaction completed) or above 1 (forward reaction completed) are clamped to [0, 1]
    compare_vectors_approx(reaction_chain.clamp_cumulative_progress({1.06, 0.5, 0.1}), {1.0, 0.5, 0.1});
    compare_vectors_approx(reaction_chain.clamp_cumulative_progress({0.8, 0.5, -0.06}), {0.8, 0.5, 0.0});

    // Non-monotonic values are clamped to enforce xi[i] <= xi[i-1]
    compare_vectors_approx(reaction_chain.clamp_cumulative_progress({0.80, 0.87, 0.1}), {0.80, 0.80, 0.1});
  }

  SECTION("With strict tolerance, out-of-range and non-monotonic values throw")
  {
    dealii::ParameterHandler strict_prm;
    ReactionChain<2> strict_reaction_chain;
    ReactionChain<2>::declare_parameters(strict_prm);

    strict_prm.enter_subsection("Reaction chain");
    {
      strict_prm.set("Kinetic models", "Interface controlled growth|Eutectoid decomposition|Interface controlled growth");
      strict_prm.set("Tolerance in reaction progress", "0.05");

      strict_prm.enter_subsection("Interface controlled growth");
      {
        strict_prm.set("Kinetic prefactors", "7.0e7|7.0e7");
        strict_prm.set("Activation enthalpies", "274e3|274e3");
        strict_prm.set("Activation volumes", "3.3e-6|3.3e-6");
      }
      strict_prm.leave_subsection();

      strict_prm.enter_subsection("Eutectoid decomposition");
      {
        strict_prm.set("Kinetic prefactors", "2.7e-16");
        strict_prm.set("Activation energies", "355e3");
      }
      strict_prm.leave_subsection();
    }
    strict_prm.leave_subsection();

    strict_reaction_chain.parse_parameters(strict_prm);

    CHECK_THROWS_AS(strict_reaction_chain.clamp_cumulative_progress({1.06, 0.5, 0.1}), dealii::ExceptionBase);
    CHECK_THROWS_AS(strict_reaction_chain.clamp_cumulative_progress({0.8, 0.5, -0.06}), dealii::ExceptionBase);
    CHECK_THROWS_AS(strict_reaction_chain.clamp_cumulative_progress({0.80, 0.87, 0.1}), dealii::ExceptionBase);
  }

  SECTION("Clamping throws on NaN or Inf input")
  {
    const double nan_val = std::numeric_limits<double>::quiet_NaN();
    const double inf_val = std::numeric_limits<double>::infinity();

    CHECK_THROWS_AS(reaction_chain.clamp_cumulative_progress({0.8, nan_val, 0.2}), dealii::ExceptionBase);
    CHECK_THROWS_AS(reaction_chain.clamp_cumulative_progress({inf_val, 0.5, 0.2}), dealii::ExceptionBase);
  }

  SECTION("Phase mass fraction extraction from cumulative progress")
  {
    SECTION("Uniform phase densities (equal volume and mass fractions)")
    {
      // xi = {1.0, 1.0, 0.4}
      // Volume fractions:
      // Phase 0: 1 - 1 = 0.0
      // Phase 1: 1 - 1 = 0.0
      // Phase 2: 1 - 0.4 = 0.6
      // Phase 3: 0.4
      const std::vector<double> xi = {1.0, 1.0, 0.4};
      const std::vector<double> densities = {3300.0, 3300.0, 3300.0, 3300.0};
      const std::vector<double> expected_X = {0.0, 0.0, 0.6, 0.4};

      const std::vector<double> computed_X = reaction_chain.compute_phase_mass_fractions(xi, densities);
      compare_vectors_approx(computed_X, expected_X);
    }

    SECTION("Variable phase densities")
    {
      // xi = {0.8, 0.5, 0.1}
      // Volume fractions: phi = {0.2, 0.3, 0.4, 0.1}
      // Densities: rho = {3200, 3500, 3700, 4100} kg/m^3
      // Mass = {640, 1050, 1480, 410} -> Bulk density = 3580 kg/m^3
      const std::vector<double> xi = {0.8, 0.5, 0.1};
      const std::vector<double> densities = {3200.0, 3500.0, 3700.0, 4100.0};

      const double bulk_density = 0.2 * 3200.0 + 0.3 * 3500.0 + 0.4 * 3700.0 + 0.1 * 4100.0;
      const std::vector<double> expected_X =
      {
        (0.2 * 3200.0) / bulk_density,
        (0.3 * 3500.0) / bulk_density,
        (0.4 * 3700.0) / bulk_density,
        (0.1 * 4100.0) / bulk_density
      };

      const std::vector<double> computed_X = reaction_chain.compute_phase_mass_fractions(xi, densities);
      compare_vectors_approx(computed_X, expected_X);
    }

    SECTION("Single phase boundary limit")
    {
      dealii::ParameterHandler single_prm;
      ReactionChain<2> single_reaction_chain;
      ReactionChain<2>::declare_parameters(single_prm);

      single_reaction_chain.parse_parameters(single_prm); // 1 default reaction

      const std::vector<double> xi = {0.3};
      const std::vector<double> densities = {3000.0, 3000.0};
      const std::vector<double> expected_X = {0.7, 0.3};

      const std::vector<double> computed_X = single_reaction_chain.compute_phase_mass_fractions(xi, densities);
      compare_vectors_approx(computed_X, expected_X);
    }

    SECTION("Asserts throw on invalid inputs")
    {
      const std::vector<double> valid_xi = {0.8, 0.5, 0.1};
      const std::vector<double> valid_densities = {3200.0, 3500.0, 3700.0, 4100.0};

      // Empty reaction progress values
      CHECK_THROWS_AS(reaction_chain.compute_phase_mass_fractions({}, valid_densities), dealii::ExceptionBase);

      // Non-positive phase density
      const std::vector<double> zero_density = {3200.0, 0.0, 3700.0, 4100.0};
      CHECK_THROWS_AS(reaction_chain.compute_phase_mass_fractions(valid_xi, zero_density), dealii::ExceptionBase);

      // Non-finite (Inf/NaN) phase density
      const std::vector<double> inf_density = {3200.0, std::numeric_limits<double>::infinity(), 3700.0, 4100.0};
      CHECK_THROWS_AS(reaction_chain.compute_phase_mass_fractions(valid_xi, inf_density), dealii::ExceptionBase);
    }
  }
}

TEST_CASE("Reaction Chain Kinetic Models Parsing")
{
  using namespace aspect::MaterialModel::ReactionModel;

  SECTION("Number of reactions matches the number of entries in 'Kinetic models'")
  {
    dealii::ParameterHandler prm;
    ReactionChain<2> reaction_chain;

    ReactionChain<2>::declare_parameters(prm);

    prm.enter_subsection("Reaction chain");
    {
      prm.set("Kinetic models", "Interface controlled growth|Eutectoid decomposition");
    }
    prm.leave_subsection();

    reaction_chain.parse_parameters(prm);

    REQUIRE(reaction_chain.n_reactions() == 2);
    REQUIRE(reaction_chain.reactions.size() == 2);
    CHECK(reaction_chain.reactions[0].kinetics != nullptr);
    CHECK(reaction_chain.reactions[1].kinetics != nullptr);
    // Distinct model names -> distinct kinetics instances, each the sole (local index 0) user of its instance.
    CHECK(reaction_chain.reactions[0].kinetics != reaction_chain.reactions[1].kinetics);
    CHECK(reaction_chain.reactions[0].local_reaction_index == 0);
    CHECK(reaction_chain.reactions[1].local_reaction_index == 0);
  }

  SECTION("Repeated model name shares one kinetics instance across local indices")
  {
    dealii::ParameterHandler prm;
    ReactionChain<2> reaction_chain;

    ReactionChain<2>::declare_parameters(prm);

    prm.enter_subsection("Reaction chain");
    {
      prm.set("Kinetic models", "Interface controlled growth|Interface controlled growth|Eutectoid decomposition");
      prm.enter_subsection("Interface controlled growth");
      {
        prm.set("Kinetic prefactors", "1.4e4|1.4e4");
        prm.set("Activation enthalpies", "274e3|274e3");
        prm.set("Activation volumes", "3.3e-6|3.3e-6");
      }
      prm.leave_subsection();
      prm.enter_subsection("Eutectoid decomposition");
      {
        prm.set("Kinetic prefactors", "1.6e-3");
        prm.set("Activation energies", "355e3");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();

    reaction_chain.parse_parameters(prm);

    REQUIRE(reaction_chain.n_reactions() == 3);
    // The two "Interface controlled growth" reactions share one instance, with consecutive local indices.
    CHECK(reaction_chain.reactions[0].kinetics == reaction_chain.reactions[1].kinetics);
    CHECK(reaction_chain.reactions[0].local_reaction_index == 0);
    CHECK(reaction_chain.reactions[1].local_reaction_index == 1);
    // The "Eutectoid decomposition" reaction gets its own instance.
    CHECK(reaction_chain.reactions[2].kinetics != reaction_chain.reactions[0].kinetics);
    CHECK(reaction_chain.reactions[2].local_reaction_index == 0);
  }

  SECTION("Default 'Kinetic models' declared in parameters")
  {
    dealii::ParameterHandler prm;
    ReactionChain<2>::declare_parameters(prm);

    prm.enter_subsection("Reaction chain");
    const std::string parsed_default = prm.get("Kinetic models");
    prm.leave_subsection();

    CHECK(parsed_default == "Interface controlled growth");
  }
}
