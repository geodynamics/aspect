// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2020-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include "common.h"
#include <aspect/particle/property/interface.h>
#include <aspect/particle/manager.h>
#include <aspect/simulator.h>
#include <deal.II/base/parameter_handler.h>

TEST_CASE("Particle Manager plugin names")
{
  // The Particle::Property::Manager needs to know about the Simulator object,
  // so we need to create a dummy one to initialize the manager. That, in
  // turn, needs to be initialized with a ParameterHandler, so we create one
  // and fill it with the minimum information needed to create a Simulator object.
  aspect::ParameterHandler sim_prm;
  aspect::Simulator<2>::declare_parameters(sim_prm, 0);
  {
    sim_prm.enter_subsection("Geometry model");
    sim_prm.set("Model name", "box");
    sim_prm.leave_subsection();
  }
  {
    sim_prm.enter_subsection("Material model");
    sim_prm.set("Model name", "simple");
    sim_prm.leave_subsection();
  }
  {
    sim_prm.enter_subsection("Gravity model");
    sim_prm.set("Model name", "vertical");
    sim_prm.leave_subsection();
  }
  aspect::Simulator<2> dummy_sim (MPI_COMM_SELF, sim_prm);

  // Now create a Particle::Property::Manager and initialize it with the dummy simulator.
  dealii::ParameterHandler prm;
  aspect::Particle::Property::Manager<2> manager;
  manager.initialize_simulator (dummy_sim);

  // The property manager needs to know about the integrator, which is declared in World
  aspect::Particle::Manager<2>::declare_parameters(prm);

  prm.enter_subsection("Particles");
  manager.declare_parameters(prm);
  prm.set("List of particle properties","composition, position");
  manager.parse_parameters(prm);
  prm.leave_subsection();

  // existing and listed pluring
  REQUIRE(manager.plugin_name_exists("composition") == true);
  REQUIRE(manager.plugin_name_exists("position") == true);

  // existing but not listed plugin
  REQUIRE(manager.plugin_name_exists("pT path") == false);

  // non-existed plugin
  REQUIRE(manager.plugin_name_exists("non-existent plugin name") == false);

  // check that one is before the other
  REQUIRE(manager.check_plugin_order("composition", "position") == true);
  REQUIRE(manager.check_plugin_order("position", "composition") == false);

  // Check the plugin indices
  REQUIRE(manager.get_plugin_index_by_name("composition") == 0);
  REQUIRE(manager.get_plugin_index_by_name("position") == 1);
}
