#include <aspect/simulator_signals.h>
#include <aspect/simulator_access.h>
#include <aspect/particle/world.h>
#include <deal.II/particles/particle_handler.h>

#include <iostream>

using namespace aspect;


template <int dim>
void post_mesh_deformation (const SimulatorAccess<dim> &simulator_access)
{
  // Get the reference location of the one particle,
  // and check that it has the correct value.                         
  const Point<dim> predicted_particle_reference_location = simulator_access.get_particle_world().get_particle_handler().begin()->get_reference_location();
  Point<dim> correct_particle_reference_location; 
  correct_particle_reference_location[0] = 0.5;
  correct_particle_reference_location[dim-1] = 1./3.;
  AssertThrow(predicted_particle_reference_location == correct_particle_reference_location,
    ExcMessage("The particle reference location is not correct."));
}


template <int dim>
void signal_connector (SimulatorSignals<dim> &signals)
{
  signals.post_mesh_deformation.connect (&post_mesh_deformation<dim>);
}


ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
