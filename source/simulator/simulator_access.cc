/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/simulator.h>

namespace aspect
{


  template <int dim>
  SimulatorAccess<dim>::~SimulatorAccess ()
  {}



  template <int dim>
  void
  SimulatorAccess<dim>::initialize (const Simulator<dim> &simulator_object)
  {
    simulator = &simulator_object;
  }

  template <int dim>
  const Simulator<dim> &
  SimulatorAccess<dim>::get_simulator() const
  {
    return *simulator;
  }



  template <int dim>
  const Introspection<dim> &
  SimulatorAccess<dim>::introspection () const
  {
    return simulator->introspection;
  }


  template <int dim>
  MPI_Comm SimulatorAccess<dim>::get_mpi_communicator () const
  {
    return simulator->mpi_communicator;
  }

  template <int dim>
  const ConditionalOStream &
  SimulatorAccess<dim>::get_pcout () const
  {
    return simulator->pcout;
  }

  template <int dim>
  double SimulatorAccess<dim>::get_time () const
  {
    return simulator->time;
  }

  template <int dim>
  double SimulatorAccess<dim>::get_timestep () const
  {
    return simulator->time_step;
  }



  template <int dim>
  unsigned int SimulatorAccess<dim>::get_timestep_number () const
  {
    return simulator->timestep_number;
  }



  template <int dim>
  const parallel::distributed::Triangulation<dim> &
  SimulatorAccess<dim>::get_triangulation () const
  {
    return simulator->triangulation;
  }



  template <int dim>
  double
  SimulatorAccess<dim>::get_volume () const
  {
    return simulator->global_volume;
  }



  template <int dim>
  const Mapping<dim> &
  SimulatorAccess<dim>::get_mapping () const
  {
    return simulator->mapping;
  }



  template <int dim>
  std::string
  SimulatorAccess<dim>::get_output_directory () const
  {
    return simulator->parameters.output_directory;
  }



  template <int dim>
  bool
  SimulatorAccess<dim>::convert_output_to_years () const
  {
    return simulator->parameters.convert_to_years;
  }


  template <int dim>
  unsigned int
  SimulatorAccess<dim>::n_compositional_fields () const
  {
    return simulator->parameters.n_compositional_fields;
  }



  template <int dim>
  bool
  SimulatorAccess<dim>::include_adiabatic_heating () const
  {
    return simulator->parameters.include_adiabatic_heating;
  }

  template <int dim>
  bool
  SimulatorAccess<dim>::include_latent_heat () const
  {
    return simulator->parameters.include_latent_heat;
  }

  template <int dim>
  int
  SimulatorAccess<dim>::get_stokes_velocity_degree () const
  {
    return simulator->parameters.stokes_velocity_degree;
  }

  template <int dim>
  double
  SimulatorAccess<dim>::get_adiabatic_surface_temperature () const
  {
    return simulator->parameters.adiabatic_surface_temperature;
  }

  template <int dim>
  double
  SimulatorAccess<dim>::get_surface_pressure () const
  {
    return simulator->parameters.surface_pressure;
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_refinement_criteria (Vector<float> &estimated_error_per_cell) const
  {
    simulator->mesh_refinement_manager.execute (estimated_error_per_cell);
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_artificial_viscosity (Vector<float> &viscosity_per_cell) const
  {
    simulator->get_artificial_viscosity(viscosity_per_cell);
  }

  template <int dim>
  const LinearAlgebra::BlockVector &
  SimulatorAccess<dim>::get_solution () const
  {
    return simulator->solution;
  }



  template <int dim>
  const LinearAlgebra::BlockVector &
  SimulatorAccess<dim>::get_old_solution () const
  {
    return simulator->old_solution;
  }



  template <int dim>
  const DoFHandler<dim> &
  SimulatorAccess<dim>::get_dof_handler () const
  {
    return simulator->dof_handler;
  }



  template <int dim>
  const FiniteElement<dim> &
  SimulatorAccess<dim>::get_fe () const
  {
    return simulator->dof_handler.get_fe();
  }



  template <int dim>
  void
  SimulatorAccess<dim>::get_depth_average_temperature(std::vector<double> &values) const
  {
    simulator->compute_depth_average_field(Simulator<dim>::AdvectionField::temperature(),
                                           values);
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_depth_average_composition(const unsigned int composition_index,
                                                      std::vector<double> &values) const
  {
    // make sure that what we get here is really an index of one of the compositional fields
    AssertIndexRange(composition_index,this->n_compositional_fields());

    simulator->compute_depth_average_field(Simulator<dim>::AdvectionField::composition(composition_index),
                                           values);
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_depth_average_viscosity(std::vector<double> &values) const
  {
    simulator->compute_depth_average_viscosity(values);
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_depth_average_velocity_magnitude(std::vector<double> &values) const
  {
    simulator->compute_depth_average_velocity_magnitude(values);
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_depth_average_sinking_velocity(std::vector<double> &values) const
  {
    simulator->compute_depth_average_sinking_velocity(values);
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_depth_average_Vs(std::vector<double> &values) const
  {
    simulator->compute_depth_average_Vs(values);
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_depth_average_Vp(std::vector<double> &values) const
  {
    simulator->compute_depth_average_Vp(values);
  }

  template <int dim>
  const MaterialModel::Interface<dim> &
  SimulatorAccess<dim>::get_material_model () const
  {
    return *simulator->material_model.get();
  }



  template <int dim>
  const BoundaryTemperature::Interface<dim> &
  SimulatorAccess<dim>::get_boundary_temperature () const
  {
    return *simulator->boundary_temperature.get();
  }


  template <int dim>
  const std::set<types::boundary_id> &
  SimulatorAccess<dim>::get_fixed_temperature_boundary_indicators () const
  {
    return simulator->parameters.fixed_temperature_boundary_indicators;
  }


  template <int dim>
  const std::set<types::boundary_id> &
  SimulatorAccess<dim>::get_fixed_composition_boundary_indicators () const
  {
    return simulator->parameters.fixed_composition_boundary_indicators;
  }


  template <int dim>
  const std::map<types::boundary_id,std_cxx1x::shared_ptr<VelocityBoundaryConditions::Interface<dim> > >
  SimulatorAccess<dim>::get_prescribed_velocity_boundary_conditions () const
  {
    return simulator->velocity_boundary_conditions;
  }


  template <int dim>
  const GeometryModel::Interface<dim> &
  SimulatorAccess<dim>::get_geometry_model () const
  {
    return *simulator->geometry_model.get();
  }

  template <int dim>
  const GravityModel::Interface<dim> &
  SimulatorAccess<dim>::get_gravity_model () const
  {
    return *simulator->gravity_model.get();
  }


  template <int dim>
  const AdiabaticConditions::Interface<dim> &
  SimulatorAccess<dim>::get_adiabatic_conditions () const
  {
    return *simulator->adiabatic_conditions.get();
  }


  template <int dim>
  const InitialConditions::Interface<dim> &
  SimulatorAccess<dim>::get_initial_conditions () const
  {
    return *simulator->initial_conditions.get();
  }


  template <int dim>
  const CompositionalInitialConditions::Interface<dim> &
  SimulatorAccess<dim>::get_compositional_initial_conditions () const
  {
    return *simulator->compositional_initial_conditions.get();
  }

  template <int dim>
  const HeatingModel::Interface<dim> &
  SimulatorAccess<dim>::get_heating_model () const
  {
    return *simulator->heating_model.get();
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_composition_values_at_q_point (const std::vector<std::vector<double> > &composition_values,
                                                           const unsigned int                      q,
                                                           std::vector<double>                    &composition_values_at_q_point)
  {
    for (unsigned int k=0; k < composition_values_at_q_point.size(); ++k)
      composition_values_at_q_point[k] = composition_values[k][q];
  }

}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template class SimulatorAccess<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
