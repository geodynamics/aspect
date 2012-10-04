/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id: interface.cc 1019 2012-05-19 11:37:00Z bangerth $  */


#include <aspect/simulator.h>

#include <typeinfo>

namespace aspect
{


  template <int dim>
  SimulatorAccess<dim>::~SimulatorAccess ()
  {}



  template <int dim>
  void
  SimulatorAccess<dim>::initialize (const Simulator<dim> &simulator_object)
  {
    simulator = SmartPointer<const Simulator<dim> > (&simulator_object, typeid(*this).name());
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
  double
  SimulatorAccess<dim>::get_adiabatic_surface_temperature () const
  {
    return simulator->parameters.adiabatic_surface_temperature;
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_refinement_criteria (Vector<float> &estimated_error_per_cell) const
  {
    simulator->compute_refinement_criterion(estimated_error_per_cell);
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
  void
  SimulatorAccess<dim>::get_depth_average_field(const unsigned int block_number,
						std::vector<double> &values) const
  {
    simulator->compute_depth_average_field(block_number, values);
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
  void
  SimulatorAccess<dim>::get_Vs_anomaly(Vector<float> &values) const
  {
    simulator->compute_Vs_anomaly(values);
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_Vp_anomaly(Vector<float> &values) const
  {
    simulator->compute_Vp_anomaly(values);
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
  const AdiabaticConditions<dim> &
  SimulatorAccess<dim>::get_adiabatic_conditions () const
  {
    return *simulator->adiabatic_conditions.get();
  }
}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template class SimulatorAccess<dim>; \

  ASPECT_INSTANTIATE(INSTANTIATE)
}
