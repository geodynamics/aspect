//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/postprocess_base.h>
#include <aspect/simulator.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    SimulatorAccess<dim>::SimulatorAccess (const Simulator<dim> &simulator_object)
      :
      simulator (&simulator_object)
    {}



    template <int dim>
    double SimulatorAccess<dim>::get_time () const
    {
      return simulator->time;
    }



    template <int dim>
    double SimulatorAccess<dim>::get_timestep_number () const
    {
      return simulator->time_step;
    }



    template <int dim>
    const TrilinosWrappers::MPI::BlockVector &
    SimulatorAccess<dim>::get_stokes_solution () const
    {
      return simulator->stokes_solution;
    }



    template <int dim>
    const TrilinosWrappers::MPI::BlockVector &
    SimulatorAccess<dim>::get_old_stokes_solution () const
    {
      return simulator->old_stokes_solution;
    }



    template <int dim>
    const DoFHandler<dim> &
    SimulatorAccess<dim>::get_stokes_dof_handler () const
    {
      return simulator->stokes_dof_handler;
    }



    template <int dim>
    const TrilinosWrappers::MPI::Vector &
    SimulatorAccess<dim>::get_temperature_solution () const
    {
      return simulator->temperature_solution;
    }



    template <int dim>
    const TrilinosWrappers::MPI::Vector &
    SimulatorAccess<dim>::get_old_temperature_solution () const
    {
      return simulator->old_temperature_solution;
    }



    template <int dim>
    const DoFHandler<dim> &
    SimulatorAccess<dim>::get_temperature_dof_handler () const
    {
      return simulator->temperature_dof_handler;
    }
  }
}


// explicit instantiation
namespace aspect
{
  namespace Postprocess
  {
    template class SimulatorAccess<deal_II_dimension>;
  }
}
