//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/postprocess_base.h>
#include <aspect/simulator.h>

#include <typeinfo>


namespace aspect
{
  namespace Postprocess
  {
// ------------------------------ Interface -----------------------------

    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void 
    Interface<dim>::declare_parameters (ParameterHandler &)
    {}



    template <int dim>
    void 
    Interface<dim>::parse_parameters (ParameterHandler &)
    {}



    template <int dim>
    void 
    Interface<dim>::save (std::map<std::string,std::string> &) const
    {}


    template <int dim>
    void
    Interface<dim>::load (const std::map<std::string,std::string> &)
    {}


// ------------------------------ SimulatorAccess -----------------------

    template <int dim>
    SimulatorAccess<dim>::SimulatorAccess (const Simulator<dim> &simulator_object)
      :
      simulator (&simulator_object, typeid(*this).name())
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
    const parallel::distributed::Triangulation<dim> &
    SimulatorAccess<dim>::get_triangulation () const
    {
      return simulator->triangulation;
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


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    template class Interface<deal_II_dimension>;
    template class SimulatorAccess<deal_II_dimension>;
  }
}
