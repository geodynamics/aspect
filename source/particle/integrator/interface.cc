/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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

#include <aspect/particle/integrator/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/std_cxx1x/tuple.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
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
      bool
      Interface<dim>::new_integration_step()
      {
        return false;
      }

      template <int dim>
      unsigned int
      Interface<dim>::get_data_size() const
      {
        return 0;
      }

      template <int dim>
      const void *
      Interface<dim>::read_data(const void *data,
                                const types::particle_index /*id*/)
      {
        return data;
      }

      template <int dim>
      void *
      Interface<dim>::write_data(void *data,
                                 const types::particle_index /*id*/) const
      {
        return data;
      }



// -------------------------------- Deal with registering models and automating
// -------------------------------- their setup and selection at run time

      namespace
      {
        std_cxx1x::tuple
        <void *,
        void *,
        aspect::internal::Plugins::PluginList<Interface<2> >,
        aspect::internal::Plugins::PluginList<Interface<3> > > registered_plugins;
      }



      template <int dim>
      void
      register_particle_integrator (const std::string &name,
                                    const std::string &description,
                                    void (*declare_parameters_function) (ParameterHandler &),
                                    Interface<dim> *(*factory_function) ())
      {
        std_cxx1x::get<dim>(registered_plugins).register_plugin (name,
                                                                 description,
                                                                 declare_parameters_function,
                                                                 factory_function);
      }


      template <int dim>
      Interface<dim> *
      create_particle_integrator (ParameterHandler &prm)
      {
        std::string name;
        prm.enter_subsection ("Postprocess");
        {
          prm.enter_subsection ("Particles");
          {
            name = prm.get ("Integration scheme");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();

        return std_cxx1x::get<dim>(registered_plugins).create_plugin (name,
                                                                      "Particle::Integrator name");
      }

      template <int dim>
      void
      declare_parameters (ParameterHandler &prm)
      {
        // declare the entry in the parameter file
        prm.enter_subsection ("Postprocess");
        {
          prm.enter_subsection ("Particles");
          {
            const std::string pattern_of_names
              = std_cxx1x::get<dim>(registered_plugins).get_pattern_of_names ();

            prm.declare_entry ("Integration scheme", "rk2",
                               Patterns::Selection (pattern_of_names),
                               "This parameter is used to decide which method to "
                               "use to solve the equation that describes the position "
                               "of particles, i.e., $\\frac{d}{dt}\\mathbf x_k(t) = "
                               "\\mathbf u(\\mathbf x_k(t),t)$, where $k$ is an index "
                               "that runs over all particles, and $\\mathbf u(\\mathbf x,t)$ "
                               "is the velocity field that results from the Stokes "
                               "equations."
                               "\n\n"
                               "In practice, the exact velocity $\\mathbf u(\\mathbf x,t)$ "
                               "is of course not available, but only a numerical "
                               "approximation $\\mathbf u_h(\\mathbf x,t)$. Furthermore, "
                               "this approximation is only available at discrete time steps, "
                               "$\\mathbf u^n(\\mathbf x)=\\mathbf u(\\mathbf x,t^n)$, and "
                               "these need to be interpolated between time steps if the "
                               "integrator for the equation above requires an evaluation at "
                               "time points between the discrete time steps. If we denote this "
                               "interpolation in time by $\\tilde{\\mathbf u}_h(\\mathbf x,t)$ "
                               "where $\\tilde{\\mathbf u}_h(\\mathbf x,t^n)="
                               "\\mathbf u^n(\\mathbf x)$, then the equation the differential "
                               "equation solver really tries to solve is "
                               "$\\frac{d}{dt}\\tilde{\\mathbf x}_k(t) = "
                               " \\tilde{\\mathbf u}_h(\\mathbf x_k(t),t)$."
                               "\n\n"
                               "As a consequence of these considerations, if you try to "
                               "assess convergence properties of an ODE integrator -- for "
                               "example to verify that the RK4 integrator converges with "
                               "fourth order --, it is important to recall that the "
                               "integrator may not solve the equation you think it "
                               "solves. If, for example, we call the numerical solution "
                               "of the ODE $\\tilde{\\mathbf x}_{k,h}(t)$, then the "
                               "error will typically satisfy a relationship like "
                               "$$"
                               "  \\| \\tilde{\\mathbf x}_k(T) - \\tilde{\\mathbf x}_{k,h}(T) \\|"
                               "  \\le"
                               "  C(T) \\Delta t^p"
                               "$$ "
                               "where $\\Delta t$ is the time step and $p$ the convergence order "
                               "of the method, and $C(T)$ is a (generally unknown) constant "
                               "that depends on the end time $T$ at which one compares the "
                               "solutions. On the other hand, an analytically computed "
                               "trajectory would likely use the \\textit{exact} velocity, "
                               "and one may be tempted to compute "
                               "$\\| \\mathbf x_k(T) - \\tilde{\\mathbf x}_{k,h}(T) \\|$, "
                               "but this quantity will, in the best case, only satisfy an "
                               "estimate of the form "
                               "$$"
                               "  \\| \\mathbf x_k(T) - \\tilde{\\mathbf x}_{k,h}(T) \\|"
                               "  \\le"
                               "  C_1(T) \\Delta t^p"
                               "  + C_2(T) \\| \\mathbf u-\\mathbf u_h \\|"
                               "  + C_3(T) \\| \\mathbf u_h-\\tilde{\\mathbf u}_h \\|"
                               "$$ "
                               "with appropriately chosen norms for the second and third "
                               "term. These second and third terms typically converge to "
                               "zero at relatively low rates (compared to the order $p$ of "
                               "the integrator, which can often be chosen relatively high) "
                               "in the mesh size $h$ and the time step size $\\\\Delta t$, "
                               "limiting the overall accuracy of the ODE integrator."
                               "\n\n"
                               "Select one of the following models:\n\n"
                               +
                               std_cxx1x::get<dim>(registered_plugins).get_description_string());
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();

        std_cxx1x::get<dim>(registered_plugins).declare_parameters (prm);
      }



      template <int dim>
      void
      write_plugin_graph (std::ostream &out)
      {
        std_cxx11::get<dim>(registered_plugins).write_plugin_graph ("Particle integrator interface",
                                                                    out);
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<Particle::Integrator::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<Particle::Integrator::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<Particle::Integrator::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<Particle::Integrator::Interface<3> >::plugins = 0;
    }
  }

  namespace Particle
  {
    namespace Integrator
    {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_particle_integrator<dim> (const std::string &, \
                                     const std::string &, \
                                     void ( *) (ParameterHandler &), \
                                     Interface<dim> *( *) ()); \
  \
  template  \
  void \
  declare_parameters<dim> (ParameterHandler &); \
  \
  template \
  void \
  write_plugin_graph<dim> (std::ostream &); \
  \
  template \
  Interface<dim> * \
  create_particle_integrator<dim> (ParameterHandler &prm);

      ASPECT_INSTANTIATE(INSTANTIATE)
    }
  }
}


