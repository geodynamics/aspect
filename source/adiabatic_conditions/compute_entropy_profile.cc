/*
  Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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


#include <aspect/adiabatic_conditions/compute_entropy_profile.h>
#include <aspect/gravity_model/interface.h>

#include <deal.II/base/signaling_nan.h>


namespace aspect
{
  namespace AdiabaticConditions
  {
    template <int dim>
    ComputeEntropyProfile<dim>::ComputeEntropyProfile()
      :
      initialized(false)
    {}



    template <int dim>
    void
    ComputeEntropyProfile<dim>::initialize()
    {
      if (initialized)
        return;

      // The simulator only keeps the initial conditions around for
      // the first time step. As a consequence, we have to save a
      // shared pointer to that object ourselves the first time we get
      // here.
      if (initial_composition_manager == nullptr)
        initial_composition_manager = this->get_initial_composition_manager_pointer();

      temperatures.resize(n_points, numbers::signaling_nan<double>());
      pressures.resize(n_points, numbers::signaling_nan<double>());
      densities.resize(n_points, numbers::signaling_nan<double>());

      delta_z = this->get_geometry_model().maximal_depth() / (n_points-1);

      MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(1, this->n_compositional_fields());
      this->get_material_model().create_additional_named_outputs (out);

      const std::shared_ptr<MaterialModel::PrescribedTemperatureOutputs<dim>> prescribed_temperature_out
        = out.template get_additional_output_object<MaterialModel::PrescribedTemperatureOutputs<dim>>();

      // check if the material model computes prescribed temperature outputs
      AssertThrow(prescribed_temperature_out != nullptr,
                  ExcMessage("The material model you use does not provide "
                             "PrescribedTemperatureOutputs, which is required "
                             "for this adiabatic conditions plugin."));

      const std::vector<unsigned int> &entropy_indices = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::entropy);

      AssertThrow(entropy_indices.size() >= 1,
                  ExcMessage("The 'compute entropy' adiabatic conditions plugin "
                             "requires at least one field of type 'entropy'."));

      // We only need the material model to compute the density
      // and prescribed temperature. Unfortunately 'additional_outputs' computes
      // a lot of other outputs as well, but we have currently no way to prevent this.
      in.requested_properties = MaterialModel::MaterialProperties::density | MaterialModel::MaterialProperties::additional_outputs;

      // No deformation on the reference profile
      in.velocity[0] = Tensor <1,dim> ();
      in.strain_rate[0] = SymmetricTensor<2,dim>();

      // Strictly speaking the temperature on the adiabat should be defined by the pressure
      // and entropy alone. However, we only compute the temperature in the call to the material
      // model below. This is a problem if we use multiple material models via a compositing
      // material model. We cannot copy the computed temperature into the MaterialModelInputs
      // until after all material models have been evaluated, and material models that do not
      // rely on entropy will access the temperature beforehand to compute their properties.
      // Therefore, we provide a reasonable temperature guess anyway. It is important to note
      // that all properties that are relevant for the equation of state will be provided by
      // the entropy material model, and will therefore not be affected by this temperature.
      in.temperature[0] = this->get_adiabatic_surface_temperature();

      // Set all chemical composition to the initial composition, except the entropies, which
      // are set to the surface entropy (since entropy is constant along an adiabat).
      // Note, that if there a multiple entropy components they could have different entropies.
      // However, since we are only interested in setting the
      // equilibrated entropy, we do not need to compute the individual entropies for all components,
      // and instead set all components to the equilibrated value.
      // TODO : provide more ways to specify compositional fields like in compute_profile.cc
      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        {
          if (this->introspection().get_composition_descriptions()[c].type == CompositionalFieldDescription::entropy)
            in.composition[0][c] = surface_entropy;
          else
            in.composition[0][c] = initial_composition_manager->initial_composition(this->get_geometry_model().representative_point(0), c);
        }

      // Check whether gravity is pointing up / out or down / in. In the normal case it should
      // point down / in and therefore gravity should be positive, leading to increasing
      // adiabatic pressures and temperatures with depth. In some cases it will point up / out
      // (e.g. for backward advection), in which case the pressures and temperatures should
      // decrease with depth and therefore gravity has to be negative in the following equations.
      const Tensor <1,dim> g = this->get_gravity_model().gravity_vector(this->get_geometry_model().representative_point(0));
      const Point<dim> point_surf = this->get_geometry_model().representative_point(0);
      const Point<dim> point_bot = this->get_geometry_model().representative_point(this->get_geometry_model().maximal_depth());
      const int gravity_direction =  (g * (point_bot - point_surf) >= 0) ?
                                     1 :
                                     -1;

      // now integrate downward using the explicit Euler method for simplicity
      //
      // note: p'(z) = rho(p,T) * |g|
      //       T(z) = look up for reference entropy and current p(z)
      for (unsigned int i=0; i<n_points; ++i)
        {
          if (i==0)
            {
              pressures[0] = this->get_surface_pressure();
            }
          else
            {
              // use material properties calculated at i-1
              const double density = out.densities[0];
              // get the magnitude of gravity. we assume
              // that gravity always points along the depth direction. this
              // may not strictly be true always but is likely a good enough
              // approximation here.
              const double gravity = gravity_direction * this->get_gravity_model().gravity_vector(in.position[0]).norm();

              pressures[i] = pressures[i-1] + density * gravity * delta_z;
            }

          const double z = static_cast<double>(i)/static_cast<double>(n_points-1)*this->get_geometry_model().maximal_depth();
          const Point<dim> representative_point = this->get_geometry_model().representative_point (z);

          in.position[0] = representative_point;
          in.pressure[0] = pressures[i];
          this->get_material_model().evaluate(in, out);

          densities[i] = out.densities[0];
          temperatures[i] = prescribed_temperature_out->prescribed_temperature_outputs[0];
          in.temperature[0] = temperatures[i];
        }

      if (gravity_direction == 1 && this->get_surface_pressure() >= 0)
        {
          Assert (*std::min_element (pressures.begin(), pressures.end()) >=
                  -std::numeric_limits<double>::epsilon() * pressures.size(),
                  ExcMessage("Adiabatic ComputeProfile encountered a negative pressure of "
                             + dealii::Utilities::to_string(*std::min_element (pressures.begin(), pressures.end()))));
        }
      else if (gravity_direction == -1 && this->get_surface_pressure() <= 0)
        {
          Assert (*std::max_element (pressures.begin(), pressures.end()) <=
                  std::numeric_limits<double>::epsilon() * pressures.size(),
                  ExcMessage("Adiabatic ComputeProfile encountered a positive pressure of "
                             + dealii::Utilities::to_string(*std::max_element (pressures.begin(), pressures.end()))));
        }

      Assert (*std::min_element (temperatures.begin(), temperatures.end()) >=
              -std::numeric_limits<double>::epsilon() * temperatures.size(),
              ExcMessage("Adiabatic ComputeProfile encountered a negative temperature."));

      initialized = true;
    }



    template <int dim>
    bool
    ComputeEntropyProfile<dim>::is_initialized() const
    {
      return initialized;
    }



    template <int dim>
    double ComputeEntropyProfile<dim>::pressure (const Point<dim> &p) const
    {
      return get_property(p,pressures);
    }



    template <int dim>
    double ComputeEntropyProfile<dim>::temperature (const Point<dim> &p) const
    {
      return get_property(p,temperatures);
    }



    template <int dim>
    double ComputeEntropyProfile<dim>::density (const Point<dim> &p) const
    {
      return get_property(p,densities);
    }



    template <int dim>
    double ComputeEntropyProfile<dim>::density_derivative (const Point<dim> &p) const
    {
      const double z = this->get_geometry_model().depth(p);

      if (z >= this->get_geometry_model().maximal_depth())
        {
          Assert (z <= this->get_geometry_model().maximal_depth() + delta_z,
                  ExcInternalError());
          return (densities.back() - densities[densities.size()-2]) / delta_z;
        }

      if (z < 0)
        {
          Assert (z >= -delta_z, ExcInternalError());
          return (densities[1] - densities.front()) / delta_z;
        }

      // if z/delta_z is within [k-eps, k+eps] of a whole number k, round it down to k-1
      const unsigned int i = static_cast<unsigned int>((z/delta_z) * (1. - 2. * std::numeric_limits<double>::epsilon()));
      Assert (i < densities.size() - 1, ExcInternalError());

      return (densities[i+1]-densities[i])/delta_z;
    }



    template <int dim>
    double ComputeEntropyProfile<dim>::get_property (const Point<dim> &p,
                                                     const std::vector<double> &property) const
    {
      const double z = this->get_geometry_model().depth(p);

      if (z >= this->get_geometry_model().maximal_depth())
        {
          Assert (z <= this->get_geometry_model().maximal_depth() + delta_z,
                  ExcInternalError());
          return property.back();
        }

      if (z <= 0)
        {
          Assert (z >= -delta_z, ExcInternalError());
          return property.front();
        }

      const double floating_index = z/delta_z;
      const unsigned int i = static_cast<unsigned int>(floating_index);

      // If p is close to an existing value use that one. This prevents
      // asking for values at i+1 while initializing i+1 (when p is at the
      // depth of index i).
      if (std::abs(floating_index-std::floor(floating_index+0.5)) < 1e-6)
        return property[i];

      Assert (i+1 < property.size(), ExcInternalError());

      // now do the linear interpolation
      const double d = floating_index - i;
      Assert ((d>=0) && (d<=1), ExcInternalError());

      return d*property[i+1] + (1.-d)*property[i];
    }




    template <int dim>
    void
    ComputeEntropyProfile<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Adiabatic conditions model");
      {
        prm.enter_subsection("Compute entropy profile");
        {
          prm.declare_entry ("Number of points", "1000",
                             Patterns::Integer (5),
                             "The number of points we use to compute the adiabatic "
                             "profile. The higher the number of points, the more accurate "
                             "the downward integration from the adiabatic surface "
                             "conditions will be.");

          prm.declare_entry ("Surface entropy", "0",
                             Patterns::Double(),
                             "The surface entropy for the profile.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    ComputeEntropyProfile<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Adiabatic conditions model");
      {
        prm.enter_subsection("Compute entropy profile");
        {
          n_points = prm.get_integer ("Number of points");
          surface_entropy = prm.get_double ("Surface entropy");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace AdiabaticConditions
  {
    ASPECT_REGISTER_ADIABATIC_CONDITIONS_MODEL(ComputeEntropyProfile,
                                               "compute entropy profile",
                                               "A model in which the adiabatic profile is "
                                               "calculated by solving the hydrostatic equations for "
                                               "pressure and entropy in depth. "
                                               "Of course the entropy along an adiabat is constant. "
                                               "This plugin requires the material model to provide an "
                                               "additional output object of type PrescribedTemperatureOutputs. "
                                               "It also requires that there is a compositional field of type "
                                               "'entropy' that represents the entropy of the material.")
  }
}
