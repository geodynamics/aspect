/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/adiabatic_conditions/compute_profile.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/initial_composition/interface.h>

#include <deal.II/base/signaling_nan.h>


namespace aspect
{
  namespace AdiabaticConditions
  {
    template <int dim>
    ComputeProfile<dim>::ComputeProfile()
      :
      initialized(false),
      surface_condition_function(2)
    {}



    template <int dim>
    void
    ComputeProfile<dim>::update()
    {
      if (use_surface_condition_function)
        {
          initialized = false;
          surface_condition_function.set_time(this->get_time());
          initialize();
        }
    }


    template <int dim>
    void
    ComputeProfile<dim>::initialize()
    {
      if (initialized)
        return;

      temperatures.resize(n_points, numbers::signaling_nan<double>());
      pressures.resize(n_points, numbers::signaling_nan<double>());
      densities.resize(n_points, numbers::signaling_nan<double>());

      delta_z = this->get_geometry_model().maximal_depth() / (n_points-1);

      MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(1, this->n_compositional_fields());

      // Constant properties on the reference profile
      in.requested_properties = MaterialModel::MaterialProperties::equation_of_state_properties;
      in.velocity[0] = Tensor <1,dim> ();

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
      //       T'(z) = alpha |g| T / C_p
      for (unsigned int i=0; i<n_points; ++i)
        {
          if (i==0)
            {
              if (!use_surface_condition_function)
                {
                  pressures[0] = this->get_surface_pressure();
                  temperatures[0] = this->get_adiabatic_surface_temperature();
                }
              else
                {
                  pressures[0] = surface_condition_function.value(Point<1>(0.0),0);
                  temperatures[0] = surface_condition_function.value(Point<1>(0.0),1);
                }
            }
          else
            {
              // use material properties calculated at i-1
              const double density = out.densities[0];
              const double alpha = out.thermal_expansion_coefficients[0];
              // Handle the case that cp is zero (happens in simple Stokes test problems like sol_cx). By setting
              // 1/cp = 0.0 we will have a constant temperature profile with depth.
              const double one_over_cp = (out.specific_heat[0]>0.0) ? 1.0/out.specific_heat[0] : 0.0;
              // get the magnitude of gravity. we assume
              // that gravity always points along the depth direction. this
              // may not strictly be true always but is likely a good enough
              // approximation here.
              const double gravity = gravity_direction * this->get_gravity_model().gravity_vector(in.position[0]).norm();

              pressures[i] = pressures[i-1] + density * gravity * delta_z;
              temperatures[i] = (this->include_adiabatic_heating())
                                ?
                                temperatures[i-1] * (1 + alpha * gravity * delta_z * one_over_cp)
                                :
                                temperatures[0];
            }

          const double z = double(i)/double(n_points-1)*this->get_geometry_model().maximal_depth();
          const Point<dim> representative_point = this->get_geometry_model().representative_point (z);
          const Tensor <1,dim> g = this->get_gravity_model().gravity_vector(representative_point);

          in.position[0] = representative_point;
          in.temperature[0] = temperatures[i];
          in.pressure[0] = pressures[i];

          // we approximate the pressure gradient by extrapolating the values
          // from the two points above
          if (i>0)
            in.pressure_gradient[0] = g/(g.norm() != 0.0 ? g.norm() : 1.0)
                                      * (pressures[i] - pressures[i-1]) / delta_z;
          else
            in.pressure_gradient[0] = Tensor <1,dim> ();

          if (reference_composition == initial_composition)
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[0][c] = this->get_initial_composition_manager().initial_composition(representative_point, c);
          else if (reference_composition == reference_function)
            {
              const double depth = this->get_geometry_model().depth(representative_point);
              const Point<1> p(depth);
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                in.composition[0][c] = composition_function->value(p, c);
            }
          else
            AssertThrow(false,ExcNotImplemented());

          this->get_material_model().evaluate(in, out);

          densities[i] = out.densities[0];
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
    ComputeProfile<dim>::is_initialized() const
    {
      return initialized;
    }



    template <int dim>
    double ComputeProfile<dim>::pressure (const Point<dim> &p) const
    {
      return get_property(p,pressures);
    }



    template <int dim>
    double ComputeProfile<dim>::temperature (const Point<dim> &p) const
    {
      return get_property(p,temperatures);
    }



    template <int dim>
    double ComputeProfile<dim>::density (const Point<dim> &p) const
    {
      return get_property(p,densities);
    }



    template <int dim>
    double ComputeProfile<dim>::density_derivative (const Point<dim> &p) const
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
    double ComputeProfile<dim>::get_property (const Point<dim> &p,
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
    ComputeProfile<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Adiabatic conditions model");
      {
        prm.enter_subsection("Compute profile");
        {
          Functions::ParsedFunction<1>::declare_parameters (prm, 1);
          prm.declare_entry("Composition reference profile","initial composition",
                            Patterns::Selection("initial composition|function"),
                            "Select how the reference profile for composition "
                            "is computed. This profile is used to evaluate the "
                            "material model, when computing the pressure and "
                            "temperature profile.");
          prm.declare_entry ("Number of points", "1000",
                             Patterns::Integer (5),
                             "The number of points we use to compute the adiabatic "
                             "profile. The higher the number of points, the more accurate "
                             "the downward integration from the adiabatic surface "
                             "temperature will be.");
          prm.declare_entry ("Use surface condition function", "false",
                             Patterns::Bool(),
                             "Whether to use the 'Surface condition function' to determine surface "
                             "conditions, or the 'Adiabatic surface temperature' and 'Surface pressure' "
                             "parameters. If this is set to true the reference profile is updated "
                             "every timestep. The function expression of the function should be "
                             "independent of space, but can depend on time 't'. The function must "
                             "return two components, the first one being reference surface pressure, "
                             "the second one being reference surface temperature.");

          prm.enter_subsection("Surface condition function");
          {
            Functions::ParsedFunction<1>::declare_parameters (prm, 2);
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    ComputeProfile<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Adiabatic conditions model");
      {
        prm.enter_subsection("Compute profile");
        {
          const std::string composition_profile = prm.get("Composition reference profile");

          if (composition_profile == "initial composition")
            reference_composition = initial_composition;
          else if (composition_profile == "function")
            reference_composition = reference_function;
          else
            AssertThrow(false, ExcNotImplemented());

          if ((this->n_compositional_fields() > 0) && (reference_composition == reference_function))
            {
              composition_function
                = std::make_unique<Functions::ParsedFunction<1>>(this->n_compositional_fields());
              try
                {
                  composition_function->parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "ERROR: FunctionParser failed to parse\n"
                            << "\t'Adiabatic conditions model.compute profile'\n"
                            << "with expression\n"
                            << "\t'" << prm.get("Function expression") << "'"
                            << "More information about the cause of the parse error \n"
                            << "is shown below.\n";
                  throw;
                }
            }

          n_points = prm.get_integer ("Number of points");
          use_surface_condition_function = prm.get_bool("Use surface condition function");
          if (use_surface_condition_function)
            {
              prm.enter_subsection("Surface condition function");
              try
                {
                  surface_condition_function.parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "ERROR: FunctionParser failed to parse\n"
                            << "\t'Adiabatic conditions model.Initial profile.Surface condition function'\n"
                            << "with expression\n"
                            << "\t'" << prm.get("Function expression") << "'"
                            << "More information about the cause of the parse error \n"
                            << "is shown below.\n";
                  throw;
                }
              prm.leave_subsection();
            }
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
    ASPECT_REGISTER_ADIABATIC_CONDITIONS_MODEL(ComputeProfile,
                                               "compute profile",
                                               "A model in which the adiabatic profile is "
                                               "calculated by solving the hydrostatic equations for "
                                               "pressure and temperature in depth. "
                                               "The gravity is assumed to be in depth direction "
                                               "and the composition is either given by the initial "
                                               "composition at reference points or computed "
                                               "as a reference depth-function. "
                                               "All material parameters are computed by the "
                                               "material model plugin. The surface conditions are "
                                               "either constant or changing over time as prescribed "
                                               "by a user-provided function.")
  }
}
