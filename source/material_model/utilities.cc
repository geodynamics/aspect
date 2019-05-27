/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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
#include <aspect/material_model/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>

#include <deal.II/base/exceptions.h>

#include <list>


namespace aspect
{
  namespace MaterialModel
  {
    namespace MaterialUtilities
    {
      std::vector<double>
      compute_volume_fractions(const std::vector<double> &compositional_fields,
                               const ComponentMask &field_mask)
      {
        std::vector<double> volume_fractions(compositional_fields.size()+1);

        // Clip the compositional fields so they are between zero and one,
        // and sum the compositional fields for normalization purposes.
        double sum_composition = 0.0;
        std::vector<double> x_comp = compositional_fields;
        for (unsigned int i=0; i < x_comp.size(); ++i)
          if (field_mask[i] == true)
            {
              x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);
              sum_composition += x_comp[i];
            }

        // Compute background material fraction
        if (sum_composition >= 1.0)
          volume_fractions[0] = 0.0;
        else
          volume_fractions[0] = 1.0 - sum_composition;

        // Compute and possibly normalize volume fractions
        for (unsigned int i=0; i < x_comp.size(); ++i)
          if (field_mask[i] == true)
            {
              if (sum_composition >= 1.0)
                volume_fractions[i+1] = x_comp[i]/sum_composition;
              else
                volume_fractions[i+1] = x_comp[i];
            }

        return volume_fractions;
      }



      CompositionalAveragingOperation
      parse_compositional_averaging_operation (const std::string &parameter_name,
                                               const ParameterHandler &prm)
      {
        CompositionalAveragingOperation averaging_operation;
        if (prm.get (parameter_name) == "harmonic")
          averaging_operation = MaterialUtilities::harmonic;
        else if (prm.get (parameter_name) == "arithmetic")
          averaging_operation = MaterialUtilities::arithmetic;
        else if (prm.get (parameter_name) == "geometric")
          averaging_operation = MaterialUtilities::geometric;
        else if (prm.get (parameter_name) == "maximum composition")
          averaging_operation = MaterialUtilities::maximum_composition;
        else
          {
            AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));

            //We will never get here, but we have to return something so the compiler does not complain
            return MaterialUtilities::harmonic;
          }

        return averaging_operation;
      }



      double
      average_value (const std::vector<double> &volume_fractions,
                     const std::vector<double> &parameter_values,
                     const enum CompositionalAveragingOperation &average_type)
      {
        Assert(volume_fractions.size() == parameter_values.size(),
               ExcMessage ("The volume fractions and parameter values vectors used for averaging "
                           "have to have the same length!"));

        double averaged_parameter = 0.0;

        switch (average_type)
          {
            case arithmetic:
            {
              for (unsigned int i=0; i<volume_fractions.size(); ++i)
                averaged_parameter += volume_fractions[i] * parameter_values[i];
              break;
            }
            case harmonic:
            {
              for (unsigned int i=0; i<volume_fractions.size(); ++i)
                {
                  AssertThrow(parameter_values[i] > 0,
                              ExcMessage ("All parameter values must be greater than 0 for harmonic averaging!"));
                  averaged_parameter += volume_fractions[i]/(parameter_values[i]);
                }
              averaged_parameter = 1.0/averaged_parameter;
              break;
            }
            case geometric:
            {
              for (unsigned int i=0; i<volume_fractions.size(); ++i)
                {
                  AssertThrow(parameter_values[i] > 0,
                              ExcMessage ("All parameter values must be greater than 0 for geometric averaging!"));
                  averaged_parameter += volume_fractions[i] * std::log(parameter_values[i]);
                }
              averaged_parameter = std::exp(averaged_parameter);
              break;
            }
            case maximum_composition:
            {
              const unsigned int idx = static_cast<unsigned int>(std::max_element( volume_fractions.begin(),
                                                                                   volume_fractions.end() )
                                                                 - volume_fractions.begin());
              averaged_parameter = parameter_values[idx];
              break;
            }
            default:
            {
              AssertThrow(false, ExcNotImplemented());
              break;
            }
          }
        return averaged_parameter;
      }



      DruckerPragerInputs::DruckerPragerInputs(const double cohesion_,
                                               const double friction_angle_,
                                               const double pressure_,
                                               const double effective_strain_rate_,
                                               const double max_yield_strength_)
        :
        cohesion(cohesion_),
        friction_angle(friction_angle_),
        pressure(pressure_),
        effective_strain_rate(effective_strain_rate_),
        max_yield_strength(max_yield_strength_)
      {}


      DruckerPragerOutputs::DruckerPragerOutputs ()
        :
        yield_strength(numbers::signaling_nan<double>()),
        plastic_viscosity(numbers::signaling_nan<double>()),
        viscosity_pressure_derivative(numbers::signaling_nan<double>())
      {}



      template <int dim>
      void
      compute_drucker_prager_yielding (const DruckerPragerInputs &in,
                                       DruckerPragerOutputs &out)
      {
        // plasticity
        const double sin_phi = std::sin(in.friction_angle);
        const double cos_phi = std::cos(in.friction_angle);
        const double strength_inv_part = 1. / (std::sqrt(3.0) * (3.0 + sin_phi));

        out.yield_strength = ( (dim==3)
                               ?
                               ( 6.0 * in.cohesion * cos_phi + 6.0 * in.pressure * sin_phi) * strength_inv_part
                               :
                               in.cohesion * cos_phi + in.pressure * sin_phi);

        out.yield_strength = std::min(out.yield_strength, in.max_yield_strength);

        // Rescale the viscosity back onto the yield surface
        const double strain_rate_effective_inv = 1./(2.*in.effective_strain_rate);
        out.plastic_viscosity = out.yield_strength * strain_rate_effective_inv;

        out.viscosity_pressure_derivative = sin_phi * strain_rate_effective_inv *
                                            (dim == 3
                                             ?
                                             (6.0 * strength_inv_part)
                                             :
                                             1);

        return;
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace MaterialUtilities
    {
#define INSTANTIATE(dim) \
  template \
  void \
  compute_drucker_prager_yielding<dim> (const DruckerPragerInputs &, \
                                        DruckerPragerOutputs &); \
   
      ASPECT_INSTANTIATE(INSTANTIATE)
    }
  }
}
