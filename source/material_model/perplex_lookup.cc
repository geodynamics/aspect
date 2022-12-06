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

#include <aspect/material_model/perplex_lookup.h>

#ifdef ASPECT_WITH_PERPLEX
extern "C" {
#include <perplex_c.h>
}
#include <deal.II/base/multithread_info.h>
#endif

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    PerpleXLookup<dim>::initialize()
    {
#ifdef ASPECT_WITH_PERPLEX
      AssertThrow(dealii::MultithreadInfo::is_running_single_threaded(),
                  ExcMessage("The PerpleXLookup MaterialModel only works in single threaded mode (do not use -j)!"));

      ini_phaseq(perplex_file_name.c_str()); // this line initializes meemum
#else
      Assert (false, ExcMessage("ASPECT has not been compiled with the PerpleX libraries"));
#endif
    }

    template <int dim>
    bool
    PerpleXLookup<dim>::
    is_compressible () const
    {
      return true;
    }

    template <int dim>
    void
    PerpleXLookup<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      /* Instead of evaluating at every quadrature point per cell,
       * we here average the P, T and X values, and evaluate once.
       * This is much quicker than evaluating at all quadrature
       * points, and if the grid is fine, it should be a reasonable
       * approximation
       */

#ifdef ASPECT_WITH_PERPLEX
      std::vector<double> wtphases(p_size_phases);
      std::vector<double> cphases(p_size_phases * p_size_components);
      std::vector<char> namephases(p_size_phases * p_pname_len);
      std::vector<double> sysprop(p_size_sysprops);

      int phaseq_dbg = 0;

      unsigned int n_quad = in.n_evaluation_points(); // number of quadrature points in cell
      unsigned int n_comp = in.composition[0].size(); // number of components in rock

      const double average_temperature = std::min(max_temperature,
                                                  std::max(min_temperature,
                                                           (accumulate( in.temperature.begin(), in.temperature.end(), 0.0) /
                                                            n_quad)));
      const double average_pressure = std::min(max_pressure,
                                               std::max(min_pressure,
                                                        (accumulate( in.pressure.begin(), in.pressure.end(), 0.0) /
                                                         n_quad)));

      std::vector<double> comp;
      comp.resize(n_comp);

      for (unsigned int c=0; c<n_comp; ++c)
        {
          for (unsigned int i=0; i<n_quad; ++i)
            {
              comp[c] += in.composition[i][c];
              out.reaction_terms[i][c] = 0.0;
            }
          comp[c] /= (double)n_quad;
        }

      // Here is the call to PerpleX/meemum
      int nphases;

      phaseq(average_pressure/1.e5, average_temperature,
             n_comp, comp.data(), &nphases, wtphases.data(), cphases.data(),
             sysprop.data(), namephases.data(), phaseq_dbg);

      AssertThrow(!isnan(sysprop[9]) && !isnan(sysprop[11]) && !isnan(sysprop[12]) && !isnan(sysprop[13]),
                  ExcMessage("PerpleX returned NaN for at least one material property at " +
                             std::to_string(average_pressure) +" bar, " +
                             std::to_string(average_temperature) + " K. Aborting. " +
                             "Please adjust the P-T bounds in the parameter file or adjust the PerpleX files."));

      for (unsigned int i=0; i<n_quad; ++i)
        {
          out.viscosities[i] = eta;
          out.thermal_conductivities[i] = k_value;
          out.densities[i] = sysprop[9];
          out.specific_heat[i] = sysprop[11]*(1000./sysprop[16]); // molar Cp * (1000/molar mass) (g)
          out.thermal_expansion_coefficients[i] = sysprop[12];
          out.compressibilities[i] = sysprop[13]*1.e5;
        }
#else
      (void)in;
      (void)out;
      Assert (false, ExcMessage("ASPECT has not been compiled with the PerpleX libraries"));
#endif

    }


    template <int dim>
    void
    PerpleXLookup<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("PerpleX lookup model");
        {

          prm.declare_entry ("PerpleX input file name", "rock.dat",
                             Patterns::Anything (),
                             "The name of the PerpleX input file (should end with .dat).");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0.),
                             "The value of the viscosity $\\eta$. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Minimum material temperature", "0.",
                             Patterns::Double (0.),
                             "The value of the minimum temperature used to query PerpleX. "
                             "Units: \\si{\\kelvin}.");
          prm.declare_entry ("Maximum material temperature", "6000.",
                             Patterns::Double (0.),
                             "The value of the maximum temperature used to query PerpleX. "
                             "Units: \\si{\\kelvin}.");
          prm.declare_entry ("Minimum material pressure", "1.e5",
                             Patterns::Double (0.),
                             "The value of the minimum pressure used to query PerpleX. "
                             "Units: \\si{\\pascal}.");
          prm.declare_entry ("Maximum material pressure", "1.e12",
                             Patterns::Double (0.),
                             "The value of the maximum pressure used to query PerpleX. "
                             "Units: \\si{\\pascal}.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    PerpleXLookup<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("PerpleX lookup model");
        {
          perplex_file_name   = prm.get ("PerpleX input file name");
          eta                 = prm.get_double ("Viscosity");
          k_value             = prm.get_double ("Thermal conductivity");
          min_temperature     = prm.get_double ("Minimum material temperature");
          max_temperature     = prm.get_double ("Maximum material temperature");
          min_pressure        = prm.get_double ("Minimum material pressure");
          max_pressure        = prm.get_double ("Maximum material pressure");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;

      this->model_dependence.density = NonlinearDependence::temperature
                                       | NonlinearDependence::pressure
                                       | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::temperature
                                               | NonlinearDependence::pressure
                                               | NonlinearDependence::compositional_fields;
      this->model_dependence.specific_heat = NonlinearDependence::temperature
                                             | NonlinearDependence::pressure
                                             | NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(PerpleXLookup,
                                   "perplex lookup",
                                   "A material model that has constant values "
                                   "for viscosity and thermal conductivity, and "
                                   "calculates other properties on-the-fly using "
                                   "PerpleX meemum. Compositional fields correspond "
                                   "to the individual components in the order given "
                                   "in the PerpleX file.")
  }
}
