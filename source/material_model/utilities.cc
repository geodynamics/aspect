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
#include <aspect/simulator_access.h>

#include <aspect/adiabatic_conditions/interface.h>

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
      namespace Lookup
      {
        double
        MaterialLookup::specific_heat(const double temperature,
                                      const double pressure) const
        {
          return value(temperature,pressure,specific_heat_values,interpolation);
        }

        double
        MaterialLookup::density(const double temperature,
                                const double pressure) const
        {
          return value(temperature,pressure,density_values,interpolation);
        }

        double
        MaterialLookup::thermal_expansivity(const double temperature,
                                            const double pressure) const
        {
          return value(temperature,pressure,thermal_expansivity_values,interpolation);
        }

        double
        MaterialLookup::seismic_Vp(const double temperature,
                                   const double pressure) const
        {
          return value(temperature,pressure,vp_values,false);
        }

        double
        MaterialLookup::seismic_Vs(const double temperature,
                                   const double pressure) const
        {
          return value(temperature,pressure,vs_values,false);
        }

        double
        MaterialLookup::enthalpy(const double temperature,
                                 const double pressure) const
        {
          return value(temperature,pressure,enthalpy_values,true);
        }

        double
        MaterialLookup::dHdT (const double temperature,
                              const double pressure) const
        {
          const double h = value(temperature,pressure,enthalpy_values,interpolation);
          const double dh = value(temperature+delta_temp,pressure,enthalpy_values,interpolation);
          return (dh - h) / delta_temp;
        }

        double
        MaterialLookup::dHdp (const double temperature,
                              const double pressure) const
        {
          const double h = value(temperature,pressure,enthalpy_values,interpolation);
          const double dh = value(temperature,pressure+delta_press,enthalpy_values,interpolation);
          return (dh - h) / delta_press;
        }

        std::array<std::pair<double, unsigned int>,2>
        MaterialLookup::enthalpy_derivatives(const std::vector<double> &temperatures,
                                             const std::vector<double> &pressures,
                                             const unsigned int n_substeps) const
        {
          Assert(temperatures.size() == pressures.size(),ExcInternalError());
          const unsigned int n_q_points = temperatures.size();
          unsigned int n_T(0), n_p(0);
          double dHdT(0.0), dHdp(0.0);

          for (unsigned int q=0; q<n_q_points; ++q)
            {
              for (unsigned int p=0; p<n_q_points; ++p)
                {
                  if (std::fabs(temperatures[q] - temperatures[p]) > 100.0 * std::numeric_limits<double>::epsilon() * std::fabs(temperatures[q] + std::numeric_limits<double>::epsilon()))
                    {
                      for (unsigned int substep = 0; substep < n_substeps; ++substep)
                        {
                          const double step_ratio = static_cast<double>(substep)/static_cast<double>(n_substeps);
                          const double step_ratio_next = static_cast<double>(substep+1)/static_cast<double>(n_substeps);

                          const double current_pressure = pressures[q]
                                                          + step_ratio
                                                          * (pressures[p]-pressures[q]);
                          const double T1_substep = temperatures[q]
                                                    + step_ratio
                                                    * (temperatures[p]-temperatures[q]);
                          const double T2_substep = temperatures[q]
                                                    + step_ratio_next
                                                    * (temperatures[p]-temperatures[q]);
                          const double enthalpy2 = enthalpy(T2_substep,current_pressure);
                          const double enthalpy1 = enthalpy(T1_substep,current_pressure);
                          dHdT += (enthalpy2-enthalpy1)/(T2_substep-T1_substep);
                          ++n_T;
                        }
                    }
                  if (std::fabs(pressures[q] - pressures[p]) > 100.0 * std::numeric_limits<double>::epsilon() * std::fabs(pressures[q] + std::numeric_limits<double>::epsilon()))
                    {
                      for (unsigned int substep = 0; substep < n_substeps; ++substep)
                        {
                          const double step_ratio = static_cast<double>(substep)/static_cast<double>(n_substeps);
                          const double step_ratio_next = static_cast<double>(substep+1)/static_cast<double>(n_substeps);
                          const double current_temperature = temperatures[q]
                                                             + step_ratio
                                                             * (temperatures[p]-temperatures[q]);
                          const double p1_substep = pressures[q]
                                                    + step_ratio
                                                    * (pressures[p]-pressures[q]);
                          const double p2_substep = pressures[q]
                                                    + step_ratio_next
                                                    * (pressures[p]-pressures[q]);
                          const double enthalpy2 = enthalpy(current_temperature,p2_substep);
                          const double enthalpy1 = enthalpy(current_temperature,p1_substep);
                          dHdp += (enthalpy2-enthalpy1)/(p2_substep-p1_substep);
                          ++n_p;
                        }
                    }
                }
            }

          if ((n_T > 0)
              && (n_p > 0))
            {
              dHdT /= n_T;
              dHdp /= n_p;
            }

          std::array<std::pair<double, unsigned int>,2> derivatives;
          derivatives[0] = std::make_pair(dHdT,n_T);
          derivatives[1] = std::make_pair(dHdp,n_p);
          return derivatives;
        }

        double
        MaterialLookup::dRhodp (const double temperature,
                                const double pressure) const
        {
          const double rho = value(temperature,pressure,density_values,interpolation);
          const double drho = value(temperature,pressure+delta_press,density_values,interpolation);
          return (drho - rho) / delta_press;
        }

        double
        MaterialLookup::value (const double temperature,
                               const double pressure,
                               const Table<2, double> &values,
                               const bool interpol) const
        {
          const double nT = get_nT(temperature);
          const unsigned int inT = static_cast<unsigned int>(nT);

          const double np = get_np(pressure);
          const unsigned int inp = static_cast<unsigned int>(np);

          Assert(inT<values.n_rows(), ExcMessage("Attempting to look up a temperature value with index greater than the number of rows."));
          Assert(inp<values.n_cols(), ExcMessage("Attempting to look up a pressure value with index greater than the number of columns."));

          if (!interpol)
            return values[inT][inp];
          else
            {
              // compute the coordinates of this point in the
              // reference cell between the data points
              const double xi = nT-inT;
              const double eta = np-inp;

              Assert ((0 <= xi) && (xi <= 1), ExcInternalError());
              Assert ((0 <= eta) && (eta <= 1), ExcInternalError());

              // use these coordinates for a bilinear interpolation
              return ((1-xi)*(1-eta)*values[inT][inp] +
                      xi    *(1-eta)*values[inT+1][inp] +
                      (1-xi)*eta    *values[inT][inp+1] +
                      xi    *eta    *values[inT+1][inp+1]);
            }
        }

        std::array<double,2>
        MaterialLookup::get_pT_steps() const
        {
          std::array<double,2> pt_steps;
          pt_steps[0] = delta_press;
          pt_steps[1] = delta_temp;
          return pt_steps;
        }

        double
        MaterialLookup::get_nT(const double temperature) const
        {
          double bounded_temperature=std::max(min_temp, temperature);
          bounded_temperature=std::min(bounded_temperature, max_temp-delta_temp);

          return (bounded_temperature-min_temp)/delta_temp;
        }

        double
        MaterialLookup::get_np(const double pressure) const
        {
          double bounded_pressure=std::max(min_press, pressure);
          bounded_pressure=std::min(bounded_pressure, max_press-delta_press);

          return (bounded_pressure-min_press)/delta_press;
        }

        HeFESToReader::HeFESToReader(const std::string &material_filename,
                                     const std::string &derivatives_filename,
                                     const bool interpol,
                                     const MPI_Comm &comm)
        {
          /* Initializing variables */
          interpolation = interpol;
          delta_press=numbers::signaling_nan<double>();
          min_press=std::numeric_limits<double>::max();
          max_press=-std::numeric_limits<double>::max();
          delta_temp=numbers::signaling_nan<double>();
          min_temp=std::numeric_limits<double>::max();
          max_temp=-std::numeric_limits<double>::max();
          n_temperature=0;
          n_pressure=0;

          std::string temp;

          // Read material data
          {
            // Read data from disk and distribute among processes
            std::istringstream in(Utilities::read_and_distribute_file_content(material_filename, comm));

            bool parsed_first_column = false;
            unsigned int i = 0;
            double current_pressure = 0.0;
            double old_pressure = -1.0;
            while (!in.eof())
              {
                in >> current_pressure;
                if (in.fail())
                  {
                    in.clear();
                  }

                if (!parsed_first_column)
                  {
                    if (current_pressure > old_pressure)
                      old_pressure = current_pressure;
                    else if (current_pressure <= old_pressure)
                      {
                        n_pressure = i;
                        parsed_first_column = true;
                      }
                  }

                std::getline(in, temp);
                if (in.eof())
                  break;
                i++;
              }

            in.clear();
            in.seekg (0, in.beg);

            n_temperature = i / n_pressure;

            Assert(i == n_temperature * n_pressure,
                   ExcMessage("Material table size not consistent."));

            density_values.reinit(n_temperature,n_pressure);
            thermal_expansivity_values.reinit(n_temperature,n_pressure);
            specific_heat_values.reinit(n_temperature,n_pressure);
            vp_values.reinit(n_temperature,n_pressure);
            vs_values.reinit(n_temperature,n_pressure);
            enthalpy_values.reinit(n_temperature,n_pressure);

            i = 0;
            while (!in.eof())
              {
                double P = 0.0;
                double depth,T;
                double rho,vb,vs,vp,vsq,vpq,h;
                std::string code;
                double alpha = 0.0;
                double cp = 0.0;

                in >> P >> depth >> T;
                if (in.fail())
                  in.clear();
                // conversion from [GPa] to [Pa]
                P *= 1e9;

                min_press=std::min(P,min_press);
                min_temp=std::min(T,min_temp);
                max_temp = std::max(T,max_temp);
                max_press = std::max(P,max_press);

                in >> rho;
                if (in.fail())
                  {
                    in.clear();
                    rho = density_values[(i-1)%n_temperature][(i-1)/n_temperature];
                  }
                else
                  rho *= 1e3; // conversion from [g/cm^3] to [kg/m^3]

                in >> vb;
                if (in.fail())
                  in.clear();

                in >> vs;
                if (in.fail())
                  {
                    in.clear();
                    vs = vs_values[(i-1)%n_temperature][(i-1)/n_temperature];
                  }
                in >> vp;
                if (in.fail())
                  {
                    in.clear();
                    vp = vp_values[(i-1)%n_temperature][(i-1)/n_temperature];
                  }
                in >> vsq >> vpq;

                in >> h;
                if (in.fail())
                  {
                    in.clear();
                    h = enthalpy_values[(i-1)%n_temperature][(i-1)/n_temperature];
                  }
                else
                  h *= 1e6; // conversion from [kJ/g] to [J/kg]

                std::getline(in, temp);
                if (in.eof())
                  break;

                density_values[i/n_pressure][i%n_pressure]=rho;
                thermal_expansivity_values[i/n_pressure][i%n_pressure]=alpha;
                specific_heat_values[i/n_pressure][i%n_pressure]=cp;
                vp_values[i/n_pressure][i%n_pressure]=vp;
                vs_values[i/n_pressure][i%n_pressure]=vs;
                enthalpy_values[i/n_pressure][i%n_pressure]=h;

                i++;
              }

            delta_temp = (max_temp - min_temp) / (n_temperature - 1);
            delta_press = (max_press - min_press) / (n_pressure - 1);

            AssertThrow(max_temp >= 0.0, ExcMessage("Read in of Material header failed (max_temp)."));
            AssertThrow(delta_temp > 0, ExcMessage("Read in of Material header failed (delta_temp)."));
            AssertThrow(n_temperature > 0, ExcMessage("Read in of Material header failed (numtemp)."));
            AssertThrow(max_press >= 0, ExcMessage("Read in of Material header failed (max_press)."));
            AssertThrow(delta_press > 0, ExcMessage("Read in of Material header failed (delta_press)."));
            AssertThrow(n_pressure > 0, ExcMessage("Read in of Material header failed (numpress)."));
          }

          // If requested read derivative data
          if (derivatives_filename != "")
            {
              std::string temp;
              // Read data from disk and distribute among processes
              std::istringstream in(Utilities::read_and_distribute_file_content(derivatives_filename, comm));

              int i = 0;
              while (!in.eof())
                {
                  double P = 0.0;
                  double depth,T;
                  double cp,alpha,alpha_eff;
                  double temp1,temp2;

                  in >> P >> depth >> T;
                  if (in.fail())
                    in.clear();


                  in >> cp;
                  if (in.fail() || (cp <= std::numeric_limits<double>::min()))
                    {
                      in.clear();
                      cp = specific_heat_values[(i-1)%n_temperature][(i-1)/n_temperature];
                    }
                  else
                    cp *= 1e3; // conversion from [J/g/K] to [J/kg/K]

                  in >> alpha >> alpha_eff;
                  if (in.fail() || (alpha_eff <= std::numeric_limits<double>::min()))
                    {
                      in.clear();
                      alpha_eff = thermal_expansivity_values[(i-1)%n_temperature][(i-1)/n_temperature];
                    }
                  else
                    {
                      alpha *= 1e-5;
                      alpha_eff *= 1e-5;
                    }

                  in >> temp1 >> temp2;
                  if (in.fail())
                    in.clear();


                  std::getline(in, temp);
                  if (in.eof())
                    break;

                  specific_heat_values[i/n_pressure][i%n_pressure]=cp;
                  thermal_expansivity_values[i/n_pressure][i%n_pressure]=alpha_eff;

                  i++;
                }
            }
        }

        PerplexReader::PerplexReader(const std::string &filename,
                                     const bool interpol,
                                     const MPI_Comm &comm)
        {
          /* Initializing variables */
          interpolation = interpol;
          delta_press=numbers::signaling_nan<double>();
          min_press=std::numeric_limits<double>::max();
          max_press=-std::numeric_limits<double>::max();
          delta_temp=numbers::signaling_nan<double>();
          min_temp=std::numeric_limits<double>::max();
          max_temp=-std::numeric_limits<double>::max();
          n_temperature=0;
          n_pressure=0;

          std::string temp;
          // Read data from disk and distribute among processes
          std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

          std::getline(in, temp); // eat first line
          std::getline(in, temp); // eat next line
          std::getline(in, temp); // eat next line
          std::getline(in, temp); // eat next line

          in >> min_temp;
          std::getline(in, temp);
          in >> delta_temp;
          std::getline(in, temp);
          in >> n_temperature;
          std::getline(in, temp);
          std::getline(in, temp);
          in >> min_press;
          min_press *= 1e5;  // conversion from [bar] to [Pa]
          std::getline(in, temp);
          in >> delta_press;
          delta_press *= 1e5; // conversion from [bar] to [Pa]
          std::getline(in, temp);
          in >> n_pressure;
          std::getline(in, temp);
          std::getline(in, temp);
          std::getline(in, temp);

          AssertThrow(min_temp >= 0.0, ExcMessage("Read in of Material header failed (mintemp)."));
          AssertThrow(delta_temp > 0, ExcMessage("Read in of Material header failed (delta_temp)."));
          AssertThrow(n_temperature > 0, ExcMessage("Read in of Material header failed (numtemp)."));
          AssertThrow(min_press >= 0, ExcMessage("Read in of Material header failed (min_press)."));
          AssertThrow(delta_press > 0, ExcMessage("Read in of Material header failed (delta_press)."));
          AssertThrow(n_pressure > 0, ExcMessage("Read in of Material header failed (numpress)."));


          max_temp = min_temp + (n_temperature-1) * delta_temp;
          max_press = min_press + (n_pressure-1) * delta_press;

          density_values.reinit(n_temperature,n_pressure);
          thermal_expansivity_values.reinit(n_temperature,n_pressure);
          specific_heat_values.reinit(n_temperature,n_pressure);
          vp_values.reinit(n_temperature,n_pressure);
          vs_values.reinit(n_temperature,n_pressure);
          enthalpy_values.reinit(n_temperature,n_pressure);

          unsigned int i = 0;
          while (!in.eof())
            {
              double temp1,temp2;
              double rho,alpha,cp,vp,vs,h;
              in >> temp1 >> temp2;
              in >> rho;
              if (in.fail())
                {
                  in.clear();
                  rho = density_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }
              in >> alpha;
              if (in.fail())
                {
                  in.clear();
                  alpha = thermal_expansivity_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }
              in >> cp;
              if (in.fail())
                {
                  in.clear();
                  cp = specific_heat_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }
              in >> vp;
              if (in.fail())
                {
                  in.clear();
                  vp = vp_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }
              in >> vs;
              if (in.fail())
                {
                  in.clear();
                  vs = vs_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }
              in >> h;
              if (in.fail())
                {
                  in.clear();
                  h = enthalpy_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }

              std::getline(in, temp);
              if (in.eof())
                break;

              density_values[i%n_temperature][i/n_temperature]=rho;
              thermal_expansivity_values[i%n_temperature][i/n_temperature]=alpha;
              specific_heat_values[i%n_temperature][i/n_temperature]=cp;
              vp_values[i%n_temperature][i/n_temperature]=vp;
              vs_values[i%n_temperature][i/n_temperature]=vs;
              enthalpy_values[i%n_temperature][i/n_temperature]=h;

              i++;
            }
          AssertThrow(i == n_temperature*n_pressure, ExcMessage("Material table size not consistent with header."));

        }
      }



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



      PhaseFunctionInputs::PhaseFunctionInputs(const double temperature_,
                                               const double pressure_,
                                               const unsigned int phase_index_,
                                               const unsigned int composition_index_)

        :
        temperature(temperature_),
        pressure(pressure_),
        phase_index(phase_index_),
        composition_index(composition_index_)
      {}


      template <int dim>
      double
      PhaseFunction<dim>::phase_function (const PhaseFunctionInputs &in) const
      {
        double transition_pressure;
        double transition_pressure_width;

        if (use_depth_instead_of_pressure)
          {
            const std::pair<double, double> phase_transition_pressure_range =
              transition_depth_to_pressure(Point<dim>(), in.phase_index);

            transition_pressure = phase_transition_pressure_range.first;
            transition_pressure_width = phase_transition_pressure_range.second;
          }
        else
          {
            transition_pressure = transition_pressures[in.phase_index];
            transition_pressure_width = transition_pressure_widths[in.phase_index];
          }

        // Calculate the deviation from the phase transition point
        // (pressure here), with the current pressure and temperature.
        double pressure_deviation = in.pressure - transition_pressure -
                                    transition_slopes[in.phase_index] * (in.temperature - transition_temperatures[in.phase_index]);

        double phase_function;

        // Calculate the percentage of material that has undergone the phase transition as a function
        // of the pressure deviation and phase transition width (defined in terms of a pressure range here).
        if (transition_pressure_width == 0)
          (pressure_deviation > 0) ? phase_function = 1 : phase_function = 0;
        else
          phase_function = 0.5*(1.0 + std::tanh(pressure_deviation / transition_pressure_width));

        return phase_function;
      }



      template <int dim>
      double
      PhaseFunction<dim>::phase_function_derivative (const PhaseFunctionInputs &in) const
      {
        double transition_pressure;
        double transition_pressure_width;

        if (use_depth_instead_of_pressure)
          {
            const std::pair<double, double> phase_transition_pressure_range =
              transition_depth_to_pressure(Point<dim>(), in.phase_index);

            transition_pressure = phase_transition_pressure_range.first;
            transition_pressure_width = phase_transition_pressure_range.second;
          }
        else
          {
            transition_pressure = transition_pressures[in.phase_index];
            transition_pressure_width = transition_pressure_widths[in.phase_index];
          }

        // Calculate the pressure deviation from the transition pressure (both in temperature
        // and in pressure)
        double pressure_deviation = in.pressure - transition_pressure -
                                    transition_slopes[in.phase_index] * (in.temperature - transition_temperatures[in.phase_index]);

        double phase_function_derivative;

        // Calculate the analytical derivative of the phase function
        if (transition_pressure_width == 0)
          {
            phase_function_derivative = 0.;
          }
        else
          {
            phase_function_derivative = 0.5 / transition_pressure_width * (1.0 - std::tanh(pressure_deviation / transition_pressure_width)
                                                                           * std::tanh(pressure_deviation / transition_pressure_width));
          }

        return phase_function_derivative;
      }



      template <int dim>
      unsigned int
      PhaseFunction<dim>::
      n_phase_transitions () const
      {
        if (use_depth_instead_of_pressure)
          return transition_depths.size();
        else
          return transition_pressures.size();
      }



      template <int dim>
      double
      PhaseFunction<dim>::
      get_transition_slope (const unsigned int phase_index) const
      {
        return transition_slopes[phase_index];
      }



      template <int dim>
      std::pair<double, double>
      PhaseFunction<dim>::
      transition_depth_to_pressure (const Point<dim> &position,
                                    const int phase) const
      {
        double transition_pressure = 0.;
        double transition_pressure_width = 0.;

        if (this->get_adiabatic_conditions().is_initialized() && this->include_latent_heat())
          {
            const Point<dim,double> representative_transition_depth = this->get_geometry_model().representative_point(transition_depths[phase]);

            const Point<dim,double> plus_width = this->get_geometry_model().representative_point(transition_depths[phase] + transition_widths[phase]);
            const Point<dim,double> minus_width = this->get_geometry_model().representative_point(transition_depths[phase] - transition_widths[phase]);

            transition_pressure = (this->get_adiabatic_conditions().pressure(representative_transition_depth));

            transition_pressure_width = 0.5 * ((this->get_adiabatic_conditions().pressure(plus_width))
                                               - (this->get_adiabatic_conditions().pressure(minus_width)));
          }
//        else
//          {
//            const double gravity = this->get_gravity_model().gravity_vector(position).norm();
//
//            transition_pressure = reference_rho * gravity * transition_depths[phase];
//
//            transition_pressure_width = 0.5 * ((reference_rho * gravity * (transition_depths[phase] + transition_widths[phase])) -
//                                               (reference_rho * gravity * (transition_depths[phase] - transition_widths[phase])));
//          }

        return std::make_pair (transition_pressure, transition_pressure_width);
      }



      template <int dim>
      void
      PhaseFunction<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Phase transition depths", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of depths where phase transitions occur. Values must "
                           "monotonically increase. "
                           "Units: $m$.");
        prm.declare_entry ("Phase transition widths", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of widths for each phase transition, in terms of depth. The phase functions "
                           "are scaled with these values, leading to a jump between phases "
                           "for a value of zero and a gradual transition for larger values. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: $m$.");
        prm.declare_entry ("Phase transition pressures", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of pressures where phase transitions occur. Values must "
                           "monotonically increase. Define transition by depth instead of "
                           "pressure must be set to false to use this parameter. "
                           "Units: $Pa$.");
        prm.declare_entry ("Phase transition pressure widths", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of widths for each phase transition, in terms of pressure. The phase functions "
                           "are scaled with these values, leading to a jump between phases "
                           "for a value of zero and a gradual transition for larger values. "
                           "List must have the same number of entries as Phase transition pressures. "
                           "Define transition by depth instead of pressure must be set to false "
                           "to use this parameter. "
                           "Units: $Pa$.");
        prm.declare_entry ("Define transition by depth instead of pressure", "true",
                           Patterns::Bool (),
                           "Whether to list phase transitions by depth or pressure. If this parameter is true, "
                           "then the input file will use Phase transitions depths and Phase transition widths "
                           "to define the phase transition. If it is false, the parameter file will read in "
                           "phase transition data from Phase transition pressures and "
                           "Phase transition pressure widths.");
        prm.declare_entry ("Phase transition temperatures", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of temperatures where phase transitions occur. Higher or lower "
                           "temperatures lead to phase transition occurring in smaller or greater "
                           "depths than given in Phase transition depths, depending on the "
                           "Clapeyron slope given in Phase transition Clapeyron slopes. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: $\\si{K}$.");
        prm.declare_entry ("Phase transition Clapeyron slopes", "",
                           Patterns::List (Patterns::Double()),
                           "A list of Clapeyron slopes for each phase transition. A positive "
                           "Clapeyron slope indicates that the phase transition will occur in "
                           "a greater depth, if the temperature is higher than the one given in "
                           "Phase transition temperatures and in a smaller depth, if the "
                           "temperature is smaller than the one given in Phase transition temperatures. "
                           "For negative slopes the other way round. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: $Pa/K$.");
      }



      template <int dim>
      void
      PhaseFunction<dim>::parse_parameters (ParameterHandler &prm)
      {
        transition_depths = Utilities::string_to_double
                            (Utilities::split_string_list(prm.get ("Phase transition depths")));
        transition_widths= Utilities::string_to_double
                           (Utilities::split_string_list(prm.get ("Phase transition widths")));
        transition_pressures = Utilities::string_to_double
                               (Utilities::split_string_list(prm.get ("Phase transition pressures")));
        transition_pressure_widths= Utilities::string_to_double
                                    (Utilities::split_string_list(prm.get ("Phase transition pressure widths")));
        use_depth_instead_of_pressure = prm.get_bool ("Define transition by depth instead of pressure");
        transition_temperatures = Utilities::string_to_double
                                  (Utilities::split_string_list(prm.get ("Phase transition temperatures")));
        transition_slopes = Utilities::string_to_double
                            (Utilities::split_string_list(prm.get ("Phase transition Clapeyron slopes")));

        // make sure to check against the depth lists for size errors, since using depth
        if (use_depth_instead_of_pressure)
          {
            if (transition_widths.size() != transition_depths.size() ||
                transition_temperatures.size() != transition_depths.size() ||
                transition_slopes.size() != transition_depths.size())
              AssertThrow(false, ExcMessage("Error: At least one list that gives input parameters for the phase "
                                            "transitions has the wrong size. Currently checking against transition depths. "
                                            "If phase transitions in terms of pressure inputs are desired, check to make sure "
                                            "'Define transition by depth instead of pressure = false'."));
          }
        // make sure to check against the pressure lists for size errors,
        // since pressure is being used instead of depth.
        else
          {
            if (transition_pressure_widths.size() != transition_pressures.size() ||
                transition_temperatures.size() != transition_pressures.size() ||
                transition_slopes.size() != transition_pressures.size())
              AssertThrow(false, ExcMessage("Error: At least one list that gives input parameters for the phase "
                                            "transitions has the wrong size. Currently checking against transition pressures. "
                                            "If phase transitions in terms of depth inputs are desired, check to make sure "
                                            "'Define transition by depth instead of pressure = true'."));
          }
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
  template class PhaseFunction<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)
    }
  }
}
