/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
#include <aspect/structured_data.h>

#include <aspect/geometry_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>

#include <aspect/material_model/interface.h>
#include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>

#include <deal.II/base/exceptions.h>
#include <boost/lexical_cast.hpp>
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

        unsigned int
        MaterialLookup::dominant_phase (const double temperature,
                                        const double pressure) const
        {
          if (!has_dominant_phase_column)
            AssertThrow(false, ExcMessage("You can not ask for the column with the dominant phase if it does not exist in the data file."));
          return value(temperature, pressure, dominant_phase_indices);
        }

        bool
        MaterialLookup::has_dominant_phase() const
        {
          return has_dominant_phase_column;
        }

        std::vector<std::string>
        MaterialLookup::phase_volume_column_names() const
        {
          return phase_column_names;
        }

        double
        MaterialLookup::phase_volume_fraction(const int phase_id,
                                              const double temperature,
                                              const double pressure) const
        {
          return value(temperature,pressure,phase_volume_fractions[phase_id],interpolation);
        }


        const std::vector<std::string> &
        MaterialLookup::get_dominant_phase_names() const
        {
          return dominant_phase_names;
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

        unsigned int
        MaterialLookup::value (const double temperature,
                               const double pressure,
                               const Table<2, unsigned int> &values) const
        {
          const double nT = get_nT(temperature);
          const unsigned int inT = static_cast<unsigned int>(nT);

          const double np = get_np(pressure);
          const unsigned int inp = static_cast<unsigned int>(np);

          Assert(inT<values.n_rows(), ExcMessage("Attempting to look up a temperature value with index greater than the number of rows."));
          Assert(inp<values.n_cols(), ExcMessage("Attempting to look up a pressure value with index greater than the number of columns."));

          return values[inT][inp];
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
                                     const MPI_Comm comm)
        {
          /* Initializing variables */
          interpolation = interpol;
          delta_press=numbers::signaling_nan<double>();
          min_press=std::numeric_limits<double>::max();
          max_press=std::numeric_limits<double>::lowest();
          delta_temp=numbers::signaling_nan<double>();
          min_temp=std::numeric_limits<double>::max();
          max_temp=std::numeric_limits<double>::lowest();
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
                ++i;
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

                ++i;
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

                  ++i;
                }
            }
        }

        PerplexReader::PerplexReader(const std::string &filename,
                                     const bool interpol,
                                     const MPI_Comm comm)
        {
          /* Initializing variables */
          interpolation = interpol;
          delta_press=numbers::signaling_nan<double>();
          min_press=std::numeric_limits<double>::max();
          max_press=std::numeric_limits<double>::lowest();
          delta_temp=numbers::signaling_nan<double>();
          min_temp=std::numeric_limits<double>::max();
          max_temp=std::numeric_limits<double>::lowest();
          n_temperature=0;
          n_pressure=0;
          has_dominant_phase_column = false;

          std::string temp;
          // Read data from disk and distribute among processes
          std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

          // The following lines read in a PerpleX tab file in standard format
          // The first 13 lines are a header in the format:
          // |<perplex version>
          // <table filename>
          // <grid dim>
          // <grid variable 1> (usually T(K) or P(bar))
          // <min grid variable 1>
          // <delta grid variable 1>
          // <n steps grid variable 1>
          // <grid variable 2> (usually T(K) or P(bar))
          // <min grid variable 2>
          // <delta grid variable 2>
          // <n steps grid variable 2>
          // Number of property columns in the table
          // Column names

          // First line is the Perplex version number
          std::getline(in, temp); // get next line, table file name

          std::getline(in, temp); // get next line, dimension of table
          unsigned int n_variables;
          in >> n_variables;
          AssertThrow (n_variables==2, ExcMessage("The PerpleX file " + filename + " must be two dimensional (P(bar)-T(K))."));

          std::getline(in, temp); // get next line, either T(K) or P(bar)

          for (unsigned int i=0; i<2; ++i)
            {
              std::string natural_variable;
              in >> natural_variable;

              if (natural_variable == "T(K)")
                {
                  std::getline(in, temp);
                  in >> min_temp;
                  std::getline(in, temp);
                  in >> delta_temp;
                  std::getline(in, temp);
                  in >> n_temperature;
                  std::getline(in, temp); // get next line, either T(K), P(bar) or number of columns
                }
              else if (natural_variable == "P(bar)")
                {
                  std::getline(in, temp);
                  in >> min_press;
                  min_press *= 1e5;  // conversion from [bar] to [Pa]
                  std::getline(in, temp);
                  in >> delta_press;
                  delta_press *= 1e5; // conversion from [bar] to [Pa]
                  std::getline(in, temp);
                  in >> n_pressure;
                  std::getline(in, temp); // get next line, either T(K), P(bar) or number of columns
                }
              else
                {
                  AssertThrow (false, ExcMessage("The start of the PerpleX file " + filename + " does not have the expected format."));
                }
            }

          in >> n_columns;
          std::getline(in, temp); // get next line, column labels

          // here we string match to assign properties to columns
          // column i in text file -> column j in properties
          // Properties are stored in the order rho, alpha, cp, vp, vs, h
          std::vector<int> prp_indices(6, -1);
          std::vector<int> phase_column_indices;
          unsigned int dominant_phase_column_index = numbers::invalid_unsigned_int;

          // First two columns should be P(bar) and T(K).
          // Here we find the order.
          std::string column_name;
          in >> column_name;

          std::string first_natural_variable;
          if (column_name == "P(bar)")
            {
              first_natural_variable = column_name;
              in >> column_name;
              AssertThrow(column_name == "T(K)", ExcMessage("The second column name in PerpleX lookup file " + filename + " should be T(K)."));
            }
          else if (column_name == "T(K)")
            {
              first_natural_variable = column_name;
              in >> column_name;
              AssertThrow(column_name == "P(bar)", ExcMessage("The second column name in PerpleX lookup file " + filename + " should be P(bar)."));
            }
          else
            {
              AssertThrow(false, ExcMessage("The first column name in the PerpleX lookup file " + filename + " should be P(bar) or T(K)."));
            }

          for (unsigned int n=2; n<n_columns; ++n)
            {
              in >> column_name;
              if (column_name == "rho,kg/m3")
                prp_indices[0] = n;
              else if (column_name == "alpha,1/K")
                prp_indices[1] = n;
              else if (column_name == "cp,J/K/kg")
                prp_indices[2] = n;
              else if (column_name == "vp,km/s")
                prp_indices[3] = n;
              else if (column_name == "vs,km/s")
                prp_indices[4] = n;
              else if (column_name == "h,J/kg")
                prp_indices[5] = n;
              else if (column_name == "phase")
                {
                  has_dominant_phase_column = true;
                  dominant_phase_column_index = n;
                }
              else if (column_name.length() > 3)
                {
                  if (column_name.substr(0,13).compare("vol_fraction_") == 0)
                    {
                      if (std::find(phase_column_names.begin(),
                                    phase_column_names.end(),
                                    column_name) != phase_column_names.end())
                        {
                          AssertThrow(false,
                                      ExcMessage("The PerpleX lookup file " + filename + " must have unique column names. "
                                                 "Sometimes, the same phase is stable with >1 composition at the same "
                                                 "pressure and temperature, so you may see several columns with the same name. "
                                                 "Either combine columns with the same name, or change the names."));
                        }
                      // Populate phase_column_names with the column name
                      // and phase_column_indices with the column index in the current lookup file.
                      phase_column_indices.push_back(n);
                      phase_column_names.push_back(column_name);
                    }
                }
            }
          AssertThrow(std::all_of(prp_indices.begin(), prp_indices.end(), [](int i)
          {
            return i>=0;
          }),
          ExcMessage("The PerpleX lookup file " + filename + " must contain columns with names "
                     "rho,kg/m3, alpha,1/K, cp,J/K/kg, vp,km/s, vs,km/s and h,J/kg."));

          std::getline(in, temp); // first data line

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

          if (has_dominant_phase_column)
            dominant_phase_indices.reinit(n_temperature,n_pressure);

          phase_volume_fractions.resize(phase_column_names.size());
          for (auto &phase_volume_fraction : phase_volume_fractions)
            phase_volume_fraction.reinit(n_temperature,n_pressure);

          unsigned int i = 0;
          std::vector<double> previous_row_values(n_columns, 0.);

          while (!in.eof())
            {
              std::vector<double> row_values(n_columns);
              std::string phase;

              for (unsigned int n=0; n<n_columns; ++n)
                {
                  if (n == dominant_phase_column_index)
                    in >> phase;
                  else
                    in >> row_values[n]; // assigned as 0 if in.fail() == True

                  // P-T grids created with PerpleX-werami sometimes contain rows
                  // filled with NaNs at extreme P-T conditions where the thermodynamic
                  // models break down. These P-T regions are typically not relevant to
                  // geodynamic modeling (they most commonly appear above
                  // mantle liquidus temperatures at low pressures).
                  // More frustratingly, PerpleX-vertex occasionally fails to find a
                  // valid mineral assemblage in small, isolated regions within the domain,
                  // and so PerpleX-werami also returns NaNs for pixels within these regions.
                  // It is recommended that the user preprocesses their input
                  // files to replace these NaNs before plugging them into ASPECT.
                  // If this lookup encounters invalid doubles it replaces them
                  // with the most recent valid double.
                  if (in.fail())
                    {
                      in.clear();
                      row_values[n] = previous_row_values[n];
                    }
                }
              previous_row_values = row_values;

              std::getline(in, temp); // read next line
              if (in.eof())
                break;

              if (std::find(dominant_phase_names.begin(), dominant_phase_names.end(), phase) == dominant_phase_names.end())
                dominant_phase_names.push_back(phase);

              // The ordering of the first two columns in the PerpleX table files
              // dictates whether the inner loop is over temperature or pressure.
              // The first column is always the inner loop.
              // The following lines populate the material property tables
              // according to that implicit loop structure.
              if (first_natural_variable == "T(K)")
                {
                  density_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[0]];
                  thermal_expansivity_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[1]];
                  specific_heat_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[2]];
                  vp_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[3]];
                  vs_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[4]];
                  enthalpy_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[5]];

                  if (has_dominant_phase_column)
                    {
                      std::vector<std::string>::iterator it = std::find(dominant_phase_names.begin(), dominant_phase_names.end(), phase);
                      dominant_phase_indices[i%n_temperature][i/n_temperature] = std::distance(dominant_phase_names.begin(), it);
                    }

                  for (unsigned int n=0; n<phase_volume_fractions.size(); ++n)
                    {
                      phase_volume_fractions[n][i%n_temperature][i/n_temperature]=row_values[phase_column_indices[n]];
                    }
                }
              else // first_natural_variable == "P(bar)"
                {
                  density_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[0]];
                  thermal_expansivity_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[1]];
                  specific_heat_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[2]];
                  vp_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[3]];
                  vs_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[4]];
                  enthalpy_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[5]];

                  if (has_dominant_phase_column)
                    {
                      std::vector<std::string>::iterator it = std::find(dominant_phase_names.begin(), dominant_phase_names.end(), phase);
                      dominant_phase_indices[i/n_pressure][i%n_pressure] = std::distance(dominant_phase_names.begin(), it);
                    }

                  for (unsigned int n=0; n<phase_volume_fractions.size(); ++n)
                    {
                      phase_volume_fractions[n][i/n_pressure][i%n_pressure]=row_values[phase_column_indices[n]];
                    }
                }
              ++i;
            }
          AssertThrow(i == n_temperature*n_pressure, ExcMessage("Material table size not consistent with header."));

        }



        void
        EntropyReader::initialize(const MPI_Comm comm,
                                  const std::string &data_directory,
                                  const std::string &material_file_name)
        {
          material_lookup = std::make_unique<Utilities::StructuredDataLookup<2>>(7,1.0);
          material_lookup->load_file(data_directory+material_file_name,
                                     comm);
        }



        double
        EntropyReader::specific_heat(const double entropy,
                                     const double pressure) const
        {
          const double specific_heat = material_lookup->get_data({entropy,pressure}, 3);
          return specific_heat;
        }



        double
        EntropyReader::density(const double entropy,
                               const double pressure) const
        {
          const double density = material_lookup->get_data({entropy,pressure}, 1);
          return density;
        }



        double
        EntropyReader::thermal_expansivity(const double entropy,
                                           const double pressure) const
        {
          const double thermal_expansivity = material_lookup->get_data({entropy,pressure}, 2);
          return thermal_expansivity;
        }



        double
        EntropyReader::temperature(const double entropy,
                                   const double pressure) const
        {
          const double temperature = material_lookup->get_data({entropy,pressure}, 0, true);
          return temperature;
        }



        double
        EntropyReader::seismic_vp(const double entropy,
                                  const double pressure) const
        {
          const double seismic_vp = material_lookup->get_data({entropy,pressure}, 4);
          return seismic_vp;
        }



        double
        EntropyReader::seismic_vs(const double entropy,
                                  const double pressure) const
        {
          const double seismic_vs = material_lookup->get_data({entropy,pressure}, 5);
          return seismic_vs;
        }



        Tensor<1, 2>
        EntropyReader::density_gradient(const double entropy,
                                        const double pressure) const
        {
          const Tensor<1, 2> density_gradient= material_lookup->get_gradients({entropy,pressure}, 1);
          return density_gradient;
        }


      }



      std::vector<double>
      compute_only_composition_fractions(const std::vector<double> &compositional_fields,
                                         const std::vector<unsigned int> &indices_to_use)
      {
        std::vector<double> composition_fractions(indices_to_use.size()+1);

        // Clip the compositional fields so they are between zero and one,
        // and sum the compositional fields for normalization purposes.
        double sum_composition = 0.0;
        std::vector<double> x_comp (indices_to_use.size());

        for (unsigned int i=0; i < x_comp.size(); ++i)
          {
            x_comp[i] = std::min(std::max(compositional_fields[indices_to_use[i]], 0.0), 1.0);
            sum_composition += x_comp[i];
          }

        // Compute background field fraction
        if (sum_composition >= 1.0)
          composition_fractions[0] = 0.0;
        else
          composition_fractions[0] = 1.0 - sum_composition;

        // Compute and possibly normalize field fractions
        for (unsigned int i=0; i < x_comp.size(); ++i)
          {
            if (sum_composition >= 1.0)
              composition_fractions[i+1] = x_comp[i]/sum_composition;
            else
              composition_fractions[i+1] = x_comp[i];
          }

        return composition_fractions;
      }



      std::vector<double>
      compute_composition_fractions(const std::vector<double> &compositional_fields,
                                    const ComponentMask &field_mask)
      {
        std::vector<double> composition_fractions(compositional_fields.size()+1);

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

        // Compute background field fraction
        if (sum_composition >= 1.0)
          composition_fractions[0] = 0.0;
        else
          composition_fractions[0] = 1.0 - sum_composition;

        // Compute and possibly normalize field fractions
        for (unsigned int i=0; i < x_comp.size(); ++i)
          if (field_mask[i] == true)
            {
              if (sum_composition >= 1.0)
                composition_fractions[i+1] = x_comp[i]/sum_composition;
              else
                composition_fractions[i+1] = x_comp[i];
            }

        return composition_fractions;
      }



      std::vector<double>
      compute_volumes_from_masses(const std::vector<double> &masses,
                                  const std::vector<double> &densities,
                                  const bool return_as_fraction)
      {
        Assert(masses.size() == densities.size(),
               ExcMessage ("The mass fractions and densities vectors used for computing "
                           "volumes from masses have to have the same length! "
                           "You have provided "
                           + Utilities::int_to_string(masses.size()) +
                           " mass fractions and "
                           + Utilities::int_to_string(densities.size()) +
                           " densities."));

        const unsigned int n_fields = masses.size();
        std::vector<double> volumes(n_fields);

        if (n_fields == 1 && return_as_fraction)
          {
            volumes[0] = 1.0;
            return volumes;
          }

        double sum_volumes = 0.0;
        for (unsigned int j=0; j < n_fields; ++j)
          {
            volumes[j] = masses[j] / densities[j];
            sum_volumes += volumes[j];
          }

        if (return_as_fraction)
          {
            for (unsigned int j=0; j < n_fields; ++j)
              volumes[j] /= sum_volumes;
          }
        return volumes;
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
                           "have to have the same length! You have provided "
                           + Utilities::int_to_string(volume_fractions.size()) +
                           " volume fractions and "
                           + Utilities::int_to_string(parameter_values.size()) +
                           " parameter values."));

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



      template <int dim>
      void
      fill_averaged_equation_of_state_outputs(const EquationOfStateOutputs<dim> &eos_outputs,
                                              const std::vector<double> &mass_fractions,
                                              const std::vector<double> &volume_fractions,
                                              const unsigned int i,
                                              MaterialModelOutputs<dim> &out)
      {
        // The density, isothermal compressibility and thermal expansivity are volume-averaged
        // The specific entropy derivatives and heat capacity are mass-averaged
        out.densities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
        out.compressibilities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);
        out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);
        out.entropy_derivative_pressure[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.entropy_derivative_pressure, MaterialUtilities::arithmetic);
        out.entropy_derivative_temperature[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.entropy_derivative_temperature, MaterialUtilities::arithmetic);
        out.specific_heat[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);
      }



      double phase_average_value (const std::vector<double> &phase_function_values,
                                  const std::vector<unsigned int> &n_phase_transitions_per_composition,
                                  const std::vector<double> &parameter_values,
                                  const unsigned int composition_index,
                                  const PhaseUtilities::PhaseAveragingOperation operation)
      {
        // Calculate base index and assign base value
        unsigned int start_phase_index = 0;
        for (unsigned int i=0; i<composition_index; ++i)
          start_phase_index += n_phase_transitions_per_composition[i] + 1;

        double averaged_parameter = parameter_values[start_phase_index];
        if (n_phase_transitions_per_composition[composition_index] > 0)
          {
            // Do averaging when there are multiple phases
            if (operation == PhaseUtilities::logarithmic)
              averaged_parameter = std::log(averaged_parameter);

            for (unsigned int i=0; i<n_phase_transitions_per_composition[composition_index]; ++i)
              {
                const unsigned int phase_index = start_phase_index + i;

                Assert(phase_index+1<parameter_values.size(), ExcInternalError());
                if (operation == PhaseUtilities::logarithmic)
                  {
                    // First average by log values and then take the exponential.
                    // This is used for averaging prefactors in flow laws.
                    averaged_parameter += phase_function_values[phase_index-composition_index] * (std::log(parameter_values[phase_index+1]) - averaged_parameter);
                  }
                else if (operation == PhaseUtilities::arithmetic)
                  averaged_parameter += phase_function_values[phase_index-composition_index] * (parameter_values[phase_index+1] - averaged_parameter);

                else
                  AssertThrow(false, ExcInternalError());
              }
            if (operation == PhaseUtilities::logarithmic)
              averaged_parameter = std::exp(averaged_parameter);
          }
        return averaged_parameter;
      }



      template <int dim>
      PhaseFunctionInputs<dim>::PhaseFunctionInputs(const double temperature_,
                                                    const double pressure_,
                                                    const double depth_,
                                                    const double pressure_depth_derivative_,
                                                    const unsigned int phase_transition_index_)

        :
        temperature(temperature_),
        pressure(pressure_),
        depth(depth_),
        pressure_depth_derivative(pressure_depth_derivative_),
        phase_transition_index(phase_transition_index_)
      {}



      template <int dim>
      void
      PhaseFunctionDiscrete<dim>::initialize()
      {
        AssertThrow (this->introspection().get_number_of_fields_of_type(CompositionalFieldDescription::chemical_composition)+1 == material_file_names.size(),
                     ExcMessage(" The 'discrete phase function' requires that the number of compositional fields of type chemical composition plus one (for a background field) "
                                "matches the number of lookup tables."));


        minimum_temperature = std::vector<double>(material_file_names.size());
        maximum_temperature = std::vector<double>(material_file_names.size());
        interval_temperature = std::vector<double>(material_file_names.size());
        minimum_pressure = std::vector<double>(material_file_names.size());
        maximum_pressure = std::vector<double>(material_file_names.size());
        interval_pressure = std::vector<double>(material_file_names.size());
        for (unsigned int i = 0; i < material_file_names.size(); ++i)
          {
            material_lookup
            .push_back(std::make_unique<Utilities::StructuredDataLookup<2>>(8,1.0));

            material_lookup[i]->load_file(data_directory+material_file_names[i],
                                          this->get_mpi_communicator());

            Assert(material_lookup[i]->has_equidistant_coordinates(), ExcMessage("The loaded lookup tables do not use equidistant coordinates. The class 'PhaseFunctionDiscrete' cannot currently handle non equidistant coordinates."));

            minimum_temperature[i] = material_lookup[i]->get_interpolation_point_coordinates(0).front();
            maximum_temperature[i] = material_lookup[i]->get_interpolation_point_coordinates(0).back();
            interval_temperature[i] = (maximum_temperature[i]-minimum_temperature[i])/(material_lookup[i]->get_number_of_coordinates(0)-1);
            minimum_pressure[i] = material_lookup[i]->get_interpolation_point_coordinates(1).front();
            maximum_pressure[i] = material_lookup[i]->get_interpolation_point_coordinates(1).back();
            interval_pressure[i] = (maximum_pressure[i]-minimum_pressure[i])/(material_lookup[i]->get_number_of_coordinates(1)-1);
          }
      }


      template <int dim>
      double
      PhaseFunctionDiscrete<dim>::compute_value (const PhaseFunctionInputs<dim> &in) const
      {
        // the percentage of material that has undergone the transition
        double function_value;

        // lookup the most abundant phase
        unsigned int start_phase_transition_index = 0;
        unsigned int n_comp = 0;
        for (unsigned int n_relevant_fields  = 0 ; n_relevant_fields  < this->introspection().n_chemical_composition_fields() + 1 ; n_relevant_fields ++)
          {
            if (in.phase_transition_index < start_phase_transition_index + n_phase_transitions_per_chemical_composition[n_relevant_fields])
              {
                n_comp = n_relevant_fields ;
                break;
              }
            start_phase_transition_index += n_phase_transitions_per_chemical_composition[n_relevant_fields];
          }

        const double pressure_in_bar = in.pressure/1.e5;

        Assert (in.temperature >= minimum_temperature[n_comp] && in.temperature < maximum_temperature[n_comp], ExcInternalError());
        Assert (pressure_in_bar >= minimum_pressure[n_comp] && pressure_in_bar < maximum_pressure[n_comp], ExcInternalError());

        const std::vector<double> &temperature_points = material_lookup[n_comp]->get_interpolation_point_coordinates(0);
        const std::vector<double> &pressure_points = material_lookup[n_comp]->get_interpolation_point_coordinates(1);

        // round T and p points to exact values in the table
        const unsigned int i_T = static_cast<unsigned int>(std::round((in.temperature - minimum_temperature[n_comp]) / interval_temperature[n_comp]));
        const unsigned int i_p = static_cast<unsigned int>(std::round((pressure_in_bar - minimum_pressure[n_comp]) / interval_pressure[n_comp]));

        Point<2> temperature_pressure(temperature_points[i_T], pressure_points[i_p]);

        // determine the dominant phase index
        const unsigned int dominant_phase_index = static_cast<unsigned int>(std::round(material_lookup[n_comp]->get_data(temperature_pressure, 7)));

        // match the dominant phase index to one of the transition indicators
        unsigned int matched_phase_transition_index = numbers::invalid_unsigned_int;
        for (unsigned int i = start_phase_transition_index;
             i < start_phase_transition_index + n_phase_transitions_per_chemical_composition[n_comp];
             i++)
          {
            if (transition_indicators[i] == dominant_phase_index)
              {
                matched_phase_transition_index = i;
              }
          }

        // determine the value of phase function to facilitate the exact transition
        if ((matched_phase_transition_index != numbers::invalid_unsigned_int) && in.phase_transition_index <= matched_phase_transition_index)
          function_value = 1.0;
        else
          function_value = 0.0;

        return function_value;
      }


      template <int dim>
      double
      PhaseFunctionDiscrete<dim>::compute_derivative () const
      {
        // raises an error to ensure that a phase derivative request is not made for this phase function.
        AssertThrow(false, ExcNotImplemented());
      }


      template <int dim>
      unsigned int
      PhaseFunctionDiscrete<dim>::
      n_phase_transitions () const
      {
        return transition_indicators.size();
      }



      template <int dim>
      unsigned int
      PhaseFunctionDiscrete<dim>::
      n_phases () const
      {
        return n_phases_total;
      }



      template <int dim>
      unsigned int
      PhaseFunctionDiscrete<dim>::
      n_phases_over_all_chemical_compositions () const
      {
        return n_phases_total_chemical_compositions;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunctionDiscrete<dim>::n_phase_transitions_for_each_composition () const
      {
        return *n_phase_transitions_per_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunctionDiscrete<dim>::n_phases_for_each_composition () const
      {
        return n_phases_per_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunctionDiscrete<dim>::n_phase_transitions_for_each_chemical_composition () const
      {
        return n_phase_transitions_per_chemical_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunctionDiscrete<dim>::n_phases_for_each_chemical_composition () const
      {
        return n_phases_per_chemical_composition;
      }



      template <int dim>
      void
      PhaseFunctionDiscrete<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Phase transition indicators", "",
                           Patterns::Anything(),
                           "A list of phase indicators in a look-up table for each phase transition. "
                           "This parameter selectively assign different rheologies to specific phases, "
                           "rather than having a unique rheology for each phase in the table. "
                           "For example, if the table has phases 0, 1, and 2, and one only want "
                           "a distinct rheology for phase 2, then only phase 2 is needed "
                           "in the list of indicator. And phases 0, 1 will just be assigned "
                           "the rheology of the base phase. ");

        prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/entropy-table/pyrtable",
                           Patterns::DirectoryName (),
                           "The path to the model data. The path may also include the special "
                           "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                           "in which the ASPECT source files were located when ASPECT was "
                           "compiled. This interpretation allows, for example, to reference "
                           "files located in the `data/' subdirectory of ASPECT. ");

        prm.declare_entry ("Material file names", "material_table_temperature_pressure_small.txt",
                           Patterns::List (Patterns::Anything()),
                           "The file names of the material data (material "
                           "data is assumed to be in order with the ordering "
                           "of the compositional fields). Note that there are "
                           "two options on how many files need to be listed "
                           "here: 1. If only one file is provided, it is used "
                           "for the whole model domain, and compositional fields "
                           "are ignored. 2. If there is one more file name than the "
                           "number of compositional fields, then the first file is "
                           "assumed to define a `background composition' that is "
                           "modified by the compositional fields. These data files need "
                           "to have the same structure as the one necessary for equation "
                           "of state plus a new column for the phase indexes, which amounts "
                           "to 8 columns in total.");
      }



      template <int dim>
      void
      PhaseFunctionDiscrete<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        data_directory               = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
        material_file_names          = Utilities::split_string_list(prm.get ("Material file names"));

        // Make options file for parsing maps to double arrays
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();
        chemical_field_names.insert(chemical_field_names.begin(),"background");

        Utilities::MapParsing::Options options(chemical_field_names, "Phase transition indicators");
        options.allow_missing_keys = true;
        options.allow_multiple_values_per_key = true;
        options.store_values_per_key = true;
        options.n_values_per_key = std::vector<unsigned int>();

        const std::vector<double> transition_indicators_double          = Utilities::MapParsing::parse_map_to_double_array (prm.get("Phase transition indicators"),
                                                                          options);

        transition_indicators.reserve(transition_indicators_double.size());
        for (const double transition_indicator: transition_indicators_double)
          transition_indicators.push_back(static_cast<unsigned int>(std::round(transition_indicator)));
        n_phase_transitions_per_composition = std::make_unique<std::vector<unsigned int>>(options.n_values_per_key);

        n_phases_total = 0;
        n_phases_per_composition.clear();
        for (unsigned int n : *n_phase_transitions_per_composition)
          {
            n_phases_per_composition.push_back(n+1);
            n_phases_total += n+1;
          }

        // The background field is always the first composition
        n_phases_per_chemical_composition = {n_phases_per_composition[0]};
        n_phase_transitions_per_chemical_composition = {n_phases_per_composition[0] - 1};
        n_phases_total_chemical_compositions = n_phases_per_composition[0];
        for (auto i : this->introspection().chemical_composition_field_indices())
          {
            n_phases_per_chemical_composition.push_back(n_phases_per_composition[i+1]);
            n_phase_transitions_per_chemical_composition.push_back(n_phases_per_composition[i+1] - 1);
            n_phases_total_chemical_compositions += n_phases_per_composition[i+1];
          }
      }


      template <int dim>
      double
      PhaseFunction<dim>::compute_value (const PhaseFunctionInputs<dim> &in) const
      {
        AssertIndexRange (in.phase_transition_index, transition_temperature_lower_limits.size());
        AssertIndexRange (in.phase_transition_index, transition_temperature_upper_limits.size());

        // the percentage of material that has undergone the transition
        double function_value;
        if (in.temperature < transition_temperature_lower_limits[in.phase_transition_index] ||
            in.temperature >= transition_temperature_upper_limits[in.phase_transition_index])
          {
            // assign 0.0 if temperature is out of range
            function_value = 0.0;
          }
        else
          {
            if (use_depth_instead_of_pressure)
              {
                AssertIndexRange (in.phase_transition_index, transition_depths.size());

                // calculate the deviation from the transition point (convert temperature to depth)
                double depth_deviation = in.depth - transition_depths[in.phase_transition_index];

                if (in.pressure_depth_derivative != 0.0)
                  {
                    AssertIndexRange (in.phase_transition_index, transition_slopes.size());
                    AssertIndexRange (in.phase_transition_index, transition_temperatures.size());

                    depth_deviation -= transition_slopes[in.phase_transition_index] / in.pressure_depth_derivative
                                       * (in.temperature - transition_temperatures[in.phase_transition_index]);
                  }

                // use delta function for width = 0
                AssertIndexRange (in.phase_transition_index, transition_widths.size());
                if (transition_widths[in.phase_transition_index] == 0)
                  function_value = (depth_deviation > 0) ? 1. : 0.;
                else
                  function_value = 0.5*(1.0 + std::tanh(depth_deviation / transition_widths[in.phase_transition_index]));
              }
            else
              {
                // calculate the deviation from the transition point (convert temperature to pressure)
                AssertIndexRange (in.phase_transition_index, transition_pressures.size());
                const double pressure_deviation = in.pressure - transition_pressures[in.phase_transition_index]
                                                  - transition_slopes[in.phase_transition_index] * (in.temperature - transition_temperatures[in.phase_transition_index]);

                // use delta function for width = 0
                AssertIndexRange (in.phase_transition_index, transition_pressure_widths.size());
                if (transition_pressure_widths[in.phase_transition_index] == 0)
                  function_value = (pressure_deviation > 0) ? 1. : 0.;
                else
                  function_value = 0.5*(1.0 + std::tanh(pressure_deviation / transition_pressure_widths[in.phase_transition_index]));
              }
          }

        return function_value;
      }



      template <int dim>
      double
      PhaseFunction<dim>::compute_derivative (const PhaseFunctionInputs<dim> &in) const
      {
        double transition_pressure;
        double pressure_width;
        double width_temp;

        // we already should have the adiabatic conditions here
        Assert (this->get_adiabatic_conditions().is_initialized(),
                ExcMessage("The adiabatic conditions need to be already initialized "
                           "to calculate the derivative of phase functions. Either call this "
                           "function after the reference conditions have been computed, or implement "
                           "a workaround for the case without reference profile."));

        // phase transition based on depth
        if (use_depth_instead_of_pressure)
          {
            const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[in.phase_transition_index]);
            const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[in.phase_transition_index] + transition_widths[in.phase_transition_index]);
            const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[in.phase_transition_index] - transition_widths[in.phase_transition_index]);
            transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
            pressure_width = 0.5 * (this->get_adiabatic_conditions().pressure(transition_plus_width)
                                    - this->get_adiabatic_conditions().pressure(transition_minus_width));
            width_temp = transition_widths[in.phase_transition_index];
          }
        // using pressure instead of depth to define the phase transition
        else
          {
            transition_pressure = transition_pressures[in.phase_transition_index];
            pressure_width = transition_pressure_widths[in.phase_transition_index];
            width_temp = transition_pressure_widths[in.phase_transition_index];
          }

        // calculate the deviation from the transition point
        const double pressure_deviation = in.pressure - transition_pressure
                                          - transition_slopes[in.phase_transition_index] * (in.temperature - transition_temperatures[in.phase_transition_index]);

        // calculate the analytical derivative of the phase function
        if (
          (in.temperature < transition_temperature_lower_limits[in.phase_transition_index]) ||
          (in.temperature >= transition_temperature_upper_limits[in.phase_transition_index])
        )
          {
            // return 0 if temperature is out of range
            return 0;
          }
        else if (width_temp == 0)
          return 0;
        else
          return 0.5 / pressure_width * (1.0 - std::tanh(pressure_deviation / pressure_width)
                                         * std::tanh(pressure_deviation / pressure_width));
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
      unsigned int
      PhaseFunction<dim>::
      n_phases () const
      {
        return n_phases_total;
      }



      template <int dim>
      unsigned int
      PhaseFunction<dim>::
      n_phases_over_all_chemical_compositions () const
      {
        return n_phases_total_chemical_compositions;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunction<dim>::n_phase_transitions_for_each_composition () const
      {
        return *n_phase_transitions_per_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunction<dim>::n_phases_for_each_composition () const
      {
        return n_phases_per_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunction<dim>::n_phase_transitions_for_each_chemical_composition () const
      {
        return n_phase_transitions_per_chemical_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunction<dim>::n_phases_for_each_chemical_composition () const
      {
        return n_phases_per_chemical_composition;
      }



      template <int dim>
      double
      PhaseFunction<dim>::
      get_transition_slope (const unsigned int phase_transition_index) const
      {
        return transition_slopes[phase_transition_index];
      }



      template <int dim>
      double
      PhaseFunction<dim>::
      get_transition_depth (const unsigned int phase_transition_index) const
      {
        return transition_depths[phase_transition_index];
      }



      template <int dim>
      void
      PhaseFunction<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Phase transition depths", "",
                           Patterns::Anything(),
                           "A list of depths where phase transitions occur. Values must "
                           "monotonically increase. "
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Phase transition widths", "",
                           Patterns::Anything(),
                           "A list of widths for each phase transition, in terms of depth. The phase functions "
                           "are scaled with these values, leading to a jump between phases "
                           "for a value of zero and a gradual transition for larger values. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Phase transition pressures", "",
                           Patterns::Anything(),
                           "A list of pressures where phase transitions occur. Values must "
                           "monotonically increase. Define transition by depth instead of "
                           "pressure must be set to false to use this parameter. "
                           "Units: \\si{\\pascal}.");
        prm.declare_entry ("Phase transition pressure widths", "",
                           Patterns::Anything(),
                           "A list of widths for each phase transition, in terms of pressure. The phase functions "
                           "are scaled with these values, leading to a jump between phases "
                           "for a value of zero and a gradual transition for larger values. "
                           "List must have the same number of entries as Phase transition pressures. "
                           "Define transition by depth instead of pressure must be set to false "
                           "to use this parameter. "
                           "Units: \\si{\\pascal}.");
        prm.declare_entry ("Define transition by depth instead of pressure", "true",
                           Patterns::Bool (),
                           "Whether to list phase transitions by depth or pressure. If this parameter is true, "
                           "then the input file will use Phase transitions depths and Phase transition widths "
                           "to define the phase transition. If it is false, the parameter file will read in "
                           "phase transition data from Phase transition pressures and "
                           "Phase transition pressure widths.");
        prm.declare_entry ("Phase transition temperatures", "",
                           Patterns::Anything(),
                           "A list of temperatures where phase transitions occur. Higher or lower "
                           "temperatures lead to phase transition occurring in smaller or greater "
                           "depths than given in Phase transition depths, depending on the "
                           "Clapeyron slope given in Phase transition Clapeyron slopes. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: \\si{\\kelvin}.");
        prm.declare_entry ("Phase transition temperature upper limits",
                           boost::lexical_cast<std::string>(std::numeric_limits<double>::max()),
                           Patterns::Anything(),
                           "A list of upper temperature limits for each phase transition. Above "
                           "this temperature the respective phase transition is deactivated. The default "
                           "value means there is no upper limit for any phase transitions. "
                           "List must have the same number of entries as Phase transition depths. "
                           "When the optional temperature limits are applied, the user has to be "
                           "careful about the consistency between adjacent phases. Phase transitions "
                           "should be continuous in pressure-temperature space. "
                           "We recommend producing a phase diagram with "
                           "simple model setups to check the implementation as a starting point."
                           "Units: \\si{\\kelvin}.");
        prm.declare_entry ("Phase transition temperature lower limits",
                           boost::lexical_cast<std::string>(std::numeric_limits<double>::lowest()),
                           Patterns::Anything(),
                           "A list of lower temperature limits for each phase transition. Below "
                           "this temperature the respective phase transition is deactivated. The default "
                           "value means there is no lower limit for any phase transition. "
                           "List must have the same number of entries as Phase transition depths. "
                           "When the optional temperature limits are applied, the user has to be "
                           "careful about the consistency between adjacent phases. Phase transitions "
                           "should be continuous in pressure-temperature space. "
                           "We recommend producing a phase diagram with "
                           "simple model setups to check the implementation as a starting point."
                           "Units: \\si{\\kelvin}.");
        prm.declare_entry ("Phase transition Clapeyron slopes", "",
                           Patterns::Anything(),
                           "A list of Clapeyron slopes for each phase transition. A positive "
                           "Clapeyron slope indicates that the phase transition will occur in "
                           "a greater depth, if the temperature is higher than the one given in "
                           "Phase transition temperatures and in a smaller depth, if the "
                           "temperature is smaller than the one given in Phase transition temperatures. "
                           "For negative slopes the other way round. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: \\si{\\pascal\\per\\kelvin}.");
      }



      template <int dim>
      void
      PhaseFunction<dim>::parse_parameters (ParameterHandler &prm)
      {
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(),"background");
        Utilities::MapParsing::Options options(compositional_field_names, "Phase transition temperatures");
        options.allow_missing_keys = true;
        options.allow_multiple_values_per_key = true;
        options.store_values_per_key = true;
        options.n_values_per_key = std::vector<unsigned int>();

        options.property_name = "Phase transition temperatures";
        transition_temperatures = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        // Now we know how many entries we have for each composition
        // ensure that the number of phase transitions is consistent
        options.store_values_per_key = false;
        options.check_values_per_key = true;

        use_depth_instead_of_pressure = prm.get_bool ("Define transition by depth instead of pressure");

        if (use_depth_instead_of_pressure)
          {
            options.property_name = "Phase transition depths";
            transition_depths          = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            options.property_name = "Phase transition widths";
            transition_widths          = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);
          }
        else
          {
            options.property_name = "Phase transition pressures";
            transition_pressures = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            options.property_name = "Phase transition pressure widths";
            transition_pressure_widths = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);
          }

        options.property_name = "Phase transition temperature upper limits";
        transition_temperature_upper_limits = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        options.property_name = "Phase transition temperature lower limits";
        transition_temperature_lower_limits = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        options.property_name = "Phase transition Clapeyron slopes";
        transition_slopes = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        n_phase_transitions_per_composition = std::make_unique<std::vector<unsigned int>>(options.n_values_per_key);

        n_phases_total = 0;
        n_phases_per_composition.clear();
        for (unsigned int n : *n_phase_transitions_per_composition)
          {
            n_phases_per_composition.push_back(n+1);
            n_phases_total += n+1;
          }

        // The background field is always the first composition
        n_phases_per_chemical_composition = {n_phases_per_composition[0]};
        n_phase_transitions_per_chemical_composition = {n_phases_per_composition[0] - 1};
        n_phases_total_chemical_compositions = n_phases_per_composition[0];
        for (auto i : this->introspection().chemical_composition_field_indices())
          {
            n_phases_per_chemical_composition.push_back(n_phases_per_composition[i+1]);
            n_phase_transitions_per_chemical_composition.push_back(n_phases_per_composition[i+1] - 1);
            n_phases_total_chemical_compositions += n_phases_per_composition[i+1];
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
  template void fill_averaged_equation_of_state_outputs<dim> (const EquationOfStateOutputs<dim> &, \
                                                              const std::vector<double> &, \
                                                              const std::vector<double> &, \
                                                              const unsigned int, \
                                                              MaterialModelOutputs<dim> &); \
  template struct PhaseFunctionInputs<dim>; \
  template class PhaseFunction<dim>; \
  template class PhaseFunctionDiscrete<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
