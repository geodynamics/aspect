/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/material_model/damage_rheology.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <iostream>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    namespace
    {
      std::vector<std::string> make_dislocation_viscosity_outputs_names()
      {
        std::vector<std::string> names;
        names.push_back("dislocation_viscosity");
        return names;
      }
    }

    template <int dim>
    DislocationViscosityOutputs<dim>::DislocationViscosityOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_dislocation_viscosity_outputs_names()),
      dislocation_viscosities(n_points, numbers::signaling_nan<double>())
    {}


    template <int dim>
    std::vector<double>
    DislocationViscosityOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 1);

      return dislocation_viscosities;
    }

    namespace Lookup
    {
      double
      MaterialLookup::specific_heat(double temperature,
                                    double pressure) const
      {
        return value(temperature,pressure,specific_heat_values,interpolation);
      }

      double
      MaterialLookup::density(double temperature,
                              double pressure) const
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
                             const dealii::Table<2,
                             double> &values,
                             bool interpol) const
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

      std_cxx1x::array<double,2>
      MaterialLookup::get_pT_steps() const
      {
        return std_cxx1x::array<double,2> {delta_press,delta_temp};
      }

      double
      MaterialLookup::get_nT(double temperature) const
      {
        temperature=std::max(min_temp, temperature);
        temperature=std::min(temperature, max_temp-delta_temp);
        Assert(temperature>=min_temp, ExcMessage("ASPECT found a temperature less than min_T."));
        Assert(temperature<=max_temp, ExcMessage("ASPECT found a temperature greater than max_T."));
        return (temperature-min_temp)/delta_temp;
      }

      double
      MaterialLookup::get_np(double pressure) const
      {
        pressure=std::max(min_press, pressure);
        pressure=std::min(pressure, max_press-delta_press);
        Assert(pressure>=min_press, ExcMessage("ASPECT found a pressure less than min_p."));
        Assert(pressure<=max_press, ExcMessage("ASPECT found a pressure greater than max_p."));
        return (pressure-min_press)/delta_press;
      }

      HeFESToReader::HeFESToReader(const std::string &material_filename,
                                   const std::string &derivatives_filename,
                                   const bool interpol,
                                   const MPI_Comm &comm)
      {
        /* Initializing variables */
        interpolation = interpol;
        delta_press=-1.0;
        min_press=1e300;
        max_press=-1e300;
        delta_temp=-1.0;
        min_temp=1e300;
        max_temp=-1e300;
        numtemp=0;
        numpress=0;

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
                      numpress = i;
                      parsed_first_column = true;
                    }
                }

              getline(in, temp);
              if (in.eof())
                break;
              i++;
            }

          in.clear();
          in.seekg (0, in.beg);

          numtemp = i / numpress;

          Assert(i == numtemp * numpress,
                 ExcMessage("Material table size not consistent."));

          density_values.reinit(numtemp,numpress);
          thermal_expansivity_values.reinit(numtemp,numpress);
          specific_heat_values.reinit(numtemp,numpress);
          vp_values.reinit(numtemp,numpress);
          vs_values.reinit(numtemp,numpress);
          enthalpy_values.reinit(numtemp,numpress);

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
                                               rho = density_values[(i-1)%numtemp][(i-1)/numtemp];
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
                                               vs = vs_values[(i-1)%numtemp][(i-1)/numtemp];
                                             }
                                           in >> vp;
                                           if (in.fail())
                                             {
                                               in.clear();
                                               vp = vp_values[(i-1)%numtemp][(i-1)/numtemp];
                                             }
                                           in >> vsq >> vpq;

                                           in >> h;
                                           if (in.fail())
                                             {
                                               in.clear();
                                               h = enthalpy_values[(i-1)%numtemp][(i-1)/numtemp];
                                             }
                                           else
                                             h *= 1e6; // conversion from [kJ/g] to [J/kg]

                                           getline(in, temp);
                                           if (in.eof())
                                             break;

                                           density_values[i/numpress][i%numpress]=rho;
                                           thermal_expansivity_values[i/numpress][i%numpress]=alpha;
                                           specific_heat_values[i/numpress][i%numpress]=cp;
                                           vp_values[i/numpress][i%numpress]=vp;
                                           vs_values[i/numpress][i%numpress]=vs;
                                           enthalpy_values[i/numpress][i%numpress]=h;

                                           i++;
            }

          delta_temp = (max_temp - min_temp) / (numtemp - 1);
          delta_press = (max_press - min_press) / (numpress - 1);

          Assert(max_temp >= 0.0, ExcMessage("Read in of Material header failed (max_temp)."));
          Assert(delta_temp > 0, ExcMessage("Read in of Material header failed (delta_temp)."));
          Assert(numtemp > 0, ExcMessage("Read in of Material header failed (numtemp)."));
          Assert(max_press >= 0, ExcMessage("Read in of Material header failed (max_press)."));
          Assert(delta_press > 0, ExcMessage("Read in of Material header failed (delta_press)."));
          Assert(numpress > 0, ExcMessage("Read in of Material header failed (numpress)."));
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
                    cp = specific_heat_values[(i-1)%numtemp][(i-1)/numtemp];
                  }
                else
                  cp *= 1e3; // conversion from [J/g/K] to [J/kg/K]

                                                            in >> alpha >> alpha_eff;
                                                            if (in.fail() || (alpha_eff <= std::numeric_limits<double>::min()))
                                                              {
                                                                in.clear();
                                                                alpha_eff = thermal_expansivity_values[(i-1)%numtemp][(i-1)/numtemp];
                                                              }
                                                            else
                                                              {
                                                                alpha *= 1e-5;
                                                                alpha_eff *= 1e-5;
                                                              }

                                                            in >> temp1 >> temp2;
                                                            if (in.fail())
                                                              in.clear();


                                                            getline(in, temp);
                                                            if (in.eof())
                                                              break;

                                                            specific_heat_values[i/numpress][i%numpress]=cp;
                                                            thermal_expansivity_values[i/numpress][i%numpress]=alpha_eff;

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
        delta_press=-1.0;
        min_press=-1.0;
        delta_temp=-1.0;
        min_temp=-1.0;
        numtemp=0;
        numpress=0;

        std::string temp;
        // Read data from disk and distribute among processes
        std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

        getline(in, temp); // eat first line
        getline(in, temp); // eat next line
        getline(in, temp); // eat next line
        getline(in, temp); // eat next line

        in >> min_temp;
        getline(in, temp);
        in >> delta_temp;
        getline(in, temp);
        in >> numtemp;
        getline(in, temp);
        getline(in, temp);
        in >> min_press;
        min_press *= 1e5;  // conversion from [bar] to [Pa]
        getline(in, temp);
        in >> delta_press;
        delta_press *= 1e5; // conversion from [bar] to [Pa]
        getline(in, temp);
        in >> numpress;
        getline(in, temp);
        getline(in, temp);
        getline(in, temp);

        Assert(min_temp >= 0.0, ExcMessage("Read in of Material header failed (mintemp)."));
        Assert(delta_temp > 0, ExcMessage("Read in of Material header failed (delta_temp)."));
        Assert(numtemp > 0, ExcMessage("Read in of Material header failed (numtemp)."));
        Assert(min_press >= 0, ExcMessage("Read in of Material header failed (min_press)."));
        Assert(delta_press > 0, ExcMessage("Read in of Material header failed (delta_press)."));
        Assert(numpress > 0, ExcMessage("Read in of Material header failed (numpress)."));


        max_temp = min_temp + (numtemp-1) * delta_temp;
        max_press = min_press + (numpress-1) * delta_press;

        density_values.reinit(numtemp,numpress);
        thermal_expansivity_values.reinit(numtemp,numpress);
        specific_heat_values.reinit(numtemp,numpress);
        vp_values.reinit(numtemp,numpress);
        vs_values.reinit(numtemp,numpress);
        enthalpy_values.reinit(numtemp,numpress);

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
                rho = density_values[(i-1)%numtemp][(i-1)/numtemp];
              }
            in >> alpha;
            if (in.fail())
              {
                in.clear();
                alpha = thermal_expansivity_values[(i-1)%numtemp][(i-1)/numtemp];
              }
            in >> cp;
            if (in.fail())
              {
                in.clear();
                cp = specific_heat_values[(i-1)%numtemp][(i-1)/numtemp];
              }
            in >> vp;
            if (in.fail())
              {
                in.clear();
                vp = vp_values[(i-1)%numtemp][(i-1)/numtemp];
              }
            in >> vs;
            if (in.fail())
              {
                in.clear();
                vs = vs_values[(i-1)%numtemp][(i-1)/numtemp];
              }
            in >> h;
            if (in.fail())
              {
                in.clear();
                h = enthalpy_values[(i-1)%numtemp][(i-1)/numtemp];
              }

            getline(in, temp);
            if (in.eof())
              break;

            density_values[i%numtemp][i/numtemp]=rho;
            thermal_expansivity_values[i%numtemp][i/numtemp]=alpha;
            specific_heat_values[i%numtemp][i/numtemp]=cp;
            vp_values[i%numtemp][i/numtemp]=vp;
            vs_values[i%numtemp][i/numtemp]=vs;
            enthalpy_values[i%numtemp][i/numtemp]=h;

            i++;
          }
        Assert(i==numtemp*numpress, ExcMessage("Material table size not consistent with header."));

      }
    }


    template <int dim>
    void
    DamageRheology<dim>::initialize()
    {
      n_material_data = material_file_names.size();
      for (unsigned i = 0; i < n_material_data; i++)
        {
          if (material_file_format == perplex)
            material_lookup.push_back(std_cxx1x::shared_ptr<Lookup::MaterialLookup>
                                      (new Lookup::PerplexReader(datadirectory+material_file_names[i],
                                                                 use_bilinear_interpolation,
                                                                 this->get_mpi_communicator())));
          else if (material_file_format == hefesto)
            material_lookup.push_back(std_cxx1x::shared_ptr<Lookup::MaterialLookup>
                                      (new Lookup::HeFESToReader(datadirectory+material_file_names[i],
                                                                 datadirectory+derivatives_file_names[i],
                                                                 use_bilinear_interpolation,
                                                                 this->get_mpi_communicator())));
          else
            AssertThrow (false, ExcNotImplemented());
        }
    }


    template <int dim>
    void
    DamageRheology<dim>::
    update()
    {
    }


    template <int dim>
    double
    DamageRheology<dim>::
    phase_function (const Point<dim> &position,
                    const double temperature,
                    const double pressure,
                    const unsigned int phase) const
    {
      Assert(phase < transition_depths.size(),
             ExcMessage("Error: Phase index is too large. This phase index does not exist!"));

      // if we already have the adiabatic conditions, we can use them
      if (this->get_adiabatic_conditions().is_initialized())
        {
          // first, get the pressure at which the phase transition occurs normally
          const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
          const double transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);

          // then calculate the deviation from the transition point (both in temperature
          // and in pressure)
          double pressure_deviation = pressure - transition_pressure
                                      - transition_slopes[phase] * (temperature - transition_temperatures[phase]);

          // last, calculate the percentage of material that has undergone the transition
          return (pressure_deviation > 0) ? 1 : 0;
        }

      // if we do not have the adiabatic conditions, we have to use the depth instead
      // this is less precise, because we do not have the exact pressure gradient, instead we use pressure/depth
      // (this is for calculating e.g. the density in the adiabatic profile)
      else
        {
          double depth = this->get_geometry_model().depth(position);
          double depth_deviation = (pressure > 0
                                    ?
                                    depth - transition_depths[phase]
                                    - transition_slopes[phase] * (depth / pressure) * (temperature - transition_temperatures[phase])
                                    :
                                    depth - transition_depths[phase]
                                    - transition_slopes[phase] / (this->get_gravity_model().gravity_vector(position).norm() * reference_rho)
                                    * (temperature - transition_temperatures[phase]));

          return (depth_deviation > 0) ? 1 : 0;
        }
    }



    template <int dim>
    unsigned int
    DamageRheology<dim>::
    get_phase_index (const Point<dim> &position,
                     const double temperature,
                     const double pressure) const
    {
      Assert(grain_growth_activation_energy.size()>0,
             ExcMessage("Error: No grain evolution parameters are given!"));

      unsigned int phase_index = 0;
      if (transition_depths.size()>0)
        if (phase_function(position, temperature, pressure, transition_depths.size()-1) == 1)
          phase_index = transition_depths.size();

      for (unsigned int j=1; j<transition_depths.size(); ++j)
        if (phase_function(position, temperature, pressure, j) != phase_function(position, temperature, pressure, j-1))
          phase_index = j;

      return phase_index;
    }

    template <int dim>
    void
    DamageRheology<dim>::
    convert_log_grain_size (const bool normal_to_log,
                            std::vector<double> &composition) const
    {
      // get grain size and limit it to a global minimum
      const std::string field_name = "olivine_grain_size";
      if (!this->introspection().compositional_name_exists(field_name))
        return;

      double grain_size = composition[this->introspection().compositional_index_for_name(field_name)];

      if (normal_to_log)
        grain_size = -std::log(std::max(grain_size,min_grain_size));
      else
        grain_size = std::max(std::exp(-grain_size),min_grain_size);

      composition[this->introspection().compositional_index_for_name(field_name)] = grain_size;

      return;
    }

    template <int dim>
    double
    DamageRheology<dim>::
    grain_size_growth_rate (const double                  temperature,
                            const double                  pressure,
                            const std::vector<double>    &compositional_fields,
                            const SymmetricTensor<2,dim> &strain_rate,
                            const Tensor<1,dim>          &/*velocity*/,
                            const Point<dim>             &position,
                            const unsigned int            field_index,
                            const int                     crossed_transition) const
    {
      // we want to iterate over the grain size evolution here, as we solve in fact an ordinary differential equation
      // and it is not correct to use the starting grain size (and introduces instabilities)
      const double original_grain_size = compositional_fields[field_index];
      if ((original_grain_size != original_grain_size) || this->get_timestep() == 0.0
          || original_grain_size < std::numeric_limits<double>::min())
        return 0.0;

      // set up the parameters for the sub-timestepping of grain size evolution
      std::vector<double> current_composition = compositional_fields;
      double grain_size = original_grain_size;
      double grain_size_change = 0.0;
      const double timestep = this->get_timestep();
      double grain_growth_timestep = 500 * 3600 * 24 * 365.25; // 500 yrs
      double time = 0;

      // find out in which phase we are
      const unsigned int ol_index = get_phase_index(position, temperature, pressure);

      // we keep the dislocation viscosity of the last iteration as guess
      // for the next one
      double current_dislocation_viscosity = 0.0;

      do
        {
          time += grain_growth_timestep;

          if (timestep - time < 0)
            {
              grain_growth_timestep = timestep - (time - grain_growth_timestep);
              time = timestep;
            }

          // grain size growth due to Ostwald ripening
          const double m = grain_growth_exponent[ol_index];
          const double grain_size_growth_rate = grain_growth_rate_constant[ol_index] / (m * pow(grain_size,m-1))
                                                * exp(- (grain_growth_activation_energy[ol_index] + pressure * grain_growth_activation_volume[ol_index])
                                                      / (gas_constant * temperature));
          const double grain_size_growth = grain_size_growth_rate * grain_growth_timestep;

          // grain size reduction in dislocation creep regime
          const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
          const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

          const double current_diffusion_viscosity   = diffusion_viscosity(temperature, pressure, current_composition, strain_rate, position);
          current_dislocation_viscosity              = dislocation_viscosity(temperature, pressure, current_composition, strain_rate, position, current_dislocation_viscosity);

          double current_viscosity;
          if (std::abs(second_strain_rate_invariant) > 1e-30)
            current_viscosity = current_dislocation_viscosity * current_diffusion_viscosity / (current_dislocation_viscosity + current_diffusion_viscosity);
          else
            current_viscosity = current_diffusion_viscosity;

          const double dislocation_strain_rate = second_strain_rate_invariant
                                                 * current_viscosity / current_dislocation_viscosity;

          double grain_size_reduction = 0.0;

          if (use_paleowattmeter)
            {
              // paleowattmeter: Austin and Evans (2007): Paleowattmeters: A scaling relation for dynamically recrystallized grain size. Geology 35, 343-346
              const double stress = 2.0 * second_strain_rate_invariant * current_viscosity;
              const double grain_size_reduction_rate = 2.0 * stress * boundary_area_change_work_fraction[ol_index] * dislocation_strain_rate * pow(grain_size,2)
                                                       / (geometric_constant[ol_index] * grain_boundary_energy[ol_index]);
              grain_size_reduction = grain_size_reduction_rate * grain_growth_timestep;
            }
          else
            {
              //TODO: check equations! are we missing a factor 2 from the conversion between strain_rate and second invariant of strain_rate?
              // paleopiezometer: Hall and Parmentier (2003): Influence of grain size evolution on convective instability. Geochem. Geophys. Geosyst., 4(3).
              grain_size_reduction = reciprocal_required_strain[ol_index] * dislocation_strain_rate * grain_size * grain_growth_timestep;
            }

          grain_size_change = grain_size_growth - grain_size_reduction;

          if ((grain_size_change / grain_size < 0.001 && grain_size_growth / grain_size < 0.1
               && grain_size_reduction / grain_size < 0.1) || grain_size == 0.0)
            grain_growth_timestep *= 2;
          else if (grain_size_change / grain_size > 0.1 || grain_size_growth / grain_size > 0.5
                   || grain_size_reduction / grain_size > 0.5)
            {
              grain_size_change = 0.0;
              time -= grain_growth_timestep;
              grain_growth_timestep /= 2.0;
            }

          grain_size += grain_size_change;
          current_composition[field_index] = grain_size;

          if (grain_size < 0)
            {
              std::cout << "Grain size smaller 0:  " << grain_size << " ," << grain_size_growth
                        << " ," << grain_size_reduction << ", timestep: " << grain_growth_timestep << "! \n ";
              break;
            }
        }
      while (time < timestep);

      // reduce grain size to recrystallized_grain_size when crossing phase transitions
      // if the distance in radial direction a grain moved compared to the last time step
      // is crossing a phase transition, reduce grain size

      // TODO: recrystallize first, and then do grain size growth/reduction for grains that crossed the transition
      // in dependence of the distance they have moved
      double phase_grain_size_reduction = 0.0;
      if (this->introspection().name_for_compositional_index(field_index) == "olivine_grain_size"
          &&
          this->get_timestep_number() > 0)
        {
          // check if material has crossed any phase transition, if yes, reset grain size
          if (crossed_transition != -1)
            if (recrystallized_grain_size[crossed_transition] > 0.0)
              phase_grain_size_reduction = grain_size - recrystallized_grain_size[crossed_transition];
        }
      else if (this->introspection().name_for_compositional_index(field_index) == "pyroxene_grain_size")
        {
          phase_grain_size_reduction = 0.0;
        }

      if (grain_size < 5.e-6)
        {
          std::cout << "Grain size is " << grain_size << "! It needs to be larger than 5e-6.\n";
          grain_size = 5e-6;
        }

      return grain_size - original_grain_size - phase_grain_size_reduction;
    }

    template <int dim>
    double
    DamageRheology<dim>::
    diffusion_viscosity (const double                  temperature,
                         const double                  pressure,
                         const std::vector<double>    &composition,
                         const SymmetricTensor<2,dim> &strain_rate,
                         const Point<dim>             &position) const
    {
      const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
      const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

      // TODO: make this more general, for more phases we have to average grain size somehow
      // TODO: default when field is not given & warning
      // limit the grain size to a global minimum
      const std::string field_name = "olivine_grain_size";
      const double grain_size = this->introspection().compositional_name_exists(field_name)
                                ?
                                composition[this->introspection().compositional_index_for_name(field_name)]
                                :
                                1.0;

      // Currently this will never be called without adiabatic_conditions initialized, but just in case
      const double adiabatic_pressure = this->get_adiabatic_conditions().is_initialized()
                                        ?
                                        this->get_adiabatic_conditions().pressure(position)
                                        :
                                        pressure;

      // find out in which phase we are
      const unsigned int ol_index = get_phase_index(position, temperature, adiabatic_pressure);

      // TODO: we use the prefactors from Behn et al., 2009 as default values, but their laws use the strain rate
      // and we use the second invariant --> check if the prefactors should be changed
      double energy_term = exp((diffusion_activation_energy[ol_index] + diffusion_activation_volume[ol_index] * adiabatic_pressure)
                               / (diffusion_creep_exponent[ol_index] * gas_constant * temperature));
      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_energy_term
            = exp((diffusion_activation_energy[ol_index] + diffusion_activation_volume[ol_index] * adiabatic_pressure)
                  / (diffusion_creep_exponent[ol_index] * gas_constant * this->get_adiabatic_conditions().temperature(position)));

          const double temperature_dependence = energy_term / adiabatic_energy_term;
          if (temperature_dependence > max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
          if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
        }

      const double strain_rate_dependence = (1.0 - diffusion_creep_exponent[ol_index]) / diffusion_creep_exponent[ol_index];

      return pow(diffusion_creep_prefactor[ol_index],-1.0/diffusion_creep_exponent[ol_index])
             * std::pow(second_strain_rate_invariant,strain_rate_dependence)
             * pow(grain_size, diffusion_creep_grain_size_exponent[ol_index]/diffusion_creep_exponent[ol_index])
             * energy_term;
    }

    template <int dim>
    double
    DamageRheology<dim>::
    dislocation_viscosity (const double      temperature,
                           const double      pressure,
                           const std::vector<double> &composition,
                           const SymmetricTensor<2,dim> &strain_rate,
                           const Point<dim> &position,
                           const double viscosity_guess) const
    {
      const double diff_viscosity = diffusion_viscosity(temperature,pressure,composition,strain_rate,position) ;

      // Start the iteration with the full strain rate
      double dis_viscosity;
      if (viscosity_guess == 0)
        dis_viscosity = dislocation_viscosity_fixed_strain_rate(temperature,pressure,std::vector<double>(),strain_rate,position);
      else
        dis_viscosity = viscosity_guess;

      double dis_viscosity_old = 0;
      unsigned int i = 0;
      while (std::abs((dis_viscosity-dis_viscosity_old) / dis_viscosity) > dislocation_viscosity_iteration_threshold && i < dislocation_viscosity_iteration_number)
        {
          SymmetricTensor<2,dim> dislocation_strain_rate = diff_viscosity
                                                           / (diff_viscosity + dis_viscosity) * strain_rate;
          dis_viscosity_old = dis_viscosity;
          dis_viscosity = dislocation_viscosity_fixed_strain_rate(temperature,
                                                                  pressure,
                                                                  std::vector<double>(),
                                                                  dislocation_strain_rate,
                                                                  position);
          i++;
        }
      return dis_viscosity;
    }

    template <int dim>
    double
    DamageRheology<dim>::
    dislocation_viscosity_fixed_strain_rate (const double      temperature,
                                             const double      pressure,
                                             const std::vector<double> &,
                                             const SymmetricTensor<2,dim> &dislocation_strain_rate,
                                             const Point<dim> &position) const
    {
      const SymmetricTensor<2,dim> shear_strain_rate = dislocation_strain_rate - 1./dim * trace(dislocation_strain_rate) * unit_symmetric_tensor<dim>();
      const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

      // Currently this will never be called without adiabatic_conditions initialized, but just in case
      const double adiabatic_pressure = this->get_adiabatic_conditions().is_initialized()
                                        ?
                                        this->get_adiabatic_conditions().pressure(position)
                                        :
                                        pressure;

      // find out in which phase we are
      const unsigned int ol_index = get_phase_index(position, temperature, adiabatic_pressure);

      double energy_term = exp((dislocation_activation_energy[ol_index] + dislocation_activation_volume[ol_index] * adiabatic_pressure)
                               / (dislocation_creep_exponent[ol_index] * gas_constant * temperature));
      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_energy_term
            = exp((dislocation_activation_energy[ol_index] + dislocation_activation_volume[ol_index] * adiabatic_pressure)
                  / (dislocation_creep_exponent[ol_index] * gas_constant * this->get_adiabatic_conditions().temperature(position)));

          const double temperature_dependence = energy_term / adiabatic_energy_term;
          if (temperature_dependence > max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
          if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
        }

      const double strain_rate_dependence = (1.0 - dislocation_creep_exponent[ol_index]) / dislocation_creep_exponent[ol_index];

      return pow(dislocation_creep_prefactor[ol_index],-1.0/dislocation_creep_exponent[ol_index])
             * std::pow(second_strain_rate_invariant,strain_rate_dependence)
             * energy_term;
    }

    template <int dim>
    double
    DamageRheology<dim>::
    viscosity_ratio (const double temperature,
                     const double pressure,
                     const std::vector<double> &composition_,
                     const SymmetricTensor<2,dim> &strain_rate,
                     const Point<dim> &position) const
    {
      // convert the grain size from log to normal
      std::vector<double> composition (composition_);
      if (advect_log_gransize)
        convert_log_grain_size(false,composition);

      return dislocation_viscosity(temperature,pressure,composition,strain_rate,position)
             / diffusion_viscosity(temperature,pressure,composition,strain_rate,position);
    }

    template <int dim>
    double
    DamageRheology<dim>::
    viscosity (const double temperature,
               const double pressure,
               const std::vector<double> &composition,
               const SymmetricTensor<2,dim> &strain_rate,
               const Point<dim> &position) const
    {
      //TODO: add assert
      /*if (this->get_timestep_number() > 0)
        Assert (grain_size >= 1.e-6, ExcMessage ("Error: The grain size should not be smaller than 1e-6 m."));*/

      const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
      const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

      const double diff_viscosity = diffusion_viscosity(temperature, pressure, composition, strain_rate, position);

      double effective_viscosity;
      if (std::abs(second_strain_rate_invariant) > 1e-30)
        {
          const double disl_viscosity = dislocation_viscosity(temperature, pressure, composition, strain_rate, position);
          effective_viscosity = disl_viscosity * diff_viscosity / (disl_viscosity + diff_viscosity);
        }
      else
        effective_viscosity = diff_viscosity;
      return effective_viscosity;
    }


    template <int dim>
    double
    DamageRheology<dim>::
    enthalpy (const double      temperature,
              const double      pressure,
              const std::vector<double> &compositional_fields,
              const Point<dim> &/*position*/) const
    {
      //this function is not called from evaluate() so we need to care about
      //corrections for temperature and pressure
      //const double corrected_temperature = get_corrected_temperature(temperature,pressure,position);
      //const double corrected_pressure = get_corrected_pressure(temperature,pressure,position);
      const double corrected_temperature = temperature;
      const double corrected_pressure = pressure;
      AssertThrow ((reference_compressibility != 0.0) || use_table_properties,
                   ExcMessage("Currently only compressible models are supported for seismic output."));

      double enthalpy = 0.0;
      if (n_material_data == 1)
        enthalpy += material_lookup[0]->enthalpy(corrected_temperature,corrected_pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            enthalpy += compositional_fields[i] * material_lookup[i]->enthalpy(corrected_temperature,corrected_pressure);
        }
      return enthalpy;
    }


    template <int dim>
    double
    DamageRheology<dim>::
    seismic_Vp (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &/*position*/) const
    {
      //this function is not called from evaluate() so we need to care about
      //corrections for temperature and pressure
      //const double corrected_temperature = get_corrected_temperature(temperature,pressure,position);
      //const double corrected_pressure = get_corrected_pressure(temperature,pressure,position);
      const double corrected_temperature = temperature;
      const double corrected_pressure = pressure;
      AssertThrow ((reference_compressibility != 0.0) || use_table_properties,
                   ExcMessage("Currently only compressible models are supported for seismic output."));

      double vp = 0.0;
      if (n_material_data == 1)
        vp += material_lookup[0]->seismic_Vp(corrected_temperature,corrected_pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            vp += compositional_fields[i] * material_lookup[i]->seismic_Vp(corrected_temperature,corrected_pressure);
        }
      return vp;
    }


    template <int dim>
    double
    DamageRheology<dim>::
    seismic_Vs (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &/*position*/) const
    {
      //this function is not called from evaluate() so we need to care about
      //corrections for temperature and pressure
      //const double corrected_temperature = get_corrected_temperature(temperature,pressure,position);
      //const double corrected_pressure = get_corrected_pressure(temperature,pressure,position);
      const double corrected_temperature = temperature;
      const double corrected_pressure = pressure;
      AssertThrow ((reference_compressibility != 0.0) || use_table_properties,
                   ExcMessage("Currently only compressible models are supported for seismic output."));

      double vs = 0.0;
      if (n_material_data == 1)
        vs += material_lookup[0]->seismic_Vs(corrected_temperature,corrected_pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            vs += compositional_fields[i] * material_lookup[i]->seismic_Vs(corrected_temperature,corrected_pressure);
        }
      return vs;
    }


    template <int dim>
    double
    DamageRheology<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    double
    DamageRheology<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &compositional_fields, /*composition*/
             const Point<dim> &) const
    {
      if (!use_table_properties)
        {
          const double composition_dependence = compositional_fields.size()>0
                                                ?
                                                compositional_delta_rho * compositional_fields[0]
                                                :
                                                0.0;

          return (reference_rho + composition_dependence) * std::exp(reference_compressibility * (pressure - this->get_surface_pressure()))
                 * (1 - thermal_alpha * (temperature - reference_T));
        }
      else
        {
          double rho = 0.0;
          if (n_material_data == 1)
            {
              rho = material_lookup[0]->density(temperature,pressure);
            }
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                rho += compositional_fields[i] * material_lookup[i]->density(temperature,pressure);
            }

          return rho;
        }
    }



    template <int dim>
    bool
    DamageRheology<dim>::
    is_compressible () const
    {
      return (reference_compressibility != 0)
             || use_table_properties;
    }

    template <int dim>
    double
    DamageRheology<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields,
                     const Point<dim> &position) const
    {
      if (!use_table_properties)
        return reference_compressibility;

      double dRhodp = 0.0;
      if (n_material_data == 1)
        dRhodp += material_lookup[0]->dRhodp(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            dRhodp += compositional_fields[i] * material_lookup[i]->dRhodp(temperature,pressure);
        }
      const double rho = density(temperature,pressure,compositional_fields,position);
      return (1/rho)*dRhodp;
    }

    template <int dim>
    double
    DamageRheology<dim>::
    thermal_expansion_coefficient (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &/*position*/) const
    {
      double alpha = 0.0;
      if (!use_table_properties)
        return thermal_alpha;
      else
        {
          if (n_material_data == 1)
            alpha = material_lookup[0]->thermal_expansivity(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                alpha += compositional_fields[i] * material_lookup[i]->thermal_expansivity(temperature,pressure);
            }
        }
      alpha = std::max(std::min(alpha,max_thermal_expansivity),min_thermal_expansivity);
      return alpha;
    }

    template <int dim>
    double
    DamageRheology<dim>::
    specific_heat (const double temperature,
                   const double pressure,
                   const std::vector<double> &compositional_fields,
                   const Point<dim> &/*position*/) const
    {
      double cp = 0.0;
      if (!use_table_properties)
        return reference_specific_heat;
      else
        {
          if (n_material_data == 1)
            cp = material_lookup[0]->specific_heat(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                cp += compositional_fields[i] * material_lookup[i]->specific_heat(temperature,pressure);
            }
        }
      cp = std::max(std::min(cp,max_specific_heat),min_specific_heat);
      return cp;
    }

    template <int dim>
    std_cxx1x::array<std::pair<double, unsigned int>,2>
    DamageRheology<dim>::
    enthalpy_derivative (const typename Interface<dim>::MaterialModelInputs &in) const
    {
      std_cxx1x::array<std::pair<double, unsigned int>,2> derivative;
      unsigned int T_points(0), p_points(0);
      double dHdT(0.0), dHdp(0.0);

      if (in.current_cell.state() == IteratorState::valid)
        {
          const QTrapez<dim> quadrature_formula;
          const unsigned int n_q_points = quadrature_formula.size();
          FEValues<dim> fe_values (this->get_mapping(),
                                   this->get_fe(),
                                   quadrature_formula,
                                   update_values);

          std::vector<double> temperatures(n_q_points), pressures(n_q_points);
          std::vector<std::vector<double> > compositions (quadrature_formula.size(),std::vector<double> (this->n_compositional_fields()));
          std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

          fe_values.reinit (in.current_cell);

          // get the various components of the solution, then
          // evaluate the material properties there
          fe_values[this->introspection().extractors.temperature]
          .get_function_values (this->get_current_linearization_point(), temperatures);
          fe_values[this->introspection().extractors.pressure]
          .get_function_values (this->get_current_linearization_point(), pressures);

          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            fe_values[this->introspection().extractors.compositional_fields[c]]
            .get_function_values(this->get_current_linearization_point(),
                                 composition_values[c]);
          for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
            {
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                compositions[q][c] = composition_values[c][q];
            }

          AssertThrow (material_lookup.size() == 1,
                       ExcMessage("This formalism is only implemented for one material "
                                  "table."));

          // We have to take into account here that the p,T spacing of the table of material properties
          // we use might be on a finer grid than our model. Because of that we average the enthalpy
          // derivatives over the whole cell, using finer temperature and pressure intervals
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              for (unsigned int p=0; p<n_q_points; ++p)
                {
                  if (std::fabs(temperatures[q] - temperatures[p]) > material_lookup[0]->get_pT_steps()[1] / 4.0)
                    {
                      for (unsigned int substep = 0; substep < max_latent_heat_substeps; ++substep)
                        {
                          const double current_pressure = pressures[q]
                                                          + ((double)(substep)/(double)(max_latent_heat_substeps))
                                                          * (pressures[p]-pressures[q]);
                          const double T1_substep = temperatures[q]
                                                    + ((double)(substep)/(double)(max_latent_heat_substeps))
                                                    * (temperatures[p]-temperatures[q]);
                          const double T2_substep = temperatures[q]
                                                    + ((double)(substep+1)/(double)(max_latent_heat_substeps))
                                                    * (temperatures[p]-temperatures[q]);
                          const double enthalpy2 = material_lookup[0]->enthalpy(T2_substep,current_pressure);
                          const double enthalpy1 = material_lookup[0]->enthalpy(T1_substep,current_pressure);
                          dHdT += (enthalpy2-enthalpy1)/(T2_substep-T1_substep);
                          T_points++;
                        }
                    }
                  if (std::fabs(pressures[q] - pressures[p]) > material_lookup[0]->get_pT_steps()[0] / 4.0)
                    {
                      for (unsigned int substep = 0; substep < max_latent_heat_substeps; ++substep)
                        {
                          const double current_temperature = temperatures[q]
                                                             + ((double)(substep)/(double)(max_latent_heat_substeps))
                                                             * (temperatures[p]-temperatures[q]);
                          const double p1_substep = pressures[q]
                                                    + ((double)(substep)/(double)(max_latent_heat_substeps))
                                                    * (pressures[p]-pressures[q]);
                          const double p2_substep = pressures[q]
                                                    + ((double)(substep+1)/(double)(max_latent_heat_substeps))
                                                    * (pressures[p]-pressures[q]);
                          const double enthalpy2 = material_lookup[0]->enthalpy(current_temperature,p2_substep);
                          const double enthalpy1 = material_lookup[0]->enthalpy(current_temperature,p1_substep);
                          dHdp += (enthalpy2-enthalpy1)/(p2_substep-p1_substep);
                          p_points++;
                        }
                    }
                }
            }

          if ((T_points > 0)
              && (p_points > 0))
            {
              dHdT /= T_points;
              dHdp /= p_points;
            }
        }

      derivative[0] = std::make_pair(dHdT,T_points);
      derivative[1] = std::make_pair(dHdp,p_points);
      return derivative;
    }

    template <int dim>
    void
    DamageRheology<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      for (unsigned int i=0; i<in.position.size(); ++i)
        {
          //Use the adiabatic pressure instead of the real one, because of oscillations
          const double pressure = (this->get_adiabatic_conditions().is_initialized())
                                  ?
                                  this->get_adiabatic_conditions().pressure(in.position[i])
                                  :
                                  in.pressure[i];

          // convert the grain size from log to normal
          std::vector<double> composition (in.composition[i]);
          if (advect_log_gransize)
            convert_log_grain_size(false,composition);
          else
            for (unsigned int c=0; c<composition.size(); ++c)
              composition[c] = std::max(min_grain_size,composition[c]);

          // set up an integer that tells us which phase transition has been crossed inside of the cell
          int crossed_transition(-1);

          if (this->get_adiabatic_conditions().is_initialized())
            for (unsigned int phase=0; phase<transition_depths.size(); ++phase)
              {
                // first, get the pressure at which the phase transition occurs normally
                const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
                const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[phase] + transition_widths[phase]);
                const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[phase] - transition_widths[phase]);
                const double transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
                const double pressure_width = 0.5 * (this->get_adiabatic_conditions().pressure(transition_plus_width)
                                                     - this->get_adiabatic_conditions().pressure(transition_minus_width));


                // then calculate the deviation from the transition point (both in temperature
                // and in pressure)
                double pressure_deviation = pressure - transition_pressure
                                            - transition_slopes[phase] * (in.temperature[i] - transition_temperatures[phase]);

                if ((std::abs(pressure_deviation) < pressure_width)
                    &&
                    ((in.velocity[i] * this->get_gravity_model().gravity_vector(in.position[i])) * pressure_deviation > 0))
                  crossed_transition = phase;
              }
          else
            for (unsigned int j=0; j<in.position.size(); ++j)
              for (unsigned int k=0; k<transition_depths.size(); ++k)
                if ((phase_function(in.position[i], in.temperature[i], pressure, k)
                     != phase_function(in.position[j], in.temperature[j], in.pressure[j], k))
                    &&
                    ((in.velocity[i] * this->get_gravity_model().gravity_vector(in.position[i]))
                     * ((in.position[i] - in.position[j]) * this->get_gravity_model().gravity_vector(in.position[i])) > 0))
                  crossed_transition = k;


          if (in.strain_rate.size() > 0)
            {
              double effective_viscosity;
              double disl_viscosity = std::numeric_limits<double>::max();

              const SymmetricTensor<2,dim> shear_strain_rate = in.strain_rate[i] - 1./dim * trace(in.strain_rate[i]) * unit_symmetric_tensor<dim>();
              const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

              const double diff_viscosity = diffusion_viscosity(in.temperature[i], pressure, composition, in.strain_rate[i], in.position[i]);

              if (std::abs(second_strain_rate_invariant) > 1e-30)
                {
                  disl_viscosity = dislocation_viscosity(in.temperature[i], pressure, composition, in.strain_rate[i], in.position[i]);
                  effective_viscosity = disl_viscosity * diff_viscosity / (disl_viscosity + diff_viscosity);
                }
              else
                effective_viscosity = diff_viscosity;

              out.viscosities[i] = std::min(std::max(min_eta,effective_viscosity),max_eta);

              if (DislocationViscosityOutputs<dim> *disl_viscosities_out = out.template get_additional_output<DislocationViscosityOutputs<dim> >())
                {
                  disl_viscosities_out->dislocation_viscosities[i] = std::min(std::max(min_eta,disl_viscosity),1e300);
                }
            }

          out.densities[i] = density(in.temperature[i], pressure, in.composition[i], in.position[i]);
          out.thermal_conductivities[i] = k_value;
          out.compressibilities[i] = compressibility(in.temperature[i], pressure, composition, in.position[i]);

          out.boundary_area_change_work_fraction[i] = boundary_area_change_work_fraction[get_phase_index(in.position[i],in.temperature[i],pressure)];

          // TODO: make this more general for not just olivine grains
          if (in.strain_rate.size() > 0)
            for (unsigned int c=0; c<composition.size(); ++c)
              {
                if (this->introspection().name_for_compositional_index(c) == "olivine_grain_size")
                  {
                    out.reaction_terms[i][c] = grain_size_growth_rate(in.temperature[i], pressure, composition,
                                                                      in.strain_rate[i], in.velocity[i], in.position[i], c, crossed_transition);
                    if (advect_log_gransize)
                      out.reaction_terms[i][c] = - out.reaction_terms[i][c] / composition[c];
                  }
                else
                  out.reaction_terms[i][c] = 0.0;
              }
        }

      /* We separate the calculation of specific heat and thermal expansivity,
       * because they depend on cell-wise averaged values that are only available
       * here
       */
      double average_temperature(0.0);
      double average_density(0.0);
      for (unsigned int i = 0; i < in.position.size(); ++i)
        {
          average_temperature += in.temperature[i];
          average_density += out.densities[i];
        }
      average_temperature /= in.position.size();
      average_density /= in.position.size();

      std_cxx1x::array<std::pair<double, unsigned int>,2> dH;

      if (use_table_properties && use_enthalpy)
        dH = enthalpy_derivative(in);

      for (unsigned int i = 0; i < in.position.size(); ++i)
        {
          //Use the adiabatic pressure instead of the real one, because of oscillations
          const double pressure = (this->get_adiabatic_conditions().is_initialized())
                                  ?
                                  this->get_adiabatic_conditions().pressure(in.position[i])
                                  :
                                  in.pressure[i];

          if (!use_table_properties)
            {
              out.thermal_expansion_coefficients[i] = thermal_alpha;
              out.specific_heat[i] = reference_specific_heat;
            }
          else if (use_enthalpy)
            {
              if (this->get_adiabatic_conditions().is_initialized()
                  && (in.current_cell.state() == IteratorState::valid)
                  && (dH[0].second > 0)
                  && (dH[1].second > 0))
                {
                  out.thermal_expansion_coefficients[i] = (1 - average_density * dH[1].first) / average_temperature;
                  out.specific_heat[i] = dH[0].first;
                }
              else
                {
                  if (material_lookup.size() == 1)
                    {
                      out.thermal_expansion_coefficients[i] = (1 - out.densities[i] * material_lookup[0]->dHdp(in.temperature[i],pressure)) / in.temperature[i];
                      out.specific_heat[i] = material_lookup[0]->dHdT(in.temperature[i],pressure);
                    }
                  else
                    {
                      ExcNotImplemented();
                    }
                }
            }
          else
            {
              out.thermal_expansion_coefficients[i] = thermal_expansion_coefficient(in.temperature[i], pressure, in.composition[i], in.position[i]);
              out.specific_heat[i] = specific_heat(in.temperature[i], pressure, in.composition[i], in.position[i]);
            }

          out.thermal_expansion_coefficients[i] = std::max(std::min(out.thermal_expansion_coefficients[i],max_thermal_expansivity),min_thermal_expansivity);
          out.specific_heat[i] = std::max(std::min(out.specific_heat[i],max_specific_heat),min_specific_heat);
        }
    }


    template <int dim>
    void
    DamageRheology<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Damage rheology model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the constant viscosity. Units: $kg/m/s$.");
          prm.declare_entry ("Composition viscosity prefactor 1", "1.0",
                             Patterns::Double (0),
                             "A linear dependency of viscosity on the first compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition.");
          prm.declare_entry ("Composition viscosity prefactor 2", "1.0",
                             Patterns::Double (0),
                             "A linear dependency of viscosity on the second compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition.");
          prm.declare_entry ("Compositional density difference", "100.0",
                             Patterns::Double (),
                             "Density excess of the first compositional field."
                             "Units: $kg/m^3$");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Reference compressibility", "4e-12",
                             Patterns::Double (0),
                             "The value of the reference compressibility. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Phase transition depths", "",
                             Patterns::List (Patterns::Double(0)),
                             "A list of depths where phase transitions occur. Values must "
                             "monotonically increase. "
                             "Units: $m$.");
          prm.declare_entry ("Phase transition temperatures", "",
                             Patterns::List (Patterns::Double(0)),
                             "A list of temperatures where phase transitions occur. Higher or lower "
                             "temperatures lead to phase transition ocurring in smaller or greater "
                             "depths than given in Phase transition depths, depending on the "
                             "Clapeyron slope given in Phase transition Clapeyron slopes. "
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: $K$.");
          prm.declare_entry ("Phase transition widths", "",
                             Patterns::List (Patterns::Double(0)),
                             "A list of widths for each phase transition. This is only use to specify "
                             "the region where the recrystallized grain size is assigned after material "
                             "has crossed a phase transition and should accordingly be chosen similar "
                             "to the maximum cell width expected at the phase transition."
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: $m$.");
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
          prm.declare_entry ("Corresponding phase for transition", "",
                             Patterns::List(Patterns::Anything()),
                             "A user-defined list of phases, which correspond to the name of the phase the "
                             "transition should occur in. "
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: $Pa/K$.");
          prm.declare_entry ("Grain growth activation energy", "3.5e5",
                             Patterns::List (Patterns::Double(0)),
                             "The activation energy for grain growth $E_g$. "
                             "Units: $J/mol$.");
          prm.declare_entry ("Grain growth activation volume", "8e-6",
                             Patterns::List (Patterns::Double(0)),
                             "The activation volume for grain growth $E_g$. "
                             "Units: $m^3/mol$.");
          prm.declare_entry ("Grain growth exponent", "3",
                             Patterns::List (Patterns::Double(0)),
                             "Exponent of the grain growth law $p_g$. This is an experimentally determined "
                             "grain growth constant. "
                             "Units: none.");
          prm.declare_entry ("Grain growth rate constant", "1.5e-5",
                             Patterns::List (Patterns::Double(0)),
                             "Prefactor of the Ostwald ripening grain growth law $G_0$. "
                             "This is dependent on water content, which is assumed to be "
                             "50 H/10^6 Si for the default value. "
                             "Units: $m^{p_g}/s$.");
          prm.declare_entry ("Reciprocal required strain", "10",
                             Patterns::List (Patterns::Double(0)),
                             "This parameters $\\lambda$ gives an estimate of the strain necessary "
                             "to achieve a new grain size. ");
          prm.declare_entry ("Recrystallized grain size", "",
                             Patterns::List (Patterns::Double(0)),
                             "The grain size $d_{ph}$ to that a phase will be reduced to when crossing a phase transition. "
                             "When set to zero, grain size will not be reduced. "
                             "Units: m.");
          prm.declare_entry ("Use paleowattmeter", "true",
                             Patterns::Bool (),
                             "A flag indicating whether the computation should be use the "
                             "paleowattmeter approach of Austin and Evans (2007) for grain size reduction "
                             "in the dislocation creep regime (if true) or the paleopiezometer aprroach "
                             "from Hall and Parmetier (2003) (if false).");
          prm.declare_entry ("Average specific grain boundary energy", "1.0",
                             Patterns::List (Patterns::Double(0)),
                             "The average specific grain boundary energy $\\gamma$. "
                             "Units: J/m^2.");
          prm.declare_entry ("Work fraction for boundary area change", "0.1",
                             Patterns::List (Patterns::Double(0)),
                             "The fraction $\\chi$ of work done by dislocation creep to change the grain boundary area. "
                             "Units: J/m^2.");
          prm.declare_entry ("Geometric constant", "3",
                             Patterns::List (Patterns::Double(0)),
                             "Geometric constant $c$ used in the paleowattmeter grain size reduction law. "
                             "Units: none.");
          prm.declare_entry ("Dislocation viscosity iteration threshold", "1e-3",
                             Patterns::Double(0),
                             "We need to perform an iteration inside the computation "
                             "of the dislocation viscosity, because it depends on the "
                             "dislocation strain rate, which depends on the dislocation "
                             "viscosity itself. This number determines the termination "
                             "accuracy, i.e. if the dislocation viscosity changes by less "
                             "than this factor we terminate the iteration.");
          prm.declare_entry ("Dislocation viscosity iteration number", "10",
                             Patterns::Integer(0),
                             "We need to perform an iteration inside the computation "
                             "of the dislocation viscosity, because it depends on the "
                             "dislocation strain rate, which depends on the dislocation "
                             "viscosity itself. This number determines the maximum "
                             "number of iterations that are performed. ");
          prm.declare_entry ("Dislocation creep exponent", "3.5",
                             Patterns::List (Patterns::Double(0)),
                             "Power-law exponent $n_{dis}$ for dislocation creep. "
                             "Units: none.");
          prm.declare_entry ("Dislocation activation energy", "4.8e5",
                             Patterns::List (Patterns::Double(0)),
                             "The activation energy for dislocation creep $E_{dis}$. "
                             "Units: $J/mol$.");
          prm.declare_entry ("Dislocation activation volume", "1.1e-5",
                             Patterns::List (Patterns::Double(0)),
                             "The activation volume for dislocation creep $V_{dis}$. "
                             "Units: $m^3/mol$.");
          prm.declare_entry ("Dislocation creep prefactor", "4.5e-15",
                             Patterns::List (Patterns::Double(0)),
                             "Prefactor for the dislocation creep law $A_{dis}$. "
                             "Units: $Pa^{-n_{dis}}/s$.");
          prm.declare_entry ("Diffusion creep exponent", "1",
                             Patterns::List (Patterns::Double(0)),
                             "Power-law exponent $n_{diff}$ for diffusion creep. "
                             "Units: none.");
          prm.declare_entry ("Diffusion activation energy", "3.35e5",
                             Patterns::List (Patterns::Double(0)),
                             "The activation energy for diffusion creep $E_{diff}$. "
                             "Units: $J/mol$.");
          prm.declare_entry ("Diffusion activation volume", "4e-6",
                             Patterns::List (Patterns::Double(0)),
                             "The activation volume for diffusion creep $V_{diff}$. "
                             "Units: $m^3/mol$.");
          prm.declare_entry ("Diffusion creep prefactor", "7.4e-15",
                             Patterns::List (Patterns::Double(0)),
                             "Prefactor for the diffusion creep law $A_{diff}$. "
                             "Units: $m^{p_{diff}} Pa^{-n_{diff}}/s$.");
          prm.declare_entry ("Diffusion creep grain size exponent", "3",
                             Patterns::List (Patterns::Double(0)),
                             "Diffusion creep grain size exponent $p_{diff}$ that determines the "
                             "dependence of vescosity on grain size. "
                             "Units: none.");
          prm.declare_entry ("Maximum temperature dependence of viscosity", "100",
                             Patterns::Double (0),
                             "The factor by which viscosity at adiabatic temperature and ambient temperature "
                             "are allowed to differ (a value of x means that the viscosity can be x times higher "
                             "or x times lower compared to the value at adiabatic temperature. This parameter "
                             "is introduced to limit local viscosity contrasts, but still allow for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: none.");
          prm.declare_entry ("Minimum viscosity", "1e18",
                             Patterns::Double (0),
                             "The minimum viscosity that is allowed in the whole model domain. This parameter "
                             "is introduced to limit global viscosity contrasts, but still allows for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: Pa s.");
          prm.declare_entry ("Maximum viscosity", "1e26",
                             Patterns::Double (0),
                             "The maximum viscosity that is allowed in the whole model domain. This parameter "
                             "is introduced to limit global viscosity contrasts, but still allows for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: Pa s.");
          prm.declare_entry ("Minimum specific heat", "500",
                             Patterns::Double (0),
                             "The minimum specific heat that is allowed in the whole model domain. This parameter "
                             "is introduced to limit global viscosity contrasts, but still allows for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: J/kg/K.");
          prm.declare_entry ("Maximum specific heat", "6000",
                             Patterns::Double (0),
                             "The maximum specific heat that is allowed in the whole model domain. This parameter "
                             "is introduced to limit global viscosity contrasts, but still allows for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: J/kg/K.");
          prm.declare_entry ("Minimum thermal expansivity", "1e-5",
                             Patterns::Double (),
                             "The minimum thermal expansivity that is allowed in the whole model domain. This parameter "
                             "is introduced to limit global viscosity contrasts, but still allows for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: 1/K.");
          prm.declare_entry ("Maximum thermal expansivity", "1e-3",
                             Patterns::Double (),
                             "The maximum thermal expansivity that is allowed in the whole model domain. "
                             "This parameter is introduced to limit global viscosity contrasts, but still "
                             "allows for a widely varying viscosity over the whole mantle range. "
                             "Units: 1/K.");
          prm.declare_entry ("Maximum latent heat substeps", "1",
                             Patterns::Integer (1),
                             "The maximum number of substeps over the temperature pressure range "
                             "to calculate the averaged enthalpy gradient over a cell.");
          prm.declare_entry ("Minimum grain size", "1e-5",
                             Patterns::Double (0),
                             "The minimum grain size that is used for the material model. This parameter "
                             "is introduced to limit local viscosity contrasts, but still allows for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: Pa s.");
          prm.declare_entry ("Lower mantle grain size scaling", "1.0",
                             Patterns::Double (0),
                             "A scaling factor for the grain size in the lower mantle. In models where the "
                             "high grain size contrast between the upper and lower mantle causes numerical "
                             "problems, the grain size in the lower mantle can be scaled to a larger value, "
                             "simultaneously scaling the viscosity prefactors and grain growth parameters "
                             "to keep the same physical behavior. Differences to the original formulation "
                             "only occur when material with a smaller grain size than the recrystallization "
                             "grain size cross the upper-lower mantle boundary. "
                             "The real grain size can be obtained by dividing the model grain size by this value. "
                             "Units: none.");
          prm.declare_entry ("Advect logarithm of grain size", "false",
                             Patterns::Bool (),
                             "Whether to advect the logarithm of the grain size or the "
                             "grain size. The equation and the physics are the same, "
                             "but for problems with high grain size gradients it might "
                             "be preferable to advect the logarithm. ");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/steinberger/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Material file names", "pyr-ringwood88.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file names of the material data. "
                             "List with as many components as active "
                             "compositional fields (material data is assumed to "
                             "be in order with the ordering of the fields). ");
          prm.declare_entry ("Derivatives file names", "",
                             Patterns::List (Patterns::Anything()),
                             "The file names of the enthalpy derivatives data. "
                             "List with as many components as active "
                             "compositional fields (material data is assumed to "
                             "be in order with the ordering of the fields). ");
          prm.declare_entry ("Use table properties", "false",
                             Patterns::Bool(),
                             "Whether to use the table properties also for "
                             "density, thermal expansivity and specific heat. "
                             "If false the properties are generated as in the "
                             "simple compressible plugin.");
          prm.declare_entry ("Material file format", "perplex",
                             Patterns::Selection ("perplex|hefesto"),
                             "The material file format to be read in the property "
                             " tables.");
          prm.declare_entry ("Use enthalpy for material properties", "true",
                             Patterns::Bool(),
                             "Whether to use the enthalpy to calculate thermal "
                             "expansivity and specific heat (if true) or use the "
                             "thermal expansivity and specific heat values from "
                             "the material properties table directly (if false).");
          prm.declare_entry ("Bilinear interpolation", "true",
                             Patterns::Bool (),
                             "Whether to use bilinear interpolation to compute "
                             "material properties (slower but more accurate). ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    DamageRheology<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Damage rheology model");
        {
          gas_constant               = 8.314462;
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Viscosity");
          composition_viscosity_prefactor_1 = prm.get_double ("Composition viscosity prefactor 1");
          composition_viscosity_prefactor_2 = prm.get_double ("Composition viscosity prefactor 2");
          compositional_delta_rho = prm.get_double("Compositional density difference");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          reference_compressibility  = prm.get_double ("Reference compressibility");


          transition_depths         = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Phase transition depths")));
          transition_temperatures   = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Phase transition temperatures")));
          transition_slopes         = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Phase transition Clapeyron slopes")));
          transition_phases         = Utilities::split_string_list (prm.get("Corresponding phase for transition"));
          recrystallized_grain_size = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Recrystallized grain size")));
          transition_widths         = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Phase transition widths")));

          if (transition_temperatures.size() != transition_depths.size() ||
              transition_slopes.size() != transition_depths.size() ||
              transition_phases.size() != transition_depths.size() ||
              transition_widths.size() != transition_depths.size() ||
              recrystallized_grain_size.size() != transition_depths.size() )
            AssertThrow(false,
                        ExcMessage("Error: At least one list that gives input parameters for the phase transitions has the wrong size."));

          if (transition_depths.size()>1)
            for (unsigned int i=0; i<transition_depths.size()-2; ++i)
              AssertThrow(transition_depths[i]<transition_depths[i+1],
                          ExcMessage("Error: Phase transition depths have to be sorted in ascending order!"));

          // grain evolution parameters
          grain_growth_activation_energy        = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth activation energy")));
          grain_growth_activation_volume        = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth activation volume")));
          grain_growth_rate_constant            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth rate constant")));
          grain_growth_exponent                 = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth exponent")));
          reciprocal_required_strain            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Reciprocal required strain")));

          use_paleowattmeter                    = prm.get_bool ("Use paleowattmeter");
          grain_boundary_energy                 = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Average specific grain boundary energy")));
          boundary_area_change_work_fraction    = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Work fraction for boundary area change")));
          geometric_constant                    = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Geometric constant")));

          // rheology parameters
          dislocation_viscosity_iteration_threshold = prm.get_double("Dislocation viscosity iteration threshold");
          dislocation_viscosity_iteration_number = prm.get_integer("Dislocation viscosity iteration number");
          dislocation_creep_exponent            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation creep exponent")));
          dislocation_activation_energy         = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation activation energy")));
          dislocation_activation_volume         = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation activation volume")));
          dislocation_creep_prefactor           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation creep prefactor")));
          diffusion_creep_exponent              = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep exponent")));
          diffusion_activation_energy           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion activation energy")));
          diffusion_activation_volume           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion activation volume")));
          diffusion_creep_prefactor             = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep prefactor")));
          diffusion_creep_grain_size_exponent   = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep grain size exponent")));
          max_temperature_dependence_of_eta     = prm.get_double ("Maximum temperature dependence of viscosity");
          min_eta                               = prm.get_double ("Minimum viscosity");
          max_eta                               = prm.get_double ("Maximum viscosity");
          min_specific_heat                     = prm.get_double ("Minimum specific heat");
          max_specific_heat                     = prm.get_double ("Maximum specific heat");
          min_thermal_expansivity               = prm.get_double ("Minimum thermal expansivity");
          max_thermal_expansivity               = prm.get_double ("Maximum thermal expansivity");
          max_latent_heat_substeps              = prm.get_integer ("Maximum latent heat substeps");
          min_grain_size                        = prm.get_double ("Minimum grain size");
          pv_grain_size_scaling                 = prm.get_double ("Lower mantle grain size scaling");

          // scale recrystallized grain size, diffusion creep and grain growth prefactor accordingly
          diffusion_creep_prefactor[diffusion_creep_prefactor.size()-1] *= pow(pv_grain_size_scaling,diffusion_creep_grain_size_exponent[diffusion_creep_grain_size_exponent.size()-1]);
          grain_growth_rate_constant[grain_growth_rate_constant.size()-1] *= pow(pv_grain_size_scaling,grain_growth_exponent[grain_growth_exponent.size()-1]);
          if (recrystallized_grain_size.size()>0)
            recrystallized_grain_size[recrystallized_grain_size.size()-1] *= pv_grain_size_scaling;

          if (use_paleowattmeter)
            boundary_area_change_work_fraction[boundary_area_change_work_fraction.size()-1] /= pv_grain_size_scaling;



          advect_log_gransize                   = prm.get_bool ("Advect logarithm of grain size");

          if (grain_growth_activation_energy.size() != grain_growth_activation_volume.size() ||
              grain_growth_activation_energy.size() != grain_growth_rate_constant.size() ||
              grain_growth_activation_energy.size() != grain_growth_exponent.size() ||
              grain_growth_activation_energy.size() != dislocation_creep_exponent.size() ||
              grain_growth_activation_energy.size() != dislocation_activation_energy.size() ||
              grain_growth_activation_energy.size() != dislocation_activation_volume.size() ||
              grain_growth_activation_energy.size() != dislocation_creep_prefactor.size() ||
              grain_growth_activation_energy.size() != diffusion_creep_exponent.size() ||
              grain_growth_activation_energy.size() != diffusion_activation_energy.size() ||
              grain_growth_activation_energy.size() != diffusion_activation_volume.size() ||
              grain_growth_activation_energy.size() != diffusion_creep_prefactor.size() ||
              grain_growth_activation_energy.size() != diffusion_creep_grain_size_exponent.size() )
            AssertThrow(false,
                        ExcMessage("Error: The lists of grain size evolution and flow law parameters "
                                   "need to have the same length!"));

          if (use_paleowattmeter)
            {
              if (grain_growth_activation_energy.size() != grain_boundary_energy.size() ||
                  grain_growth_activation_energy.size() != boundary_area_change_work_fraction.size() ||
                  grain_growth_activation_energy.size() != geometric_constant.size() )
                AssertThrow(false,
                            ExcMessage("Error: One of the lists of grain size evolution parameters "
                                       "given for the paleowattmeter does not have the correct length!"));
            }
          else
            AssertThrow(grain_growth_activation_energy.size() == reciprocal_required_strain.size(),
                        ExcMessage("Error: The list of grain size evolution parameters in the "
                                   "paleopiezometer does not have the correct length!"));

          AssertThrow(grain_growth_activation_energy.size() == transition_depths.size()+1,
                      ExcMessage("Error: The lists of grain size evolution and flow law parameters need to "
                                 "have exactly one more entry than the number of phase transitions "
                                 "(which is defined by the length of the lists of phase transition depths, ...)!"));

          // parameters for reading in tables with material properties
          datadirectory        = prm.get ("Data directory");
          {
            const std::string      subst_text = "$ASPECT_SOURCE_DIR";
            std::string::size_type position;
            while (position = datadirectory.find (subst_text),  position!=std::string::npos)
              datadirectory.replace (datadirectory.begin()+position,
                                     datadirectory.begin()+position+subst_text.size(),
                                     ASPECT_SOURCE_DIR);
          }
          material_file_names  = Utilities::split_string_list
                                 (prm.get ("Material file names"));
          derivatives_file_names = Utilities::split_string_list
                                   (prm.get ("Derivatives file names"));
          use_table_properties = prm.get_bool ("Use table properties");
          use_enthalpy = prm.get_bool ("Use enthalpy for material properties");

          if (prm.get ("Material file format") == "perplex")
            material_file_format = perplex;
          else if (prm.get ("Material file format") == "hefesto")
            material_file_format = hefesto;
          else
            AssertThrow (false, ExcNotImplemented());

          use_bilinear_interpolation = prm.get_bool ("Bilinear interpolation");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();


      // Declare dependencies on solution variables
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;

      this->model_dependence.viscosity = NonlinearDependence::temperature
                                         | NonlinearDependence::pressure
                                         | NonlinearDependence::strain_rate
                                         | NonlinearDependence::compositional_fields;

      this->model_dependence.density = NonlinearDependence::none;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;

      if (use_table_properties)
        {
          this->model_dependence.density |= NonlinearDependence::temperature
                                            | NonlinearDependence::pressure
                                            | NonlinearDependence::compositional_fields;
          this->model_dependence.compressibility = NonlinearDependence::temperature
                                                   | NonlinearDependence::pressure
                                                   | NonlinearDependence::compositional_fields;
          this->model_dependence.specific_heat = NonlinearDependence::temperature
                                                 | NonlinearDependence::pressure
                                                 | NonlinearDependence::compositional_fields;
        }
      else
        {
          if (thermal_alpha != 0)
            this->model_dependence.density |=NonlinearDependence::temperature;
          if (reference_compressibility != 0)
            this->model_dependence.density |=NonlinearDependence::pressure;
          if (compositional_delta_rho != 0)
            this->model_dependence.density |=NonlinearDependence::compositional_fields;
        }
    }


    template <int dim>
    void
    DamageRheology<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<DislocationViscosityOutputs<dim> >() == NULL)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std_cxx11::shared_ptr<MaterialModel::AdditionalMaterialOutputs<dim> >
            (new MaterialModel::DislocationViscosityOutputs<dim> (n_points)));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(DamageRheology,
                                   "damage rheology",
                                   "A material model that behaves in the same way as "
                                   "the simple material model, but includes compositional "
                                   "fields that stand for average grain sizes of a mineral "
                                   "phase and source terms for them that determine the grain "
                                   "size evolution in dependence of the strain rate, "
                                   "temperature, phase transitions, and the creep regime. "
                                   "In the diffusion creep regime, the viscosity depends "
                                   "on this grain size."
                                   "We use the grain size evolution laws described in Behn "
                                   "et al., 2009. Implications of grain size evolution on the "
                                   "seismic structure of the oceanic upper mantle, "
                                   "Earth Planet. Sci. Letters, 282, 178–189.")
  }
}
