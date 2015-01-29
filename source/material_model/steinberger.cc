/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/material_model/steinberger.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {


    namespace internal
    {

      class MaterialLookup
      {
        public:
          MaterialLookup(const std::string &filename,
                         const bool interpol)
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
            std::ifstream in(filename.c_str(), std::ios::in);
            AssertThrow (in,
                         ExcMessage (std::string("Couldn't open file <") + filename));

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

          double
          specific_heat(double temperature,
                        double pressure) const
          {
            return value(temperature,pressure,specific_heat_values,interpolation);
          }

          double
          density(double temperature,
                  double pressure) const
          {
            return value(temperature,pressure,density_values,interpolation);
          }

          double
          thermal_expansivity(const double temperature,
                              const double pressure) const
          {
            return value(temperature,pressure,thermal_expansivity_values,interpolation);
          }

          double
          seismic_Vp(const double temperature,
                     const double pressure) const
          {
            return value(temperature,pressure,vp_values,false);
          }

          double
          seismic_Vs(const double temperature,
                     const double pressure) const
          {
            return value(temperature,pressure,vs_values,false);
          }

          double
          dHdT (const double temperature,
                const double pressure) const
          {
            const double h = value(temperature,pressure,enthalpy_values,interpolation);
            const double dh = value(temperature+delta_temp,pressure,enthalpy_values,interpolation);
            return (dh - h) / delta_temp;
          }

          double
          dHdp (const double temperature,
                const double pressure) const
          {
            const double h = value(temperature,pressure,enthalpy_values,interpolation);
            const double dh = value(temperature,pressure+delta_press,enthalpy_values,interpolation);
            return (dh - h) / delta_press;
          }

          double
          dRhodp (const double temperature,
                  const double pressure) const
          {
            const double rho = value(temperature,pressure,density_values,interpolation);
            const double drho = value(temperature,pressure+delta_press,density_values,interpolation);
            return (drho - rho) / delta_press;
          }

          double
          value (const double temperature,
                 const double pressure,
                 const dealii::Table<2,
                 double> &values,
                 bool interpol) const
          {
            const double nT = get_nT(temperature);
            const unsigned int inT = static_cast<unsigned int>(nT);

            const double np = get_np(pressure);
            const unsigned int inp = static_cast<unsigned int>(np);

            Assert(inT<values.n_rows(), ExcMessage("not in range"));
            Assert(inp<values.n_cols(), ExcMessage("not in range"));

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



        private:


          double get_nT(double temperature) const
          {
            temperature=std::max(min_temp, temperature);
            temperature=std::min(temperature, max_temp-delta_temp);
            Assert(temperature>=min_temp, ExcMessage("not in range"));
            Assert(temperature<=max_temp, ExcMessage("not in range"));
            return (temperature-min_temp)/delta_temp;
          }

          double get_np(double pressure) const
          {
            pressure=std::max(min_press, pressure);
            pressure=std::min(pressure, max_press-delta_press);
            Assert(pressure>=min_press, ExcMessage("not in range"));
            Assert(pressure<=max_press, ExcMessage("not in range"));
            return (pressure-min_press)/delta_press;
          }


          dealii::Table<2,double> density_values;
          dealii::Table<2,double> thermal_expansivity_values;
          dealii::Table<2,double> specific_heat_values;
          dealii::Table<2,double> vp_values;
          dealii::Table<2,double> vs_values;
          dealii::Table<2,double> enthalpy_values;


          double delta_press;
          double min_press;
          double max_press;
          double delta_temp;
          double min_temp;
          double max_temp;
          unsigned int numtemp;
          unsigned int numpress;
          bool interpolation;
      };

      class LateralViscosityLookup
      {
        public:
          LateralViscosityLookup(const std::string &filename)
          {
            std::string temp;
            std::ifstream in(filename.c_str(), std::ios::in);
            AssertThrow (in,
                         ExcMessage (std::string("Couldn't open file <") + filename));

            getline(in, temp); // eat first line

            min_depth=1e20;
            max_depth=-1;

            while (!in.eof())
              {
                double visc, depth;
                in >> visc;;
                if (in.eof())
                  break;
                in >> depth;
                depth *=1000.0;
                getline(in, temp);

                min_depth = std::min(depth, min_depth);
                max_depth = std::max(depth, max_depth);

                values.push_back(visc);
              }
            delta_depth = (max_depth-min_depth)/(values.size()-1);
          }

          double lateral_viscosity(double depth)
          {
            depth=std::max(min_depth, depth);
            depth=std::min(depth, max_depth);

            Assert(depth>=min_depth, ExcMessage("not in range"));
            Assert(depth<=max_depth, ExcMessage("not in range"));
            const unsigned int idx = static_cast<unsigned int>((depth-min_depth)/delta_depth);
            Assert(idx<values.size(), ExcMessage("not in range"));
            return values[idx];
          }

          int get_nslices() const
          {
            return values.size();
          }

        private:
          std::vector<double> values;
          double min_depth;
          double delta_depth;
          double max_depth;

      };

      class RadialViscosityLookup
      {
        public:
          RadialViscosityLookup(const std::string &filename)
          {
            std::string temp;
            std::ifstream in(filename.c_str(), std::ios::in);
            AssertThrow (in,
                         ExcMessage (std::string("Couldn't open file <") + filename));

            min_depth=1e20;
            max_depth=-1;

            while (!in.eof())
              {
                double visc, depth;
                in >> visc;;
                if (in.eof())
                  break;
                in >> depth;
                depth *=1000.0;
                getline(in, temp);

                min_depth = std::min(depth, min_depth);
                max_depth = std::max(depth, max_depth);

                values.push_back(visc);
              }
            delta_depth = (max_depth-min_depth)/(values.size()-1);
          }

          double radial_viscosity(double depth)
          {
            depth=std::max(min_depth, depth);
            depth=std::min(depth, max_depth);

            Assert(depth>=min_depth, ExcMessage("not in range"));
            Assert(depth<=max_depth, ExcMessage("not in range"));
            const unsigned int idx = static_cast<unsigned int>((depth-min_depth)/delta_depth);
            Assert(idx<values.size(), ExcMessage("not in range"));
            return values[idx];
          }

        private:
          std::vector<double> values;
          double min_depth;
          double delta_depth;
          double max_depth;

      };
    }



    template <int dim>
    void
    Steinberger<dim>::initialize()
    {

      n_material_data = material_file_names.size();
      for (unsigned i = 0; i < n_material_data; i++)
        material_lookup.push_back(std_cxx1x::shared_ptr<internal::MaterialLookup>
                                  (new internal::MaterialLookup(datadirectory+material_file_names[i],interpolation)));
      lateral_viscosity_lookup.reset(new internal::LateralViscosityLookup(datadirectory+lateral_viscosity_file_name));
      radial_viscosity_lookup.reset(new internal::RadialViscosityLookup(datadirectory+radial_viscosity_file_name));
      avg_temp.resize(lateral_viscosity_lookup->get_nslices());
    }



    template <int dim>
    void
    Steinberger<dim>::
    update()
    {
      if (use_lateral_average_temperature)
        this->get_depth_average_temperature(avg_temp);
    }



    template <int dim>
    double
    Steinberger<dim>::
    viscosity (const double temperature,
               const double /*pressure*/,
               const std::vector<double> &compositional_fields,
               const SymmetricTensor<2,dim> &,
               const Point<dim> &position) const
    {
      const double depth = this->get_geometry_model().depth(position);
      const double adiabatic_temperature = this->get_adiabatic_conditions().temperature(position);

      double delta_temperature;
      if (use_lateral_average_temperature)
        {
          const unsigned int idx = static_cast<unsigned int>((avg_temp.size()-1) * depth / this->get_geometry_model().maximal_depth());
          delta_temperature = temperature-avg_temp[idx];
        }
      else
        delta_temperature = temperature-adiabatic_temperature;

      // For an explanation on this formula see the Steinberger & Calderwood 2006 paper
      const double vis_lateral_exp = -1.0*lateral_viscosity_lookup->lateral_viscosity(depth)*delta_temperature/(temperature*adiabatic_temperature);
      // Limit the lateral viscosity variation to a reasonable interval
      const double vis_lateral = std::max(std::min(std::exp(vis_lateral_exp),max_lateral_eta_variation),1/max_lateral_eta_variation);

      const double vis_radial = radial_viscosity_lookup->radial_viscosity(depth);

      return std::max(std::min(vis_lateral * vis_radial,max_eta),min_eta);
    }



    template <int dim>
    double
    Steinberger<dim>::
    get_corrected_temperature (const double temperature,
                               const double pressure,
                               const Point<dim> &position) const
    {
      if (!(this->get_adiabatic_conditions().is_initialized())
          || this->include_adiabatic_heating()
          || compressible)
        return temperature;

      return temperature
             + this->get_adiabatic_conditions().temperature(position)
             - this->get_adiabatic_surface_temperature();
    }



    template <int dim>
    double
    Steinberger<dim>::
    get_corrected_pressure (const double temperature,
                            const double pressure,
                            const Point<dim> &position) const
    {
      if (!(this->get_adiabatic_conditions().is_initialized())
          || compressible)
        return pressure;

      return this->get_adiabatic_conditions().pressure(position);
    }

    template <int dim>
    double
    Steinberger<dim>::
    get_corrected_density (const double temperature,
                           const double pressure,
                           const std::vector<double> &compositional_fields,
                           const Point<dim> &position) const
    {
      const double rho = get_compressible_density(temperature,pressure,compositional_fields,position);

      const double adiabatic_temperature = this->get_adiabatic_conditions().temperature(position);
      const double adiabatic_rho = get_compressible_density(adiabatic_temperature,
                                                            pressure,
                                                            compositional_fields,
                                                            position);

      const Point<dim> surface_point = this->get_geometry_model().representative_point(0.0);
      const double surface_temperature = this->get_adiabatic_surface_temperature();
      const double surface_pressure = this->get_surface_pressure();
      const double surface_rho = get_compressible_density(surface_temperature,
                                                          surface_pressure,
                                                          compositional_fields,
                                                          surface_point);

      //Return the density scaled to an incompressible profile
      const double scaled_density = (rho / adiabatic_rho) * surface_rho;
      return scaled_density;
    }



    template <int dim>
    double
    Steinberger<dim>::
    reference_viscosity () const
    {
      return reference_eta;
    }



    template <int dim>
    double
    Steinberger<dim>::
    reference_density () const
    {
      const double reference_density    = 3300e0;
      return reference_density;
    }



    template <int dim>
    double
    Steinberger<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return 0;
    }


    template <int dim>
    double
    Steinberger<dim>::
    specific_heat (const double temperature,
                   const double pressure,
                   const std::vector<double> &compositional_fields,
                   const Point<dim> &position) const
    {
      double cp = 0.0;
      if (!latent_heat)
        {
          if (n_material_data == 1)
            cp = material_lookup[0]->specific_heat(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                cp += compositional_fields[i] * material_lookup[i]->specific_heat(temperature,pressure);
            }
        }
      else
        {
          if (n_material_data == 1)
            cp = material_lookup[0]->dHdT(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                cp += compositional_fields[i] * material_lookup[i]->dHdT(temperature,pressure);
              cp = std::max(std::min(cp,6000.0),500.0);
            }
        }
      return cp;
    }



    template <int dim>
    double
    Steinberger<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &,
                          const Point<dim> &) const
    {
      return 4.7;
    }



    template <int dim>
    double
    Steinberger<dim>::
    get_compressible_density (const double temperature,
                              const double pressure,
                              const std::vector<double> &compositional_fields,
                              const Point<dim> &position) const
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

    template <int dim>
    double
    Steinberger<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &compositional_fields,
             const Point<dim> &position) const
    {
      if (compressible
          || !(this->get_adiabatic_conditions().is_initialized()))
        return get_compressible_density(temperature,pressure,compositional_fields,position);
      else
        return get_corrected_density(temperature,pressure,compositional_fields,position);
    }



    template <int dim>
    double
    Steinberger<dim>::
    thermal_expansion_coefficient (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &position) const
    {
      double alpha = 0.0;
      if (!latent_heat)
        {
          if (n_material_data == 1)
            alpha = material_lookup[0]->thermal_expansivity(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                alpha += compositional_fields[i] * material_lookup[i]->thermal_expansivity(temperature,pressure);
            }
        }
      else
        {
          double dHdp = 0.0;
          if (n_material_data == 1)
            dHdp += material_lookup[0]->dHdp(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                dHdp += compositional_fields[i] * material_lookup[i]->dHdp(temperature,pressure);
            }
          alpha = (1 - density(temperature,pressure,compositional_fields,position) * dHdp) / temperature;
          alpha = std::max(std::min(alpha,1e-3),1e-5);
        }
      return alpha;
    }



    template <int dim>
    double
    Steinberger<dim>::
    seismic_Vp (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &position) const
    {
      //this function is not called from evaluate() so we need to care about
      //corrections for temperature and pressure
      const double corrected_temperature = get_corrected_temperature(temperature,pressure,position);
      const double corrected_pressure = get_corrected_pressure(temperature,pressure,position);

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
    Steinberger<dim>::
    seismic_Vs (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &position) const
    {
      //this function is not called from evaluate() so we need to care about
      //corrections for temperature and pressure
      const double corrected_temperature = get_corrected_temperature(temperature,pressure,position);
      const double corrected_pressure = get_corrected_pressure(temperature,pressure,position);


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
    Steinberger<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields,
                     const Point<dim> &position) const
    {
      if (!compressible)
        return 0.0;

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
    bool
    Steinberger<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else
        return false;
    }



    template <int dim>
    bool
    Steinberger<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::pressure) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }



    template <int dim>
    bool
    Steinberger<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::pressure) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }



    template <int dim>
    bool
    Steinberger<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::pressure) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }



    template <int dim>
    bool
    Steinberger<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    Steinberger<dim>::
    is_compressible () const
    {
      return compressible;
    }

    template <int dim>
    void
    Steinberger<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                               typename Interface<dim>::MaterialModelOutputs &out) const
    {

      Assert ((n_material_data <= in.composition[0].size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const double temperature = get_corrected_temperature(in.temperature[i],
                                                               in.pressure[i],
                                                               in.position[i]);
          const double pressure    = get_corrected_pressure(in.temperature[i],
                                                            in.pressure[i],
                                                            in.position[i]);

          /* We are only asked to give viscosities if strain_rate.size() > 0
           * and we can only calculate it if adiabatic_conditions are available.
           * Note that the used viscosity formulation needs the not
           * corrected temperatures in case we compare it to the lateral
           * temperature average.
           */
          if (this->get_adiabatic_conditions().is_initialized() && in.strain_rate.size())
            {
              if (use_lateral_average_temperature)
                {
                  out.viscosities[i]            = viscosity                     (in.temperature[i], in.pressure[i], in.composition[i], in.strain_rate[i], in.position[i]);
                }
              else
                {
                  out.viscosities[i]            = viscosity                     (temperature, pressure, in.composition[i], in.strain_rate[i], in.position[i]);
                }
            }
          out.densities[i]                      = density                       (temperature, pressure, in.composition[i], in.position[i]);
          out.thermal_expansion_coefficients[i] = thermal_expansion_coefficient (temperature, pressure, in.composition[i], in.position[i]);
          out.specific_heat[i]                  = specific_heat                 (temperature, pressure, in.composition[i], in.position[i]);
          out.thermal_conductivities[i]         = thermal_conductivity          (temperature, pressure, in.composition[i], in.position[i]);
          out.compressibilities[i]              = compressibility               (temperature, pressure, in.composition[i], in.position[i]);
          out.entropy_derivative_pressure[i]    = 0;
          out.entropy_derivative_temperature[i] = 0;
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c]            = 0;
        }
    }


    template <int dim>
    void
    Steinberger<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Steinberger model");
        {
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
          prm.declare_entry ("Radial viscosity file name", "radial-visc.txt",
                             Patterns::Anything (),
                             "The file name of the radial viscosity data. ");
          prm.declare_entry ("Lateral viscosity file name", "temp-viscosity-prefactor.txt",
                             Patterns::Anything (),
                             "The file name of the lateral viscosity data. ");
          prm.declare_entry ("Use lateral average temperature for viscosity", "true",
                             Patterns::Bool (),
                             "Whether to use to use the laterally averaged temperature "
                             "instead of the adiabatic temperature for the viscosity "
                             "calculation. This ensures that the laterally averaged "
                             "viscosities remain more or less constant over the model "
                             "runtime. This behaviour might or might not be desired.");
          prm.declare_entry ("Bilinear interpolation", "true",
                             Patterns::Bool (),
                             "Whether to use bilinear interpolation to compute "
                             "material properties (slower but more accurate). ");
          prm.declare_entry ("Latent heat", "false",
                             Patterns::Bool (),
                             "Whether to include latent heat effects in the "
                             "calculation of thermal expansivity and specific heat. "
                             "Following the approach of Nakagawa et al. 2009. ");
          prm.declare_entry ("Compressible", "false",
                             Patterns::Bool (),
                             "Whether to include a compressible material description."
                             "For a description see the manual section. ");
          prm.declare_entry ("Reference viscosity", "1e23",
                             Patterns::Double(0),
                             "The reference viscosity that is used for pressure scaling. ");
          prm.declare_entry ("Minimum viscosity", "1e19",
                             Patterns::Double(0),
                             "The minimum viscosity that is allowed in the viscosity "
                             "calculation. Smaller values will be cut off.");
          prm.declare_entry ("Maximum viscosity", "1e23",
                             Patterns::Double(0),
                             "The maximum viscosity that is allowed in the viscosity "
                             "calculation. Larger values will be cut off.");
          prm.declare_entry ("Maximum lateral viscosity variation", "1e2",
                             Patterns::Double(0),
                             "The relative cutoff value for lateral viscosity variations "
                             "caused by temperature deviations. The viscosity may vary "
                             "laterally by this factor squared.");
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }



    template <int dim>
    void
    Steinberger<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Steinberger model");
        {
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
          radial_viscosity_file_name   = prm.get ("Radial viscosity file name");
          lateral_viscosity_file_name  = prm.get ("Lateral viscosity file name");
          use_lateral_average_temperature = prm.get_bool ("Use lateral average temperature for viscosity");
          interpolation        = prm.get_bool ("Bilinear interpolation");
          latent_heat          = prm.get_bool ("Latent heat");
          compressible         = prm.get_bool ("Compressible");
          reference_eta        = prm.get_double ("Reference viscosity");
          min_eta              = prm.get_double ("Minimum viscosity");
          max_eta              = prm.get_double ("Maximum viscosity");
          max_lateral_eta_variation    = prm.get_double ("Maximum lateral viscosity variation");

          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Steinberger,
                                   "Steinberger",
                                   "This material model looks up the viscosity from the tables that "
                                   "correspond to the paper of Steinberger and Calderwood "
                                   "2006 (``Models of large-scale viscous flow in the Earth's "
                                   "mantle with constraints from mineral physics and surface observations'', "
                                   "Geophys. J. Int., 167, 1461-1481, "
                                   "\\url{http://dx.doi.org/10.1111/j.1365-246X.2006.03131.x}) and material "
                                   "data from a database generated by the thermodynamics code \\texttt{Perplex}, "
                                   "see \\url{http://www.perplex.ethz.ch/}. "
                                   "The default example data builds upon the thermodynamic "
                                   "database by Stixrude 2011 and assumes a pyrolitic composition by "
                                   "Ringwood 1988 but is easily replaceable by other data files. ")
  }
}
