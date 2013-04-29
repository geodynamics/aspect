/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */


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
            max_press = min_press + (numtemp-1) * delta_press;

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
          value (const double temperature,
                 const double pressure,
                 const dealii::Table<2,
                 double>& values,
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
            temperature=std::max(min_temp+delta_temp, temperature);
            temperature=std::min(temperature, max_temp-delta_temp);
            Assert(temperature>=min_temp, ExcMessage("not in range"));
            Assert(temperature<=max_temp, ExcMessage("not in range"));
            return (temperature-min_temp)/delta_temp;
          }

          double get_np(double pressure) const
          {
            pressure=std::max(min_press+delta_press, pressure);
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
    }



    template <int dim>
    void
    Steinberger<dim>::
    update()
    {
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
      const unsigned int idx = static_cast<unsigned int>(avg_temp.size() * depth / this->get_geometry_model().maximal_depth());
      const double delta_temp = temperature-avg_temp[idx];
      const double adia_temp = this->get_adiabatic_conditions().temperature(position);

      const double vis_lateral_exp = -1.0*lateral_viscosity_lookup->lateral_viscosity(depth)*delta_temp/(temperature*adia_temp);

      const double vis_lateral = std::max(std::min(std::exp(vis_lateral_exp),1e2),1e-2);
      const double vis_radial = radial_viscosity_lookup->radial_viscosity(depth);

      return std::max(std::min(vis_lateral * vis_radial,1e23),1e19);
    }



    template <int dim>
    double
    Steinberger<dim>::
    get_deltat (const Point<dim> &position) const
    {
      if (!(&this->get_adiabatic_conditions()))
        return 0.0;
      static const bool a = this->include_adiabatic_heating();
      return a ? 0.0 : (this->get_adiabatic_conditions().temperature(position)
                        - this->get_adiabatic_surface_temperature());
    }



    template <int dim>
    double
    Steinberger<dim>::
    reference_viscosity () const
    {
      const double reference_eta    = 1e23;
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
      Assert ((n_material_data <= compositional_fields.size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));
      double cp = 0.0;
      if (!latent_heat)
        {
          if (n_material_data == 1)
            cp = material_lookup[0]->specific_heat(temperature+get_deltat(position),pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                cp += compositional_fields[i] * material_lookup[i]->specific_heat(temperature+get_deltat(position),pressure);
            }
        }
      else
        {
          if (n_material_data == 1)
            cp = material_lookup[0]->dHdT(temperature+get_deltat(position),pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                cp += compositional_fields[i] * material_lookup[i]->dHdT(temperature+get_deltat(position),pressure);
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
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return 4.7;
    }



    template <int dim>
    double
    Steinberger<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &compositional_fields,
             const Point<dim> &position) const
    {
      Assert ((n_material_data <= compositional_fields.size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));
      double rho = 0.0;
      if (n_material_data == 1)
        rho = material_lookup[0]->density(temperature+get_deltat(position),pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            rho += compositional_fields[i] * material_lookup[i]->density(temperature+get_deltat(position),pressure);
        }
      return rho;
    }



    template <int dim>
    double
    Steinberger<dim>::
    thermal_expansion_coefficient (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &position) const
    {
      Assert ((n_material_data <= compositional_fields.size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));
      double alpha = 0.0;
      if (!latent_heat)
        {
          if (n_material_data == 1)
            alpha = material_lookup[0]->thermal_expansivity(temperature+get_deltat(position),pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                alpha += compositional_fields[i] * material_lookup[i]->thermal_expansivity(temperature+get_deltat(position),pressure);
            }
        }
      else
        {
          double dHdp = 0.0;
          if (n_material_data == 1)
            dHdp += material_lookup[0]->dHdp(temperature+get_deltat(position),pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                dHdp += compositional_fields[i] * material_lookup[i]->dHdp(temperature+get_deltat(position),pressure);
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
      Assert ((n_material_data <= compositional_fields.size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));
      double vp = 0.0;
      if (n_material_data == 1)
        vp += material_lookup[0]->seismic_Vp(temperature+get_deltat(position),pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            vp += compositional_fields[i] * material_lookup[i]->seismic_Vp(temperature+get_deltat(position),pressure);
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
      Assert ((n_material_data <= compositional_fields.size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));
      double vs = 0.0;
      if (n_material_data == 1)
        vs += material_lookup[0]->seismic_Vs(temperature+get_deltat(position),pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            vs += compositional_fields[i] * material_lookup[i]->seismic_Vs(temperature+get_deltat(position),pressure);
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
      return 0.0;
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
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
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
      return false;
    }



    template <int dim>
    void
    Steinberger<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Steinberger model");
        {
          prm.declare_entry ("Data directory", "data/material-model/steinberger/",
                             Patterns::DirectoryName (),
                             "The path to the model data. ");
          prm.declare_entry ("Material file names", "pyr-ringwood88.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file names of the material data. "
                             "List with as many components as active"
                             "compositional fields (material data is assumed to"
                             "be in order with the ordering of the fields). ");
          prm.declare_entry ("Radial viscosity file name", "radial-visc.txt",
                             Patterns::Anything (),
                             "The file name of the radial viscosity data. ");
          prm.declare_entry ("Lateral viscosity file name", "temp-viscosity-prefactor.txt",
                             Patterns::Anything (),
                             "The file name of the lateral viscosity data. ");
          prm.declare_entry ("Bilinear interpolation", "true",
                             Patterns::Bool (),
                             "whether to use bilinear interpolation to compute "
                             "material properties (slower but more accurate).");
          prm.declare_entry ("Latent heat", "false",
                             Patterns::Bool (),
                             "whether to include latent heat effects in the"
                             "calculation of thermal expansivity and specific heat."
                             "Following the approach of Nakagawa et al. 2009.");
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
          material_file_names  = Utilities::split_string_list
                                 (prm.get ("Material file names"));
          radial_viscosity_file_name   = prm.get ("Radial viscosity file name");
          lateral_viscosity_file_name   = prm.get ("Lateral viscosity file name");
          interpolation        = prm.get_bool ("Bilinear interpolation");
          latent_heat          = prm.get_bool ("Latent heat");

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
                                   "lookup viscosity from the paper of Steinberger/Calderwood"
                                   "2006 and material data from a database generated by Perplex. "
                                   "The database builds upon the thermodynamic database by "
                                   "Stixrude 2011 and assumes a pyrolitic composition by "
                                   "Ringwood 1988. ")
  }
}
