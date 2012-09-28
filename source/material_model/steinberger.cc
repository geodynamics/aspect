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
#include <aspect/simulator.h>
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

      class material_lookup
         {
           public:
             material_lookup(const std::string &filename)
             {

            	/* Initializing variables */

                delta_press=-1.0;
                min_press=-1.0;
                delta_temp=-1.0;
                min_temp=-1.0;
                numtemp=0;
                numpress=0;

               std::string temp;
               std::ifstream in(filename, std::ios::in);
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

               int i = 0;
               while (!in.eof())
                 {
            	   double temp1,temp2;
                   double rho,alpha,cp;
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
                   getline(in, temp);
                   if (in.eof())
                     break;

                   density_values[i%numtemp][i/numtemp]=rho;
                   thermal_expansivity_values[i%numtemp][i/numtemp]=alpha;
                   specific_heat_values[i%numtemp][i/numtemp]=cp;
                   i++;
                 }
               Assert(i==numtemp*numpress, ExcMessage("Material table size not consistent with header."));

             }

             double specific_heat(double temperature,double pressure)
             {
            	return specific_heat_values[get_idT(temperature)][get_idp(pressure)];
             }

             double density(double temperature,double pressure)
             {
            	return density_values[get_idT(temperature)][get_idp(pressure)];
             }

             double thermal_expansivity(double temperature,double pressure)
             {
            	return thermal_expansivity_values[get_idT(temperature)][get_idp(pressure)];
             }


           private:


             int get_idT(double temperature)
             {
            	temperature=std::max(min_temp, temperature);
            	temperature=std::min(temperature, max_temp);
            	Assert(temperature>=min_temp, ExcMessage("not in range"));
            	Assert(temperature<=max_temp, ExcMessage("not in range"));
            	const unsigned int idT = static_cast<unsigned int>((temperature-min_temp)/delta_temp);
            	// This check only works if all 3 parameter tables are equally sized
            	Assert(idT<density_values.n_rows(), ExcMessage("not in range"));
            	return idT;
             }

             int get_idp(double pressure)
             {
            	pressure=std::max(min_press, pressure);
            	pressure=std::min(pressure, max_press);
            	Assert(pressure>=min_press, ExcMessage("not in range"));
            	Assert(pressure<=max_press, ExcMessage("not in range"));
            	const unsigned int idp = static_cast<unsigned int>((pressure-min_press)/delta_press);
            	// This check only works if all 3 parameter tables are equally sized
            	Assert(idp<density_values.n_cols(), ExcMessage("not in range"));
            	return idp;
             }


             dealii::Table<2,double> density_values;
             dealii::Table<2,double> thermal_expansivity_values;
             dealii::Table<2,double> specific_heat_values;

             double delta_press;
             double min_press;
             double max_press;
             double delta_temp;
             double min_temp;
             double max_temp;
             int numtemp;
             int numpress;
         };

      class lateral_viscosity_lookup
          {
            public:
              lateral_viscosity_lookup(const std::string &filename)
              {
                std::string temp;
                std::ifstream in(filename, std::ios::in);
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

      class radial_viscosity_lookup
      {
        public:
          radial_viscosity_lookup(const std::string &filename)
          {
            std::string temp;
            std::ifstream in(filename, std::ios::in);
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
    Steinberger<dim>::
    update()
    {
      this->get_depth_average_temperature(avg_temp);
    }

    template <int dim>
    double
    Steinberger<dim>::
    viscosity (const double temperature,
               const double,
               const SymmetricTensor<2,dim> &,
               const Point<dim> &position) const
    {
      static internal::radial_viscosity_lookup table(datadirectory+radial_viscosity_file_name);
      static internal::lateral_viscosity_lookup lat_table(datadirectory+lateral_viscosity_file_name);

      const double depth = this->get_geometry_model().depth(position);
      const unsigned int idx = 100 * depth / this->get_geometry_model().maximal_depth();
      const double delta_temp = temperature-avg_temp[idx];
      const double adia_temp = this->get_adiabatic_conditions().temperature(position);

      const double vis_lateral_exp = -1.0*lat_table.lateral_viscosity(depth)*delta_temp/(temperature*adia_temp);

      const double vis_lateral = std::max(std::min(std::exp(vis_lateral_exp),1e4),1e-4);
      const double vis_radial = table.radial_viscosity(depth);

      return std::max(std::min(vis_lateral * vis_radial,1e23),1e19);
    }

    internal::material_lookup&
    get_material_data (std::string &datapath)
    {
      static internal::material_lookup mat(datapath);
      return mat;
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
      const double reference_density    = 4.89e-4;
      return 0;
    }


    template <int dim>
    double
    Steinberger<dim>::
    specific_heat (const double temperature,
                   const double pressure,
                   const Point<dim> &position) const
    {
      static std::string datapath = datadirectory+material_file_name;
      return get_material_data(datapath).specific_heat(temperature+get_deltat(position),pressure);
    }



    template <int dim>
    double
    Steinberger<dim>::
    thermal_conductivity (const double,
                          const double,
                          const Point<dim> &) const
    {
      return 4.7;
    }


    template <int dim>
    double
    Steinberger<dim>::
    density (const double temperature,
             const double pressure,
             const Point<dim> &position) const
    {
      static std::string datapath = datadirectory+material_file_name;
      return get_material_data(datapath).density(temperature+get_deltat(position),pressure);
    }

    template <int dim>
    double
    Steinberger<dim>::
    thermal_expansion_coefficient (const double      temperature,
                                                  const double      pressure,
                                                  const Point<dim> &position) const
    {
      static std::string datapath = datadirectory+material_file_name;
      return get_material_data(datapath).thermal_expansivity(temperature+get_deltat(position),pressure);
    }


    template <int dim>
    double
    Steinberger<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const Point<dim> &position) const
    {
      return 0.0;
    }


    template <int dim>
    bool
    Steinberger<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence) const
    {
      //Assert (false, ExcMessage("Need to go through this model and figure out the correct answer."));
      return false;
    }


    template <int dim>
    bool
    Steinberger<dim>::
    density_depends_on (const NonlinearDependence::Dependence) const
    {
      //Assert (false, ExcMessage("Need to go through this model and figure out the correct answer."));
      return false;
    }

    template <int dim>
    bool
    Steinberger<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      //Assert (false, ExcMessage("Need to go through this model and figure out the correct answer."));
      return false;
    }

    template <int dim>
    bool
    Steinberger<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      //Assert (false, ExcMessage("Need to go through this model and figure out the correct answer."));
      return false;
    }

    template <int dim>
    bool
    Steinberger<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      //Assert (false, ExcMessage("Need to go through this model and figure out the correct answer."));
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
		   prm.declare_entry ("Material file name", "pyr-ringwood88.txt",
							  Patterns::Anything (),
							  "The file name of the material data. ");
		   prm.declare_entry ("Radial viscosity file name", "radial_visc.txt",
							  Patterns::Anything (),
							  "The file name of the radial viscosity data. ");
		   prm.declare_entry ("Lateral viscosity file name", "temp_viscosity_prefactor.txt",
							  Patterns::Anything (),
							  "The file name of the lateral viscosity data. ");
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
			 material_file_name   = prm.get ("Material file name");
			 radial_viscosity_file_name   = prm.get ("Radial viscosity file name");
			 lateral_viscosity_file_name   = prm.get ("Lateral viscosity file name");
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
