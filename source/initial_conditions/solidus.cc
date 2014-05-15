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
/*  $Id: function.cc 1389 2012-11-28 14:08:16Z bangerth $  */

//#define _USE_MATH_DEFINES
#include <cmath>
#include <aspect/initial_conditions/solidus.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/boundary_temperature/interface.h>

namespace aspect
{
	namespace melting
	{
		Melting_curve::Melting_curve(const std::string &filename)
		{
			read(filename);
		}
		void Melting_curve::read(const std::string &filename)
		{
			std::ifstream in(filename.c_str(), std::ios::in);
			char temp[256];
			std::string T_Unit,P_Unit;
          	        Num_points=0;
			if(in.fail())return;
			in.getline(temp,256);
			in>>T_Unit>>P_Unit;
			in.getline(temp,256);
          	while(!in.eof())
          	{
            	double T,p;
                in>>T>>p;
                if(!in.fail())
                {
					//Unit switching
					if(T_Unit=="C")    T+=273.15; // Degree C to K
					if(P_Unit=="kbar") p*=1.e8;   // kbar to Pa
					if(P_Unit=="GPa")  p*=1.e9;   // GPa to Pa
					is_radius=false;              // Second column is pressure
					if(P_Unit=="km")
					{
						is_radius=true;           // Second column in radius instead of pressure
						p*=1.e3;                  // km to meters
					}
                	T_array.push_back(T);
                	P_array.push_back(p);
					Num_points++;
              	}
				in.getline(temp,256);
          	}
      	}
		
		double Melting_curve::T(const double p, const double radius) const
		{
			double T_value,P_value=is_radius?radius:p;
			if(T_array.size()==0)return(0.);
			for(unsigned i=1;i<Num_points;i++)
			{
				if(     (i==Num_points-1) || 
						(is_radius && P_value>P_array[i]) ||
						(!is_radius && P_value<P_array[i]) )
				{
					T_value=T_array[i-1]+(T_array[i]-T_array[i-1])/(P_array[i]-P_array[i-1])*(P_value-P_array[i-1]);
					break;
				}
			}
			return(T_value);
		}
		
	}

  namespace InitialConditions
  {
    template <int dim>
    Solidus<dim>::Solidus ()
    {}

    template <int dim>
    double
    Solidus<dim>::
    initial_temperature (const Point<dim> &position) const
    {
		double R1,T_min,T_litho;
		double Theta=atan2(position(0),position(1));
		double T_solidus,T_perturbation;
		double litho_thick_theta;
		double Depth=this->geometry_model->depth(position);
		static melting::Melting_curve Solidus_curve(solidus_filename);

		AssertThrow(Solidus_curve.is_radius==true,ExcMessage("The solidus curve has to be depth dependent."));
		AssertThrow(Solidus_curve.Num_points!=0,ExcMessage("Error eading solidus file."));
                AssertThrow(dynamic_cast< const GeometryModel::SphericalShell<dim> *>( &this->get_geometry_model() )!=0,
				ExcMessage("This initial condition can only be work with sphereical shell geometry model."));
        R1=(dynamic_cast< const GeometryModel::SphericalShell<dim> &>(this->get_geometry_model())).R1;
		T_min=(this->get_boundary_temperature()).minimal_temperature();

		litho_thick_theta=litho_thick-Magnitude_lith*sin(n*Theta);
		T_litho=Solidus_curve.T(0,R1-litho_thick_theta)+deltaT;
		
		if(litho_thick_theta>0 && Depth<litho_thick_theta)
			T_solidus=T_min+(T_litho-T_min)*(Depth/litho_thick_theta);
		else
			T_solidus=Solidus_curve.T(0,sqrt(position.square()))+deltaT;
		
		//Currently only work in 2D
		T_perturbation=Depth/( this->geometry_model->maximal_depth() )*Magnitude_T*sin(n*Theta);
        return T_solidus+T_perturbation;
    }

    template <int dim>
    void
    Solidus<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
		prm.declare_entry ("Supersolidus","0e0",
				           Patterns::Double (),
						   "The difference from solidus.");
        prm.declare_entry ("Lithosphere thickness","0",
                           Patterns::Double (0),
                           "The thickness of lithosphere thickness. Unit: m");        
		prm.enter_subsection("Perturbation");
        {
			prm.declare_entry ("Temperature amplitude", "0e0",
					           Patterns::Double (0),
							   "The amplitude of the initial spherical temperature perturbation in (K)");
			prm.declare_entry ("Lithosphere thickness amplitude", "0e0",
					           Patterns::Double (),
							   "The amplitude of the initial lithosphere thickness perturbation in (m)");
			prm.declare_entry ("Degree","1",
					           Patterns::Integer(0),
							   "The degree of the perturbation.");
        }
        prm.leave_subsection();
		prm.enter_subsection ("Data");
		{
        	prm.declare_entry ("Solidus filename", "",
                               Patterns::Anything(),
                               "The solidus filename.");
      	}
        prm.leave_subsection();
      }
	  prm.leave_subsection();
	}


    template <int dim>
    void
    Solidus<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
		deltaT=prm.get_double("Supersolidus");
		litho_thick=prm.get_double("Lithosphere thickness");
		prm.enter_subsection("Perturbation");
		{
			Magnitude_T    = prm.get_double("Temperature amplitude");
			Magnitude_lith = prm.get_double("Lithosphere thickness amplitude");
			n              = prm.get_integer("Degree");

		}
	    prm.leave_subsection();
        prm.enter_subsection("Data");
        {
			solidus_filename=prm.get ("Solidus filename");
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
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(Solidus,
                                       "solidus",
                                       "Temperature initial condition as solidus," 
									   "with perturbation as sin() funciton.");
  }
}
