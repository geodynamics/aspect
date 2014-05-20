/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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
/*  $Id  */

//#define _USE_MATH_DEFINES
#include <cmath>
#include <aspect/initial_conditions/solidus.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/boundary_temperature/interface.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace aspect
{
    namespace InitialConditions
    {
        void Melting_curve::read(const std::string &filename)
        {
            data_filename=filename;
            std::ifstream in(data_filename.c_str(), std::ios::in);
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
                    if(T_Unit=="C")
                    {
                        T+=273.15;                      // Degree C to K
                    }
                    else if(T_Unit!="K")
                    {
                        AssertThrow(false,ExcMessage ("Unit of the first column of melting curve data "
                                                            "has to be one of the following: C/K."))
                    }
                    
                    is_radius=false;                    // Second column is pressure
                    if(P_Unit=="kbar")
                        p*=1.e8;                        // kbar to Pa
                    else if(P_Unit=="GPa")
                        p*=1.e9;                        // GPa to Pa
                    else 
                    {
                        is_radius=true;                 // Second column in radius instead of pressure
                        if(P_Unit=="km")
                            p*=1.e3;                    // km to meters
                        else if(P_Unit!="m")
                            AssertThrow(false,ExcMessage ("Unit of the second column of melting curve data "
                                                            "has to be one of the following: Pa/Gpa/km/m."))
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
            if(T_array.size()==0)return(0);
            for(unsigned i=1;i<Num_points;i++)
            {
                if(     (i==Num_points-1) || 
                        (is_radius && P_value>P_array[i]) ||
                        (!is_radius && P_value<P_array[i]) )
                {
                    T_value=T_array[i-1]+(T_array[i]-T_array[i-1])/(P_array[i]-P_array[i-1])*(P_value-P_array[i-1]);
                    return(T_value);
                }
            }
            AssertThrow(false,ExcMessage(std::string("Something wrong with the melting curve data ")+ data_filename ));
            return(-1);
        }
        
    }

  namespace InitialConditions
  {
    template <int dim>
    Solidus<dim>::Solidus ():
    solidus_curve()
    {}

    template <int dim>
    double
    Solidus<dim>::
    initial_temperature (const Point<dim> &position) const
    {
        double T_min,T_litho;
        double T_solidus,T_perturbation;
        double litho_thick_theta;
        double lateral_perturbation;
        double Depth=this->geometry_model->depth(position);

        AssertThrow(solidus_curve.is_radius==true,ExcMessage("The solidus curve has to be radius dependent."));
        AssertThrow(solidus_curve.Num_points!=0,ExcMessage("Error eading solidus file."));
        const GeometryModel::SphericalShell<dim> *spherical_geometry_model=
                dynamic_cast< const GeometryModel::SphericalShell<dim> *>(this->geometry_model);
        AssertThrow(spherical_geometry_model!=0,
                ExcMessage("This initial condition can only be work with spherical shell geometry model."));
        T_min=(this->get_boundary_temperature()).minimal_temperature();

        // In case of spherical shell calculate spherical coordinates
        const Tensor<1,dim> scoord = spherical_surface_coordinates(position);
        if(dim==2)
        {
           // Use a sine as lateral perturbation that is scaled to the opening angle of the geometry.
           // This way the perturbation is alway 0 at the model boundaries.
           const double opening_angle = spherical_geometry_model->opening_angle()*numbers::PI/180.0;
           lateral_perturbation = std::sin(lateral_wave_number_1*scoord[1]*numbers::PI/opening_angle);	
        }
        else if(dim==3)
        {
            // Spherical harmonics are only defined for order <= degree
            // and degree >= 0. Verify that it is indeed.
            Assert ( std::abs(lateral_wave_number_2) <= lateral_wave_number_1,
                       ExcMessage ("Spherical harmonics can only be computed for "
                                   "order <= degree."));
            Assert ( lateral_wave_number_1 >= 0,
                       ExcMessage ("Spherical harmonics can only be computed for "
                                   "degree >= 0."));
            // use a spherical harmonic function as lateral perturbation
            lateral_perturbation = boost::math::spherical_harmonic_r(lateral_wave_number_1,lateral_wave_number_2,scoord[2],scoord[1]);
        }
        litho_thick_theta=litho_thick-magnitude_lith*lateral_perturbation;
        T_litho=solidus_curve.T(0,spherical_geometry_model->R1-litho_thick_theta)+deltaT;
        
        if(litho_thick_theta>0 && Depth<litho_thick_theta)
            T_solidus=T_min+(T_litho-T_min)*(Depth/litho_thick_theta);
        else
            T_solidus=solidus_curve.T(0,sqrt(position.square()))+deltaT;

        T_perturbation=Depth/( this->geometry_model->maximal_depth() )*magnitude_T*lateral_perturbation;
        return T_solidus+T_perturbation;
    }
    
    
    
    template <int dim>
    const Tensor<1,dim>
    Solidus<dim>::
    spherical_surface_coordinates(const Tensor<1,dim> &position) const
    {
      Tensor<1,dim> scoord;

      scoord[0] = std::sqrt(position.norm_square());                // R
      scoord[1] = std::atan2(position[1],position[0]);              // Phi
      if (scoord[1] < 0.0) scoord[1] = 2*numbers::PI + scoord[1];   // correct phi to [0,2*pi]
      if (dim==3)
        scoord[2] = std::acos(position[2]/std::sqrt(position.norm_square())); // Theta

      return scoord;
    }
    
    

    template <int dim>
    void
    Solidus<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Solidus");
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
                prm.declare_entry ("Lateral wave number one","3",
                                   Patterns::Integer(),
                                   "Doubled first lateral wave number of the harmonic perturbation. "
                                   "Equals the spherical harmonic degree in 3D spherical shells. "
                                   "In all other cases one equals half of a sine period over "
                                   "the model domain. This allows for single up-/downswings. "
                                   "Negative numbers reverse the sign of the perturbation but are "
                                   "not allowed for the spherical harmonic case.");
                prm.declare_entry ("Lateral wave number two", "2",
                                   Patterns::Integer (),
                                   "Doubled second lateral wave number of the harmonic perturbation. "
                                   "Equals the spherical harmonic order in 3D spherical shells. "
                                   "In all other cases one equals half of a sine period over "
                                   "the model domain. This allows for single up-/downswings. "
                                   "Negative numbers reverse the sign of the perturbation.");
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
      prm.leave_subsection();
    }



    template <int dim>
    void
    Solidus<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Solidus");
        {
            deltaT=prm.get_double("Supersolidus");
            litho_thick=prm.get_double("Lithosphere thickness");
            prm.enter_subsection("Perturbation");
            {
                magnitude_T    = prm.get_double("Temperature amplitude");
                magnitude_lith = prm.get_double("Lithosphere thickness amplitude");
                lateral_wave_number_1 = prm.get_integer ("Lateral wave number one");
                lateral_wave_number_2 = prm.get_integer ("Lateral wave number two");
            }
            prm.leave_subsection();
            prm.enter_subsection("Data");
            {
                solidus_filename=prm.get ("Solidus filename");
                solidus_curve.read(solidus_filename);
            }
            prm.leave_subsection();
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
                                       "with lithosphere thickness and perturbation "
                                       "in shpere harmonic functions.");
  }
}
