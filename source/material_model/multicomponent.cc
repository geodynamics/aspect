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


#include <aspect/material_model/multicomponent.h>
#include <deal.II/base/parameter_handler.h>

#include <numeric>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    Multicomponent<dim>::
    compute_volume_fractions( const std::vector<double> &compositional_fields, 
                                    std::vector<double> &fractions) const
    {
      Assert( compositional_fields.size()+1 == fractions.size(), 
              ExcMessage("Size mismatch") );

      //clip the compositional fields so they are between zero and one
      std::vector<double> x_comp = compositional_fields;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);

      //sum the compositional fields for normalization purposes
      double sum_composition = 0.0;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        sum_composition += x_comp[i]; 
      
      if(sum_composition >= 1.0)
        {
          fractions[0] = 0.0;  //background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            fractions[i] = x_comp[i-1]/sum_composition;
        }
      else
        {
          fractions[0] = 1.0 - sum_composition; //background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            fractions[i] = x_comp[i-1];
        }
    }    
          
    

    template <int dim>
    double
    Multicomponent<dim>::
    viscosity (const double temperature,
               const double,
               const std::vector<double> &composition,       /*composition*/
               const SymmetricTensor<2,dim> &,
               const Point<dim> &p) const
    {
      double visc = 0.0;
      std::vector<double> volume_fractions(this->n_compositional_fields() + 1);
      compute_volume_fractions( composition, volume_fractions);

      switch (ViscosityAveraging)
        {
          case Arithmetic:
            {
              for(unsigned int i=0; i< volume_fractions.size(); ++i)
                  visc += volume_fractions[i]*viscosities[i];
              break;
            }
          case Harmonic:
            {
              for(unsigned int i=0; i< volume_fractions.size(); ++i)
                visc += volume_fractions[i]/(viscosities[i]);
              visc = 1.0/visc;
              break;
            }
          case Geometric:
            {
              double geometric_mean = 0.0;
              for(unsigned int i=0; i < volume_fractions.size(); ++i)
                visc += volume_fractions[i]*std::log(viscosities[i]);
              visc = std::exp(visc);
              break;
            }
          case MaximumComposition:
            {
              unsigned int i = (unsigned int)(std::max_element( volume_fractions.begin(), 
                                                              volume_fractions.end() )
                                            - volume_fractions.begin());
              visc = viscosities[i];
              break;
            }
        }
      return visc;
    }


    template <int dim>
    double
    Multicomponent<dim>::
    reference_viscosity () const
    {
      return viscosities[0]; //background
    }

    template <int dim>
    double
    Multicomponent<dim>::
    reference_density () const
    { 
      return densities[0];  //background
    }

    template <int dim>
    double
    Multicomponent<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_expansivities[0]; //background
    }

    template <int dim>
    double
    Multicomponent<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &composition,
                   const Point<dim> &) const
    {
      double cp = 0.0;

      //Arithmetic averaging of specific heats
      std::vector<double> volume_fractions(this->n_compositional_fields() + 1);
      compute_volume_fractions( composition, volume_fractions);
      for(unsigned int i=0; i< volume_fractions.size(); ++i)
        cp += volume_fractions[i]*specific_heats[i];
  
      return cp;
    }

    template <int dim>
    double
    Multicomponent<dim>::
    reference_cp () const
    {
      return specific_heats[0]; //background
    }

    template <int dim>
    double
    Multicomponent<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &composition,
                          const Point<dim> &) const
    {
      double k = 0.0;

      //Arithmetic averaging of thermal conductivities
      //This may not be strictly the most reasonable thing,
      //but for most Earth materials we hope that they do
      //not vary so much that it is a big problem.
      std::vector<double> volume_fractions(this->n_compositional_fields() + 1);
      compute_volume_fractions( composition, volume_fractions);
      for(unsigned int i=0; i< volume_fractions.size(); ++i)
        k += volume_fractions[i]*thermal_conductivities[i];

      return k;     
    }

    template <int dim>
    double
    Multicomponent<dim>::
    reference_thermal_diffusivity () const
    {
      return thermal_conductivities[0] /( densities[0]* specific_heats[0] ); //background
    }

    template <int dim>
    double
    Multicomponent<dim>::
    density (const double temperature,
             const double,
             const std::vector<double> &composition,
             const Point<dim> &) const
    {
      double rho = 0.0;

      //Arithmetic averaging of densities
      std::vector<double> volume_fractions(this->n_compositional_fields() + 1);
      compute_volume_fractions( composition, volume_fractions);
      for(unsigned int i=0; i< volume_fractions.size(); ++i)
      {
        //TODO not strictly correct if thermal expansivities are different, since we are interpreting
        //these compositions as volume fractions, but the error introduced should not be too bad.
        const double temperature_factor= (1.0 - thermal_expansivities[i] * (temperature - reference_T));
        rho += volume_fractions[i]*densities[i]*temperature_factor;
      }
 
      return rho;
    }


    template <int dim>
    double
    Multicomponent<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const std::vector<double> &composition,
                                   const Point<dim> &) const
    {
      double alpha = 0.0;

      //Arithmetic averaging of thermal expansivities
      std::vector<double> volume_fractions(this->n_compositional_fields() + 1);
      compute_volume_fractions( composition, volume_fractions);
      for(unsigned int i=0; i< volume_fractions.size(); ++i)
        alpha += volume_fractions[i]*thermal_expansivities[i];

      return alpha;     
    }


    template <int dim>
    double
    Multicomponent<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0; //incompressible
    }

    template <int dim>
    bool
    Multicomponent<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none))
        return true;
      else
        return false;
    }


    template <int dim>
    bool
    Multicomponent<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if (((dependence & NonlinearDependence::temperature) != NonlinearDependence::none))
        return true;
      else if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none))
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    Multicomponent<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Multicomponent<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none))
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    Multicomponent<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none))
        return true;
      else
        return false;
    }


    template <int dim>
    bool
    Multicomponent<dim>::
    is_compressible () const
    {
      return false;
    }



    template <int dim>
    void
    Multicomponent<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent");
        {
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If one value is given, then all use the same value");
          prm.declare_entry ("Viscosities", "1.e21",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If one value is given, then all use the same value");
          prm.declare_entry ("Thermal expansivities", "4.e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If one value is given, then all use the same value");
          prm.declare_entry ("Specific heats", "1250.",
                             Patterns::List(Patterns::Double(0)),
                             "List of specific heats for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If one value is given, then all use the same value");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If one value is given, then all use the same value");
          prm.declare_entry("Viscosity averaging scheme", "Harmonic",
                             Patterns::Selection("Arithmetic|Harmonic|Geometric|Maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Multicomponent<dim>::parse_parameters (ParameterHandler &prm)
    {
      //not pretty, but we need to get the number of compositional fields before
      //simulatoraccess has been initialized here...
      int n_fields;
      prm.enter_subsection ("Compositional fields");
        {
          n_fields = prm.get_integer ("Number of fields");
        }
      prm.leave_subsection();
      n_fields++; //increment for background

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Multicomponent");
        {
          reference_T                = prm.get_double ("Reference temperature");
 
          if (prm.get ("Viscosity averaging scheme") == "Harmonic")
            ViscosityAveraging = Harmonic; 
          else if (prm.get ("Viscosity averaging scheme") == "Arithmetic")
            ViscosityAveraging = Arithmetic; 
          else if (prm.get ("Viscosity averaging scheme") == "Geometric")
            ViscosityAveraging = Geometric; 
          else if (prm.get ("Viscosity averaging scheme") == "Maximum composition")
            ViscosityAveraging = MaximumComposition; 

          std::vector<double> x_values;

          //Parse densities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Densities")));
          AssertThrow(x_values.size() == 1 || (x_values.size() == n_fields), 
                      ExcMessage("Length of list must be either one, or n_compositional_fields+1"));
          if(x_values.size() == 1) 
            densities.assign( n_fields , x_values[0]);
          else
            densities.assign(x_values.begin(), x_values.end());

          //Parse viscosities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Viscosities")));
          AssertThrow(x_values.size() == 1 || (x_values.size() == n_fields), 
                      ExcMessage("Length of list must be either one, or n_compositional_fields+1"));
          if(x_values.size() == 1) 
            viscosities.assign( n_fields , x_values[0]);
          else
            viscosities.assign(x_values.begin(), x_values.end());

          //Parse thermal conductivities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Thermal conductivities")));
          AssertThrow(x_values.size() == 1 || (x_values.size() == n_fields), 
                      ExcMessage("Length of list must be either one, or n_compositional_fields+1"));
          if(x_values.size() == 1) 
            thermal_conductivities.assign( n_fields , x_values[0]);
          else
            thermal_conductivities.assign(x_values.begin(), x_values.end());
 
          //Parse thermal expansivities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Thermal expansivities")));
          AssertThrow(x_values.size() == 1 || (x_values.size() == n_fields), 
                      ExcMessage("Length of list must be either one, or n_compositional_fields+1"));
          if(x_values.size() == 1) 
            thermal_expansivities.assign( n_fields , x_values[0]);
          else
            thermal_expansivities.assign(x_values.begin(), x_values.end());

          //Parse specific heats
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Specific heats")));
          AssertThrow(x_values.size() == 1 || (x_values.size() == n_fields), 
                      ExcMessage("Length of list must be either one, or n_compositional_fields+1"));
          if(x_values.size() == 1) 
            specific_heats.assign( n_fields , x_values[0]);
          else
            specific_heats.assign(x_values.begin(), x_values.end());
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
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Multicomponent,
                                   "multicomponent",
                                   "'Multicomponent model'.  This model is for use with an arbitrary number of"
                                   "compositional fields, where each field may have completely different material"
                                   "parameters.  Within each field, however, the material behaviour is very simple,"
                                   "with constant material coefficients.  When more than one field is present, the"
                                   "material properties are averaged arithmetically.  An exception is the viscosity,"
                                   "where the averaging should make more of a difference.  For this, the user selects"
                                   "between arithmetic, harmonic, geometric, or maximum composition averaging.")
  }
}
