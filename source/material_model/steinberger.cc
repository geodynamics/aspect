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



      class thermal_exp_lookup
      {
        public:
          thermal_exp_lookup(const char *filename)
          {
            std::string temp;
            std::ifstream in(filename, std::ios::in);
            AssertThrow (in,
                         ExcMessage (std::string("Couldn't open file <") + filename));

            getline(in, temp); // eat first line

            delta_depth=1000;
            min_depth=-1;

            while (!in.eof())
              {
                double exps_real, exps_simple, depth;
                in >> exps_real;
                if (in.eof())
                  break;
                in >> exps_simple >> depth;
                getline(in, temp);
                depth*=-1000.0; //[km] to [m] and positive sign
                if (min_depth<0)
                  min_depth= depth;
                max_depth=depth;
                values.push_back(exps_simple);
              }


          }

          double thermal_expansivity(double depth)
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
          double delta_depth;
          double min_depth;
          double max_depth;

      };

      class lateral_viscosity_lookup
          {
            public:
              lateral_viscosity_lookup(const char *filename)
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
          radial_viscosity_lookup(const char *filename)
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
      static internal::radial_viscosity_lookup table("data/material-model/steinberger/radial_visc.txt");
      static internal::lateral_viscosity_lookup lat_table("data/material-model/steinberger/temp_viscosity_prefactor.txt");

      const double depth = this->get_geometry_model().depth(position);
      const unsigned int idx = 100 * depth / this->get_geometry_model().maximal_depth();
      const double delta_temp = temperature-avg_temp[idx];
      const double adia_temp = this->get_adiabatic_conditions().temperature(position);

      const double vis_lateral_exp = -1.0*lat_table.lateral_viscosity(depth)*delta_temp/(temperature*adia_temp);

      const double vis_lateral = std::max(std::min(std::exp(vis_lateral_exp),1e2),1e-2);
      const double vis_radial = table.radial_viscosity(depth);

      return std::max(std::min(vis_lateral * vis_radial,1e23),1e19);
    }

    template <int dim>
    double
    Steinberger<dim>::
    reference_viscosity () const
    {
      const double reference_eta    = 5e24;
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
      Assert (false, ExcNotImplemented());
      return 0;
    }


    template <int dim>
    double
    Steinberger<dim>::
    specific_heat (const double temperature,
                   const double pressure,
                   const Point<dim> &p) const
    {
      return 1200;

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
      static internal::thermal_exp_lookup table("data/material-model/steinberger/thermal_expansivity.d");

      const double depth = this->get_geometry_model().depth(position);

      const double reference_density = 3300;
      const double thermal_expansion_coefficient = table.thermal_expansivity(depth);
      const double reference_temperature = 293;

      return (reference_density *
              (1 - thermal_expansion_coefficient * (temperature -
                                                    reference_temperature)));
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
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Steinberger,
                                   "Steinberger",
                                   "lookup from the paper of Steinberger/Calderwood")
  }
}
