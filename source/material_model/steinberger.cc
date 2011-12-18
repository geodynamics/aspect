//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/material_model/steinberger.h>
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
            const unsigned int idx = (depth-min_depth)/delta_depth;
            Assert(idx<values.size(), ExcMessage("not in range"));
            return values[idx];
          }

        private:
          std::vector<double> values;
          double delta_depth;
          double min_depth;
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
            const unsigned int idx = (depth-min_depth)/delta_depth;
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
    double
    Steinberger<dim>::
    viscosity (const double temperature,
               const double pressure,
               const Point<dim> &position) const
    {
      static internal::radial_viscosity_lookup table("data/material-model/steinberger/radial_visc.txt");
      const double vis_lateral = 1.0; // TODO

      const double R1=6371000.0; //TODO
      const double depth = R1-std::sqrt(position.square());
      const double vis_radial = table.radial_viscosity(depth);

      return std::max(1e10,vis_lateral * vis_radial);
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


      const double R1=6371000.0; //TODO
      const double depth = R1-std::sqrt(position.square());

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
    template class Steinberger<deal_II_dimension>;

    ASPECT_REGISTER_MATERIAL_MODEL(Steinberger,
                                   "Steinberger",
                                   "lookup from the paper of Steinberger/Calderwood")
  }
}
