//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/material_model/table.h>
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
      /**
       * A class that is used to read and and evaluate the pressure and temperature
       * dependent density, thermal expansivity and c_p values.
       **/
      class P_T_LookupFunction
      {
        public:
          /**
           * @brief Constructor
           *
           * @param filename The name of the file in which the values the variable
           * represented by this object are stored.
           **/
          P_T_LookupFunction (const std::string &filename);

          /**
           * @brief Evaluate the table for a given value of pressure
           * and temperature.
           **/
          double value (const double T,
                        const double p) const;

          /**
           * @brief Evaluate the table for the derivative with respect to
           * pressure at a given value of pressure and temperature.
           **/
          double d_by_dp (const double T,
                          const double p) const;
        private:
          /**
           * Number of data points in p and T directions.
           */
          unsigned int n_p, n_T;

          /**
           * Minimal and maximal value for the pressure and temperature
           * for which data exists in the table.
           */
          double min_p, max_p;
          double min_T, max_T;
          /**
           * Step sizes in p and T directions.
           */
          double delta_p, delta_T;

          dealii::Table<2,double> values;
      };

      inline
      P_T_LookupFunction::
      P_T_LookupFunction (const std::string &filename)
      {
        {
          // read in definitions from data file
          std::string temp;
          std::ifstream in("DataDir/tabeldatastruct.txt", std::ios::in);
          AssertThrow (in,
                       ExcMessage (std::string("Couldn't open file <") +
                                   "DataDir/tabeldatastruct.txt>."));

          in >> n_p >> n_T;
          getline(in, temp); // eat remainder of the line

          in >> min_p;
          getline(in, temp); // eat remainder of the line

          in >> min_T;
          getline(in, temp); // eat remainder of the line

          in >> delta_p;
          getline(in, temp); // eat remainder of the line

          in >> delta_T;
          getline(in, temp); // eat remainder of the line

          // now note that in these files pressures are always given in GPa
          // whereas in the rest of the program we use SI (meter-kilogram-seconds)
          // units. so multiply all pressure related quantities by 1e9
          min_p *= 1e9;
          delta_p *= 1e9;

          max_T = min_T + (n_T-1)*delta_T;
          max_p = min_p + (n_p-1)*delta_p;

          values.reinit(n_p, n_T);
        }

        std::ifstream in (filename.c_str(), std::ios::binary);
        AssertThrow (in,
                     ExcMessage (std::string("Couldn't open file <") +
                                 filename + ">."));

        // allocate the following on the heap so as not to bust
        // stack size limits
        double *array = new double[n_p*n_T];
        in.read (reinterpret_cast<char *>(&(array[0])),
                 n_p*n_T*sizeof(double));

        for (unsigned int i=0; i<n_p; ++i)
          for (unsigned int j=0; j<n_T; ++j)
            values[i][j] = array[i*n_T+j];

        delete[] array;
      }


      inline
      double
      P_T_LookupFunction::value (const double T,
                                 const double p) const
      {
// TODO: clamping into the valid range in all cases okay?
        const double pressure = std::max(min_p, std::min(p, max_p-delta_p));

        Assert (pressure >= min_p, ExcMessage ("Not in range"));
        Assert (pressure <= max_p, ExcMessage ("Not in range"));

// TODO: clamping into the valid range in all cases okay?
        const double temperature = std::max(min_T, std::min(T, max_T-delta_T));

        const unsigned int i = (temperature-min_T) / delta_T;
        const unsigned int j = (pressure-min_p) / delta_p;
        Assert (i < n_T-1, ExcInternalError());
        Assert (j < n_p-1, ExcInternalError());

        // compute the coordinates of this point in the
        // reference cell between the data points
        const double xi  = ((temperature-min_T) / delta_T - i);
        const double eta = ((pressure-min_p) / delta_p - j);
        Assert ((0 <= xi) && (xi <= 1), ExcInternalError());
        Assert ((0 <= eta) && (eta <= 1), ExcInternalError());

        // use these co-ordinates for a bilinear interpolation
        return ((1-xi)*(1-eta)*values[i][j] +
                xi    *(1-eta)*values[i+1][j] +
                (1-xi)*eta    *values[i][j+1] +
                xi    *eta    *values[i+1][j+1]);
      }


      inline
      double
      P_T_LookupFunction::d_by_dp (const double T,
                                   const double p) const
      {
// TODO: clamping into the valid range in all cases okay?
        const double pressure = std::max(min_p, std::min(p, max_p-delta_p));

        Assert (pressure >= min_p, ExcMessage ("Not in range"));
        Assert (pressure <= max_p, ExcMessage ("Not in range"));
        Assert (T >= min_T, ExcMessage ("Not in range"));
        Assert (T <= max_T, ExcMessage ("Not in range"));

        const unsigned int i = (T-min_T) / delta_T;
        const unsigned int j = (pressure-min_p) / delta_p;
        Assert (i < n_T-1, ExcInternalError());
        Assert (j < n_p-1, ExcInternalError());

        // compute the coordinates of this point in the
        // reference cell between the data points
        //
        // since the derivative in p-direction (eta direction)
        // is constant for the bilinear interpolation, we really
        // only need xi here
        const double xi  = ((T-min_T) / delta_T - i);
        Assert ((0 <= xi) && (xi <= 1), ExcInternalError());

        // use these co-ordinates for a bilinear interpolation
        return ((1-xi)*(values[i][j+1] - values[i][j]) +
                xi    *(values[i+1][j+1] - values[i+1][j])) / delta_p;
      }
    }



    template <int dim>
    double
    Table<dim>::
    viscosity (const double temperature,
               const double pressure,
               const Point<dim> &position) const
    {
      const double reference_eta    = 5e24;
      return reference_eta;
    }



    template <int dim>
    double
    Table<dim>::
    reference_viscosity () const
    {
      const double reference_eta    = 5e24;
      return reference_eta;
    }



    template <int dim>
    double
    Table<dim>::
    specific_heat (const double temperature,
                   const double pressure,
                   const Point<dim> &) const
    {
//    const double reference_specific_heat = 1250;    /* J / K / kg */  //??
//      if (!IsCompressible) return reference_specific_heat; TODO
      static internal::P_T_LookupFunction cp("DataDir/cp_bin");
      return cp.value(temperature, pressure);
    }



    template <int dim>
    double
    Table<dim>::
    thermal_conductivity (const double,
                          const double,
                          const Point<dim> &) const
    {
      // this model assumes that the thermal conductivity is in fact constant
      return 4.7;
    }



    template <int dim>
    double
    Table<dim>::
    density (const double temperature,
             const double pressure,
             const Point<dim> &position) const
    {
      static internal::P_T_LookupFunction rho("DataDir/rho_bin");
      return rho.value(temperature, pressure);
    }


    template <int dim>
    double
    Table<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const Point<dim> &position) const
    {
      static internal::P_T_LookupFunction rho("DataDir/rho_bin");
      return rho.d_by_dp(temperature, pressure) / rho.value(temperature,pressure);
    }


    template <int dim>
    bool
    Table<dim>::
    is_compressible () const
    {
      return true;
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    template class Table<deal_II_dimension>;

    ASPECT_REGISTER_MATERIAL_MODEL(Table,
                                   "table",
                                   "A material model that reads tables of pressure and temperature "
                                   "dependent material coefficients from files.")
  }
}
