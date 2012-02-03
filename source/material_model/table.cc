//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/material_model/table.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>
#include <cstring>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {


    namespace internal
    {
        /**
         * A class that is used to read and and evaluate the pressure and temperature
         * dependent Phase.
        **/
        template <int dim>
        class PhaseLookupFunction
          {
            public:
              /**
               * @brief Constructor
               *
               * @param filename The name of the file in which the values the variable
               * represented by this object are stored.
               **/

              PhaseLookupFunction<dim>(const std::string &filename);

              /**
               * @brief Evaluate the table for a given value of pressure
               * and temperature.
               **/
              double value (const double T,
                            const double p) const;

            private:

             int  nPhaseFields;
             int  nIsotherms;

             std::list<std::string> PhasefieldLabels;

          };
          template <int dim>
          inline
          PhaseLookupFunction<dim>::
          PhaseLookupFunction (const std::string &filename)
          {
            {
              // read in definitions from data file
              std::string temp;
              std::string path = filename;
              int ipos = filename.find("table");
              ipos = filename.find("/", ipos+2);
              ipos = filename.find("/", ipos+2);
              path.replace(ipos+1,19,".Phases.lab"); // TODO remove this ugly hack with something more elegant!
              std::ifstream in(path.c_str(), std::ios::in);
              AssertThrow (in,
                           ExcMessage (std::string("Couldn't open file <") +
                                       path));

              /* skip the first four lines*/
              getline(in, temp); // eat remainder of the line
              getline(in, temp); // eat remainder of the line
              getline(in, temp); // eat remainder of the line
              getline(in, temp); // eat remainder of the line

              in >> nPhaseFields >> nIsotherms;
              getline(in, temp); // eat remainder of the line

              int i, j;

              /* read phase labels */
              for (i=0;i<nPhaseFields; i++){
                 in >> j  >> PhasefieldLabels.push_back;
                 getline(in, temp); // eat remainder of the line
              }
              /* skip line*/
              getline(in, temp); // eat remainder of the line

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

            for (unsigned int j=0; j<n_T; ++j)
              for (unsigned int i=0; i<n_p; ++i)
                values[i][j] = array[j*n_p+i];

            delete[] array;
       }
      /**
       * A class that is used to read and and evaluate the pressure and temperature
       * dependent Phase.
      **/
      template <int dim>
      class PhaseLookupFunction
      {
        public:
          /**
           * @brief Constructor
           *
           * @param filename The name of the file in which the values the variable
           * represented by this object are stored.
           **/

          PhaseLookupFunction<dim>(const std::string &filename);

          /**
           * @brief Evaluate the table for a given value of pressure
           * and temperature.
           **/
          int value (const double T,
                     const double p) const;

        private:

          unsigned int  nPhaseFields;
          unsigned int  nIsotherms;
          unsigned int  nTransMax;

          std::list<std::string> PhasefieldLabels;

          dealii::Table<2,double> PressureTransitions;
          dealii::Table<2,int>    Phase;
          dealii::Table<1,double> IsothermValue;
          dealii::Table<1,double> nTransIsotherm;

          double IsothermMin;
          double IsothermMax;
          double IsothermdT;

      };
      template <int dim>
      inline
      PhaseLookupFunction<dim>::
      PhaseLookupFunction (const std::string &filename)
      {
        // read in definitions from data file
        std::string temp;
        std::ifstream in(filename.c_str(), std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Couldn't open file <") +
                                 filename));

        /* skip the first four lines*/
        getline(in, temp); // eat remainder of the line
        getline(in, temp); // eat remainder of the line
        getline(in, temp); // eat remainder of the line
        getline(in, temp); // eat remainder of the line

        in >> nPhaseFields;
        getline(in, temp); // eat remainder of the line
        in >> nIsotherms;
        getline(in, temp); // eat remainder of the line

        unsigned int i, j;

        /* read phase labels */
        for (i=0; i<nPhaseFields; i++)
          {
            in >> j  >> temp;
            PhasefieldLabels.push_back(temp);
            getline(in, temp); // eat remainder of the line
          }
        /* skip line*/
        getline(in, temp); // eat remainder of the line

        nTransMax = 20;

        PressureTransitions.reinit(nTransMax, nIsotherms);
        Phase.reinit(nTransMax , nIsotherms);
        IsothermValue.reinit(nIsotherms);
        nTransIsotherm.reinit(nIsotherms);


        /* read table of phase changes as function of pressure per isotherm */
        for (i=0; i<nIsotherms; i++)
          {
            in >>  IsothermValue[i] >> nTransIsotherm[i];
            getline(in, temp); // eat remainder of the line
            Assert (nTransIsotherm[i] < nTransMax, ExcInternalError());

            for (j=0; j<nTransIsotherm[i]; j++)
              {
                in >> PressureTransitions[j][i]  >> Phase[j][i];
              }
          }
        IsothermMin = IsothermValue[0];
        IsothermMax = IsothermValue[nIsotherms-1];
        IsothermdT  = (IsothermMax-IsothermMin)/(nIsotherms-1) ;

      }

      template <int dim>
      inline
      int
      PhaseLookupFunction<dim>::value (const double T,
                                       const double p) const
      {
        unsigned int Isotherm;
        double P = p/1e9;

        Isotherm = (int) ( (T - IsothermMin)/IsothermdT + 0.5);
        Assert (Isotherm > 0, ExcInternalError());
        Assert (Isotherm < nIsotherms, ExcInternalError());

        for (int i=1; i<nTransIsotherm[Isotherm]; i++)
          {
            if (P > PressureTransitions[i-1][Isotherm] && P < PressureTransitions[i][Isotherm])
              {
                return Phase[i-1][Isotherm];
              }
          }
        return 0;
      }
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
          std::string path = filename;
          int ipos = filename.find("table");
          ipos = filename.find("/", ipos+2);
          ipos = filename.find("/", ipos+2);
          path.replace(ipos+1,19,"tabledatastruct.txt"); // TODO remove this ugly hack with something more elegant!
          std::ifstream in(path.c_str(), std::ios::in);
          AssertThrow (in,
                       ExcMessage (std::string("Couldn't open file <") +
                                   path));

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

        for (unsigned int j=0; j<n_T; ++j)
          for (unsigned int i=0; i<n_p; ++i)
            values[i][j] = array[j*n_p+i];

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

        const unsigned int i = static_cast<unsigned int>((pressure-min_p) / delta_p);
        const unsigned int j = static_cast<unsigned int>((temperature-min_T) / delta_T);
        Assert (i < n_p-1, ExcInternalError());
        Assert (j < n_T-1, ExcInternalError());

        // compute the coordinates of this point in the
        // reference cell between the data points
        const double xi  = ((temperature-min_T) / delta_T - j);
        const double eta = ((pressure-min_p) / delta_p - i);
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

        const unsigned int i = static_cast<unsigned int>((pressure-min_p) / delta_p);
        const unsigned int j = static_cast<unsigned int>((T-min_T) / delta_T);
        Assert (i < n_p-1, ExcInternalError());
        Assert (j < n_T-1, ExcInternalError());

        // compute the coordinates of this point in the
        // reference cell between the data points
        //
        // since the derivative in p-direction (eta direction)
        // is constant for the bilinear interpolation, we really
        // only need xi here
        const double xi  = ((T-min_T) / delta_T - j);
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
      double viscosity;
      if (!strcmp(ViscosityModel.c_str(),"Exponential"))
        {
          const double R0=  3591e3; //TODO
          const double R1=  6591e3; //TODO
          const double T1=  375; //TODO
          const double dT=  3498; //TODO
          const double depth = (1e0 - (position.norm()-R0)/(R1-R0));
          const double T = (temperature-T1/dT);
          viscosity = reference_eta*std::exp(- std::log(ExponentialT)*T +
                               std::log(ExponentialP)*depth);
        }
      else
        {
          viscosity = reference_eta;
        }
      return viscosity;
    }



    template <int dim>
    double
    Table<dim>::
    reference_viscosity () const
    {
      return reference_eta;
    }

    template <int dim>
    double
    Table<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    Table<dim>::
    reference_thermal_alpha () const
    {
      return reference_alpha;
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
      std::string path= data_directory;
      path +="cp_bin";
      static internal::P_T_LookupFunction cp(path);
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
      return k_value;
    }


    template <int dim>
    double
    Table<dim>::
    reference_thermal_diffusivity () const
    {
//TODO: the right thing to do would be to compute this value as
// kappa=k/rho/c_P with reference values for the three quantities
// on the right hand side.
      return reference_kappa;
    }


    template <int dim>
    double
    Table<dim>::
    density (const double temperature,
             const double pressure,
             const Point<dim> &position) const
    {
      static internal::P_T_LookupFunction rho(data_directory+"rho_bin");
      return rho.value(temperature, pressure);
    }


    template <int dim>
    double
    Table<dim>::
    seismic_Vp (const double temperature,
                const double pressure) const
    {
      static internal::P_T_LookupFunction vp(data_directory+"vseis_p_bin");
      return vp.value(temperature, pressure);
    }


    template <int dim>
    double
    Table<dim>::
    seismic_Vs (const double temperature,
                const double pressure) const
    {
      static internal::P_T_LookupFunction vs(data_directory+"vseis_s_bin");
      return vs.value(temperature, pressure);
    }


    template <int dim>
    unsigned int
    Table<dim>::
    thermodynamic_phase (const double temperature,
                         const double pressure) const
    {
      if (!ComputePhases) return 0;
      static internal::PhaseLookupFunction<dim> phase(data_directory+"Phases.lab");
      return phase.value(temperature, pressure);
    }


    template <int dim>
    double
    Table<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const Point<dim> &position) const
    {
      static internal::P_T_LookupFunction rho(data_directory+"rho_bin");
      return rho.d_by_dp(temperature, pressure) / rho.value(temperature,pressure);
    }


    template <int dim>
    bool
    Table<dim>::
    is_compressible () const
    {
      return true;
    }


    template <int dim>
    void
    Table<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Table model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Gravity", "30",
                             Patterns::Double (0),
                             "The value of the Gravity$. "
                             "Units: $m/s^2$.");
          prm.declare_entry ("Composition", "olixene",
                             Patterns::Anything (),
                             "The Composition of the model. ");
          prm.declare_entry ("Path to model data", "datadir",
                             Patterns::DirectoryName (),
                             "The path to the model data. ");
          prm.declare_entry ("ComputePhases", "false",
                             Patterns::Bool (),
                             "whether to compute phases. ");
          prm.enter_subsection ("Viscosity");
          {
            prm.declare_entry ("ViscosityModel", "Exponential",
                               Patterns::Anything (),
                               "Viscosity Model");
            prm.declare_entry ("ReferenceViscosity", "5e24",
                               Patterns::Double (0),
                               "The value of the constant viscosity. Units: $kg/m/s$.");
            prm.declare_entry ("ExponentialT", "1",
                               Patterns::Double (0),
                               "multiplication factor or Temperature exponent");
            prm.declare_entry ("ExponentialP", "1",
                               Patterns::Double (0),
                               "multiplication factor or Pressure exponent");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Table<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Table model");
        {
          reference_rho     = prm.get_double ("Reference density");
          reference_T       = prm.get_double ("Reference temperature");
          k_value           = prm.get_double ("Thermal conductivity");
          reference_alpha   = prm.get_double ("Thermal expansion coefficient");
          composition       = prm.get ("Composition");
          data_directory    = prm.get ("Path to model data") + "/" + composition +"/";
          ComputePhases     = prm.get_bool ("ComputePhases");
          prm.enter_subsection ("Viscosity");
          {
            ViscosityModel = prm.get ("ViscosityModel");
            reference_eta  = prm.get_double ("ReferenceViscosity");
            ExponentialT   = prm.get_double ("ExponentialT");
            ExponentialP   = prm.get_double ("ExponentialP");
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
  namespace MaterialModel
  {
    template class Table<deal_II_dimension>;

    ASPECT_REGISTER_MATERIAL_MODEL(Table,
                                   "table",
                                   "A material model that reads tables of pressure and temperature "
                                   "dependent material coefficients from files.")
  }
}
