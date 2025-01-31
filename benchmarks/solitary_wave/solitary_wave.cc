/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/
#include <aspect/melt.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/simulator.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>


/**
 * MPI operator used by MPI_max_and_data
 */
void myop_func(void *invec,
               void *inoutvec,
               int *len,
               MPI_Datatype *datatype)
{
  AssertThrow(*len > 1, dealii::ExcNotImplemented());
  AssertThrow(*datatype == MPI_DOUBLE, dealii::ExcNotImplemented());

  double *indata = static_cast<double *>(invec);
  double *inoutdata = static_cast<double *>(inoutvec);

  if (indata[0]>inoutdata[0])
    {
      for (int i=0; i<*len; ++i)
        inoutdata[i] = indata[i];
    }
}


/**
 * Computes MPI_MAX of @p local_max like Allreduce, but also transmits the
 * @p local_data from the rank with the largest @local_max to every rank
 * (returned in @p global_data).
 */
void MPI_max_and_data(const double &local_max,
                      const double &local_data,
                      double &global_max,
                      double &global_data)
{
  MPI_Op myop;
  MPI_Op_create(&myop_func, /* commutes? */ 1 ,&myop);

  double local[] = {local_max, local_data};
  double global[2];

  MPI_Allreduce(local, global, 2, MPI_DOUBLE, myop, MPI_COMM_WORLD);
  global_max = global[0];
  global_data = global[1];
  MPI_Op_free(&myop);
}


namespace aspect
{
  /**
   * This is the "Solitary wave" benchmark defined in the following paper:
   * @code
   *  @Article{DMGT11,
   *    author =       {T. Keller and D. A. May and B. J. P. Kaus},
   *    title =        {Numerical modelling of magma dynamics coupled
   *                    to tectonic deformation of lithosphere and crust},
   *    journal =      {Geophysical Journal International},
   *    year =         2013,
   *    volume =       195(3),
   *    pages =        {1406-1442}
   * @endcode
   *
   * To calculate the initial condition, which is a solitary wave solution of
   * the magma dynamics equations, we use the equation for the one-dinemsional
   * case and the non-dimensionalization as it is described in
   * @code
   *  @Article{SS11,
   *    author =       {G. Simpson and M. Spiegelman},
   *    title =        {Solitary Wave Benchmarks in Magma Dynamics},
   *    journal =      {Journal of Scientific Computing},
   *    year =         2011,
   *    volume =       49(3),
   *    pages =        {268-290}
   * @endcode
   *
   * Specifically, this means that we scale the porosity with the background
   * porosity, and the coordinates with the compaction length $\delta_0$, which is
   * defined as $\sqrt \frac{k(\phi_0) \xi^{*}+4/3 \eta^{*}}{\eta_f}$.  $k(\phi_0)$ is the
   * permeability at background porosity, $\xi^{*}$ is the compaction viscosity,
   * $\eta^{*}$ is the shear viscosity of the fluid and $\eta_f$ is the shear viscosity
   * of the melt.
   */
  namespace SolitaryWaveBenchmark
  {
    namespace AnalyticSolutions
    {
      // vectors to store the porosity field and the corresponding coordinate in
      const unsigned int max_points = 2e8;
      std::vector<double> porosity, coordinate;

      /**
       * @note The solitary wave solution only exists as a function x = func(phi)
       * and not phi = func(x), which is what we would like to have for describing
       * the shape of the wave. Thus, we calculate x = func(phi) for a range of phis
       * between the background porosity and the amplitude of the wave. In a next
       * step, we interpolate these values to the grid.
       *
       * @param phi The characteristic shape of the wave, with phi --> 1
       * for x --> +- infinity
       *
       * @param amplitude The amplitude of the solitary wave, which is always
       * greater than 1.
       */
      double solitary_wave_solution (const double phi, const double amplitude)
      {
        AssertThrow(phi > 1.0 && phi <= amplitude,
                    ExcMessage("The solitary wave solution can only be computed "
                               "for porosities larger than the background porosity of 1 "
                               "and smaller than or equal to the amplitude of the wave."));
        AssertThrow(amplitude > 1,
                    ExcMessage("Amplitude of the solitary wave must be larger than 1!"));
        const double A_1   = std::sqrt(amplitude - 1.0);
        const double A_phi = std::sqrt(amplitude - phi);
        return std::sqrt(amplitude + 0.5)
               * (2 * A_phi - 1.0/A_1 * std::log((A_1 - A_phi)/(A_1 + A_phi)));
      }


      /**
       * This function reads the coordinate and the porosity of the solitary wave
       * from an input file.
       *
       * @param filename Name of the input file.
       */
      void read_solitary_wave_solution (const std::string &filename,
                                        const MPI_Comm comm)
      {
        std::string temp;
        std::stringstream in(Utilities::read_and_distribute_file_content(filename, comm));

        while (!in.eof())
          {
            double x, f;
            in >> x >> f;
            if (in.eof())
              break;
            getline(in, temp);

            coordinate.insert(coordinate.begin(),x);
            porosity.insert(porosity.begin(),f);
          }
      }


      /**
       * This function gets the coordinate as an input parameters and gives
       * back the porosity of the solitary wave. As this function is only defined
       * implicitly, we have to interpolate from the coordinates where we have the
       * porosity to our mesh.
       *
       * @param amplitude The amplitude of the solitary wave, which is always
       * greater than 1.
       * @param offset The offset of the center of the solitary wave from the
       * boundary of the domain.
       */
      void compute_porosity (const double amplitude,
                             const double background_porosity,
                             const double /*offset*/,
                             const double compaction_length,
                             const bool read_solution,
                             const std::string file_name,
                             const MPI_Comm comm)
      {
        // non-dimensionalize the amplitude
        const double non_dim_amplitude = amplitude / background_porosity;

        if (read_solution)
          read_solitary_wave_solution(file_name, comm);
        else
          {
            porosity.resize(max_points);
            coordinate.resize(max_points);

            // get the coordinates where we have the solution
            for (unsigned int i=0; i<max_points; ++i)
              {
                porosity[i] = 1.0 + 1e-10*non_dim_amplitude
                              + double(i)/double(max_points-1) * (non_dim_amplitude * (1.0 - 1e-10) - 1.0);
                coordinate[i] = solitary_wave_solution(porosity[i], non_dim_amplitude);
              }
          }

        for (unsigned int i=0; i<coordinate.size(); ++i)
          {
            // re-scale porosity and position
            porosity[i] *= background_porosity;
            coordinate[i] *= compaction_length;
          }
      }


      double interpolate (const double position,
                          const double offset)
      {
        // interpolate from the solution grid to the mesh used in the simulation
        // solitary wave is a monotonically decreasing function, so the coordinates
        // should be in descending order

        // we only have the solution of the solitary wave for
        // coordinates larger than 0 (one half of the wave)
        const double x = (position > offset
                          ?
                          position - offset
                          :
                          offset - position);

        if (x > coordinate[0])
          return porosity[0];

        unsigned int j= coordinate.size()-2;
        unsigned int i = j/2;
        while (!(x < coordinate[j] && x >= coordinate[j+1]))
          {
            if (x < coordinate[j])
              j += i;
            else
              j -= i;
            if (i>1)
              i /= 2;
          }

        const double distance = (x - coordinate[j+1])
                                /(coordinate[j] - coordinate[j+1]);
        return porosity[j+1] + distance * (porosity[j] - porosity[j+1]);
      }

      /**
       * The exact solution for the Solitary wave benchmark.
       */
      template <int dim>
      class FunctionSolitaryWave : public Function<dim>
      {
        public:
          FunctionSolitaryWave (const double offset, const double delta, const std::vector<double> &initial_pressure, const double max_z, const unsigned int n_components)
            :
            Function<dim>(n_components),
            offset_(offset),
            delta_(delta),
            initial_pressure_(initial_pressure),
            max_z_(max_z)
          {}

          void set_delta(const double delta)
          {
            delta_ = delta;
          }

          void vector_value (const Point<dim> &p,
                             Vector<double>   &values) const override
          {
            unsigned int index = static_cast<int>((p[dim-1]-delta_)/max_z_ * (initial_pressure_.size()-1));
            if (p[dim-1]-delta_ < 0)
              index = 0;
            else if (p[dim-1]-delta_ > max_z_)
              index = initial_pressure_.size()-1;
            AssertThrow(index < initial_pressure_.size(), ExcMessage("not in range"));
            const double z_coordinate1 = static_cast<double>(index)/static_cast<double>(initial_pressure_.size()-1) * max_z_;
            const double z_coordinate2 = static_cast<double>(index+1)/static_cast<double>(initial_pressure_.size()-1) * max_z_;
            const double interpolated_pressure = (index == initial_pressure_.size()-1)
                                                 ?
                                                 initial_pressure_[index]
                                                 :
                                                 initial_pressure_[index] + (initial_pressure_[index+1] - initial_pressure_[index])
                                                 * (p[dim-1]-delta_ - z_coordinate1) / (z_coordinate2 - z_coordinate1);

            values[dim+2+dim+2] = AnalyticSolutions::interpolate(p[dim-1]-delta_,offset_); //porosity
            values[dim+1] = interpolated_pressure;                                   //compaction pressure
          }

        private:
          const double offset_;
          double delta_;
          const std::vector<double> initial_pressure_;
          const double max_z_;
      };
    }


    /**
     * An initial conditions model for the solitary waves benchmark.
     */
    template <int dim>
    class SolitaryWaveInitialCondition : public InitialComposition::Interface<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Initialization function. Take references to the material model and
         * get the compaction length, so that it can be used subsequently to
         * compute the analytical solution for the shape of the solitary wave.
         */
        void
        initialize () override;

        /**
         * Return the boundary velocity as a function of position.
         */
        double
        initial_composition (const Point<dim> &position, const unsigned int n_comp) const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

        double
        get_amplitude () const;

        double
        get_background_porosity () const;

        double
        get_offset () const;

      private:
        double amplitude;
        double background_porosity;
        double offset;
        double compaction_length;
        bool read_solution;
        std::string file_name;
    };


    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class SolitaryWaveMaterial : public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        bool is_compressible () const override
        {
          return false;
        }

        double reference_darcy_coefficient () const override
        {
          // Make sure we keep track of the initial composition manager and
          // that it continues to live beyond the time when the simulator
          // class releases its pointer to it.
          if (initial_composition_manager == nullptr)
            const_cast<std::shared_ptr<const aspect::InitialComposition::Manager<dim>>&>(initial_composition_manager)
              = this->get_initial_composition_manager_pointer();

          // Note that this number is based on the background porosity in the
          // solitary wave initial condition.
          const SolitaryWaveInitialCondition<dim> &initial_composition =
            initial_composition_manager->template
            get_matching_active_plugin<SolitaryWaveInitialCondition<dim>>();

          return reference_permeability * std::pow(initial_composition.get_background_porosity(), 3.0) / eta_f;

        }

        double length_scaling (const double porosity) const
        {
          return std::sqrt(reference_permeability * Utilities::fixed_power<3>(porosity) * (xi_0 + 4.0/3.0 * eta_0) / eta_f);
        }

        double velocity_scaling (const double porosity) const
        {
          const Point<dim> surface_point = this->get_geometry_model().representative_point(0.0);
          return reference_permeability * Utilities::fixed_power<2>(porosity) * (reference_rho_s - reference_rho_f)
                 * this->get_gravity_model().gravity_vector(surface_point).norm() / eta_f;
        }

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                      typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const override
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              double porosity = in.composition[i][porosity_idx];

              out.viscosities[i] = eta_0 * (1.0 - porosity);
              out.densities[i] = reference_rho_s;
              out.thermal_expansion_coefficients[i] = 0.0;
              out.specific_heat[i] = 1.0;
              out.thermal_conductivities[i] = 0.0;
              out.compressibilities[i] = 0.0;
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                out.reaction_terms[i][c] = 0.0;
            }

          // fill melt outputs if they exist
          aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim>>();

          if (melt_out != nullptr)
            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                double porosity = in.composition[i][porosity_idx];

                melt_out->compaction_viscosities[i] = xi_0 * (1.0 - porosity);
                melt_out->fluid_viscosities[i]= eta_f;
                melt_out->permeabilities[i]= reference_permeability * Utilities::fixed_power<3>(porosity);
                melt_out->fluid_densities[i]= reference_rho_f;
                melt_out->fluid_density_gradients[i] = 0.0;
              }

        }

      private:
        double reference_rho_s;
        double reference_rho_f;
        double eta_0;
        double xi_0;
        double eta_f;
        double reference_permeability;

        /**
         * A shared pointer to the initial composition object
         * that ensures that the current object can continue
         * to access the initial composition object beyond the
         * first time step.
         */
        std::shared_ptr<const aspect::InitialComposition::Manager<dim>> initial_composition_manager;
    };

    template <int dim>
    void
    SolitaryWaveMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Solitary wave");
        {
          prm.declare_entry ("Reference solid density", "3000",
                             Patterns::Double (0),
                             "Reference density of the solid $\\rho_{s,0}$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference melt density", "2500",
                             Patterns::Double (0),
                             "Reference density of the melt/fluid$\\rho_{f,0}$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference shear viscosity", "1e20",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference compaction viscosity", "1e20",
                             Patterns::Double (0),
                             "The value of the constant volumetric viscosity $\\xi_0$ of the solid matrix. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference melt viscosity", "100.0",
                             Patterns::Double (0),
                             "The value of the constant melt viscosity $\\eta_f$. "
                             "Units: \\si{\\pascal\\second}.");

          prm.declare_entry ("Reference permeability", "5e-9",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: \\si{\\meter\\squared}.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SolitaryWaveMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Solitary wave");
        {
          reference_rho_s            = prm.get_double ("Reference solid density");
          reference_rho_f            = prm.get_double ("Reference melt density");
          eta_0                      = prm.get_double ("Reference shear viscosity");
          xi_0                       = prm.get_double ("Reference compaction viscosity");
          eta_f                      = prm.get_double ("Reference melt viscosity");
          reference_permeability     = prm.get_double ("Reference permeability");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    double
    SolitaryWaveInitialCondition<dim>::get_amplitude () const
    {
      return amplitude;
    }

    template <int dim>
    double
    SolitaryWaveInitialCondition<dim>::get_background_porosity () const
    {
      return background_porosity;
    }

    template <int dim>
    double
    SolitaryWaveInitialCondition<dim>::get_offset () const
    {
      return offset;
    }

    template <int dim>
    void
    SolitaryWaveInitialCondition<dim>::initialize ()
    {
      std::cout << "Initialize solitary wave solution"
                << std::endl;

      AssertThrow(Plugins::plugin_type_matches<const SolitaryWaveMaterial<dim>>(this->get_material_model()),
                  ExcMessage("Initial condition Solitary Wave only works with the material model Solitary wave."));

      const SolitaryWaveMaterial<dim> &
      material_model
        = Plugins::get_plugin_as_type<const SolitaryWaveMaterial<dim>>(this->get_material_model());

      compaction_length = material_model.length_scaling(background_porosity);

      AnalyticSolutions::compute_porosity(amplitude,
                                          background_porosity,
                                          offset,
                                          compaction_length,
                                          read_solution,
                                          file_name,
                                          this->get_mpi_communicator());
    }


    template <int dim>
    double
    SolitaryWaveInitialCondition<dim>::
    initial_composition (const Point<dim> &position, const unsigned int /*n_comp*/) const
    {
      return AnalyticSolutions::interpolate(position[dim-1],
                                            offset);
    }


    template <int dim>
    void
    SolitaryWaveInitialCondition<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Solitary wave initial condition");
        {
          prm.declare_entry ("Amplitude", "0.01",
                             Patterns::Double (0),
                             "Amplitude of the solitary wave. Units: none.");
          prm.declare_entry ("Background porosity", "0.001",
                             Patterns::Double (0),
                             "Background porosity of the solitary wave. Units: none.");
          prm.declare_entry ("Offset", "150",
                             Patterns::Double (0),
                             "Offset of the center of the solitary wave from the boundary"
                             "of the domain. "
                             "Units: \\si{\\meter}.");
          prm.declare_entry ("Read solution from file", "false",
                             Patterns::Bool (),
                             "Whether to read the porosity initial condition from "
                             "a file or to compute it.");
          prm.declare_entry ("File name", "solitary_wave.txt",
                             Patterns::Anything (),
                             "The file name of the porosity initial condition data. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SolitaryWaveInitialCondition<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Solitary wave initial condition");
        {
          amplitude            = prm.get_double ("Amplitude");
          background_porosity  = prm.get_double ("Background porosity");
          offset               = prm.get_double ("Offset");
          read_solution        = prm.get_bool ("Read solution from file");
          file_name            = prm.get ("File name");

          AssertThrow(amplitude > background_porosity,
                      ExcMessage("Amplitude of the solitary wave must be larger "
                                 "than the background porosity."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    /**
      * A postprocessor that evaluates the accuracy of the solution.
      *
      * The implementation of error evaluators that correspond to the
      * benchmarks defined in the paper Keller et al. reference above.
      */
    template <int dim>
    class SolitaryWavePostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Initialization function. Take references to the material model and
         * initial conditions model to get the parameters necessary for computing
         * the analytical solution for the shape of the solitary wave and store them.
         */
        void
        initialize () override;

        void
        store_initial_pressure ();

        double
        compute_phase_shift ();

      private:
        double amplitude;
        double background_porosity;
        double offset;
        double compaction_length;
        double velocity_scaling;
        double boundary_velocity;
        unsigned int max_points;
        std::vector<double> initial_pressure;
        double maximum_pressure;
        std::shared_ptr<AnalyticSolutions::FunctionSolitaryWave<dim>> ref_func;

    };

    template <int dim>
    void
    SolitaryWavePostprocessor<dim>::initialize ()
    {
      // verify that we are using the "Solitary wave" initial conditions and material model,
      // then get the parameters we need

      const SolitaryWaveInitialCondition<dim> &initial_composition
        = this->get_initial_composition_manager().template get_matching_active_plugin<SolitaryWaveInitialCondition<dim>> ();

      amplitude           = initial_composition.get_amplitude();
      background_porosity = initial_composition.get_background_porosity();
      offset              = initial_composition.get_offset();

      AssertThrow(Plugins::plugin_type_matches<const SolitaryWaveMaterial<dim>>(this->get_material_model()),
                  ExcMessage("Postprocessor Solitary Wave only works with the material model Solitary wave."));

      const SolitaryWaveMaterial<dim> &material_model
        = Plugins::get_plugin_as_type<const SolitaryWaveMaterial<dim>>(this->get_material_model());

      compaction_length = material_model.length_scaling(background_porosity);
      velocity_scaling = material_model.velocity_scaling(background_porosity);


      // we also need the boundary velocity, but we can not get it from simulator access
      // TODO: write solitary wave boundary condition where the phase speed is calculated!

      max_points = 1e6;
      initial_pressure.resize(max_points);
      maximum_pressure = 0.0;
    }

    template <int dim>
    void
    SolitaryWavePostprocessor<dim>::store_initial_pressure ()
    {
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.pressure).degree);
      const unsigned int n_q_points = quadrature_formula.size();
      const double max_depth = this->get_geometry_model().maximal_depth();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      // do the same stuff we do in depth average
      std::vector<double> volume(max_points,0.0);
      std::vector<double> pressure(max_points,0.0);
      std::vector<double> p_c(n_q_points);
      double local_max_pressure = 0.0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().variable("compaction pressure").extractor_scalar()]
            .get_function_values (this->get_solution(),
                                  p_c);

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                double z = fe_values.quadrature_point(q)[dim-1];
                const unsigned int idx = static_cast<unsigned int>((z*(max_points-1))/max_depth);
                AssertThrow(idx < max_points, ExcInternalError());

                pressure[idx] += p_c[q] * fe_values.JxW(q);
                volume[idx] += fe_values.JxW(q);

                local_max_pressure = std::max (local_max_pressure, std::abs(p_c[q]));
              }
          }

      std::vector<double> volume_all(max_points, 0.0);
      Utilities::MPI::sum(volume, this->get_mpi_communicator(), volume_all);
      Utilities::MPI::sum(pressure, this->get_mpi_communicator(), initial_pressure);
      maximum_pressure = Utilities::MPI::max (local_max_pressure, this->get_mpi_communicator());

      for (unsigned int i=0; i<initial_pressure.size(); ++i)
        {
          initial_pressure[i] = initial_pressure[i] / (static_cast<double>(volume_all[i])+1e-20);
        }

      // fill the first and last element of the initial_pressure vector if they are empty
      // this makes sure they can be used for the interpolation later on
      if (initial_pressure[0] == 0.0)
        {
          unsigned int j = 1;
          while (initial_pressure[j] == 0.0)
            j++;
          initial_pressure[0] = initial_pressure[j];
        }

      if (initial_pressure[max_points-1] == 0.0)
        {
          unsigned int k = max_points-2;
          while (initial_pressure[k] == 0.0)
            k--;
          initial_pressure[max_points-1] = initial_pressure[k];
        }

      // interpolate between the non-zero elements to fill the elements that are 0
      for (unsigned int i=1; i<max_points-1; ++i)
        {
          if (initial_pressure[i] == 0.0)
            {
              // interpolate between the values we have
              unsigned int k = i-1;
              while (initial_pressure[k] == 0.0)
                {
                  Assert(k > 0, ExcInternalError());
                  k--;
                }

              unsigned int j = i+1;
              while (initial_pressure[j] == 0.0)
                j++;
              Assert(j < max_points, ExcInternalError());
              initial_pressure[i] = initial_pressure[k] + (initial_pressure[j]-initial_pressure[k]) * static_cast<double>(i-k)/static_cast<double>(j-k);
            }
        }
    }

    template <int dim>
    double
    SolitaryWavePostprocessor<dim>::compute_phase_shift ()
    {
      AssertThrow(this->introspection().compositional_name_exists("porosity"),
                  ExcMessage("Postprocessor Solitary Wave only works if there is a compositional field called porosity."));
      const unsigned int porosity_index = this->introspection().compositional_index_for_name("porosity");
      const typename Simulator<dim>::AdvectionField porosity = Simulator<dim>::AdvectionField::composition(porosity_index);

      // create a quadrature formula based on the compositional element alone.
      AssertThrow (this->introspection().n_compositional_fields > 0,
                   ExcMessage("This postprocessor cannot be used without compositional fields."));

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(porosity.base_element(this->introspection())).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      std::vector<double> compositional_values(n_q_points);

      // The idea here is to first find the maximum, and then use the analytical solution of the
      // solitary wave to calculate a phase shift for every point.
      // This has to be done separately for points left and right of the maximum, as the analytical
      // solution is only defined for coordinates > 0.
      // In the end, these values for the phase shift are averaged.

      // compute the maximum composition by quadrature (because we also need the coordinate)
      double z_max_porosity = std::numeric_limits<double>::quiet_NaN();
      {
        double local_max_porosity = -std::numeric_limits<double>::max();
        double local_max_z_location = std::numeric_limits<double>::quiet_NaN();

        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values[this->introspection().extractors.compositional_fields[porosity_index]].get_function_values (this->get_solution(),
                  compositional_values);
              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double composition = compositional_values[q];

                  if (composition > local_max_porosity)
                    {
                      local_max_porosity = composition;
                      local_max_z_location = fe_values.quadrature_point(q)[dim-1];
                    }
                }
            }

        double max_porosity = 0.0;
        MPI_max_and_data(local_max_porosity, local_max_z_location, max_porosity, z_max_porosity);
      }


      // iterate over all points and calculate the phase shift
      double phase_shift_integral = 0.0;
      unsigned int number_of_points = 0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.compositional_fields[porosity_index]].get_function_values (this->get_solution(),
                compositional_values);

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                const double composition = compositional_values[q];

                // we do not want to include the constant-porosity background in the calculation
                // nor the peak of the wave where we maximum of the composition is not a good indicator
                // if we are left or right of the maximum of the analytical solution
                if (composition > background_porosity + (amplitude - background_porosity)*0.05  && composition <= amplitude*0.9)
                  {
                    double z_analytical = compaction_length
                                          * AnalyticSolutions::solitary_wave_solution(composition/background_porosity,
                                                                                      amplitude/background_porosity);
                    double z = fe_values.quadrature_point(q)[dim-1];

                    if (z > z_max_porosity)
                      {
                        z -= offset;
                        phase_shift_integral += (z - z_analytical);
                      }
                    else
                      {
                        z = offset - z;
                        phase_shift_integral -= (z - z_analytical);
                      }

                    number_of_points += 1;
                  }
              }
          }

      double integral = Utilities::MPI::sum (phase_shift_integral, this->get_mpi_communicator())
                        / Utilities::MPI::sum (static_cast<double>(number_of_points), this->get_mpi_communicator());

      // TODO: different case for moving wave (with zero boundary velocity)
      // const double phase_speed = velocity_scaling * (2.0 * amplitude / background_porosity + 1);
      return integral; // + phase_speed * this->get_time();
    }

    template <int dim>
    std::pair<std::string,std::string>
    SolitaryWavePostprocessor<dim>::execute (TableHandler & /*statistics*/)
    {
      // as we do not have an analytical solution for the pressure, we store the initial solution
      if (this->get_timestep_number()==0)
        {
          store_initial_pressure();
          ref_func = std::make_unique<AnalyticSolutions::FunctionSolitaryWave<dim>>(offset,0.0,initial_pressure,
                                                                                     this->get_geometry_model().maximal_depth(), this->introspection().n_components);
        }

      double delta=0;

      delta = compute_phase_shift();
      // reset the phase shift of the analytical solution so we can compare the shape of the wave
      ref_func->set_delta(delta);

      // what we want to compare:
      // (1) error of the numerical phase speed c:
      // c_numerical = c_analytical - Delta / time;
      const double c_analytical = velocity_scaling * (2.0 * amplitude / background_porosity + 1);
      const double c_numerical = c_analytical + (this->get_time() > 0 ? delta / this->get_time() : 0.0);
      const double error_c = std::abs (c_numerical / c_analytical - 1);

      // (3) preservation of shape of melt fraction
      // (4) preservation of the shape of compaction pressure

      Vector<float> cellwise_errors_f (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());

      // get correct components for porosity and compaction pressure
      const unsigned int n_total_comp = this->introspection().n_components;
      ComponentSelectFunction<dim> comp_f(dim+2+dim+2, n_total_comp);
      ComponentSelectFunction<dim> comp_p(dim+1, n_total_comp);

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.pressure).degree+1);

      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_f,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_f);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_p,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_p);

      const double e_f = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_f, VectorTools::L2_norm);
      const double e_p = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p, VectorTools::L2_norm);

      std::ostringstream os;
      os << std::scientific << e_f / amplitude
         << ", " << e_p / maximum_pressure
         << ", " << error_c
         << ", " << std::abs(delta);


      return std::make_pair("Errors e_f, e_p_c_bar, e_c, delta:", os.str());
    }

  }
}



// explicit instantiations
namespace aspect
{
  namespace SolitaryWaveBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SolitaryWaveMaterial,
                                   "Solitary Wave",
                                   "A material model that corresponds to the 'SolitaryWave' benchmark "
                                   "defined in Keller et al., JGI, 2013.")

    ASPECT_REGISTER_POSTPROCESSOR(SolitaryWavePostprocessor,
                                  "solitary wave statistics",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "the Keller et al., JGI, 2013, paper with the one computed by ASPECT "
                                  "and reports the error.")

    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(SolitaryWaveInitialCondition,
                                              "Solitary wave initial condition",
                                              "Composition is set to a solitary wave function.")
  }
}
