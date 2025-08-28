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


#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/interface.h>
#include <aspect/utilities.h>
#include <aspect/newton.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <tuple>
#include <list>

#ifdef DEBUG
#ifdef ASPECT_USE_FP_EXCEPTIONS
#include <cfenv>
#endif
#endif

namespace aspect
{
  namespace MaterialModel
  {
    namespace NonlinearDependence
    {
      bool
      identifies_single_variable(const Dependence dependence)
      {
        Assert (dependence != uninitialized,
                ExcMessage ("You cannot call this function on an uninitialized dependence value!"));
        return ((dependence == temperature)
                ||
                (dependence == pressure)
                ||
                (dependence == strain_rate)
                ||
                (dependence == compositional_fields));
      }



      ModelDependence::ModelDependence ()
        :
        viscosity (uninitialized),
        density (uninitialized),
        compressibility (uninitialized),
        specific_heat (uninitialized),
        thermal_conductivity (uninitialized)
      {}
    }



// -------------------------------- Deal with registering material models and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      std::tuple
      <aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::PluginList<Interface<2>>,
      aspect::internal::Plugins::PluginList<Interface<3>>> registered_plugins;
    }



    template <int dim>
    void
    register_material_model (const std::string &name,
                             const std::string &description,
                             void (*declare_parameters_function) (ParameterHandler &),
                             std::unique_ptr<Interface<dim>> (*factory_function) ())
    {
      std::get<dim>(registered_plugins).register_plugin (name,
                                                         description,
                                                         declare_parameters_function,
                                                         factory_function);
    }


    template <int dim>
    std::unique_ptr<Interface<dim>>
    create_material_model (const std::string &model_name)
    {
      return std::get<dim>(registered_plugins).create_plugin (model_name, "Material model::Model name");
    }



    template <int dim>
    std::unique_ptr<Interface<dim>>
    create_material_model (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Material model");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      // If one sets the model name to an empty string in the input file,
      // ParameterHandler produces an error while reading the file. However,
      // if one omits specifying any model name at all (not even setting it to
      // the empty string) then the value we get here is the empty string. If
      // we don't catch this case here, we end up with awkward downstream
      // errors because the value obviously does not conform to the Pattern.
      AssertThrow(model_name != "unspecified",
                  ExcMessage("You need to select a material model "
                             "(`set Model name' in `subsection Material model')."));

      return create_material_model<dim> (model_name);
    }



    template <int dim>
    const NonlinearDependence::ModelDependence &
    Interface<dim>::
    get_model_dependence() const
    {
      return model_dependence;
    }



    template <int dim>
    void
    Interface<dim>::
    create_additional_named_outputs (MaterialModelOutputs &/*outputs*/) const
    {
      // by default we do nothing!
    }



    template <int dim>
    void
    Interface<dim>::
    fill_additional_material_model_inputs(MaterialModel::MaterialModelInputs<dim> &input,
                                          const LinearAlgebra::BlockVector        &solution,
                                          const FEValuesBase<dim>                 &fe_values,
                                          const Introspection<dim>                &introspection) const
    {
      // go through the list of additional inputs and fill them
      for (unsigned int i=0; i<input.additional_inputs.size(); ++i)
        input.additional_inputs[i]->fill(solution, fe_values, introspection);
    }



    template <int dim>
    std::string
    get_valid_model_names_pattern ()
    {
      return std::get<dim>(registered_plugins).get_pattern_of_names ();
    }


    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
      // declare the actual entry in the parameter file
      prm.enter_subsection ("Material model");
      {
        const std::string pattern_of_names = get_valid_model_names_pattern<dim>();
        prm.declare_entry ("Model name", "unspecified",
                           Patterns::Selection (pattern_of_names+"|unspecified"),
                           "The name of the material model to be used in "
                           "this simulation. There are many material models "
                           "you can choose from, as listed below. They generally "
                           "fall into two category: (i) models that implement "
                           "a particular case of material behavior, (ii) models "
                           "that modify other models in some way. We sometimes "
                           "call the latter ``compositing models''. An example "
                           "of a compositing model is the ``depth dependent'' model "
                           "below in that it takes another, freely choosable "
                           "model as its base and then modifies that model's "
                           "output in some way."
                           "\n\n"
                           "You can select one of the following models:\n\n"
                           +
                           std::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Material model interface",
                                                            out);
    }


    template <int dim>
    MaterialModelInputs<dim>::MaterialModelInputs(const unsigned int n_points,
                                                  const unsigned int n_comp)
      :
      position(n_points, Point<dim>(numbers::signaling_nan<Tensor<1,dim>>())),
      temperature(n_points, numbers::signaling_nan<double>()),
      pressure(n_points, numbers::signaling_nan<double>()),
      pressure_gradient(n_points, numbers::signaling_nan<Tensor<1,dim>>()),
      velocity(n_points, numbers::signaling_nan<Tensor<1,dim>>()),
      composition(n_points, std::vector<double>(n_comp, numbers::signaling_nan<double>())),
      strain_rate(n_points, numbers::signaling_nan<SymmetricTensor<2,dim>>()),
      current_cell(),
      requested_properties(MaterialProperties::all_properties)
    {}



    template <int dim>
    MaterialModelInputs<dim>::MaterialModelInputs(const DataPostprocessorInputs::Vector<dim> &input_data,
                                                  const Introspection<dim> &introspection,
                                                  const bool compute_strain_rate)
      :
      position(input_data.evaluation_points),
      temperature(input_data.solution_values.size(), numbers::signaling_nan<double>()),
      pressure(input_data.solution_values.size(), numbers::signaling_nan<double>()),
      pressure_gradient(input_data.solution_values.size(), numbers::signaling_nan<Tensor<1,dim>>()),
      velocity(input_data.solution_values.size(), numbers::signaling_nan<Tensor<1,dim>>()),
      composition(input_data.solution_values.size(), std::vector<double>(introspection.n_compositional_fields, numbers::signaling_nan<double>())),
      strain_rate(input_data.solution_values.size(), numbers::signaling_nan<SymmetricTensor<2,dim>>()),
      current_cell(input_data.template get_cell<dim>()),
      requested_properties(MaterialProperties::all_properties)
    {
      AssertThrow (compute_strain_rate == true,
                   ExcMessage ("The option to not compute the strain rate is no longer supported."));

      for (unsigned int q=0; q<input_data.solution_values.size(); ++q)
        {
          Tensor<2,dim> grad_u;
          for (unsigned int d=0; d<dim; ++d)
            {
              grad_u[d] = input_data.solution_gradients[q][introspection.component_indices.velocities[d]];
              this->velocity[q][d] = input_data.solution_values[q][introspection.component_indices.velocities[d]];
              this->pressure_gradient[q][d] = input_data.solution_gradients[q][introspection.component_indices.pressure][d];
            }

          this->strain_rate[q] = symmetrize (grad_u);
          this->pressure[q] = input_data.solution_values[q][introspection.component_indices.pressure];
          this->temperature[q] = input_data.solution_values[q][introspection.component_indices.temperature];

          for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
            this->composition[q][c] = input_data.solution_values[q][introspection.component_indices.compositional_fields[c]];
        }
    }



    template <int dim>
    MaterialModelInputs<dim>::MaterialModelInputs(const FEValuesBase<dim,dim> &fe_values,
                                                  const typename DoFHandler<dim>::active_cell_iterator &cell_x,
                                                  const Introspection<dim> &introspection,
                                                  const LinearAlgebra::BlockVector &solution_vector,
                                                  const bool compute_strain_rate)
      :
      position(fe_values.get_quadrature_points()),
      temperature(fe_values.n_quadrature_points, numbers::signaling_nan<double>()),
      pressure(fe_values.n_quadrature_points, numbers::signaling_nan<double>()),
      pressure_gradient(fe_values.n_quadrature_points, numbers::signaling_nan<Tensor<1,dim>>()),
      velocity(fe_values.n_quadrature_points, numbers::signaling_nan<Tensor<1,dim>>()),
      composition(fe_values.n_quadrature_points, std::vector<double>(introspection.n_compositional_fields, numbers::signaling_nan<double>())),
      strain_rate(fe_values.n_quadrature_points, numbers::signaling_nan<SymmetricTensor<2,dim>>()),
      current_cell (cell_x),
      requested_properties(MaterialProperties::all_properties)
    {
      // Call the function reinit to populate the new arrays.
      this->reinit(fe_values, current_cell, introspection, solution_vector, compute_strain_rate);
    }



    template <int dim>
    MaterialModelInputs<dim>::MaterialModelInputs(const MaterialModelInputs &source)
      :
      position(source.position),
      temperature(source.temperature),
      pressure(source.pressure),
      pressure_gradient(source.pressure_gradient),
      velocity(source.velocity),
      composition(source.composition),
      strain_rate(source.strain_rate),
      current_cell(source.current_cell),
      requested_properties(source.requested_properties)
    {
      Assert (source.additional_inputs.size() == 0,
              ExcMessage ("You can not copy MaterialModelInputs objects that have "
                          "additional input objects attached"));
    }



    template <int dim>
    void
    MaterialModelInputs<dim>::reinit(const FEValuesBase<dim,dim> &fe_values,
                                     const typename DoFHandler<dim>::active_cell_iterator &cell_x,
                                     const Introspection<dim> &introspection,
                                     const LinearAlgebra::BlockVector &solution_vector,
                                     const bool compute_strain_rate)
    {
      AssertThrow (compute_strain_rate == true,
                   ExcMessage ("The option to not compute the strain rate is no longer supported."));

      // Populate the arrays that hold solution values and gradients
      fe_values[introspection.extractors.temperature].get_function_values (solution_vector, this->temperature);
      fe_values[introspection.extractors.velocities].get_function_values (solution_vector, this->velocity);
      fe_values[introspection.extractors.pressure].get_function_values (solution_vector, this->pressure);
      fe_values[introspection.extractors.pressure].get_function_gradients (solution_vector, this->pressure_gradient);
      fe_values[introspection.extractors.velocities].get_function_symmetric_gradients (solution_vector, this->strain_rate);

      // Vectors for evaluating the compositional field parts of the finite element solution
      std::vector<std::vector<double>> composition_values (introspection.n_compositional_fields,
                                                            std::vector<double> (fe_values.n_quadrature_points));
      for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
        fe_values[introspection.extractors.compositional_fields[c]]
        .get_function_values(solution_vector,composition_values[c]);

      // Then copy these values to exchange the inner and outer vector, because for the material
      // model we need a vector with values of all the compositional fields for every quadrature point
      for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
        {
          for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
            this->composition[q][c] = composition_values[c][q];
        }

      // Finally also record quadrature point positions and the cell
      this->position = fe_values.get_quadrature_points();
      this->current_cell = cell_x;
    }



    template <int dim>
    unsigned int
    MaterialModelInputs<dim>::n_evaluation_points() const
    {
      return position.size();
    }



    template <int dim>
    bool
    MaterialModelInputs<dim>::requests_property(const MaterialProperties::Property &property) const
    {
      // Note that this means 'requested_properties' can include other properties than
      // just 'property', but in any case it at least requests 'property'.
      return (requested_properties & property) != 0;
    }



    template <int dim>
    MaterialModelOutputs<dim>::MaterialModelOutputs(const unsigned int n_points,
                                                    const unsigned int n_comp)
      :
      viscosities(n_points, numbers::signaling_nan<double>()),
      densities(n_points, numbers::signaling_nan<double>()),
      thermal_expansion_coefficients(n_points, numbers::signaling_nan<double>()),
      specific_heat(n_points, numbers::signaling_nan<double>()),
      thermal_conductivities(n_points, numbers::signaling_nan<double>()),
      compressibilities(n_points, numbers::signaling_nan<double>()),
      entropy_derivative_pressure(n_points, numbers::signaling_nan<double>()),
      entropy_derivative_temperature(n_points, numbers::signaling_nan<double>()),
      reaction_terms(n_points, std::vector<double>(n_comp, numbers::signaling_nan<double>()))
    {}


    template <int dim>
    MaterialModelOutputs<dim>::MaterialModelOutputs(const MaterialModelOutputs<dim> &source)
      :
      viscosities(source.viscosities),
      densities(source.densities),
      thermal_expansion_coefficients(source.thermal_expansion_coefficients),
      specific_heat(source.specific_heat),
      thermal_conductivities(source.thermal_conductivities),
      compressibilities(source.compressibilities),
      entropy_derivative_pressure(source.entropy_derivative_pressure),
      entropy_derivative_temperature(source.entropy_derivative_temperature),
      reaction_terms(source.reaction_terms),
      additional_outputs()
    {
      Assert (source.additional_outputs.size() == 0,
              ExcMessage ("You can not copy MaterialModelOutputs objects that have "
                          "additional output objects attached"));
    }



    template <int dim>
    unsigned int
    MaterialModelOutputs<dim>::n_evaluation_points() const
    {
      return densities.size();
    }



    namespace MaterialAveraging
    {
      std::string get_averaging_operation_names ()
      {
        return "none|default averaging|arithmetic average|harmonic average|geometric average|pick largest|project to Q1|log average|harmonic average only viscosity|geometric average only viscosity|project to Q1 only viscosity";
      }


      AveragingOperation parse_averaging_operation_name (const std::string &s)
      {
        if (s == "none")
          return none;
        else if (s == "arithmetic average")
          return arithmetic_average;
        else if (s == "harmonic average")
          return harmonic_average;
        else if (s == "geometric average")
          return geometric_average;
        else if (s == "pick largest")
          return pick_largest;
        else if (s == "project to Q1")
          return project_to_Q1;
        else if (s == "log average")
          return log_average;
        else if (s == "harmonic average only viscosity")
          return harmonic_average_only_viscosity;
        else if (s == "geometric average only viscosity")
          return geometric_average_only_viscosity;
        else if (s == "project to Q1 only viscosity")
          return project_to_Q1_only_viscosity;
        else if (s == "default averaging")
          return default_averaging;
        else
          AssertThrow (false,
                       ExcMessage ("The value <" + s + "> for a material "
                                   "averaging operation is not one of the "
                                   "valid values."));

        return none;
      }



      namespace
      {
        bool
        all_entries_NaN (const std::vector<double> &values)
        {
          for (const double value : values)
            if (std::isnan(value) == false)
              return false;

          return true;
        }
      }



      // Do the requested averaging operation for one array. The
      // projection matrix argument is only used if the operation
      // chosen is project_to_Q1.
      void average_property (const AveragingOperation  operation,
                             const FullMatrix<double>      &projection_matrix,
                             const FullMatrix<double>      &expansion_matrix,
                             std::vector<double>           &values_out)
      {
#ifdef DEBUG
#ifdef ASPECT_USE_FP_EXCEPTIONS
        // disable floating point exceptions while averaging. Errors will be reported
        // as soon as somebody will try to use the averaged values later.
        fedisableexcept(FE_DIVBYZERO|FE_INVALID);
#endif
#endif

        // if an output field has not been filled (because it was
        // not requested), then simply do nothing -- no harm no foul
        // note that it is still an error if only some entries are NaN
        if (values_out.size() == 0 || all_entries_NaN(values_out) == true)
          return;

        const unsigned int N = values_out.size();
        const unsigned int P = expansion_matrix.n();
        Assert ((P==0) || (/*dim=2*/ P==4) || (/*dim=3*/ P==8),
                ExcInternalError());
        Assert (((operation == project_to_Q1) &&
                 (projection_matrix.m() == P) &&
                 (projection_matrix.n() == N) &&
                 (expansion_matrix.m() == N) &&
                 (expansion_matrix.n() == P))
                ||
                ((projection_matrix.m() == 0) &&
                 (projection_matrix.n() == 0)),
                ExcInternalError());

        // otherwise do as instructed
        switch (operation)
          {
            case none:
            {
              break;
            }

            case arithmetic_average:
            {
              double sum = 0;
              for (unsigned int i=0; i<N; ++i)
                sum += values_out[i];

              const double average = sum/N;
              for (unsigned int i=0; i<N; ++i)
                values_out[i] = average;
              break;
            }

            case harmonic_average:
            {
              // if one of the values is zero, the average is 0.0
              for (unsigned int i=0; i<N; ++i)
                if (values_out[i] == 0.0)
                  {
                    for (unsigned int j=0; j<N; ++j)
                      values_out[j] = 0.0;
                    return;
                  }

              double sum = 0;
              for (unsigned int i=0; i<N; ++i)
                sum += 1./values_out[i];

              const double average = 1./(sum/N);
              for (unsigned int i=0; i<N; ++i)
                values_out[i] = average;
              break;
            }

            case geometric_average:
            {
              double average = 1;
              for (unsigned int i=0; i<N; ++i)
                {
                  Assert (values_out[i] >= 0,
                          ExcMessage ("Computing the geometric average "
                                      "only makes sense for non-negative "
                                      "quantities."));
                  average *= std::pow (values_out[i], 1./N);
                }

              for (unsigned int i=0; i<N; ++i)
                values_out[i] = average;
              break;
            }

            case pick_largest:
            {
              double max = std::numeric_limits<double>::lowest();
              for (unsigned int i=0; i<N; ++i)
                max = std::max(max, values_out[i]);

              for (unsigned int i=0; i<N; ++i)
                values_out[i] = max;
              break;
            }

            case project_to_Q1:
            {
              // we will need the min/max values below, for use
              // after the projection operation
              double min = std::numeric_limits<double>::max();
              for (unsigned int i=0; i<N; ++i)
                min = std::min(min, values_out[i]);

              double max = std::numeric_limits<double>::lowest();
              for (unsigned int i=0; i<N; ++i)
                max = std::max(max, values_out[i]);

              // take the projection matrix and apply it to the
              // values. as explained in the documentation of the
              // compute_projection_matrix, this performs the operation
              // we want in the current context
              Vector<double> x (N), z(P), y(N);
              for (unsigned int i=0; i<N; ++i)
                y(i) = values_out[i];
              projection_matrix.vmult (z, y);

              // now that we have the Q1 values, restrict them to
              // the min/max range of the original data
              for (unsigned int i=0; i<P; ++i)
                z[i] = std::max (min,
                                 std::min (max,
                                           z[i]));

              // then expand back to the quadrature points
              expansion_matrix.vmult (x, z);
              for (unsigned int i=0; i<N; ++i)
                values_out[i] = x(i);

              break;
            }

            case log_average:
            {
              double sum = 0;
              for (unsigned int i=0; i<N; ++i)
                {
                  if (values_out[i] == 0.0)
                    {
                      sum = -std::numeric_limits<double>::infinity();
                      break;
                    }
                  Assert (values_out[i] > 0.0,
                          ExcMessage ("Computing the log average "
                                      "only makes sense for positive "
                                      "quantities."));
                  sum += std::log10(values_out[i]);
                }
              const double log_value_average = std::pow (10., sum/N);
              for (unsigned int i=0; i<N; ++i)
                values_out[i] = log_value_average;
              break;
            }

            default:
            {
              AssertThrow (false,
                           ExcMessage ("This averaging operation is not implemented."));
            }
          }

#ifdef DEBUG
#ifdef ASPECT_USE_FP_EXCEPTIONS
        // enable floating point exceptions again:
        feenableexcept(FE_DIVBYZERO|FE_INVALID);
#endif
#endif
      }


      /**
       * Given a quadrature formula, compute a matrices $E, M^{-1}F$
       * representing a linear operator in the following way: Let
       * there be a vector $F$ with $N$ elements where the elements
       * are data stored at each of the $N$ quadrature points of the
       * given quadrature formula. (For this quadrature formula, we only
       * care about the locations of the quadrature points on the
       * reference cell, but not about quadrature weights.)
       * Then project this data into a $Q_1$ space and evaluate this
       * projection at the same quadrature points.
       *
       * This operator can be expressed in the following way where $P=2^{dim}$
       * is the number of degrees of freedom of the $Q_1$ element:
       *
       * Let $y$ be the input vector with $N$ elements. We would like to find
       * a function $u_h(x)$ that is a $Q_1$ function that is closest to
       * the data points, i.e., we seek
       * @f{align*}{
       *   min_{u_h \in Q_1(K)} \frac 12 \sum_q |u_h(x_q) - y_q|^2 |J(x_q)| w_q.
       * @f}
       * where $x_q$ are the evaluation points, $J$ is the Jacobian matrix of the
       * mapping from reference to real cell, and $w_q$ are the
       * quadrature weights.
       * Expanding $u_h = \sum_j U_j \varphi_j$, this leads to the following
       * minimization problem for the vector of coefficients $U$:
       * @f{align*}{
       *   min_{U} \frac 12 \sum_q |\sum_j U_j \varphi_j(x_q) - y_q|^2 w_q,
       * @f}
       * which can be rewritten as
       * @f{align*}{
       *   min_{U} \frac 12 \sum_q
       *       [ \sum_i \sum_j U_i U_j \varphi_i(x_q)\varphi_j(x_q)
       *        - 2 \sum_j U_j \varphi_j(x_q) y_q
       *        + y_q^2                                             ].
       * @f}
       * This optimization problem has the following solution that one can arrive at
       * simply by setting the derivative of the objective function with regard to
       * $U$ to zero:
       * @f{align*}{
       *       \sum q \sum_i \sum_j \varphi_i(x_q)\varphi_j(x_q) |J(x_q)| w_q U_j
       *       =
       *       \sum_q \varphi_i(x_q) |J(x_q)| w_q y_q
       * @f}
       * If we let $F$ be the $P \times N$ matrix so that
       * @f{align*}{
       *   F_{iq} = \varphi_i(x_q) |J(x_q)| w_q
       * @f}
       * Next, let $M_{ij}$ be the $P\times P$ mass matrix on the
       * $Q_1$ space with regard to the given evaluation points and corresponding
       * weights. Then $U = M^{-1}Fy$ corresponds to the projection
       * of $y$ into the $Q_1$ space (or, more correctly, it
       * corresponds to vector of the nodal values of this projection).
       *
       * Finally, let $E_{qi} = \varphi_i(x_q)$ be the evaluation
       * operation of shape functions at quadrature points.
       *
       * Then, the operation $X=EM^{-1}F$ is the operation we seek.
       * This function computes the matrices $E$ and $M^{-1}F$.
       */
      template <int dim>
      void compute_projection_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                      const Quadrature<dim>   &quadrature_formula,
                                      const Mapping<dim>      &mapping,
                                      FullMatrix<double>      &projection_matrix,
                                      FullMatrix<double>      &expansion_matrix)
      {
        static const FE_Q<dim> fe(1);
        FEValues<dim> fe_values (mapping, fe, quadrature_formula,
                                 update_values | update_JxW_values);

        const unsigned int P = fe.dofs_per_cell;
        const unsigned int N = quadrature_formula.size();

        FullMatrix<double> F (P, N);
        FullMatrix<double> M (P, P);

        projection_matrix.reinit (P, N);
        expansion_matrix.reinit (N, P);

        // reinitialize the fe_values object with the current cell. we get a
        // DoFHandler cell, but we are not going to use it with the
        // finite element associated with that DoFHandler, so cast it back
        // to just a tria iterator (all we need anyway is the geometry)
        fe_values.reinit (typename Triangulation<dim>::active_cell_iterator(cell));

        // compute the matrices F, M, E
        for (unsigned int i=0; i<P; ++i)
          for (unsigned int q=0; q<N; ++q)
            F(i,q) = fe_values.shape_value(i,q) *
                     fe_values.JxW(q);

        for (unsigned int i=0; i<P; ++i)
          for (unsigned int j=0; j<P; ++j)
            for (unsigned int q=0; q<N; ++q)
              M(i,j) += fe_values.shape_value(i,q) *
                        fe_values.shape_value(j,q) *
                        fe_values.JxW(q);

        for (unsigned int q=0; q<N; ++q)
          for (unsigned int i=0; i<P; ++i)
            expansion_matrix(q,i) = fe_values.shape_value(i,q);

        // replace M by M^{-1}
        M.gauss_jordan();

        // form M^{-1} F
        M.mmult (projection_matrix, F);
      }


      /**
       * Calculate the weight for viscosity derivative, which depends on
       * the material averaging scheme. Currently the newton method is
       * only compatible with arithmetic average, harmonic average and
       * geometric/log average for viscosity.
       */
      double
      compute_viscosity_derivative_averaging_weight(const AveragingOperation operation,
                                                    const double average_viscosity,
                                                    const double viscosity_before_averaging,
                                                    const double one_over_Nq)
      {
        switch (operation)
          {
            case none:
              // should never reach here:
              return numbers::signaling_nan<double>();

            case arithmetic_average:
              return one_over_Nq;

            case harmonic_average:
            case harmonic_average_only_viscosity:
              return Utilities::fixed_power<2,double>(average_viscosity / viscosity_before_averaging)
                     * one_over_Nq;

            case geometric_average:
            case geometric_average_only_viscosity:
            case log_average:
              return (average_viscosity / viscosity_before_averaging)
                     * one_over_Nq;

            default:
              AssertThrow(false,
                          ExcMessage("The Newton method currently only works if the material "
                                     "averaging scheme is ``none'', ``arithmetic average'', "
                                     "``harmonic average (only viscosity)'', ``geometric "
                                     "average (only viscosity)'' or ``log average''."));
          }
      }


      template <int dim>
      void average (const AveragingOperation operation,
                    const typename DoFHandler<dim>::active_cell_iterator &cell,
                    const Quadrature<dim>         &quadrature_formula,
                    const Mapping<dim>            &mapping,
                    const MaterialProperties::Property &requested_properties,
                    MaterialModelOutputs<dim>     &values_out)
      {
        if (operation == none)
          return;

        FullMatrix<double> projection_matrix;
        FullMatrix<double> expansion_matrix;

        const bool average_viscosity = requested_properties & MaterialProperties::Property::viscosity;

        if (operation == project_to_Q1
            ||
            operation == project_to_Q1_only_viscosity)
          {
            Assert (quadrature_formula.size() == values_out.n_evaluation_points(),
                    ExcMessage("When asking for a Q1-type averaging operation, "
                               "this function requires to know the locations of "
                               "the evaluation points."));
            projection_matrix.reinit (quadrature_formula.size(),
                                      quadrature_formula.size());
            compute_projection_matrix (cell,
                                       quadrature_formula,
                                       mapping,
                                       projection_matrix,
                                       expansion_matrix);
          }

        // store the original viscosities if we need to compute the
        // system jacobian later on
        const std::shared_ptr<MaterialModelDerivatives<dim>> derivatives =
          values_out.template get_additional_output_object<MaterialModelDerivatives<dim>>();

        std::vector<double> viscosity_before_averaging;
        if (derivatives != nullptr)
          viscosity_before_averaging = values_out.viscosities;

        // compute the average of viscosity
        if (average_viscosity)
          {
            if (operation == harmonic_average_only_viscosity)
              average_property (harmonic_average, projection_matrix, expansion_matrix,
                                values_out.viscosities);

            else if (operation == geometric_average_only_viscosity)
              average_property (geometric_average, projection_matrix, expansion_matrix,
                                values_out.viscosities);

            else if (operation == project_to_Q1_only_viscosity)
              average_property (project_to_Q1, projection_matrix, expansion_matrix,
                                values_out.viscosities);

            else
              average_property (operation, projection_matrix, expansion_matrix,
                                values_out.viscosities);

            // calculate the weight of viscosity derivative at each
            // quadrature point
            if (derivatives != nullptr)
              {
                for (unsigned int q = 0; q < values_out.n_evaluation_points(); ++q)
                  derivatives->viscosity_derivative_averaging_weights[q] =
                    compute_viscosity_derivative_averaging_weight(
                      operation, values_out.viscosities[q], viscosity_before_averaging[q],
                      1. / values_out.n_evaluation_points());
              }
          }

        if (operation == harmonic_average_only_viscosity ||
            operation == geometric_average_only_viscosity ||
            operation == project_to_Q1_only_viscosity)
          return;

        average_property (operation, projection_matrix, expansion_matrix,
                          values_out.densities);
        average_property (operation, projection_matrix, expansion_matrix,
                          values_out.thermal_expansion_coefficients);
        average_property (operation, projection_matrix, expansion_matrix,
                          values_out.specific_heat);
        average_property (operation, projection_matrix, expansion_matrix,
                          values_out.thermal_conductivities);
        average_property (operation, projection_matrix, expansion_matrix,
                          values_out.compressibilities);
        average_property (operation, projection_matrix, expansion_matrix,
                          values_out.entropy_derivative_pressure);
        average_property (operation, projection_matrix, expansion_matrix,
                          values_out.entropy_derivative_temperature);

        // The reaction terms are unfortunately stored in reverse
        // indexing. It's also not quite clear whether these should
        // really be averaged, so avoid this for now.

        // average all additional outputs
        for (unsigned int i=0; i<values_out.additional_outputs.size(); ++i)
          values_out.additional_outputs[i]->average (operation, projection_matrix, expansion_matrix);
      }



      AveragingOperation
      get_averaging_operation_for_viscosity(const AveragingOperation operation)
      {
        AveragingOperation operation_for_viscosity = operation;
        switch (operation)
          {
            case harmonic_average:
              operation_for_viscosity = harmonic_average_only_viscosity;
              break;

            case geometric_average:
              operation_for_viscosity = geometric_average_only_viscosity;
              break;

            case project_to_Q1:
              operation_for_viscosity = project_to_Q1_only_viscosity;
              break;

            default:
              operation_for_viscosity = operation;
          }

        return operation_for_viscosity;
      }
    }



    template <int dim>
    NamedAdditionalMaterialOutputs<dim>::
    NamedAdditionalMaterialOutputs(const std::vector<std::string> &output_names)
      :
      names(output_names)
    {}



    template <int dim>
    NamedAdditionalMaterialOutputs<dim>::
    NamedAdditionalMaterialOutputs(const std::vector<std::string> &output_names,
                                   const unsigned int n_points)
      :
      output_values(output_names.size(), std::vector<double>(n_points, numbers::signaling_nan<double>())),
      names(output_names)
    {}



    template <int dim>
    NamedAdditionalMaterialOutputs<dim>::
    ~NamedAdditionalMaterialOutputs()
      = default;



    template <int dim>
    const std::vector<std::string> &
    NamedAdditionalMaterialOutputs<dim>::get_names() const
    {
      return names;
    }



    template <int dim>
    std::vector<double>
    NamedAdditionalMaterialOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      // In this function we extract the values for the nth output
      // The number of outputs is the outer vector
      Assert (output_values.size() > idx,
              ExcMessage ("The requested output index is out of range for output_values."));
      Assert (output_values[idx].size() > 0,
              ExcMessage ("There must be one or more points for the nth output."));
      return output_values[idx];
    }



    namespace
    {
      std::vector<std::string> make_prescribed_dilation_outputs_names()
      {
        std::vector<std::string> names;
        names.emplace_back("dilation_lhs_term");
        names.emplace_back("dilation_rhs_term");
        return names;
      }
    }



    template <int dim>
    PrescribedPlasticDilation<dim>::PrescribedPlasticDilation (const unsigned int n_points)
      : NamedAdditionalMaterialOutputs<dim>(make_prescribed_dilation_outputs_names()),
        dilation_lhs_term(n_points, numbers::signaling_nan<double>()),
        dilation_rhs_term(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double> PrescribedPlasticDilation<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 2);
      switch (idx)
        {
          case 0:
            return dilation_lhs_term;

          case 1:
            return dilation_rhs_term;

          default:
            AssertThrow(false, ExcInternalError());
        }
      // we will never get here, so just return something
      return std::vector<double>();
    }



    namespace
    {
      std::vector<std::string> make_seismic_additional_outputs_names()
      {
        std::vector<std::string> names;
        names.emplace_back("seismic_Vs");
        names.emplace_back("seismic_Vp");
        return names;
      }
    }



    template <int dim>
    SeismicAdditionalOutputs<dim>::SeismicAdditionalOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_seismic_additional_outputs_names()),
      vs(n_points, numbers::signaling_nan<double>()),
      vp(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
    SeismicAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 2);
      switch (idx)
        {
          case 0:
            return vs;

          case 1:
            return vp;

          default:
            AssertThrow(false, ExcInternalError());
        }
      // we will never get here, so just return something
      return vs;
    }



    namespace
    {
      std::vector<std::string> make_reaction_rate_outputs_names(const unsigned int n_comp)
      {
        std::vector<std::string> names;
        for (unsigned int c=0; c<n_comp; ++c)
          names.push_back("reaction_rate_C" + Utilities::int_to_string(c));

        return names;
      }



      std::vector<std::string> make_prescribed_field_output_names(const unsigned int n_comp)
      {
        std::vector<std::string> names;
        for (unsigned int c=0; c<n_comp; ++c)
          names.push_back("prescribed_field_output_C" + Utilities::int_to_string(c));
        return names;
      }
    }



    template <int dim>
    ReactionRateOutputs<dim>::ReactionRateOutputs (const unsigned int n_points,
                                                   const unsigned int n_comp)
      :
      NamedAdditionalMaterialOutputs<dim>(make_reaction_rate_outputs_names(n_comp)),
      reaction_rates(n_points, std::vector<double>(n_comp, std::numeric_limits<double>::quiet_NaN()))
    {}



    template <int dim>
    std::vector<double>
    ReactionRateOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      // we have to extract the reaction rate outputs for one particular compositional
      // field, but the vector in the material model outputs is sorted so that the
      // number of evaluation points (and not the compositional fields) is the outer
      // vector
      std::vector<double> cth_reaction_rates(reaction_rates.size());
      for (unsigned int q=0; q<reaction_rates.size(); ++q)
        cth_reaction_rates[q] = reaction_rates[q][idx];

      return cth_reaction_rates;
    }



    namespace
    {
      std::vector<std::string> make_phase_outputs_names()
      {
        std::vector<std::string> names;
        names.emplace_back("phase");
        return names;
      }
    }



    template <int dim>
    PhaseOutputs<dim>::PhaseOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_phase_outputs_names(), n_points)
    {}



    template <int dim>
    PrescribedFieldOutputs<dim>::PrescribedFieldOutputs (const unsigned int n_points,
                                                         const unsigned int n_comp)
      :
      NamedAdditionalMaterialOutputs<dim>(make_prescribed_field_output_names(n_comp)),
      prescribed_field_outputs(n_points, std::vector<double>(n_comp, std::numeric_limits<double>::quiet_NaN()))
    {}



    template <int dim>
    std::vector<double>
    PrescribedFieldOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      // we have to extract the prescribed field outputs for one particular compositional
      // field, but the vector in the material model outputs is sorted so that the
      // number of evaluation points (and not the compositional fields) is the outer
      // vector
      std::vector<double> nth_prescribed_field_output(prescribed_field_outputs.size());
      for (unsigned int q=0; q<prescribed_field_outputs.size(); ++q)
        nth_prescribed_field_output[q] = prescribed_field_outputs[q][idx];
      return nth_prescribed_field_output;
    }



    template <int dim>
    PrescribedTemperatureOutputs<dim>::PrescribedTemperatureOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(std::vector<std::string>(1,"prescribed_temperature")),
      prescribed_temperature_outputs(n_points, std::numeric_limits<double>::quiet_NaN())
    {}



    template <int dim>
    std::vector<double>
    PrescribedTemperatureOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      (void)idx;
      AssertIndexRange (idx, 1);
      return prescribed_temperature_outputs;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<MaterialModel::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::Interface<2>>::plugins = nullptr;

      template <>
      std::list<internal::Plugins::PluginList<MaterialModel::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::Interface<3>>::plugins = nullptr;
    }
  }

  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_material_model<dim> (const std::string &, \
                                const std::string &, \
                                void ( *) (ParameterHandler &), \
                                std::unique_ptr<Interface<dim>>( *) ()); \
  \
  template \
  std::string \
  get_valid_model_names_pattern<dim> (); \
  \
  template  \
  void \
  declare_parameters<dim> (ParameterHandler &); \
  \
  template \
  std::unique_ptr<Interface<dim>> \
  create_material_model<dim> (const std::string &model_name); \
  \
  template \
  void \
  write_plugin_graph<dim> (std::ostream &); \
  \
  template \
  std::unique_ptr<Interface<dim>> \
  create_material_model<dim> (ParameterHandler &prm); \
  \
  template class MaterialModelInputs<dim>; \
  \
  template class MaterialModelOutputs<dim>; \
  \
  template class AdditionalMaterialOutputs<dim>; \
  \
  template class NamedAdditionalMaterialOutputs<dim>; \
  \
  template class SeismicAdditionalOutputs<dim>; \
  \
  template class ReactionRateOutputs<dim>; \
  \
  template class PhaseOutputs<dim>; \
  \
  template class PrescribedPlasticDilation<dim>; \
  \
  template class PrescribedFieldOutputs<dim>; \
  \
  template class PrescribedTemperatureOutputs<dim>; \
  \
  namespace MaterialAveraging \
  { \
    template                \
    void average (const AveragingOperation operation, \
                  const DoFHandler<dim>::active_cell_iterator &cell, \
                  const Quadrature<dim>     &quadrature_formula, \
                  const Mapping<dim>        &mapping, \
                  const MaterialProperties::Property &requested_properties, \
                  MaterialModelOutputs<dim> &values_out); \
  }


    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
