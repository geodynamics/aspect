/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/interface.h>
#include <aspect/utilities.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx11/tuple.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <list>


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


    template <int dim>
    Interface<dim>::~Interface ()
    {}

    template <int dim>
    void
    Interface<dim>::initialize ()
    {}

    template <int dim>
    void
    Interface<dim>::update ()
    {}



    template <int dim>
    double
    Interface<dim>::reference_thermal_expansion_coefficient () const
    {
      Assert(false, ExcMessage("Implement individual functions or evaluate() in material model."));
      return 1.0;
    }

    template <int dim>
    void
    Interface<dim>::
    declare_parameters (dealii::ParameterHandler &)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &)
    {}


// -------------------------------- Deal with registering material models and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      std_cxx11::tuple
      <void *,
      void *,
      internal::Plugins::PluginList<Interface<2> >,
      internal::Plugins::PluginList<Interface<3> > > registered_plugins;
    }



    template <int dim>
    void
    register_material_model (const std::string &name,
                             const std::string &description,
                             void (*declare_parameters_function) (ParameterHandler &),
                             Interface<dim> *(*factory_function) ())
    {
      std_cxx11::get<dim>(registered_plugins).register_plugin (name,
                                                               description,
                                                               declare_parameters_function,
                                                               factory_function);
    }


    template <int dim>
    Interface<dim> *
    create_material_model (const std::string &model_name)
    {
      Interface<dim> *plugin = std_cxx11::get<dim>(registered_plugins).create_plugin (model_name,
                                                                                      "Material model::Model name");
      return plugin;
    }


    template <int dim>
    Interface<dim> *
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
                             "('set Model name' in 'subsection Material model')."));

      return create_material_model<dim> (model_name);
    }


    template <int dim>
    double
    Interface<dim>::
    viscosity_ratio (const double,
                     const double,
                     const std::vector<double> &,
                     const SymmetricTensor<2,dim> &,
                     const Point<dim> &) const
    {
      return 1.0;
    }


    template <int dim>
    const NonlinearDependence::ModelDependence &
    Interface<dim>::
    get_model_dependence() const
    {
      return model_dependence;
    }


    template <int dim>
    double
    Interface<dim>::
    seismic_Vp (const double,
                const double,
                const std::vector<double> &, /*composition*/
                const Point<dim> &) const
    {
      return -1.0;
    }


    template <int dim>
    double
    Interface<dim>::
    seismic_Vs (const double,
                const double,
                const std::vector<double> &, /*composition*/
                const Point<dim> &) const
    {
      return -1.0;
    }


    template <int dim>
    unsigned int
    Interface<dim>::
    thermodynamic_phase (const double,
                         const double,
                         const std::vector<double> & /*composition*/) const
    {
      return 0;
    }


    template <int dim>
    std::string
    get_valid_model_names_pattern ()
    {
      return std_cxx11::get<dim>(registered_plugins).get_pattern_of_names ();
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
                           std_cxx11::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection ();

      std_cxx11::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    MaterialModelInputs<dim>::MaterialModelInputs(const unsigned int n_points,
                                                  const unsigned int n_comp)
      :
      position(n_points, Point<dim>(aspect::Utilities::signaling_nan<Tensor<1,dim> >())),
      temperature(n_points, aspect::Utilities::signaling_nan<double>()),
      pressure(n_points, aspect::Utilities::signaling_nan<double>()),
      pressure_gradient(n_points, aspect::Utilities::signaling_nan<Tensor<1,dim> >()),
      velocity(n_points, aspect::Utilities::signaling_nan<Tensor<1,dim> >()),
      composition(n_points, std::vector<double>(n_comp, aspect::Utilities::signaling_nan<double>())),
      strain_rate(n_points, aspect::Utilities::signaling_nan<SymmetricTensor<2,dim> >()),
      cell (NULL)
    {}



    template <int dim>
    MaterialModelOutputs<dim>::MaterialModelOutputs(const unsigned int n_points,
                                                    const unsigned int n_comp)
      :
      viscosities(n_points, aspect::Utilities::signaling_nan<double>()),
      stress_strain_directors(n_points, dealii::identity_tensor<dim> ()),
      densities(n_points, aspect::Utilities::signaling_nan<double>()),
      thermal_expansion_coefficients(n_points, aspect::Utilities::signaling_nan<double>()),
      specific_heat(n_points, aspect::Utilities::signaling_nan<double>()),
      thermal_conductivities(n_points, aspect::Utilities::signaling_nan<double>()),
      compressibilities(n_points, aspect::Utilities::signaling_nan<double>()),
      entropy_derivative_pressure(n_points, aspect::Utilities::signaling_nan<double>()),
      entropy_derivative_temperature(n_points, aspect::Utilities::signaling_nan<double>()),
      reaction_terms(n_points, std::vector<double>(n_comp, aspect::Utilities::signaling_nan<double>()))
    {}


    template <int dim>
    double
    InterfaceCompatibility<dim>::
    entropy_derivative (const double,
                        const double,
                        const std::vector<double> &,
                        const Point<dim> &,
                        const NonlinearDependence::Dependence) const
    {
      return 0.0;
    }


    template <int dim>
    double
    InterfaceCompatibility<dim>::
    reaction_term (const double,
                   const double,
                   const std::vector<double> &,
                   const Point<dim> &,
                   const unsigned int) const
    {
      return 0.0;
    }


    template <int dim>
    void
    InterfaceCompatibility<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                                          MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          // as documented, if the strain rate array is empty, then do not compute the
          // viscosities
          if (in.strain_rate.size() > 0)
            out.viscosities[i]                  = viscosity                     (in.temperature[i], in.pressure[i], in.composition[i], in.strain_rate[i], in.position[i]);

          out.densities[i]                      = density                       (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.thermal_expansion_coefficients[i] = thermal_expansion_coefficient (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.specific_heat[i]                  = specific_heat                 (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.thermal_conductivities[i]         = thermal_conductivity          (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.compressibilities[i]              = compressibility               (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.entropy_derivative_pressure[i]    = entropy_derivative            (in.temperature[i], in.pressure[i], in.composition[i], in.position[i], NonlinearDependence::pressure);
          out.entropy_derivative_temperature[i] = entropy_derivative            (in.temperature[i], in.pressure[i], in.composition[i], in.position[i], NonlinearDependence::temperature);
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c]            = reaction_term                 (in.temperature[i], in.pressure[i], in.composition[i], in.position[i], c);
        }
    }


    namespace MaterialAveraging
    {
      std::string get_averaging_operation_names ()
      {
        return "none|arithmetic average|harmonic average|geometric average|pick largest|project to Q1|log average";
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
        else
          AssertThrow (false,
                       ExcMessage ("The value <" + s + "> for a material "
                                   "averaging operation is not one of the "
                                   "valid values."));

        return none;
      }

      namespace
      {
        // Do the requested averaging operation for one array. The
        // projection matrix argument is only used if the operation
        // chosen is project_to_Q1
        void average (const AveragingOperation  operation,
                      const FullMatrix<double>      &projection_matrix,
                      const FullMatrix<double>      &expansion_matrix,
                      std::vector<double>           &values_out)
        {
          // if an output field has not been filled (because it was
          // not requested), then simply do nothing -- no harm no foul
          if (values_out.size() == 0)
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
                double prod = 1;
                for (unsigned int i=0; i<N; ++i)
                  {
                    Assert (values_out[i] >= 0,
                            ExcMessage ("Computing the geometric average "
                                        "only makes sense for non-negative "
                                        "quantities."));
                    prod *= values_out[i];
                  }

                const double average = std::pow (prod, 1./N);
                for (unsigned int i=0; i<N; ++i)
                  values_out[i] = average;
                break;
              }

              case pick_largest:
              {
                double max = -std::numeric_limits<double>::max();
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

                double max = -std::numeric_limits<double>::max();
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
        }


        /**
         * Given a quadrature formula, compute a matrices $E, M^{-1}F$
         * representing a linear operator in the following way: Let
         * there be a vector $F$ with $N$ elements where the elements
         * are data stored at each of the $N$ quadrature points. Then
         * project this data into a Q1 space and evaluate this
         * projection at the quadrature points. This operator can be
         * expressed in the following way where $P=2^{dim}$ is the
         * number of degrees of freedom of the $Q_1$ element:
         *
         * Let $y$ be the input vector with $N$ elements. Then let
         * $F$ be the $P \times N$ matrix so that
         * @f{align}
         *   F_{iq} = \varphi_i(x_q) |J(x_q)| w_q
         * @f}
         * where $\varphi_i$ are the $Q_1$ shape functions, $x_q$ are
         * the quadrature points, $J$ is the Jacobian matrix of the
         * mapping from reference to real cell, and $w_q$ are the
         * quadrature weights.
         *
         * Next, let $M_{ij}$ be the $P\times P$ mass matrix on the
         * $Q_1$ space. Then $M^{-1}Fy$ corresponds to the projection
         * of $y$ into the $Q_1$ space (or, more correctly, it
         * corresponds to the nodal values of this projection).
         *
         * Finally, let $E_{qi} = \varphi_i(x_q)$ be the evaluation
         * operation of shape functions at quadrature points.
         *
         * Then, the operation $X=EM^{-1}F$ is the operation we seek.
         * This function computes the matrices E and M^{-1}F.
         */
        template <int dim>
        void compute_projection_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                        const Quadrature<dim>   &quadrature_formula,
                                        const Mapping<dim>      &mapping,
                                        FullMatrix<double>      &projection_matrix,
                                        FullMatrix<double>      &expansion_matrix)
        {
          static FE_Q<dim> fe(1);
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
      }


      template <int dim>
      void average (const AveragingOperation operation,
                    const typename DoFHandler<dim>::active_cell_iterator &cell,
                    const Quadrature<dim>         &quadrature_formula,
                    const Mapping<dim>            &mapping,
                    MaterialModelOutputs<dim>          &values_out)
      {
        FullMatrix<double> projection_matrix;
        FullMatrix<double> expansion_matrix;

        if (operation == project_to_Q1)
          {
            projection_matrix.reinit (quadrature_formula.size(),
                                      quadrature_formula.size());
            compute_projection_matrix (cell,
                                       quadrature_formula,
                                       mapping,
                                       projection_matrix,
                                       expansion_matrix);
          }

        average (operation, projection_matrix, expansion_matrix, values_out.viscosities);
        average (operation, projection_matrix, expansion_matrix,
                 values_out.densities);
        average (operation, projection_matrix, expansion_matrix,
                 values_out.thermal_expansion_coefficients);
        average (operation, projection_matrix, expansion_matrix,
                 values_out.specific_heat);
        average (operation, projection_matrix, expansion_matrix,
                 values_out.compressibilities);
        average (operation, projection_matrix, expansion_matrix,
                 values_out.entropy_derivative_pressure);
        average (operation, projection_matrix, expansion_matrix,
                 values_out.entropy_derivative_temperature);

        // the reaction terms are unfortunately stored in reverse
        // indexing. it's also not quite clear whether these should
        // really be averaged, so avoid this for now
      }

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
      std::list<internal::Plugins::PluginList<MaterialModel::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::Interface<2> >::plugins = 0;

      template <>
      std::list<internal::Plugins::PluginList<MaterialModel::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::Interface<3> >::plugins = 0;
    }
  }

  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template class InterfaceCompatibility<dim>; \
  \
  template \
  void \
  register_material_model<dim> (const std::string &, \
                                const std::string &, \
                                void ( *) (ParameterHandler &), \
                                Interface<dim> *( *) ()); \
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
  Interface<dim> * \
  create_material_model<dim> (const std::string &model_name); \
  \
  template \
  Interface<dim> * \
  create_material_model<dim> (ParameterHandler &prm); \
  \
  template struct MaterialModelInputs<dim>; \
  \
  template struct MaterialModelOutputs<dim>; \
  \
  template class AdditionalMaterialOutputs<dim>; \
  \
  namespace MaterialAveraging \
  { \
    template                \
    void average (const AveragingOperation operation, \
                  const DoFHandler<dim>::active_cell_iterator &cell, \
                  const Quadrature<dim>     &quadrature_formula, \
                  const Mapping<dim>        &mapping, \
                  MaterialModelOutputs<dim>      &values_out); \
  }


    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
