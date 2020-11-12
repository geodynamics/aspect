/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#include <aspect/lateral_averaging.h>
#include <aspect/material_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/two_merged_boxes.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>



namespace aspect
{
  /**
   * This namespace contains all the implemented functors. They are used to
   * compute various properties of the solution that will be laterally averaged.
   */
  namespace
  {
    template <int dim>
    class FunctorDepthAverageField: public internal::FunctorBase<dim>
    {
      public:
        FunctorDepthAverageField(const FEValuesExtractors::Scalar &field)
          : field_(field)
        {}

        void operator()(const MaterialModel::MaterialModelInputs<dim> &,
                        const MaterialModel::MaterialModelOutputs<dim> &,
                        const FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution,
                        std::vector<double> &output) override
        {
          fe_values[field_].get_function_values (solution, output);
        }

        const FEValuesExtractors::Scalar field_;
    };



    template <int dim>
    class FunctorDepthAverageViscosity: public internal::FunctorBase<dim>
    {
      public:
        bool need_material_properties() const override
        {
          return true;
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &,
                        const MaterialModel::MaterialModelOutputs<dim> &out,
                        const FEValues<dim> &,
                        const LinearAlgebra::BlockVector &,
                        std::vector<double> &output) override
        {
          output = out.viscosities;
        }
    };



    template <int dim>
    class FunctorDepthAverageVelocityMagnitude: public internal::FunctorBase<dim>
    {
      public:
        FunctorDepthAverageVelocityMagnitude(const FEValuesExtractors::Vector &field,
                                             bool convert_to_years)
          : field_(field), convert_to_years_(convert_to_years)
        {}

        void setup(const unsigned int q_points) override
        {
          velocity_values.resize(q_points);
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &,
                        const MaterialModel::MaterialModelOutputs<dim> &,
                        const FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution,
                        std::vector<double> &output) override
        {
          fe_values[field_].get_function_values (solution, velocity_values);
          for (unsigned int q=0; q<output.size(); ++q)
            output[q] = std::sqrt( velocity_values[q] * velocity_values[q] ) *
                        (convert_to_years_ ? year_in_seconds : 1.0);
        }

        std::vector<Tensor<1,dim> > velocity_values;
        const FEValuesExtractors::Vector field_;
        const bool convert_to_years_;
    };



    template <int dim>
    class FunctorDepthAverageSinkingVelocity: public internal::FunctorBase<dim>
    {
      public:
        FunctorDepthAverageSinkingVelocity(const FEValuesExtractors::Vector &field,
                                           const GravityModel::Interface<dim> *gravity,
                                           bool convert_to_years)
          : field_(field),
            gravity_(gravity),
            convert_to_years_(convert_to_years)
        {}

        bool need_material_properties() const override
        {
          // this is needed because we want to access in.position in operator()
          return true;
        }

        void setup(const unsigned int q_points) override
        {
          velocity_values.resize(q_points);
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &in,
                        const MaterialModel::MaterialModelOutputs<dim> &,
                        const FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution,
                        std::vector<double> &output) override
        {
          fe_values[field_].get_function_values (solution, velocity_values);
          for (unsigned int q=0; q<output.size(); ++q)
            {
              const Tensor<1,dim> g = gravity_->gravity_vector(in.position[q]);
              const Tensor<1,dim> vertical = (g.norm() > 0 ? g/g.norm() : Tensor<1,dim>());

              output[q] = std::fabs(std::min(0.0, velocity_values[q] * vertical))
                          * (convert_to_years_ ? year_in_seconds : 1.0);
            }
        }

        std::vector<Tensor<1,dim> > velocity_values;
        const FEValuesExtractors::Vector field_;
        const GravityModel::Interface<dim> *gravity_;
        const bool convert_to_years_;
    };



    template <int dim>
    class FunctorDepthAverageVsVp: public internal::FunctorBase<dim>
    {
      public:
        FunctorDepthAverageVsVp(bool vs)
          : vs_(vs)
        {}

        bool need_material_properties() const override
        {
          return true;
        }

        void
        create_additional_material_model_outputs (const unsigned int n_points,
                                                  MaterialModel::MaterialModelOutputs<dim> &outputs) const override
        {
          outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &,
                        const MaterialModel::MaterialModelOutputs<dim> &out,
                        const FEValues<dim> &,
                        const LinearAlgebra::BlockVector &,
                        std::vector<double> &output) override
        {
          const MaterialModel::SeismicAdditionalOutputs<dim> *seismic_outputs
            = out.template get_additional_output<const MaterialModel::SeismicAdditionalOutputs<dim> >();

          Assert(seismic_outputs != nullptr,ExcInternalError());

          if (vs_)
            for (unsigned int q=0; q<output.size(); ++q)
              output[q] = seismic_outputs->vs[q];
          else
            for (unsigned int q=0; q<output.size(); ++q)
              output[q] = seismic_outputs->vp[q];
        }

        bool vs_;
    };



    template <int dim>
    class FunctorDepthAverageVerticalHeatFlux: public internal::FunctorBase<dim>
    {
      public:
        FunctorDepthAverageVerticalHeatFlux(const FEValuesExtractors::Vector &velocity_field,
                                            const FEValuesExtractors::Scalar &temperature_field,
                                            const GravityModel::Interface<dim> *gm)
          : velocity_field_(velocity_field),
            temperature_field_(temperature_field),
            gravity_model(gm)
        {}

        bool need_material_properties() const override
        {
          return true;
        }

        void setup(const unsigned int q_points) override
        {
          velocity_values.resize(q_points);
          temperature_values.resize(q_points);
          temperature_gradients.resize(q_points);
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &in,
                        const MaterialModel::MaterialModelOutputs<dim> &out,
                        const FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution,
                        std::vector<double> &output) override
        {
          fe_values[velocity_field_].get_function_values (solution, velocity_values);
          fe_values[temperature_field_].get_function_values (solution, temperature_values);
          fe_values[temperature_field_].get_function_gradients (solution, temperature_gradients);

          for (unsigned int q=0; q<output.size(); ++q)
            {
              const Tensor<1,dim> gravity = gravity_model->gravity_vector(in.position[q]);
              const Tensor<1,dim> vertical = -gravity/( gravity.norm() != 0.0 ?
                                                        gravity.norm() : 1.0 );
              const double advective_flux = (velocity_values[q] * vertical) * in.temperature[q] *
                                            out.densities[q]*out.specific_heat[q];
              const double conductive_flux = -(temperature_gradients[q]*vertical) *
                                             out.thermal_conductivities[q];
              output[q] = advective_flux + conductive_flux;
            }
        }

        const FEValuesExtractors::Vector velocity_field_;
        const FEValuesExtractors::Scalar temperature_field_;
        const GravityModel::Interface<dim> *gravity_model;
        std::vector<Tensor<1,dim> > velocity_values;
        std::vector<Tensor<1,dim> > temperature_gradients;
        std::vector<double> temperature_values;
    };



    template <int dim>
    class FunctorDepthAverageVerticalMassFlux: public internal::FunctorBase<dim>
    {
      public:
        FunctorDepthAverageVerticalMassFlux(const FEValuesExtractors::Vector &velocity_field,
                                            const GravityModel::Interface<dim> *gm)
          : velocity_field_(velocity_field),
            gravity_model(gm)
        {}

        bool need_material_properties() const override
        {
          return true;
        }

        void setup(const unsigned int q_points) override
        {
          velocity_values.resize(q_points);
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &in,
                        const MaterialModel::MaterialModelOutputs<dim> &out,
                        const FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution,
                        std::vector<double> &output) override
        {
          fe_values[velocity_field_].get_function_values (solution, velocity_values);

          for (unsigned int q=0; q<output.size(); ++q)
            {
              const Tensor<1,dim> gravity = gravity_model->gravity_vector(in.position[q]);
              const Tensor<1,dim> vertical = -gravity/( gravity.norm() != 0.0 ?
                                                        gravity.norm() : 1.0 );

              output[q] = std::fabs(velocity_values[q] * vertical) * out.densities[q];
            }
        }

        const FEValuesExtractors::Vector velocity_field_;
        const GravityModel::Interface<dim> *gravity_model;
        std::vector<Tensor<1,dim> > velocity_values;
    };



    template <int dim>
    class FunctorDepthAverageFieldMass: public internal::FunctorBase<dim>
    {
      public:
        FunctorDepthAverageFieldMass(const FEValuesExtractors::Scalar &field)
          : field_(field)
        {}

        bool need_material_properties() const override
        {
          return true;
        }

        void setup(const unsigned int q_points) override
        {
          field_values.resize(q_points);
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &,
                        const MaterialModel::MaterialModelOutputs<dim> &out,
                        const FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution,
                        std::vector<double> &output) override
        {
          fe_values[field_].get_function_values (solution, field_values);

          for (unsigned int q=0; q<output.size(); ++q)
            output[q] = field_values[q] * out.densities[q];
        }

        const FEValuesExtractors::Scalar field_;
        std::vector<double> field_values;
    };
  }

  namespace internal
  {
    template <int dim>
    bool
    FunctorBase<dim>::need_material_properties() const
    {
      return false;
    }



    template <int dim>
    FunctorBase<dim>::~FunctorBase()
    {}



    template <int dim>
    void
    FunctorBase<dim>::create_additional_material_model_outputs (const unsigned int /*n_points*/,
                                                                MaterialModel::MaterialModelOutputs<dim> &/*outputs*/) const
    {}



    template <int dim>
    void
    FunctorBase<dim>::setup(const unsigned int /*q_points*/)
    {}

    template <int dim>
    QAnisotropic<dim>
    get_quadrature_formula(const unsigned int lateral_quadrature_degree,
                           const unsigned int depth_dimension);

    template <>
    QAnisotropic<3>
    get_quadrature_formula(const unsigned int lateral_quadrature_degree,
                           const unsigned int depth_dimension)
    {
      const unsigned int dim = 3;

      if (depth_dimension == dim)
        return QAnisotropic<dim> (QGauss<1>(lateral_quadrature_degree),
                                  QGauss<1>(lateral_quadrature_degree),
                                  QIterated<1>(QMidpoint<1>(),10));
      else if (depth_dimension == 1)
        return QAnisotropic<dim> (QIterated<1>(QMidpoint<1>(),10),
                                  QGauss<1>(lateral_quadrature_degree),
                                  QGauss<1>(lateral_quadrature_degree));

      return QAnisotropic<dim>(QGauss<1>(lateral_quadrature_degree),
                               QGauss<1>(lateral_quadrature_degree),
                               QIterated<1>(QMidpoint<1>(),10));
    }

    template <>
    QAnisotropic<2>
    get_quadrature_formula(const unsigned int lateral_quadrature_degree,
                           const unsigned int depth_dimension)
    {
      const unsigned int dim = 2;

      if (depth_dimension == dim)
        return QAnisotropic<dim> (QGauss<1>(lateral_quadrature_degree),
                                  QIterated<1>(QMidpoint<1>(),10));
      else if (depth_dimension == 1)
        return QAnisotropic<dim> (QIterated<1>(QMidpoint<1>(),10),
                                  QGauss<1>(lateral_quadrature_degree));

      return QAnisotropic<dim>(QGauss<1>(lateral_quadrature_degree),
                               QIterated<1>(QMidpoint<1>(),10));
    }
  }



  template <int dim>
  std::vector<std::vector<double> >
  LateralAveraging<dim>::compute_lateral_averages(const std::vector<double> &depth_bounds,
                                                  std::vector<std::unique_ptr<internal::FunctorBase<dim> > > &functors) const
  {
    Assert (functors.size() > 0,
            ExcMessage ("To call this function, you need to request a positive "
                        "number of properties to compute."));
    Assert (depth_bounds.size() > 1,
            ExcMessage ("To call this function, you need to request at least two "
                        "depth boundaries."));
    Assert(std::is_sorted(depth_bounds.begin(),depth_bounds.end()),
           ExcMessage ("To call this function the depth boundaries need to be ordered "
                       "with increasing depth."));

    const unsigned int n_properties = functors.size();
    const unsigned int n_slices = depth_bounds.size()-1;

    std::vector<std::vector<double> > values(n_properties,
                                             std::vector<double>(n_slices,0.0));
    std::vector<double> volume(n_slices,0.0);

    // We would like to use a quadrature formula that is appropriately accurate laterally,
    // but has a higher resolution in depth (to have values for all depth slices, even if there
    // are adaptively coarsened cells that are much bigger than individual slices).
    // For that we need to know the depth direction in the unit cell coordinate system, which
    // is only unique (= the same for all cells) in some geometries. In these geometries we
    // can optimize the quadrature, otherwise we need to use a high-resolution quadrature in
    // all directions, which is more expensive.
    // The Chunk geometry model has depth as first dimension (radius, lon, lat),
    // all others with unique direction have it last; however, the Chunk and
    // EllipsoidalChunk geometry have not been successfully tested with
    // the lower quadrature, so we leave it at the conservative quadrature for now.

    unsigned int geometry_unique_depth_direction;
    if (Plugins::plugin_type_matches<GeometryModel::Box<dim> >(this->get_geometry_model()) ||
        Plugins::plugin_type_matches<GeometryModel::SphericalShell<dim> >(this->get_geometry_model()) ||
        Plugins::plugin_type_matches<GeometryModel::TwoMergedBoxes<dim> >(this->get_geometry_model()))
      geometry_unique_depth_direction = dim;
    else if (Plugins::plugin_type_matches<GeometryModel::Chunk<dim> >(this->get_geometry_model()) ||
             Plugins::plugin_type_matches<GeometryModel::EllipsoidalChunk<dim> >(this->get_geometry_model()))
      geometry_unique_depth_direction = numbers::invalid_unsigned_int;
    else
      geometry_unique_depth_direction = numbers::invalid_unsigned_int;

    const unsigned int max_fe_degree = std::max(this->introspection().polynomial_degree.velocities,
                                                std::max(this->introspection().polynomial_degree.temperature,
                                                         this->introspection().polynomial_degree.compositional_fields));

    // We want to integrate over a polynomial of degree p = max_fe_degree, for which we
    // need a quadrature of at least q, with p <= 2q-1 --> q >= (p+1)/2
    const unsigned int lateral_quadrature_degree = static_cast<unsigned int>(std::ceil((max_fe_degree+1.0)/2.0));

    std::unique_ptr<Quadrature<dim> > quadrature_formula;
    if (geometry_unique_depth_direction != numbers::invalid_unsigned_int)
      quadrature_formula = std_cxx14::make_unique<Quadrature<dim> >(internal::get_quadrature_formula<dim>(lateral_quadrature_degree,
                                                                    geometry_unique_depth_direction));
    else
      quadrature_formula = std_cxx14::make_unique<Quadrature<dim> >(QIterated<dim>(QMidpoint<1>(),10));

    const unsigned int n_q_points = quadrature_formula->size();

    FEValues<dim> fe_values (this->get_mapping(),
                             this->get_fe(),
                             *quadrature_formula,
                             update_values | update_gradients | update_quadrature_points | update_JxW_values);

    std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),
                                                          std::vector<double> (n_q_points));
    std::vector<std::vector<double> > output_values(n_properties,
                                                    std::vector<double>(n_q_points));

    MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                               this->n_compositional_fields());
    MaterialModel::MaterialModelOutputs<dim> out(n_q_points,
                                                 this->n_compositional_fields());

    bool functors_need_material_output = false;
    for (unsigned int i=0; i<n_properties; ++i)
      {
        functors[i]->setup(n_q_points);
        if (functors[i]->need_material_properties())
          functors_need_material_output = true;

        functors[i]->create_additional_material_model_outputs(n_q_points,out);
      }

    for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);

          if (functors_need_material_output)
            {
              // get the material properties at each quadrature point if necessary
              in.reinit(fe_values,
                        cell,
                        this->introspection(),
                        this->get_solution());
              this->get_material_model().evaluate(in, out);
            }

          for (unsigned int i = 0; i < n_properties; ++i)
            (*functors[i])(in, out, fe_values, this->get_solution(), output_values[i]);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double depth = this->get_geometry_model().depth(fe_values.quadrature_point(q));

              if (depth < depth_bounds.front() || depth > depth_bounds.back())
                continue;

              // This makes sure depth == front() and depth == back() are handled correctly.
              // lower_bound returns the first layer boundary larger than depth, the correct
              // layer index is then one less than this (except for depth == depth_bounds[0],
              // in which case the depth_bounds index is also the layer_index, namely 0).
              unsigned int layer_index = std::distance(depth_bounds.begin(),
                                                       std::lower_bound(depth_bounds.begin(),depth_bounds.end(),depth));
              if (layer_index > 0)
                layer_index -= 1;

              Assert(layer_index<n_slices, ExcInternalError());

              for (unsigned int i = 0; i < n_properties; ++i)
                values[i][layer_index] += output_values[i][q] * fe_values.JxW(q);

              volume[layer_index] += fe_values.JxW(q);
            }
        }

    std::vector<double> volume_all(n_slices);
    Utilities::MPI::sum(volume, this->get_mpi_communicator(), volume_all);

    bool print_under_res_warning=false;
    for (unsigned int property=0; property<n_properties; ++property)
      {
        std::vector<double> values_all(n_slices);
        Utilities::MPI::sum(values[property], this->get_mpi_communicator(), values_all);

        for (unsigned int i=0; i<n_slices; ++i)
          {
            if (volume_all[i] > 0.0)
              {
                values[property][i] = values_all[i] / (static_cast<double>(volume_all[i]));
              }
            else
              {
                print_under_res_warning = true;
                // Output nan if no quadrature points in depth block
                values[property][i] = std::numeric_limits<double>::quiet_NaN();
              }
          }
      }

    if (print_under_res_warning)
      {
        this->get_pcout() << std::endl
                          << "**** Warning: When computing depth averages, there is at least one depth band"
                          << std::endl
                          << "     that does not have any quadrature points in it."
                          << std::endl
                          << "     Consider reducing the number of depth layers for averaging."
                          << std::endl << std::endl;
      }

    return values;
  }



  template <int dim>
  void LateralAveraging<dim>::get_temperature_averages(std::vector<double> &values) const
  {
    values = compute_lateral_averages(values.size(),
                                      std::vector<std::string>(1,"temperature"))[0];
  }



  template <int dim>
  void LateralAveraging<dim>::get_composition_averages(const unsigned int c,
                                                       std::vector<double> &values) const
  {
    values = compute_lateral_averages(values.size(),
                                      std::vector<std::string>(1,"C_" + Utilities::int_to_string(c)))[0];
  }



  template <int dim>
  void LateralAveraging<dim>::get_viscosity_averages(std::vector<double> &values) const
  {
    values = compute_lateral_averages(values.size(),
                                      std::vector<std::string>(1,"viscosity"))[0];
  }



  template <int dim>
  void LateralAveraging<dim>::get_velocity_magnitude_averages(std::vector<double> &values) const
  {
    values = compute_lateral_averages(values.size(),
                                      std::vector<std::string>(1,"velocity_magnitude"))[0];
  }



  template <int dim>
  void LateralAveraging<dim>::get_sinking_velocity_averages(std::vector<double> &values) const
  {
    values = compute_lateral_averages(values.size(),
                                      std::vector<std::string>(1,"sinking_velocity"))[0];
  }



  template <int dim>
  void LateralAveraging<dim>::get_Vs_averages(std::vector<double> &values) const
  {
    values = compute_lateral_averages(values.size(),
                                      std::vector<std::string>(1,"Vs"))[0];
  }



  template <int dim>
  void LateralAveraging<dim>::get_Vp_averages(std::vector<double> &values) const
  {
    values = compute_lateral_averages(values.size(),
                                      std::vector<std::string>(1,"Vp"))[0];
  }



  template <int dim>
  void LateralAveraging<dim>::get_vertical_heat_flux_averages(std::vector<double> &values) const
  {
    values = compute_lateral_averages(values.size(),
                                      std::vector<std::string>(1,"vertical_heat_flux"))[0];
  }



  template <int dim>
  void LateralAveraging<dim>::get_vertical_mass_flux_averages(std::vector<double> &values) const
  {
    values = compute_lateral_averages(values.size(),
                                      std::vector<std::string>(1,"vertical_mass_flux"))[0];
  }



  template <int dim>
  std::vector<std::vector<double> >
  LateralAveraging<dim>::get_averages(const unsigned int n_slices,
                                      const std::vector<std::string> &property_names) const
  {
    return compute_lateral_averages(n_slices, property_names);
  }



  template <int dim>
  std::vector<std::vector<double> >
  LateralAveraging<dim>::compute_lateral_averages(const unsigned int n_slices,
                                                  const std::vector<std::string> &property_names) const
  {
    const double maximal_depth = this->get_geometry_model().maximal_depth();
    std::vector<double> depth_bounds(n_slices+1, 0.0);

    // Leave index 0 at 0.0, and generate an increasing range of equidistant depth bounds
    for (unsigned int i=1; i<depth_bounds.size(); ++i)
      depth_bounds[i] = depth_bounds[i-1] + maximal_depth / n_slices;

    return compute_lateral_averages(depth_bounds, property_names);
  }



  template <int dim>
  std::vector<std::vector<double> >
  LateralAveraging<dim>::compute_lateral_averages(const std::vector<double> &depth_thresholds,
                                                  const std::vector<std::string> &property_names) const
  {
    std::vector<std::unique_ptr<internal::FunctorBase<dim> > > functors;
    for (unsigned int property_index=0; property_index<property_names.size(); ++property_index)
      {
        if (property_names[property_index] == "temperature")
          {
            functors.push_back(std_cxx14::make_unique<FunctorDepthAverageField<dim>>
                               (this->introspection().extractors.temperature));
          }
        else if (this->introspection().compositional_name_exists(property_names[property_index]))
          {
            const unsigned int c =
              this->introspection().compositional_index_for_name(property_names[property_index]);

            functors.push_back(std_cxx14::make_unique<FunctorDepthAverageField<dim>> (
                                 this->introspection().extractors.compositional_fields[c]));
          }
        else if (property_names[property_index] == "velocity_magnitude")
          {
            functors.push_back(std_cxx14::make_unique<FunctorDepthAverageVelocityMagnitude<dim>>
                               (this->introspection().extractors.velocities,
                                this->convert_output_to_years()));
          }
        else if (property_names[property_index] == "sinking_velocity")
          {
            functors.push_back(std_cxx14::make_unique<FunctorDepthAverageSinkingVelocity<dim>>
                               (this->introspection().extractors.velocities,
                                &this->get_gravity_model(),
                                this->convert_output_to_years()));
          }
        else if (property_names[property_index] == "Vs")
          {
            functors.push_back(std_cxx14::make_unique<FunctorDepthAverageVsVp<dim>> (true /* Vs */));
          }
        else if (property_names[property_index] == "Vp")
          {
            functors.push_back(std_cxx14::make_unique<FunctorDepthAverageVsVp<dim>> (false /* Vp */));
          }
        else if (property_names[property_index] == "viscosity")
          {
            functors.push_back(std_cxx14::make_unique<FunctorDepthAverageViscosity<dim>>());
          }
        else if (property_names[property_index] == "vertical_heat_flux")
          {
            functors.push_back(std_cxx14::make_unique<FunctorDepthAverageVerticalHeatFlux<dim>>
                               (this->introspection().extractors.velocities,
                                this->introspection().extractors.temperature,
                                &this->get_gravity_model()));
          }
        else if (property_names[property_index] == "vertical_mass_flux")
          {
            functors.push_back(std_cxx14::make_unique<FunctorDepthAverageVerticalMassFlux<dim>>
                               (this->introspection().extractors.velocities,
                                &this->get_gravity_model()));
          }
        else if (this->introspection().compositional_name_exists(property_names[property_index].substr(0, property_names[property_index].size()-5)) &&
                 property_names[property_index].substr(property_names[property_index].size()-5) == "_mass")
          {
            const unsigned int c =
              this->introspection().compositional_index_for_name(property_names[property_index].substr(0, property_names[property_index].size()-5));

            functors.push_back(std_cxx14::make_unique<FunctorDepthAverageFieldMass<dim>> (
                                 this->introspection().extractors.compositional_fields[c]));
          }
        else
          {
            AssertThrow(false,
                        ExcMessage("The lateral averaging scheme was asked to average the property "
                                   "named <" + property_names[property_index] + ">, but it does not know how "
                                   "to do that. There is no functor implemented that computes this property."));
          }
      }

    // Now compute values for all selected properties.
    return compute_lateral_averages(depth_thresholds, functors);
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template class LateralAveraging<dim>; \
  namespace internal \
  { \
    template class FunctorBase<dim>; \
  }

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
