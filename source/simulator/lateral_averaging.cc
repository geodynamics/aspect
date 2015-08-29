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

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <aspect/lateral_averaging.h>
#include <aspect/material_model/interface.h>
#include <aspect/gravity_model/interface.h>


namespace aspect
{

  template <int dim>
  template<class FUNCTOR>
  void LateralAveraging<dim>::compute_lateral_average(std::vector<double> &values,
                                                      FUNCTOR &fctr) const
  {
    Assert (values.size() > 0,
            ExcMessage ("To call this function, you need to request a positive "
                        "number of depth slices."));
    const unsigned int num_slices = values.size();
    std::vector<double> volume(num_slices);

    // this yields 10^dim quadrature points evenly distributed in the interior of the cell.
    // We avoid points on the faces, as they would be counted more than once.
    const QIterated<dim> quadrature_formula (QMidpoint<1>(),
                                             10);
    const unsigned int n_q_points = quadrature_formula.size();
    const double max_depth = this->get_geometry_model().maximal_depth();

    FEValues<dim> fe_values (this->get_mapping(),
                             this->get_fe(),
                             quadrature_formula,
                             update_values | update_gradients | update_quadrature_points);

    std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (n_q_points));

    std::vector<double> output_values(quadrature_formula.size());

    typename DoFHandler<dim>::active_cell_iterator
    cell = this->get_dof_handler().begin_active(),
    endc = this->get_dof_handler().end();

    MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                               this->n_compositional_fields());
    MaterialModel::MaterialModelOutputs<dim> out(n_q_points,
                                                 this->n_compositional_fields());

    fctr.setup(quadrature_formula.size());

    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          if (fctr.need_material_properties())
            {
              fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                        in.pressure);
              fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                           in.temperature);
              fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                          in.velocity);
              fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(),
                  in.strain_rate);
              fe_values[this->introspection().extractors.pressure].get_function_gradients (this->get_solution(),
                                                                                           in.pressure_gradient);
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                    composition_values[c]);

              for (unsigned int i=0; i<n_q_points; ++i)
                {
                  in.position[i] = fe_values.quadrature_point(i);
                  for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                    in.composition[i][c] = composition_values[c][i];
                }
              in.cell = &cell;

              this->get_material_model().evaluate(in, out);
            }

          fctr(in, out, fe_values, this->get_solution(), output_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double depth = this->get_geometry_model().depth(fe_values.quadrature_point(q));
              const unsigned int idx = static_cast<unsigned int>((depth*num_slices)/max_depth);
              Assert(idx<num_slices, ExcInternalError());

              values[idx] += output_values[q] * fe_values.JxW(q);
              volume[idx] += fe_values.JxW(q);
            }
        }

    std::vector<double> values_all(num_slices);
    std::vector<double> volume_all(num_slices);
    Utilities::MPI::sum(volume, this->get_mpi_communicator(), volume_all);
    Utilities::MPI::sum(values, this->get_mpi_communicator(), values_all);

    for (unsigned int i=0; i<num_slices; ++i)
      values[i] = values_all[i] / (static_cast<double>(volume_all[i])+1e-20);
  }

  namespace
  {
    template <int dim>
    class FunctorDepthAverageField
    {
      public:
        FunctorDepthAverageField(const FEValuesExtractors::Scalar &field)
          : field_(field) {}

        bool need_material_properties() const
        {
          return false;
        }

        void setup(const unsigned int)
        {
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &,
                        const MaterialModel::MaterialModelOutputs<dim> &,
                        FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution,
                        std::vector<double> &output)
        {
          fe_values[field_].get_function_values (solution, output);
        }

        const FEValuesExtractors::Scalar field_;
    };
  }

  template <int dim>
  void LateralAveraging<dim>::get_temperature_averages(std::vector<double> &values) const
  {
    FEValuesExtractors::Scalar field = this->introspection().extractors.temperature;
    get_field_averages( field, values );
  }

  template <int dim>
  void LateralAveraging<dim>::get_composition_averages( const unsigned int c, std::vector<double> &values) const
  {
    FEValuesExtractors::Scalar field = this->introspection().extractors.compositional_fields[c];
    get_field_averages( field, values );
  }

  template <int dim>
  void LateralAveraging<dim>::get_field_averages(const FEValuesExtractors::Scalar &field,
                                                 std::vector<double> &values) const
  {
    FunctorDepthAverageField<dim> f(field);
    compute_lateral_average(values, f);
  }

  namespace
  {
    template <int dim>
    class FunctorDepthAverageViscosity
    {
      public:
        bool need_material_properties() const
        {
          return true;
        }

        void setup(const unsigned int)
        {}

        void operator()(const MaterialModel::MaterialModelInputs<dim> &,
                        const MaterialModel::MaterialModelOutputs<dim> &out,
                        FEValues<dim> &,
                        const LinearAlgebra::BlockVector &,
                        std::vector<double> &output)
        {
          output = out.viscosities;
        }
    };
  }

  template <int dim>
  void LateralAveraging<dim>::get_viscosity_averages(std::vector<double> &values) const
  {
    FunctorDepthAverageViscosity<dim> f;
    compute_lateral_average(values, f);
  }


  namespace
  {
    template <int dim>
    class FunctorDepthAverageVelocityMagnitude
    {
      public:
        FunctorDepthAverageVelocityMagnitude(const FEValuesExtractors::Vector &field,
                                             bool convert_to_years)
          : field_(field), convert_to_years_(convert_to_years) {}

        bool need_material_properties() const
        {
          return false;
        }

        void setup(const unsigned int q_points)
        {
          velocity_values.resize(q_points);
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &,
                        const MaterialModel::MaterialModelOutputs<dim> &,
                        FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution,
                        std::vector<double> &output)
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
  }


  template <int dim>
  void LateralAveraging<dim>::get_velocity_magnitude_averages(std::vector<double> &values) const
  {
    FunctorDepthAverageVelocityMagnitude<dim> f(this->introspection().extractors.velocities,
                                                this->convert_output_to_years());

    compute_lateral_average(values, f);
  }

  namespace
  {
    template <int dim>
    class FunctorDepthAverageSinkingVelocity
    {
      public:
        FunctorDepthAverageSinkingVelocity(const FEValuesExtractors::Vector &field,
                                           const GravityModel::Interface<dim> *gravity,
                                           bool convert_to_years)
          : field_(field),
            gravity_(gravity),
            convert_to_years_(convert_to_years) {}

        bool need_material_properties() const
        {
          // this is needed because we want to access in.position in operator()
          return true;
        }

        void setup(const unsigned int q_points)
        {
          velocity_values.resize(q_points);
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &in,
                        const MaterialModel::MaterialModelOutputs<dim> &,
                        FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution,
                        std::vector<double> &output)
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
  }


  template <int dim>
  void LateralAveraging<dim>::get_sinking_velocity_averages(std::vector<double> &values) const
  {
    FunctorDepthAverageSinkingVelocity<dim> f(this->introspection().extractors.velocities,
                                              &this->get_gravity_model(),
                                              this->convert_output_to_years());

    compute_lateral_average(values, f);
  }

  namespace
  {
    template <int dim>
    class FunctorDepthAverageVsVp
    {
      public:
        FunctorDepthAverageVsVp(const MaterialModel::Interface<dim> *mm, bool vs)
          : material_model(mm), vs_(vs)
        {}

        bool need_material_properties() const
        {
          return true;
        }

        void setup(const unsigned int)
        {}

        void operator()(const MaterialModel::MaterialModelInputs<dim> &in,
                        const MaterialModel::MaterialModelOutputs<dim> &,
                        FEValues<dim> &,
                        const LinearAlgebra::BlockVector &,
                        std::vector<double> &output)
        {
          if (vs_)
            for (unsigned int q=0; q<output.size(); ++q)
              output[q] = material_model->seismic_Vs(
                            in.temperature[q], in.pressure[q], in.composition[q],
                            in.position[q]);
          else
            for (unsigned int q=0; q<output.size(); ++q)
              output[q] = material_model->seismic_Vp(
                            in.temperature[q], in.pressure[q], in.composition[q],
                            in.position[q]);
        }

        const MaterialModel::Interface<dim> *material_model;
        bool vs_;
    };

  }

  template <int dim>
  void LateralAveraging<dim>::get_Vs_averages(std::vector<double> &values) const
  {
    FunctorDepthAverageVsVp<dim> f(&this->get_material_model(), true /* Vs */);

    compute_lateral_average(values, f);
  }

  template <int dim>
  void LateralAveraging<dim>::get_Vp_averages(std::vector<double> &values) const
  {
    FunctorDepthAverageVsVp<dim> f(&this->get_material_model(), false /* Vp */);

    compute_lateral_average(values, f);
  }

  namespace
  {
    template <int dim>
    class FunctorDepthAverageVerticalHeatFlux
    {
      public:
        FunctorDepthAverageVerticalHeatFlux(const FEValuesExtractors::Vector &velocity_field,
                                            const FEValuesExtractors::Scalar &temperature_field,
                                            const GravityModel::Interface<dim> *gm)
          : velocity_field_(velocity_field),
            temperature_field_(temperature_field),
            gravity_model(gm)
        {}

        bool need_material_properties() const
        {
          return true;
        }

        void setup(const unsigned int q_points)
        {
          velocity_values.resize(q_points);
          temperature_values.resize(q_points);
          temperature_gradients.resize(q_points);
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &in,
                        const MaterialModel::MaterialModelOutputs<dim> &out,
                        FEValues<dim> &fe_values,
                        const LinearAlgebra::BlockVector &solution,
                        std::vector<double> &output)
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

  }

  template <int dim>
  void LateralAveraging<dim>::get_vertical_heat_flux_averages(std::vector<double> &values) const
  {
    FunctorDepthAverageVerticalHeatFlux<dim> f(this->introspection().extractors.velocities,
                                               this->introspection().extractors.temperature,
                                               &this->get_gravity_model() );

    compute_lateral_average(values, f);
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template class LateralAveraging<dim>;
  ASPECT_INSTANTIATE(INSTANTIATE)
}
