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

#include <deal.II/base/std_cxx11/array.h>
#include <aspect/material_model/averaging.h>
#include <utility>
#include <limits>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    AveragingOperation
    Averaging<dim>::parse_averaging_operation_name (const std::string &s)
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
      else if (s == "log average")
        return log_average;
      else if (s == "nwd arithmetic average")
        return nwd_arithmetic_average;
      else if (s == "nwd harmonic average")
        return nwd_harmonic_average;
      else if (s == "nwd geometric average")
        return nwd_geometric_average;
      else
        AssertThrow (false,
                     ExcMessage ("The value <" + s + "> for a material "
                                 "averaging operation is not one of the "
                                 "valid values."));

      return none;
    }

    // Do the requested averaging operation for one array.
    template <int dim>
    void
    Averaging<dim>::average (const AveragingOperation averaging_operation,
                             const std::vector<Point<dim> >    &position,
                             std::vector<double>           &values_out) const
    {
      // if an output field has not been filled (because it was
      // not requested), then simply do nothing -- no harm no foul
      if (values_out.size() == 0)
        return;

      const unsigned int N = values_out.size();

      // alfad is a constant which is dependent on the dimension and is used to define the shape of the bell shape.
      const double alfad = (dim == 2 ? 5/(numbers::PI * bell_shape_limit * bell_shape_limit) : 106/(numbers::PI * bell_shape_limit * bell_shape_limit * bell_shape_limit));

      // perform the requested averaging
      switch (averaging_operation)
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
          case log_average:
          {
            double sum = 0;
            for (unsigned int i=0; i<N; ++i)
              {
                Assert (values_out[i] >= 0,
                        ExcMessage ("Computing the log average "
                                    "only makes sense for non-negative "
                                    "quantities."));
                sum += std::log10(values_out[i]);
              }
            const double log_value_average = std::pow (10.,sum/N);
            for (unsigned int i=0; i<N; ++i)
              values_out[i] = log_value_average;
            break;
          }
          case nwd_arithmetic_average:
          {
            //initialize variables
            double sum_value = 0;
            double sum_weights = 0;
            std::vector<double> temp_values(N,0);

            // determine the maximum distance between all the points
            double max_distance = 0;
            for (unsigned int i=0; i<N; ++i)
              {
                for (unsigned int j=0; j<N; ++j)
                  {
                    max_distance = std::max (max_distance, position[i].distance(position[j]));
                  }
              }

            // apply the averaging to the values
            for (unsigned int i=0; i<N; ++i)
              {
                sum_value = 0;
                sum_weights = 0;

                for (unsigned int j=0; j<N; ++j)
                  {
                    const double distance = position[i].distance(position[j])/max_distance;
                    double weight = alfad * ((1 + 3 * (distance / bell_shape_limit)) * (1 - (distance / bell_shape_limit)) * (1 - (distance / bell_shape_limit)) * (1 - (distance / bell_shape_limit)));
                    if (distance / bell_shape_limit > 1)
                      weight = 0;
                    sum_value += weight * values_out[j];
                    sum_weights += weight;
                  }

                const double average = sum_value / sum_weights;

                temp_values[i] = average;
              };
            for (unsigned int i = 0; i<N; ++i)
              {
                values_out[i] = temp_values[i];
              }
            break;
          }
          case nwd_harmonic_average:
          {
            // if one of the values is zero, the average is 0.0
            for (unsigned int i=0; i<N; ++i)
              if (values_out[i] == 0.0)
                {
                  for (unsigned int j=0; j<N; ++j)
                    values_out[j] = 0.0;
                  return;
                }

            //initialize variables
            double sum_value = 0;
            double sum_weights = 0;
            std::vector<double> temp_values(N,0);

            // determine the maximum distance between all the points
            double max_distance = 0;
            for (unsigned int i=0; i<N; ++i)
              {
                for (unsigned int j=0; j<N; ++j)
                  {
                    max_distance = std::max (max_distance, position[i].distance(position[j]));
                  }
              }

            // apply the averaging to the values
            for (unsigned int i=0; i<N; ++i)
              {
                sum_value = 0;
                sum_weights = 0;

                for (unsigned int j=0; j<N; ++j)
                  {
                    const double distance = position[i].distance(position[j])/max_distance;
                    double weight = alfad * ((1 + 3 * (distance / bell_shape_limit)) * (1 - (distance / bell_shape_limit)) * (1 - (distance / bell_shape_limit)) * (1 - (distance / bell_shape_limit)));
                    if (distance > bell_shape_limit )
                      weight = 0;
                    if (values_out[j] != 0)
                      {
                        sum_value += weight / values_out[j];
                      }
                    sum_weights += weight;
                  }
                const double average = sum_weights / sum_value;

                temp_values[i] = average;
              };
            for (unsigned int i = 0; i<N; ++i)
              {
                values_out[i] = temp_values[i];
              }
            break;
          }
          case nwd_geometric_average:
          {
            //initialize variables
            double sum_value = 0;
            double sum_weights = 0;
            std::vector<double> temp_values(N,0);

            // determine the maximum distance between all the points
            double max_distance = 0;
            for (unsigned int i=0; i<N; ++i)
              {
                for (unsigned int j=0; j<N; ++j)
                  {
                    max_distance = std::max (max_distance, position[i].distance(position[j]));
                  }
              }

            // apply the averaging to the values
            for (unsigned int i=0; i<N; ++i)
              {
                sum_value = 0;
                sum_weights = 0;

                Assert (values_out[i] >= 0,
                        ExcMessage ("Computing the geometric average "
                                    "only makes sense for non-negative "
                                    "quantities."));

                for (unsigned int j=0; j<N; ++j)
                  {

                    const double distance = position[i].distance(position[j])/max_distance;
                    double weight = alfad * ((1 + 3 * (distance / bell_shape_limit)) * (1 - (distance / bell_shape_limit)) * (1 - (distance / bell_shape_limit)) * (1 - (distance / bell_shape_limit)));

                    // the weigth beyond the bell shape limit should always be zero.
                    if (distance > bell_shape_limit )
                      weight = 0;

                    /**
                     * If the value is zero to begin with the log of that value will return nan.
                     * To prevent this from happening nothing is added to sum_value in this case.
                     */
                    if (values_out[j] != 0)
                      {
                        sum_value += weight*log(values_out[j]);
                      }

                    sum_weights += weight;

                  }
                const double average = std::exp (sum_value/sum_weights);


                temp_values[i] = average;
              };
            for (unsigned int i = 0; i<N; ++i)
              {
                values_out[i] = temp_values[i];
              }
            break;
          }
          default:
          {
            AssertThrow (false,
                         ExcMessage ("This averaging operation is not implemented."));
          }
        }
    }

    template <int dim>
    void
    Averaging<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                             typename Interface<dim>::MaterialModelOutputs &out) const
    {
      // fill variable out with the results form the base material model
      base_model -> evaluate(in,out);

      /**
       * Check if the size of the viscosities (and thereby all the other vectors) is larger
       * than one. Averaging over one or zero points does not make a difference anyway,
       * and the normalized weighted distance averaging schemes need the distance between
       * the points and can not handle a distance of zero.
       */
      if (out.viscosities.size() > 1)
        {
          /* Average the base model values based on the chosen average */
          average (averaging_operation,in.position,out.viscosities);
          average (averaging_operation,in.position,out.densities);
          average (averaging_operation,in.position,out.thermal_expansion_coefficients);
          average (averaging_operation,in.position,out.specific_heat);
          average (averaging_operation,in.position,out.compressibilities);
          average (averaging_operation,in.position,out.entropy_derivative_pressure);
          average (averaging_operation,in.position,out.entropy_derivative_temperature);
        }
    }

    template <int dim>
    void
    Averaging<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Averaging");
        {
          prm.declare_entry("Base model","simple",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by an"
                            "averaging operation. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "that for more information.");
          prm.declare_entry ("Averaging operation", "none",
                             Patterns::Selection ("none|arithmetic average|harmonic average|geometric average|pick largest|log average|nwd arithmetic average|nwd harmonic average|nwd geometric average"),
                             "Chose the averaging operation to use.");
          prm.declare_entry ("Bell shape limit", "1",
                             Patterns::Double(0),
                             "The limit normalized distance between 0 and 1 where the bell shape becomes zero. See the manual for a more information.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Averaging<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Averaging");
        {
          Assert( prm.get("Base model") != "averaging",
                  ExcMessage("You may not use ``averaging'' as the base model for "
                             "a averaging model.") );

          // create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model.reset(create_material_model<dim>(prm.get("Base model")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            sim->initialize_simulator (this->get_simulator());

          averaging_operation = Averaging<dim>::parse_averaging_operation_name(prm.get ("Averaging operation"));
          bell_shape_limit = prm.get_double ("Bell shape limit");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      /* After parsing the parameters for averaging, it is essential to parse
      parameters related to the base model. */
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();
    }

    template <int dim>
    bool
    Averaging<dim>::
    is_compressible () const
    {
      return base_model->is_compressible();
    }

    template <int dim>
    double
    Averaging<dim>::
    reference_viscosity() const
    {
      return base_model->reference_viscosity();
    }

    template <int dim>
    double
    Averaging<dim>::
    reference_density() const
    {
      return base_model->reference_density();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Averaging,
                                   "averaging",
                                   "The ``averaging'' Material model applies an averaging of the quadrature points "
                                   "within a cell. The values to average are supplied by any of the other available "
                                   "material models. In other words, it is a ``compositing material model''. "
                                   "Parameters related to the average model are read from a subsection "
                                   "``Material model/Averaging''. "
                                   "\n\n"
                                   "The user must specify a ``Base model'' from which material properties are "
                                   "derived. Furthermore an averaging operation must be selected, where the "
                                   "Choice should be from the list none|arithmetic average|harmonic average|"
                                   "geometric average|pick largest|log average|NWD arithmetic average|NWD harmonic average"
                                   "|NWD geometric average. "
                                   "\n\n"
                                   "NWD stands for Normalized Weighed Distance. The models with this in front "
                                   "of their name work with a weighed average, which means each quadrature point "
                                   "requires an individual weight. The weight is determined by the distance, where "
                                   "the exact relation is determined by a bell shaped curve. A bell shaped curve is "
                                   "a continuous function which is one at it's maximum and exactly zero at and beyond "
                                   "it's limit. This bell shaped curve is spanned around each quadrature point to "
                                   "determine the weighting map for each quadrature point. The used bell shape comes "
                                   "from Lucy (1977). The distance is normalized so the largest distance becomes one. "
                                   "This means that if variable ''Bell shape limit'' is exactly one, the farthest "
                                   "quadrature point is just on the limit and it's weight will be exactly zero. In "
                                   "this plugin it is not implemented as larger and equal than the limit, but larger "
                                   "than, to ensure the the quadrature point at distance zero is always included."
                                  )
  }
}
