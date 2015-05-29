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



#include <aspect/postprocess/heating_statistics.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    HeatingStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the temperature element alone.
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());

      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      const std::list<std_cxx11::shared_ptr<HeatingModel::Interface<dim> > > &heating_model_objects = this->get_heating_model_manager().get_active_heating_models();
      const std::vector<std::string> &heating_model_names = this->get_heating_model_manager().get_active_heating_model_names();

      HeatingModel::HeatingModelOutputs heating_model_outputs(n_q_points, this->n_compositional_fields());

      std::ostringstream output;
      output.precision(4);

      double average_heating_integral = 0.0;
      double total_heating_integral = 0.0;

      std::vector<double> local_heating_integrals (heating_model_objects.size());
      double local_mass = 0.0;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      // compute the integral quantities by quadrature
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {

            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.temperature]
            .get_function_values (this->get_solution(),
                                  in.temperature);
            fe_values[this->introspection().extractors.pressure]
            .get_function_values (this->get_solution(),
                                  in.pressure);
            fe_values[this->introspection().extractors.velocities]
            .get_function_values (this->get_solution(),
                                  in.velocity);
            fe_values[this->introspection().extractors.pressure]
            .get_function_gradients (this->get_solution(),
                                     in.pressure_gradient);
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              fe_values[this->introspection().extractors.compositional_fields[c]]
              .get_function_values(this->get_solution(),
                                   composition_values[c]);
            for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
              {
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  in.composition[i][c] = composition_values[c][i];
              }

            fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(),
                in.strain_rate);
            in.position = fe_values.get_quadrature_points();

            this->get_material_model().evaluate(in, out);

            for (unsigned int q=0; q<n_q_points; ++q)
              local_mass += out.densities[q] * fe_values.JxW(q);

            unsigned int index = 0;
            for (typename std::list<std_cxx11::shared_ptr<HeatingModel::Interface<dim> > >::const_iterator
                 heating_model = heating_model_objects.begin();
                 heating_model != heating_model_objects.end(); ++heating_model, ++index)
              {
                (*heating_model)->evaluate(in, out, heating_model_outputs);

                for (unsigned int q=0; q<n_q_points; ++q)
                  local_heating_integrals[index] += heating_model_outputs.heating_source_terms[q]
                                                    * fe_values.JxW(q);
              }
          }

      std::vector<double> global_heating_integrals (heating_model_objects.size());
      double global_mass = 0.0;

      // compute the sum over all processors
      Utilities::MPI::sum (local_heating_integrals,
                           this->get_mpi_communicator(),
                           global_heating_integrals);
      global_mass = Utilities::MPI::sum (local_mass, this->get_mpi_communicator());

      unsigned int index = 0;
      for (typename std::list<std_cxx11::shared_ptr<HeatingModel::Interface<dim> > >::const_iterator
           heating_model = heating_model_objects.begin();
           heating_model != heating_model_objects.end(); ++heating_model, ++index)
        {
          // finally produce something for the statistics file
          const std::string name1("Average " + heating_model_names[index] + " rate (W/kg) ");
          statistics.add_value (name1, global_heating_integrals[index]/global_mass);
          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision (name1, 8);
          statistics.set_scientific (name1, true);

          const std::string name2("Total " + heating_model_names[index] + " rate (W) ");
          statistics.add_value (name2, global_heating_integrals[index]);
          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision (name2, 8);
          statistics.set_scientific (name2, true);

          total_heating_integral += global_heating_integrals[index];
        }

      average_heating_integral = total_heating_integral/global_mass;

      output << average_heating_integral << " W/kg, "
             << total_heating_integral << " W";

      return std::pair<std::string, std::string> ("Heating rate (average/total): ",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(HeatingStatistics,
                                  "heating statistics",
                                  "A postprocessor that computes some statistics about "
                                  "heating, averaged by volume. ")
  }
}
