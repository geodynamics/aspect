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


#include <aspect/postprocess/table_heat_flux_statistics.h>
#include <aspect/material_model/table.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/box.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    TableHeatfluxStatistics<dim>::execute (TableHandler &statistics)
    {
      Assert (dynamic_cast<const MaterialModel::Table<dim> *>(&this->get_material_model())
              != 0,
              ExcMessage ("This postprocessor can only be used when using "
                          "the MaterialModel::Table implementation of the "
                          "material model interface."));

      // create a quadrature formula based on the temperature element alone.
      const QGauss<dim-1> quadrature (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature,
                                        update_gradients      | update_values |
                                        update_normal_vectors |
                                        update_q_points       | update_JxW_values);

      std::vector<Tensor<1,dim> > temperature_gradients (quadrature.size());
      std::vector<double>         temperature_values (quadrature.size());
      std::vector<double>         pressure_values (quadrature.size());
      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature.size()));
      std::vector<double> composition_values_at_q_point (this->n_compositional_fields());


      // find out which boundary indicators are related to Dirichlet temperature boundaries.
      // it only makes sense to compute heat fluxes on these boundaries.
      const std::set<types::boundary_id>
      boundary_indicators
        = this->get_geometry_model().get_used_boundary_indicators ();
      std::map<types::boundary_id, double> local_boundary_fluxes;

      typename MaterialModel::Interface<dim>::MaterialModelInputs in(quadrature.size(),
                                                                     this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(quadrature.size(),
                                                                       this->n_compositional_fields());

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      // for every surface face on which it makes sense to compute a
      // heat flux and that is owned by this processor,
      // integrate the normal heat flux given by the formula
      //   j =  - k * n . grad T
      //
      // for the spherical shell geometry, note that for the inner boundary,
      // the normal vector points *into* the core, i.e. we compute the flux
      // *out* of the mantle, not into it. we fix this when we add the local
      // contribution to the global flux
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->at_boundary(f))
              {
                fe_face_values.reinit (cell, f);
                fe_face_values[this->introspection().extractors.temperature].get_function_gradients (this->get_solution(),
                    temperature_gradients);
                fe_face_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                    in.temperature);
                fe_face_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                               in.pressure);
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  fe_face_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                      composition_values[c]);

                in.position = fe_face_values.get_quadrature_points();
                in.strain_rate.resize(0); // we do not need the viscosity
                for (unsigned int i=0; i<quadrature.size(); ++i)
                  {
                    for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                      in.composition[i][c] = composition_values[c][i];
                  }

                this->get_material_model().evaluate(in, out);

                double local_normal_flux = 0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    local_normal_flux
                    +=
                      -out.thermal_conductivities[q] *
                      (temperature_gradients[q] *
                       fe_face_values.normal_vector(q)) *
                      fe_face_values.JxW(q);
                  }

                local_boundary_fluxes[cell->face(f)->boundary_indicator()]
                += local_normal_flux;
              }

      // now communicate to get the global values
      std::map<types::boundary_id, double> global_boundary_fluxes;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        std::vector<double> local_values;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p)
          local_values.push_back (local_boundary_fluxes[*p]);

        // then collect contributions from all processors
        std::vector<double> global_values;
        Utilities::MPI::sum (local_values, this->get_mpi_communicator(), global_values);

        // and now take them apart into the global map again
        unsigned int index = 0;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p, ++index)
          global_boundary_fluxes[*p] = global_values[index];
      }


//TODO: This doesn't scale to more geometry models. simply output the data as is,
// i.e. with association from boundary id to heat flux.
      // record results and have something for the output. this depends
      // on the interpretation of what boundary is which, which we can
      // only do knowing what the geometry is
      if (dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&this->get_geometry_model())
          != 0)
        {
          // we have a spherical shell. note that in that case we want
          // to invert the sign of the core-mantle flux because we
          // have computed the flux with a normal pointing from
          // the mantle into the core, and not the other way around

          statistics.add_value ("Core-mantle heat flux (W)", -global_boundary_fluxes[0]);
          statistics.add_value ("Surface heat flux (W)",     global_boundary_fluxes[1]);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision ("Core-mantle heat flux (W)", 8);
          statistics.set_scientific ("Core-mantle heat flux (W)", true);
          statistics.set_precision ("Surface heat flux (W)", 8);
          statistics.set_scientific ("Surface heat flux (W)", true);

          /*-------------------------------------------------
             output Inner/Outer Nu number to file per timestep*/
          const GeometryModel::SphericalShell<dim> *geometry = dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&this->get_geometry_model());
          const double phi = (*geometry).opening_angle();
          const double R0 = (*geometry).inner_radius();
          const double R1 = (*geometry).outer_radius();
          const double h = R1-R0;

          // dT is only meaningful if boundary temperatures are prescribed, otherwise it is 0
          const double dT = (&this->get_boundary_temperature())
                            ?
                            this->get_boundary_temperature().maximal_temperature(this->get_fixed_temperature_boundary_indicators())
                            - this->get_boundary_temperature().minimal_temperature(this->get_fixed_temperature_boundary_indicators())
                            :
                            0;

          const double conductive_heatflux = dT/h;
          const double nusselt_outer = global_boundary_fluxes[0]/conductive_heatflux;
          const double boundary_curveLength_outer = R0*phi;

          statistics.add_value ("Outer Nusselt number", nusselt_outer/boundary_curveLength_outer);
          statistics.set_precision ("Outer Nusselt number", 4);
          statistics.set_scientific ("Outer Nusselt number", true);

          const double nusselt_inner= global_boundary_fluxes[1]/conductive_heatflux;
          const double boundary_curveLength_inner = R1*phi;
          statistics.add_value ("Inner Nusselt number", nusselt_inner/boundary_curveLength_inner);
          statistics.set_precision ("Inner Nusselt number", 4);
          statistics.set_scientific ("Inner Nusselt number", true);

          // finally have something for the screen
          std::ostringstream output;
          output.precision(4);
          output << -global_boundary_fluxes[0] << " W, "
                 << global_boundary_fluxes[1] << " W,"
                 << nusselt_outer/boundary_curveLength_outer << " W,"
                 << nusselt_inner/boundary_curveLength_inner << " W,";
          return std::pair<std::string, std::string> ("Inner/outer heat fluxes Nusselt number:",
                                                      output.str());
        }
      else if (dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model())
               != 0)
        {
          // for the box geometry we can associate boundary indicators 0 and 1
          // to left and right boundaries
          statistics.add_value ("Left boundary heat flux (W)", global_boundary_fluxes[0]);
          statistics.add_value ("Right boundary heat flux (W)", global_boundary_fluxes[1]);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision ("Left boundary heat flux (W)", 8);
          statistics.set_scientific ("Left boundary heat flux (W)", true);
          statistics.set_precision ("Right boundary heat flux (W)", 8);
          statistics.set_scientific ("Right boundary heat flux (W)", true);

          // finally have something for the screen
          std::ostringstream output;
          output.precision(4);
          output << global_boundary_fluxes[0] << " W, "
                 << global_boundary_fluxes[1] << " W";

          return std::pair<std::string, std::string> ("Left/right heat fluxes:",
                                                      output.str());
        }
      else
        AssertThrow (false, ExcNotImplemented());

      return std::pair<std::string, std::string> ("",
                                                  "");
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(TableHeatfluxStatistics,
                                  "heat flux statistics for the table model",
                                  "A postprocessor that computes some statistics about "
                                  "the heat flux across boundaries.")
  }
}
