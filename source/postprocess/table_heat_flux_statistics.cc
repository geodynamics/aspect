//-------------------------------------------------------------
//    $Id: heat_flux_statistics.cc 572 2012-01-03 08:24:22Z geenen $
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/postprocess/table_heat_flux_statistics.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/box.h>
#include <aspect/simulator.h>

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
      const QGauss<dim-1> quadrature_formula (this->get_temperature_dof_handler().get_fe().degree+1);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_temperature_dof_handler().get_fe(),
                                        quadrature_formula,
                                        update_gradients      | update_values |
                                        update_normal_vectors |
                                        update_q_points       | update_JxW_values);

      FEFaceValues<dim> stokes_fe_face_values (this->get_mapping(),
                                               this->get_stokes_dof_handler().get_fe(),
                                               quadrature_formula,
                                               update_values);
      const FEValuesExtractors::Scalar pressure (dim);

      std::vector<Tensor<1,dim> > temperature_gradients (quadrature_formula.size());
      std::vector<double>         temperature_values (quadrature_formula.size());
      std::vector<double>         pressure_values (quadrature_formula.size());

      // find out which boundary indicators are related to Dirichlet temperature boundaries.
      // it only makes sense to compute heat fluxes on these boundaries.
      const std::set<unsigned char>
      boundary_indicators
        =
          this->get_geometry_model().get_used_boundary_indicators ();
      std::map<unsigned char, double> local_boundary_fluxes;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_temperature_dof_handler().begin_active(),
      endc = this->get_temperature_dof_handler().end();
      typename DoFHandler<dim>::active_cell_iterator
      stokes_cell = this->get_stokes_dof_handler().begin_active();

      // for every surface face on which it makes sense to compute a
      // heat flux and that is owned by this processor,
      // integrate the normal heat flux given by the formula
      //   j =  - k * n . grad T
      //
      // for the spherical shell geometry, note that for the inner boundary,
      // the normal vector points *into* the core, i.e. we compute the flux
      // *out* of the mantle, not into it. we fix this when we add the local
      // contribution to the global flux
      for (; cell!=endc; ++cell, ++stokes_cell)
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            // check if the face is at the boundary and has either boundary indicator
            // zero (inner boundary) or one (outer boundary)
            if (cell->at_boundary(f))
              {
                fe_face_values.reinit (cell, f);
                fe_face_values.get_function_gradients (this->get_temperature_solution(),
                                                       temperature_gradients);
                fe_face_values.get_function_values (this->get_temperature_solution(),
                                                    temperature_values);

                stokes_fe_face_values.reinit (stokes_cell, f);
                stokes_fe_face_values[pressure].get_function_values (this->get_stokes_solution(),
                                                                     pressure_values);

                double local_normal_flux = 0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    const double thermal_conductivity
                      = this->get_material_model().thermal_conductivity(temperature_values[q],
                                                                        pressure_values[q],
                                                                        fe_face_values.quadrature_point(q));

                    local_normal_flux
                    +=
                      -thermal_conductivity *
                      (temperature_gradients[q] *
                       fe_face_values.normal_vector(q)) *
                      fe_face_values.JxW(q);
                  }

                local_boundary_fluxes[cell->face(f)->boundary_indicator()]
                += local_normal_flux;
              }

      // now communicate to get the global values
      std::map<unsigned char, double> global_boundary_fluxes;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        std::vector<double> local_values;
        for (std::set<unsigned char>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p)
          local_values.push_back (local_boundary_fluxes[*p]);

        // then collect contributions from all processors
        std::vector<double> global_values;
        Utilities::MPI::sum (local_values, MPI_COMM_WORLD, global_values);

        // and now take them apart into the global map again
        unsigned int index = 0;
        for (std::set<unsigned char>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p, ++index)
          global_boundary_fluxes[*p] = local_values[index];
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
          double phi = (*geometry).opening_angle();
          double R0 = (*geometry).inner_radius();
          double R1 = (*geometry).outer_radius();
          double h = R1-R0;

          double dT = this->get_boundary_temperature().maximal_temperature() - this->get_boundary_temperature().minimal_temperature();
          double ConductiveHeatflux = dT/h;
          double NusseltOuter = global_boundary_fluxes[0]/ConductiveHeatflux;
          double BoundaryCurveLengthOuter = R0*phi;

          statistics.add_value ("Outer Nusselt number", NusseltOuter/BoundaryCurveLengthOuter);
          statistics.set_precision ("Outer Nusselt number", 4);
          statistics.set_scientific ("Outer Nusselt number", true);

          double NusseltInner= global_boundary_fluxes[1]/ConductiveHeatflux;
          double BoundaryCurveLengthInner = R1*phi;
          statistics.add_value ("Inner Nusselt number", NusseltInner/BoundaryCurveLengthInner);
          statistics.set_precision ("Inner Nusselt number", 4);
          statistics.set_scientific ("Inner Nusselt number", true);

          // finally have something for the screen
          std::ostringstream output;
          output.precision(4);
          output << -global_boundary_fluxes[0] << " W, "
                 << global_boundary_fluxes[1] << " W,"
                 << NusseltOuter/BoundaryCurveLengthOuter << " W,"
                 << NusseltInner/BoundaryCurveLengthInner << " W,";
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
    template class TableHeatfluxStatistics<deal_II_dimension>;

    ASPECT_REGISTER_POSTPROCESSOR(TableHeatfluxStatistics,
                                  "heat flux statistics for the table model",
                                  "A postprocessor that computes some statistics about "
                                  "the heat flux across boundaries.")
  }
}
