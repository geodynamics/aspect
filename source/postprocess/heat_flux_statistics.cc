//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/postprocess/heat_flux_statistics.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/equation_data.h>
#include <aspect/simulator.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    HeatFluxStatistics<dim>::execute (TableHandler &statistics)
    {
      Assert (dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&this->get_geometry_model()) != 0,
              ExcNotImplemented());

      const QGauss<dim-1> quadrature_formula (this->get_temperature_dof_handler().get_fe().degree+1);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_temperature_dof_handler().get_fe(),
                                        quadrature_formula,
                                        update_gradients | update_normal_vectors | update_JxW_values);

      std::vector<Tensor<1,dim> > temperature_gradients (quadrature_formula.size());

      const double thermal_conductivity = this->get_material_model().thermal_conductivity();

      double local_inner_boundary_flux = 0;
      double local_outer_boundary_flux = 0;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_temperature_dof_handler().begin_active(),
      endc = this->get_temperature_dof_handler().end();

      // for every core-mantle or surface face owned by this processor,
      // integrate the normal heat flux given by the formula
      //   j =  - k * n . grad T
      //
      // note however that for the inner boundary, the normal vector
      // points *into* the core, i.e. we compute the flux *out* of the
      // mantle, not into it. we fix this when we add the local contribution
      // to the global flux
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            // check if the face is at the boundary and has either boundary indicator
            // zero (inner boundary) or one (outer boundary)
            if (cell->at_boundary(f)
                &&
                (cell->face(f)->boundary_indicator() < 2))
              {
                fe_face_values.reinit (cell, f);
                fe_face_values.get_function_gradients (this->get_temperature_solution(),
                                                       temperature_gradients);

                double local_normal_flux = 0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  local_normal_flux
                  +=
                    -thermal_conductivity *
                    (temperature_gradients[q] *
                     fe_face_values.normal_vector(q)) *
                    fe_face_values.JxW(q);

                if (cell->face(f)->boundary_indicator() == 0)
                  local_inner_boundary_flux -= local_normal_flux;
                else if (cell->face(f)->boundary_indicator() == 1)
                  local_outer_boundary_flux += local_normal_flux;
                else
                  Assert (false, ExcInternalError());
              }

      // now communicate to get the global values
      double global_inner_boundary_flux = 0;
      double global_outer_boundary_flux = 0;
      {
        double local_values[2] = { local_inner_boundary_flux, local_outer_boundary_flux };
        double global_values[2];

        Utilities::MPI::sum (local_values, MPI_COMM_WORLD, global_values);

        global_inner_boundary_flux = global_values[0];
        global_outer_boundary_flux = global_values[1];
      }

      // record results and have something for the output
      statistics.add_value ("Core-mantle heat flux (W)", global_inner_boundary_flux);
      statistics.add_value ("Surface heat flux (W)",     global_outer_boundary_flux);

      // now also compute Inner/Outer Nu number and put that into the statistics
      // object as well
      const double Scaling = 1. /(this->get_boundary_temperature().maximal_temperature()
                                  -
                                  this->get_boundary_temperature().minimal_temperature());

      // check that the following formula makes sense in 3d as well.
      Assert (dim==2, ExcNotImplemented());

      const GeometryModel::SphericalShell<dim> &geometry
        = dynamic_cast<const GeometryModel::SphericalShell<dim> &>(this->get_geometry_model());

      const double inner_boundary_curve_length = geometry.inner_radius()/(geometry.outer_radius() - geometry.inner_radius()) *
                                                 (geometry.opening_angle() / 360 * 2 * numbers::PI);
      const double outer_boundary_curve_length = geometry.outer_radius()/(geometry.outer_radius() - geometry.inner_radius()) *
                                                 (geometry.opening_angle() / 360 * 2 * numbers::PI);

      const double dTdr_inner = (global_inner_boundary_flux / thermal_conductivity)*Scaling;
      statistics.add_value ("Inner Nu number", dTdr_inner / inner_boundary_curve_length);
      const double dTdr_outer = (global_outer_boundary_flux / thermal_conductivity)*Scaling;
      statistics.add_value ("Outer Nu number", dTdr_outer / outer_boundary_curve_length);

      // finally have something for the screen
      std::ostringstream output;
      output.precision(4);
      output << global_inner_boundary_flux << " W, "
             << global_outer_boundary_flux << " W";

      return std::pair<std::string, std::string> ("Inner/outer heat fluxes:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    template class HeatFluxStatistics<deal_II_dimension>;

    ASPECT_REGISTER_POSTPROCESSOR("heat flux statistics", HeatFluxStatistics)
  }
}
