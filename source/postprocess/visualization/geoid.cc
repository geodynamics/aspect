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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/postprocess/visualization/geoid.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      std::pair<std::string, Vector<float> *>
      Geoid<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("geoid",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        Postprocess::Geoid<dim> *geoid_postprocessor = this->template find_postprocessor<Postprocess::Geoid<dim> >();
        AssertThrow(geoid_postprocessor != NULL,
                    ExcMessage("Could not find the Geoid postprocessor"
                               "Perhaps you forgot to include it in the Postprocessors list?"));
        Postprocess::internal::InternalMultipoleExpansion<dim> cmb_expansion = geoid_postprocessor->get_cmb_potential_expansion();
        Postprocess::internal::ExternalMultipoleExpansion<dim> surface_expansion = geoid_postprocessor->get_surface_potential_expansion();


        const double gravity_at_surface = this->get_gravity_model().gravity_vector(
                                            ( this->get_geometry_model().representative_point( 0.0 ) ) ).norm();
        const double gravity_at_bottom = this->get_gravity_model().gravity_vector(
                                           this->get_geometry_model().representative_point(
                                             this->get_geometry_model().maximal_depth() )  ).norm();
        const GeometryModel::SphericalShell<dim> *geometry_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>
                                                                   (&this->get_geometry_model());
        AssertThrow (geometry_model != 0,
                     ExcMessage("The geoid postprocessor is currently only implemented for "
                                "the spherical shell geometry model."));

        const double inner_radius = geometry_model->inner_radius();
        const double outer_radius = geometry_model->outer_radius();

        // loop over all of the surface cells and if one less than h/3 away from
        // the top or bottom surface
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned())
            {
              (*return_value.second)(cell_index) = 0.0;
              if (cell->at_boundary())
                {
                  // see if the cell is at the *top* or *bottom* boundary
                  bool surface_cell = false;
                  bool bottom_cell = false;

                  unsigned int f = 0;
                  for (; f<GeometryInfo<dim>::faces_per_cell; ++f)
                    {
                      if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                        {
                          surface_cell = true;
                          break;
                        }
                      if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) > (outer_radius - inner_radius - cell->face(f)->minimum_vertex_distance()/3))
                        {
                          bottom_cell = true;
                          break;
                        }
                    }


                  //Get the location of the cell for the expansion
                  const Point<dim> p = cell->face(f)->center();

                  double potential_value = 0.;
                  //Do the outer expansion at the point if it is a surface cell
                  if (surface_cell)
                    {
                      potential_value = surface_expansion.evaluate( p );
                      (*return_value.second)(cell_index)  = potential_value / gravity_at_surface;
                    }
                  //Do the inner expansion at the point if it is a surface cell
                  if (bottom_cell)
                    {
                      potential_value = cmb_expansion.evaluate( p );
                      (*return_value.second)(cell_index)  = potential_value / gravity_at_bottom;
                    }
                }
            }

        return return_value;
      }

      template <int dim>
      std::list<std::string>
      Geoid<dim>::required_other_postprocessors() const
      {
        std::list<std::string> deps;
        deps.push_back("geoid");
        return deps;
      }

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Geoid,
                                                  "geoid",
                                                  "")
    }
  }
}
