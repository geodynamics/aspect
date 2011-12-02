//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/geometry_model/box.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>


namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    Box<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      GridGenerator::hyper_cube (coarse_grid, 0, 1);
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        coarse_grid.begin_active()->face(f)->set_boundary_indicator(f);
    }


    template <int dim>
    std::set<unsigned char>
    Box<dim>::
    get_temperature_dirichlet_boundary_indicators () const
    {
      // see the documentation of this class and how we set boundary
      // indicators above to understand the following code
      const unsigned char s[] = { 0, 1 };
      return std::set<unsigned char>(&s[0],
                                     &s[sizeof(s)/sizeof(s[0])]);
    }


    template <int dim>
    std::set<unsigned char>
    Box<dim>::
    get_zero_velocity_boundary_indicators () const
    {
      // see the documentation of this class and how we set boundary
      // indicators above to understand the following code
      std::set<unsigned char> s;
      for (unsigned int i=2; i<2*dim; ++i)
        s.insert (i);
      return s;
    }



    template <int dim>
    std::set<unsigned char>
    Box<dim>::
    get_tangential_velocity_boundary_indicators () const
    {
      // see the documentation of this class and how we set boundary
      // indicators above to understand the following code
      const unsigned char s[] = { 0, 1 };
      return std::set<unsigned char>(&s[0],
                                     &s[sizeof(s)/sizeof(s[0])]);
    }


    template <int dim>
    double
    Box<dim>::
    length_scale () const
    {
      return 0.01;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace GeometryModel
  {
    template class Box<deal_II_dimension>;

    ASPECT_REGISTER_GEOMETRY_MODEL(Box,
                                   "box",
                                   "A box geometry with fixed length 1 in each coordinate direction.")
  }
}
