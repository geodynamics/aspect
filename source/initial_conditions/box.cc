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


#include <aspect/initial_conditions/box.h>
#include <aspect/geometry_model/box.h>


namespace aspect
{
  namespace InitialConditions
  {
    template <int dim>
    double
    PerturbedBox<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // this initial condition only makes sense if the geometry is a
      // Box. verify that it is indeed
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());
      AssertThrow (geometry != 0,
                   ExcMessage ("This initial condition can only be used if the geometry "
                               "is a box."));

      double perturbation = 1;
      for (unsigned int d=0; d<dim; ++d)
        perturbation *= std::sin(numbers::PI*(position[d]-geometry->get_origin()[d])/geometry->get_extents()[d]);
      return 1 + perturbation/10;
    }

    template <int dim>
    double
    PolarBox<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // this initial condition only makes sense if the geometry is a
      // Box. verify that it is indeed
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());
      AssertThrow (geometry != 0,
                   ExcMessage ("This initial condition can only be used if the geometry "
                               "is a box."));

      Point<dim> temporary1, temporary2;
      for (int d=0; d<dim; ++d)
        {
          temporary1[d]=geometry->get_extents()[d]*0.625+geometry->get_origin()[d];
          temporary2[d]=geometry->get_extents()[d]*0.375+geometry->get_origin()[d];
        }

      return 1+(1/exp(position.distance(temporary2)) - 1/exp(position.distance(temporary1)));
    }

    template<int dim>
    double
    MandelBox<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // this initial condition only makes sense if the geometry is a
      // box. verify that it is indeed
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());
      AssertThrow (geometry != 0,
                   ExcMessage ("This initial condition can only be used if the geometry "
                               "is a box."));

      double perturbation, ratio;
      Point<dim> center;
      ratio = center[0] = geometry->get_extents()[0]*0.66;
      center[1] = geometry->get_extents()[1]*0.5;
      if (center[1] < ratio)
        ratio = center[1];

      double zx = (position[0] - geometry->get_origin()[0] - center[0])/ratio;
      double zy = (position[1] - geometry->get_origin()[1] - center[1])/ratio;
      double x = zx;
      double y = zy;

      for (perturbation = 0; perturbation < 50 && (Point<2>(x,y)).norm() <= 2; ++perturbation)
        {
          x = x*x - y*y + zx;
          y = 2 * x*y + zy;
        }
      return perturbation / 50;
    }



    template<int dim>
    double
    InclusionShapeBox<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // this initial condition only makes sense if the geometry is a
      // box. verify that it is indeed
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());
      AssertThrow (geometry != 0,
                   ExcMessage ("This initial condition can only be used if the geometry "
                               "is a box."));

      double perturbation = 0;
      Point<dim> center;
      for (int d=0; d<dim; ++d)
        center[d] = (Point<3>(center_x,center_y,center_z))[d];

      if (inclusion_shape == "square")
        {
          if (inclusion_gradient == "gaussian")
            {
              perturbation = inclusion_temperature - ambient_temperature;
              for (int d=0; d<dim; ++d)
                perturbation *= exp(-8*pow(position.distance(center)/radius, 2));
            }
          else if (inclusion_gradient == "linear")
            {
              double x = position[0] - center[0];
              double y = position[1] - center[1];
              if ( x <= y && x >= -y)
                perturbation = radius - y;
              else if (x <= y && x <= -y)
                perturbation = radius + x;
              else if (x >= y && x >= -y)
                perturbation = radius - x;
              else if (x >= y && x <= -y)
                perturbation = radius + y;
              perturbation *= (inclusion_temperature - ambient_temperature) / radius;
            }
          else if (inclusion_gradient == "constant")
            {
              perturbation = inclusion_temperature - ambient_temperature;
            }
          for (int d = 0; d < dim; ++d)
            if (abs (center[d] - position[d]) > radius)
              perturbation = 0;
        }
      else if (inclusion_shape == "circle")
        {
          if (inclusion_gradient == "gaussian")
            {
              perturbation = inclusion_temperature - ambient_temperature;
              perturbation *= exp(-pow(position.distance(center),2) / (2 * pow((radius / 4), 2))) / (2 * radius);
            }
          else if (inclusion_gradient == "linear")
            {
              perturbation = ((radius - position.distance(center)) / radius) * (inclusion_temperature - ambient_temperature);
            }
          else if (inclusion_gradient == "constant")
            {
              perturbation = inclusion_temperature - ambient_temperature;
            }
          if (position.distance(center) > radius)
            perturbation = 0;
        }

      return ambient_temperature + perturbation;
    }

    template <int dim>
    void
    InclusionShapeBox<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial conditions");
      {
        prm.enter_subsection ("Inclusion shape perturbation");
        {
          prm.declare_entry("Inclusion shape", "circle",
                            Patterns::Selection("square|circle"),
                            "The shape of the inclusion to be generated.");
          prm.declare_entry("Inclusion gradient", "constant",
                            Patterns::Selection("gaussian|linear|constant"),
                            "The gradient of the inclusion to be generated.");
          prm.declare_entry("Shape radius", "1.0",
                            Patterns::Double (0),
                            "The radius of the inclusion to be generated. For "
                            "shapes with no radius (e.g. square), this will "
                            "be the width, and for shapes with no width, this "
                            "gives a general guideline for the size of the shape.");
          prm.declare_entry("Ambient temperature", "1.0",
                            Patterns::Double (),
                            "The background temperature for the temperature field.");
          prm.declare_entry("Inclusion temperature", "0.0",
                            Patterns::Double (),
                            "The temperature of the inclusion shape. This is only "
                            "the true temperature in the case of the constant "
                            "gradient. In all other cases, it gives one endpoint "
                            "of the temperature gradient for the shape.");
          prm.declare_entry("Center X", "0.5",
                            Patterns::Double (),
                            "The X coordinate for the center of the shape.");
          prm.declare_entry("Center Y", "0.5",
                            Patterns::Double (),
                            "The Y coordinate for the center of the shape.");
          prm.declare_entry("Center Z", "0.5",
                            Patterns::Double (),
                            "The Z coordinate for the center of the shape. This "
                            "is only necessary for three-dimensional fields.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void
    InclusionShapeBox<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Inclusion shape perturbation");
        {
          inclusion_shape = prm.get ("Inclusion shape");
          inclusion_gradient = prm.get ("Inclusion gradient");
          radius = prm.get_double ("Shape radius");
          ambient_temperature = prm.get_double ("Ambient temperature");
          inclusion_temperature = prm.get_double ("Inclusion temperature");
          center_x = prm.get_double ("Center X");
          center_y = prm.get_double ("Center Y");
          center_z = prm.get_double ("Center Z");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(PerturbedBox,
                                       "perturbed box",
                                       "An initial temperature field in which the temperature "
                                       "is perturbed slightly from an otherwise constant value "
                                       "equal to one. The perturbation is chosen in such a way "
                                       "that the initial temperature is constant to one along "
                                       "the entire boundary.")
    ASPECT_REGISTER_INITIAL_CONDITIONS(PolarBox,
                                       "polar box",
                                       "An initial temperature field in which the temperature "
                                       "is perturbed slightly from an otherwise constant value "
                                       "equal to one. The perturbation is such that there are "
                                       "two poles on opposing corners of the box. ")
    ASPECT_REGISTER_INITIAL_CONDITIONS(InclusionShapeBox,
                                       "inclusion shape perturbation",
                                       "An initial temperature field in which there is an "
                                       "inclusion in a constant-temperature box field. The size, "
                                       "shape, gradient, position, and temperature of the "
                                       "inclusion are defined by parameters.")
    ASPECT_REGISTER_INITIAL_CONDITIONS(MandelBox,
                                       "mandelbox",
                                       "Fractal-shaped temperature field.")
  }
}
