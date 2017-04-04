/*
  Copyright (C) 2012 by the authors of the ASPECT code.

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
/*  $Id: mckenzie.cc 1088 2013-10-22 13:45 glerum $  */


#include <aspect/initial_conditions/polygon_depth.h>
#include <aspect/geometry_model/interface.h>
#include <boost/lexical_cast.hpp>
#include <aspect/utilities.h>


namespace aspect
{
namespace InitialConditions
{

template <int dim>
double
PolygonDepth<dim>::
initial_temperature (const Point<dim> &position) const
{
	const GeometryModel::Interface<dim>* geometry_model= &this->get_geometry_model();
	const ChartManifold<dim,dim>* manifold = geometry_model->get_manifold();

	const dealii::Point<dim> internal_position = manifold != NULL ? manifold->pull_back(position) : position; // radius, longitude, lattitude or x,y and z;

	const Point<2> internal_position_surface(internal_position(1), dim == 3 ? internal_position(2) : 0);

	double temperature;
	const double depth = geometry_model->depth(position);

	if(Utilities::polygon_contains_point<dim>(coordinate_list, internal_position_surface) && depth <= bottom_depth)
	{
		temperature = reference_temperature;
	}
	else
	{
		temperature = 5;
	}
	return temperature;
}


template <int dim>
void
PolygonDepth<dim>::declare_parameters (ParameterHandler &prm)
{
	prm.enter_subsection ("Initial conditions");
	{
		prm.enter_subsection ("Polygon depth");
		{
			prm.declare_entry("Reference temperature","273", Patterns::Double(),
					          "Set the minimum parts per degree, which determines how many extra coordinate points will be generated.");
			prm.declare_entry("Bottom depth","600e3", Patterns::Double(),
					          "Set the minimum parts per degree, which determines how many extra coordinate points will be generated.");
			prm.declare_entry("List of polygon coordinates","", Patterns::Anything (),
					          "Set the minimum parts per degree, which determines how many extra coordinate points will be generated.");
		}
		prm.leave_subsection ();
	}
	prm.leave_subsection ();

}


template <int dim>
void
PolygonDepth<dim>::parse_parameters (ParameterHandler &prm)
{
	prm.enter_subsection ("Initial conditions");
	{
		prm.enter_subsection ("Polygon depth");
		{
			reference_temperature = prm.get_double ("Reference temperature");
			bottom_depth = prm.get_double ("Bottom depth");

			const GeometryModel::Interface<dim>* geometry_model= &this->get_geometry_model();
			const ChartManifold<dim,dim>* manifold = geometry_model->get_manifold();

			/**
			 * If there is no manifold, we assume we don't have to transform the coordinates to radians,
			 * otherwise it is some kind of ellipsoid/sphere and we do convert it from degree to radians.
			 */
			const bool use_degree_to_rad = (manifold != NULL ? true : false);

			/**
			 * Retrieve the coordinates
			 */
			const std::vector<std::string> coordinate_list_strings = Utilities::split_string_list(prm.get ("List of polygon coordinates"),';');
			const unsigned int number_of_coordinates = coordinate_list_strings.size();
			coordinate_list.resize(number_of_coordinates);

			AssertThrow (number_of_coordinates >= 2,ExcMessage ("Need at least two coordinate points per object."));

			for(unsigned int i=0;i<number_of_coordinates;++i)
			{
				const std::vector<double> coordinate_list_temp = Utilities::string_to_double(Utilities::split_string_list(coordinate_list_strings[i],','));
				AssertThrow (coordinate_list_temp.size() == 2,
						ExcMessage("Surface objects are defined at the surface in a 3d world, and should always 2 values in 2d and 3d. Coordinate " +
								boost::lexical_cast<std::string>(i+1) + " does not meet this requirement."));

				coordinate_list[i](0) =  coordinate_list_temp[0] * (use_degree_to_rad ? degree_to_rad : 1);
				coordinate_list[i](1) =  coordinate_list_temp[1] * (use_degree_to_rad ? degree_to_rad : 1);
			}

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
ASPECT_REGISTER_INITIAL_CONDITIONS(PolygonDepth,
		"polygon depth",
		"Temperature is prescribed as a linear profile in the lithosphere, "
		"adiabat in the mantle and according to McKenzie 1970 in the slab. "
		"Slab properties (e.g. dip) can be specified in the input file.")
}
}
