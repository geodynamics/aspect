/*
  Copyright (C) 2022 - 2023 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/geometry_model/two_merged_boxes.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>
#include <aspect/simulator_access.h>
#include <limits>

int f()
{
  using namespace aspect;
  // test whether the functions return the inverse of each other.

  // test point
  Point<2> point_2d(25000,50000);
  Point<3> point_3d(25000,25000,50000);

  std::array<double,2> inter_point_2d;
  std::array<double,3> inter_point_3d;
  Point<2> new_point_2d = Point<2>();
  Point<3> new_point_3d= Point<3>();
  // get the supported geometries
  std::cout << std::endl << "2D:" << std::endl;
  // Box 2D
  {
    GeometryModel::Box<2> box;
    ParameterHandler prm_box;
    box.declare_parameters(prm_box);
    prm_box.enter_subsection("Geometry model");
    prm_box.enter_subsection("Box");
    prm_box.set("X extent", "100e3");
    prm_box.set("Y extent", "100e3");
    prm_box.leave_subsection();
    prm_box.leave_subsection();
    box.parse_parameters(prm_box);

    inter_point_2d = box.cartesian_to_natural_coordinates(point_2d);
    new_point_2d = box.natural_to_cartesian_coordinates(inter_point_2d);

    if (std::fabs(new_point_2d(0) - point_2d(0)) < 1e-8 &&
        std::fabs(new_point_2d(1) - point_2d(1)) < 1e-8)
      std::cout << "Box 2D ok:                   point = " << point_2d << ", inter point = " << inter_point_2d[0] << ", " << inter_point_2d[1] << ", new_point = " << new_point_2d << std::endl;
    else
      std::cout << "Box 2D failed:               point = " << point_2d << ", inter point = " << inter_point_2d[0] << ", " << inter_point_2d[1] << ", new_point = " << new_point_2d << std::endl;
  }

  // Two merged boxes 2d
  {
    GeometryModel::TwoMergedBoxes<2> two_merged_boxes;
    ParameterHandler prm_two_merged_boxes;
    two_merged_boxes.declare_parameters(prm_two_merged_boxes);
    prm_two_merged_boxes.enter_subsection("Geometry model");
    prm_two_merged_boxes.enter_subsection("Box with lithosphere boundary indicators");
    prm_two_merged_boxes.set("X extent", "100e3");
    prm_two_merged_boxes.set("Y extent", "100e3");
    prm_two_merged_boxes.leave_subsection();
    prm_two_merged_boxes.leave_subsection();
    two_merged_boxes.parse_parameters(prm_two_merged_boxes);
    inter_point_2d = two_merged_boxes.cartesian_to_natural_coordinates(point_2d);
    new_point_2d = two_merged_boxes.natural_to_cartesian_coordinates(inter_point_2d);

    if (std::fabs(new_point_2d(0) - point_2d(0)) < 1e-8 &&
        std::fabs(new_point_2d(1) - point_2d(1)) < 1e-8)
      std::cout << "Two merged boxes 2d ok:      point = " << point_2d << ", inter point = " << inter_point_2d[0] << ", " << inter_point_2d[1] << ", new_point = " << new_point_2d << std::endl;
    else
      std::cout << "Two merged boxes 2d failed:  point = " << point_2d << ", inter point = " << inter_point_2d[0] << ", " << inter_point_2d[1] << ", new_point = " << new_point_2d << std::endl;
  }

  // Chunk 2d
  {
    GeometryModel::Chunk<2> chunk;
    ParameterHandler prm_chunk;
    chunk.declare_parameters(prm_chunk);
    prm_chunk.enter_subsection("Geometry model");
    prm_chunk.enter_subsection("Chunk");
    prm_chunk.set("Chunk inner radius", "6128137.0");
    prm_chunk.set("Chunk outer radius", "6378137.0");
    prm_chunk.leave_subsection();
    prm_chunk.leave_subsection();
    chunk.parse_parameters(prm_chunk);

    inter_point_2d = chunk.cartesian_to_natural_coordinates(point_2d);
    new_point_2d = chunk.natural_to_cartesian_coordinates(inter_point_2d);

    if (std::fabs(new_point_2d(0) - point_2d(0)) < 1e-8 &&
        std::fabs(new_point_2d(1) - point_2d(1)) < 1e-8)
      std::cout << "Chunk 2D ok:                 point = " << point_2d << ", inter point = " << inter_point_2d[0] << ", " << inter_point_2d[1] << ", new_point = " << new_point_2d << std::endl;
    else
      std::cout << "Chunk 2D failed:             point = " << point_2d << ", inter point = " << inter_point_2d[0] << ", " << inter_point_2d[1] << ", new_point = " << new_point_2d << std::endl;

  }
  std::cout << std::endl << "3D:" << std::endl;
  // Box 3D
  {
    GeometryModel::Box<3> box;
    //box.initialize(topo);
    ParameterHandler prm_box;
    box.declare_parameters(prm_box);
    prm_box.enter_subsection("Geometry model");
    prm_box.enter_subsection("Box");
    prm_box.set("X extent", "100e3");
    prm_box.set("Y extent", "100e3");
    prm_box.set("Z extent", "200e3");
    prm_box.leave_subsection();
    prm_box.leave_subsection();
    box.parse_parameters(prm_box);

    inter_point_3d = box.cartesian_to_natural_coordinates(point_3d);
    new_point_3d = box.natural_to_cartesian_coordinates(inter_point_3d);

    if (std::fabs(new_point_3d(0) - point_3d(0)) < 1e-8 &&
        std::fabs(new_point_3d(1) - point_3d(1)) < 1e-8 &&
        std::fabs(new_point_3d(2) - point_3d(2)) < 1e-8)
      std::cout << "Box 3D ok:                   point = " << point_3d << ", inter point = " << inter_point_3d[0] << ", " << inter_point_3d[1] << ", " << inter_point_3d[2]<< ", new_point = " << new_point_3d << std::endl;
    else
      std::cout << "Box 3D failed:               point = " << point_3d << ", inter point = " << inter_point_3d[0] << ", " << inter_point_3d[1] << ", " << inter_point_3d[2]<< ", new_point = " << new_point_3d << std::endl;
  }

  // Two merged boxes 3d
  {
    GeometryModel::TwoMergedBoxes<3> two_merged_boxes;
    ParameterHandler prm_two_merged_boxes;
    two_merged_boxes.declare_parameters(prm_two_merged_boxes);
    prm_two_merged_boxes.enter_subsection("Geometry model");
    prm_two_merged_boxes.enter_subsection("Box with lithosphere boundary indicators");
    prm_two_merged_boxes.set("X extent", "100e3");
    prm_two_merged_boxes.set("Y extent", "100e3");
    prm_two_merged_boxes.set("Z extent", "200e3");
    prm_two_merged_boxes.leave_subsection();
    prm_two_merged_boxes.leave_subsection();
    two_merged_boxes.parse_parameters(prm_two_merged_boxes);
    inter_point_3d = two_merged_boxes.cartesian_to_natural_coordinates(point_3d);
    new_point_3d = two_merged_boxes.natural_to_cartesian_coordinates(inter_point_3d);

    if (std::fabs(new_point_3d(0) - point_3d(0)) < 1e-8 &&
        std::fabs(new_point_3d(1) - point_3d(1)) < 1e-8 &&
        std::fabs(new_point_3d(2) - point_3d(2)) < 1e-8)
      std::cout << "Two merged boxes 3d ok:      point = " << point_3d << ", inter point = " << inter_point_3d[0] << ", " << inter_point_3d[1] << ", " << inter_point_3d[2]<< ", new_point = " << new_point_3d << std::endl;
    else
      std::cout << "Two merged boxes 3d failed:  point = " << point_3d << ", inter point = " << inter_point_3d[0] << ", " << inter_point_3d[1] << ", " << inter_point_3d[2]<< ", new_point = " << new_point_3d << std::endl;
  }

  // Chunk 3d
  {
    GeometryModel::Chunk<3> chunk;
    ParameterHandler prm_chunk;
    chunk.declare_parameters(prm_chunk);
    prm_chunk.enter_subsection("Geometry model");
    prm_chunk.enter_subsection("Chunk");
    prm_chunk.set("Chunk inner radius", "6128137.0");
    prm_chunk.set("Chunk outer radius", "6378137.0");
    prm_chunk.leave_subsection();
    prm_chunk.leave_subsection();
    chunk.parse_parameters(prm_chunk);

    inter_point_3d = chunk.cartesian_to_natural_coordinates(point_3d);
    new_point_3d = chunk.natural_to_cartesian_coordinates(inter_point_3d);

    if (std::fabs(new_point_3d(0) - point_3d(0)) < 1e-8 &&
        std::fabs(new_point_3d(1) - point_3d(1)) < 1e-8 &&
        std::fabs(new_point_3d(2) - point_3d(2)) < 1e-8)
      std::cout << "Chunk 3D ok:                 point = " << point_3d << ", inter point = " << inter_point_3d[0] << ", " << inter_point_3d[1] << ", new_point = " << new_point_3d << std::endl;
    else
      std::cout << "Chunk 3D failed:             point = " << point_3d << ", inter point = " << inter_point_3d[0] << ", " << inter_point_3d[1] << ", new_point = " << new_point_3d << std::endl;
  }

  // Ellipsoidal chunk 3d
  {
    GeometryModel::EllipsoidalChunk<3> ellipsoidal_chunk;
    ParameterHandler prm_ellispoidal_chunk;
    ellipsoidal_chunk.declare_parameters(prm_ellispoidal_chunk);
    prm_ellispoidal_chunk.enter_subsection("Geometry model");
    prm_ellispoidal_chunk.enter_subsection("Ellipsoidal chunk");
    prm_ellispoidal_chunk.set("NE corner", "10:10");
    prm_ellispoidal_chunk.set("SW corner", "0:0");
    prm_ellispoidal_chunk.leave_subsection();
    prm_ellispoidal_chunk.leave_subsection();
    ellipsoidal_chunk.parse_parameters(prm_ellispoidal_chunk);
    inter_point_3d = ellipsoidal_chunk.cartesian_to_natural_coordinates(point_3d);
    new_point_3d = ellipsoidal_chunk.natural_to_cartesian_coordinates(inter_point_3d);

    if (std::fabs(new_point_3d(0) - point_3d(0)) < 1e-8 &&
        std::fabs(new_point_3d(1) - point_3d(1)) < 1e-8 &&
        std::fabs(new_point_3d(2) - point_3d(2)) < 1e-8)
      std::cout << "Ellipsoidal 3d chunk ok:     point = " << point_3d << ", inter point = " << inter_point_3d[0] << ", " << inter_point_3d[1] << ", " << inter_point_3d[2]<< ", new_point = " << new_point_3d << std::endl;
    else
      std::cout << "Ellipsoidal 3d chunk failed: point = " << point_3d << ", inter point = " << inter_point_3d[0] << ", " << inter_point_3d[1] << ", " << inter_point_3d[2]<< ", new_point = " << new_point_3d << std::endl;
  }


  exit(0);
  return 42;
}
// run this function by initializing a global variable by it
int i = f();
