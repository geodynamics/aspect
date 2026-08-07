/*
  Copyright (C) 2025 - 2026 by the authors of the ASPECT code.

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


#include <aspect/mesh_deformation/landlab.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/patterns.h>
#if DEAL_II_VERSION_GTE(9,8,0)
#include <deal.II/numerics/data_out_points.h>
#endif

#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <deal.II/base/array_view.h>

#include <fstream>

#include <cfenv>

#include <aspect/python_helper.h>

#if defined(ASPECT_WITH_PYTHON) && defined(ASPECT_WITH_LANDLAB)

using namespace dealii;

/**
 * Call a Python function from module @p pModule with name @p func_name
 * Returns the result (caller must Py_DECREF). Throws on error.
 */
PyObject *call_python_function(PyObject *pModule, const char *func_name, PyObject *pArgs = nullptr)
{
  PyObject *pFunc = PyObject_GetAttrString(pModule, func_name);
  if (!pFunc || !PyCallable_Check(pFunc))
    AssertThrow(false, ExcMessage(std::string("Failed to load function: ") + func_name));

  PyObject *pValue = PyObject_CallObject(pFunc, pArgs);
  Py_DECREF(pFunc);

  if (pValue == nullptr)
    {
      if (PyErr_Occurred())
        PyErr_Print();
      AssertThrow(false, ExcMessage(std::string(func_name) + " returned NULL"));
    }
  return pValue;
}

#endif



namespace aspect
{
#if defined(ASPECT_WITH_PYTHON) && defined(ASPECT_WITH_LANDLAB)
  namespace MeshDeformation
  {
    template <int dim>
    void
    Landlab<dim>::initialize ()
    {
      // Determine whether we are in a spherical geometry or not, which affects how we interpret the coordinates of the evaluation points.
      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::CoordinateSystem::spherical)
        is_spherical = true;
      else if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::CoordinateSystem::cartesian)
        is_spherical = false;
      else
        AssertThrow(false, ExcMessage("The Landlab mesh deformation plugin only supports Cartesian and spherical geometries."));

      unsigned int rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());

      this_rank_runs_landlab = (rank < n_landlab_ranks);
      const int color = this_rank_runs_landlab?1:0;
      const int ierr = MPI_Comm_split(this->get_mpi_communicator(), color, rank, &landlab_communicator);
      AssertThrow(ierr == MPI_SUCCESS, ExcMessage("Failed to split MPI communicator for Landlab simulation"));

      if (this_rank_runs_landlab)
        {
          // Append script dirs so env packages (venv site-packages, PYTHONPATH) are found first
          // for "import landlab":
          PyRun_SimpleString("import sys");
          PyRun_SimpleString("sys.path.append(\"" ASPECT_SOURCE_DIR "/tests\")");
          PyRun_SimpleString("sys.path.append(\".\")");
          PyRun_SimpleString("sys.path.append(\"" ASPECT_SOURCE_DIR "/contrib/python/scripts\")");

          // avoid floating point exceptions in Landlab Python code:
#ifdef ASPECT_USE_FP_EXCEPTIONS
          fedisableexcept(FE_DIVBYZERO|FE_INVALID);
#endif

          std::cout << "importing '" << script_module_name << "' ..." << std::endl;
          pModule = PyImport_ImportModule(script_module_name.c_str());
          if (PyErr_Occurred())
            PyErr_Print();
          AssertThrow(pModule, ExcMessage("Failed to load Python module"));

          // Copy the landlab script to the output directory for reproducibility.
          const std::string module_filename = (script_path.empty() ? "" : script_path + "/") + script_module_name + ".py";
          std::ifstream python_source(module_filename, std::ios::binary);

          const std::string copied_module_filename = this->get_output_directory() + "original_landlab.py";
          std::ofstream python_copy(copied_module_filename, std::ios::binary);
          AssertThrow(python_copy,
                      ExcMessage("Failed to open destination file for writing: " + copied_module_filename));

          python_copy << python_source.rdbuf();

          // Call Python initialize() function with communicator handle
          PyObject *pArgs;
          if (n_landlab_ranks == 1)
            pArgs = PyTuple_Pack(1, Py_None);
          else
            pArgs = PyTuple_Pack(1, PyLong_FromLong(MPI_Comm_c2f(landlab_communicator)));
          PyObject *pValue = call_python_function(pModule, "initialize", pArgs);

          Py_DECREF(pArgs);
          Py_DECREF(pValue);
        }
    }



    template <int dim>
    void
    Landlab<dim>::update ()
    {
      if (!this->remote_point_evaluator)
        {
          if (!this_rank_runs_landlab)
            {
              // This rank does not participate, so we don't own any evaluation points:
              std::vector<Point<dim>> surface_points;
              this->set_evaluation_points(surface_points);
              return;
            }

          {
            // set_mesh_information: call with None
            PyObject *pArgs = PyTuple_Pack(1, Py_None);
            PyObject *pValue = call_python_function(pModule, "set_mesh_information", pArgs);
            Py_DECREF(pArgs);
            Py_DECREF(pValue);
          }

          {
            // get grid:
            PyObject *pArgs = PyTuple_Pack(1, PyLong_FromLong(dim));
            PyObject *pgrid_x = call_python_function(pModule, "get_grid_x", pArgs);

            // Depending on the ASPECT model geometry and the dimension, we need to
            // include the y coordinates (3D Cartesian), and the z coordinates
            // (3D spherical) of the Landlab grid.
            PyObject *pgrid_y = nullptr;
            PyObject *pgrid_z = nullptr;
            if (dim == 3)
              pgrid_y = call_python_function(pModule, "get_grid_y", pArgs);
            if (dim == 3 && is_spherical)
              pgrid_z = call_python_function(pModule, "get_grid_z", pArgs);
            Py_DECREF(pArgs);

            // Create a C++ view of the numpy arrays
            const ArrayView<double> data_x = PythonHelper::numpy_to_array_view(pgrid_x);
            const ArrayView<double> data_y = (dim == 3)
                                             ? PythonHelper::numpy_to_array_view(pgrid_y)
                                             : ArrayView<double>(nullptr, 0);
            const ArrayView<double> data_z = (dim == 3 && is_spherical)
                                             ? PythonHelper::numpy_to_array_view(pgrid_z)
                                             : ArrayView<double>(nullptr, 0);

            if (dim == 3)
              AssertThrow(data_x.size() == data_y.size(), ExcMessage("get_grid_x and get_grid_y returned different sizes"));
            if (dim == 3 && is_spherical)
              AssertThrow(data_x.size() == data_z.size(), ExcMessage("get_grid_x and get_grid_z returned different sizes"));

            // Loop over the ArrayViews and store them in a vector of Points. These
            // are the 'evaluation points'.
            std::vector<Point<dim>> surface_points(data_x.size());
            for (size_t i = 0; i < data_x.size(); i++)
              {
                const double depth = this->get_geometry_model().representative_point(0.0)[dim-1];
                if (dim == 2)
                  surface_points[i] = Point<dim>(data_x[i], depth);
                else if (dim == 3)
                  {
                    const double z = (is_spherical) ? data_z[i] : depth;
                    surface_points[i] = Point<dim>(data_x[i], data_y[i], z);
                  }
              }

            // Clean up Python objects
            Py_DECREF(pgrid_x);
            if (pgrid_y)
              Py_DECREF(pgrid_y);
            if (pgrid_z)
              Py_DECREF(pgrid_z);

            this->set_evaluation_points(surface_points);
          }
        }
    }



    template <int dim>
    std::vector<Tensor<1,dim>>
    Landlab<dim>::compute_updated_velocities_at_points (const std::vector<std::vector<double>> &current_solution_at_points) const
    {
      Assert(current_solution_at_points.size() == this->evaluation_points.size(), ExcInternalError());

      // Initialize the vector that will compute the velocities in ASPECT from the information
      // sent from Landlab.
      std::vector<Tensor<1,dim>> velocities(current_solution_at_points.size(), Tensor<1,dim>());

      if (this_rank_runs_landlab)
        {
          // Build a dictionary with solution values for each variable to pass to Python.
          // This is the x and y velocity in 2D, and the z velocity in 3D, as well as the
          // pressure, temperature, and compositional fields.

          // Create dictionary to hold variable names and their corresponding data
          PyObject *pDict_solution  = PyDict_New();
          PyObject *pDict_auxiliary = PyDict_New();

          // Add velocities
          std::vector<std::string> variable_names = {"x velocity", "y velocity"};
          if (dim == 3)
            variable_names.push_back("z velocity");

          // Add pressure and temperature
          variable_names.push_back("pressure");
          variable_names.push_back("temperature");

          // Add compositional fields
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            variable_names.push_back(this->introspection().name_for_compositional_index(c));

          // Loop over all solution variables at each Landlab evaluation point and store the
          // ASPECT solution in a vector.
          std::vector<std::vector<double>> variable_data(variable_names.size(),  std::vector<double>(current_solution_at_points.size(), 0.0));
          for (unsigned int i=0; i<variable_names.size(); ++i)
            {
              for (unsigned int j=0; j<current_solution_at_points.size(); ++j)
                {
                  variable_data[i][j] = current_solution_at_points[j][i];
                }
            }
          // Store the solution vector for each variable in a python dictionary to send to Landlab.
          for (unsigned int i=0; i<variable_names.size(); ++i)
            {
              auto pValue = PythonHelper::vector_to_numpy_object(variable_data[i]);
              PyDict_SetItemString(pDict_solution, variable_names[i].c_str(), pValue.get());
            }

          // Create a second dictionary which holds other information that is useful for Landlab to know about the ASPECT model.
          // This is used for keeping landlab and ASPECT in sync while running, and also for checkpointing/restarting and
          // postprocessing.
          PyDict_SetItemString(pDict_auxiliary, "ASPECT dimension", PyLong_FromLong(dim));
          PyDict_SetItemString(pDict_auxiliary, "ASPECT model time", PyFloat_FromDouble(this->get_time()));
          PyDict_SetItemString(pDict_auxiliary, "ASPECT timestep size", PyFloat_FromDouble(this->get_timestep()));
          PyDict_SetItemString(pDict_auxiliary, "ASPECT timestep number", PyFloat_FromDouble(this->get_timestep_number()));
          PyDict_SetItemString(pDict_auxiliary, "ASPECT output directory", PyUnicode_FromString(this->get_output_directory().c_str()));

          // Call update_until(), which is the main loop in Landlab that evolves the topography.
          // update_until() returns the change in the topography, which we convert to a mesh
          // velocity in ASPECT.
          PyObject *pArgs  = PyTuple_Pack(2, pDict_solution, pDict_auxiliary);
          PyObject *pValue = call_python_function(pModule, "update_until", pArgs);

          // Remove these python objects from memory.
          Py_DECREF(pDict_solution);
          Py_DECREF(pDict_auxiliary);
          Py_DECREF(pArgs);

          // Convert the returned numpy array to a C++ view and compute the mesh velocities in ASPECT.
          const ArrayView<const double> data = PythonHelper::numpy_to_array_view(pValue);
          const double one_over_dt = 1.0 / ((this->get_timestep() > 0.0) ? this->get_timestep() : 1.0);

          // The velocity is calculated as the change in topography divided by the time step, and is assumed
          // to move the mesh in the direction opposite to the gravity vector.
          for (size_t i=0; i<data.size(); ++i)
            {
              const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(this->evaluation_points[i]);
              Tensor<1,dim> topography_direction;
              if (gravity.norm() > 0.0)
                topography_direction = -gravity / gravity.norm();

              velocities[i] = topography_direction * data[i] * one_over_dt;
            }
          // Remove the python object from memory.
          Py_DECREF(pValue);
        }

      // Produce debug output as a vtu file
#if DEAL_II_VERSION_GTE(9,8,0)
      {
        static unsigned int output_no = 0;

        DataOutPoints<dim, dim> out;
        std::vector<Point<dim>> real_evaluation_points(this->evaluation_points.size());
        std::vector<std::vector<double>> data(this->evaluation_points.size(), std::vector<double>(dim, 0.0));
        const double one_over_dt = 1.0 / ((this->get_timestep() > 0.0) ? this->get_timestep() : 1.0);

        for (unsigned int i=0; i<this->evaluation_points.size(); ++i)
          {
            real_evaluation_points[i] = this->evaluation_points[i];  // TODO: use mapping to compute real position
            for (unsigned int c=0; c<dim; ++c)
              data[i][c] = velocities[i][c] * one_over_dt;
          }

        const std::vector<std::string> data_component_names(dim, "velocity");
        const std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretations(dim, DataComponentInterpretation::component_is_part_of_vector);

        out.build_patches(real_evaluation_points, 0, data, data_component_names, data_component_interpretations);

        out.write_vtu_with_pvtu_record(this->get_output_directory(), "surf_points", output_no, this->get_mpi_communicator(), 4, 0);

        ++output_no;
      }
#endif

      return velocities;
    }



    template <int dim>
    void Landlab<dim>::
    compute_initial_deformation_as_constraints(const Mapping<dim> &/*mapping*/,
                                               const DoFHandler<dim> &mesh_deformation_dof_handler,
                                               const types::boundary_id boundary_indicator,
                                               AffineConstraints<double> &constraints) const
    {
      // We need to initialize the evaluation points in order to extract the initial topography
      // from Landlab and apply it as constraints on the initial mesh in ASPECT. This means that
      // we need to call update(), which determines the evaluation points, which requires that
      // this function is not actually a const function. We get around this by casting away the
      // const.
      const_cast<Landlab<dim>*>(this)->update();

      // Grab the initial topography from Landlab and convert it to a mesh deformation in ASPECT.
      // This is done in three steps:
      // 1. Receive from Landlab the initial topography at the evaluation points. The topography
      //    is assumed to be in the direction opposite to the gravity vector.
      std::vector<Tensor<1,dim>> initial_deformation(this->evaluation_points.size(), Tensor<1,dim>());
      if (this_rank_runs_landlab)
        {
          PyObject *pArgs  = PyTuple_Pack(1, PyLong_FromLong(dim));
          PyObject *pValue = call_python_function(pModule, "get_initial_topography", pArgs);
          Py_DECREF(pArgs);
          ArrayView<double> data = PythonHelper::numpy_to_array_view(pValue);

          for (size_t i=0; i<data.size(); ++i)
            {
              const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(this->evaluation_points[i]);
              Tensor<1,dim> topography_direction;
              if (gravity.norm() > 0.0)
                topography_direction = -gravity / gravity.norm();

              initial_deformation[i] = data[i] * topography_direction;
            }
          Py_DECREF(pValue);
        }

      // 2. Interpolate deformation into a DoF vector:
      LinearAlgebra::Vector initial_deformation_dof_vector = this->interpolate_external_velocities_to_surface_support_points(initial_deformation);
      const DoFHandler<dim> &mesh_dof_handler = this->get_mesh_deformation_handler().get_mesh_deformation_dof_handler();
      const IndexSet mesh_locally_relevant = DoFTools::extract_locally_relevant_dofs (mesh_dof_handler);
      LinearAlgebra::Vector initial_deformation_ghosted(mesh_dof_handler.locally_owned_dofs(),
                                                        mesh_locally_relevant,
                                                        this->get_mpi_communicator());
      initial_deformation_ghosted = initial_deformation_dof_vector;

      const IndexSet constrained_dofs = DoFTools::extract_boundary_dofs(mesh_deformation_dof_handler,
                                                                        ComponentMask(dim, true),
      {boundary_indicator});

      // 3. Add constraints from DoF values:
      for (const types::global_dof_index index : constrained_dofs)
        {
          if (constraints.can_store_line(index))
            if (constraints.is_constrained(index)==false)
              {
                constraints.add_constraint(index,
                                           {},
                                           initial_deformation_ghosted(index));
              }
        }
    }



    template <int dim>
    void Landlab<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh deformation");
      {
        prm.enter_subsection("Landlab");
        {
          prm.declare_entry("MPI ranks for Landlab", "1",
                            Patterns::Integer(1),
                            "Number of ranks to use for the Landlab simulation. If set to 1, the Landlab simulation will run sequentially "
                            "without MPI. If set to -1, the Landlab simulation will run on all ranks.");
          prm.declare_entry("Script path", "",
                            Patterns::Anything(),
                            "Path to the Python script to execute. Relative paths and the placeholders "
                            "ASPECT_SOURCE_DIR and ASPECT_BINARY_DIR are allowed.");
          prm.declare_entry("Script name", "",
                            Patterns::Anything(),
                            "Name of the Python module to load (without .py extension).");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void Landlab<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Landlab");
        {
          n_landlab_ranks = prm.get_integer("MPI ranks for Landlab");
          AssertThrow(n_landlab_ranks == 1,
                      ExcMessage("The Landlab mesh deformation model currently only supports running on a single rank. "
                                 "Please set 'MPI ranks for Landlab' to 1 in the parameter file."));

          script_path        = prm.get("Script path");
          script_module_name = prm.get("Script name");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(Landlab,
                                           "Landlab",
                                           "A mesh deformation plugin that lets a Python script control the "
                                           "deformation of the surface. It is meant for coupling with the landscape evolution "
                                           "code Landlab, but any other script that provides the necessary functions can be used. "
                                           "It is necessary to have Python and numpy with their C APIs installed and that "
                                           "ASPECT_WITH_PYTHON and ASPECT_WITH_LANDLAB are enabled when ASPECT is configured with "
                                           "CMake. ")
  }
#endif
}
