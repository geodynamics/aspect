// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code
// and contributions can be found in the folder LICENSES
// and in CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/mesh_deformation/landlab.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/patterns.h>
#if DEAL_II_VERSION_GTE(9,8,0)
#include <deal.II/numerics/data_out_points.h>
#endif

#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <deal.II/base/array_view.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>
#include <numeric>

#include <cfenv>

#include <aspect/python_helper.h>

using namespace dealii;
namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    void
    Landlab<dim>::initialize ()
    {
#ifdef ASPECT_WITH_LANDLAB
      // Determine whether we are in a spherical geometry or not, which affects how we interpret the coordinates of the evaluation points.
      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::CoordinateSystem::spherical)
        is_spherical = true;
      else if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::CoordinateSystem::cartesian)
        is_spherical = false;
      else
        AssertThrow(false, ExcMessage("The Landlab mesh deformation plugin only supports Cartesian and spherical geometries."));

      const unsigned int rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());

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
          PyRun_SimpleString(("sys.path.append(\"" + script_path + "\")").c_str());

          // disable floating point exceptions in Landlab Python code during the module import
          // ("import landlab" crashes otherwise)
#ifdef ASPECT_USE_FP_EXCEPTIONS
          fedisableexcept(FE_DIVBYZERO|FE_INVALID);
#endif

          std::cout << "importing '" << script_module_name << "' ..." << std::endl;
          pModule = PyImport_ImportModule(script_module_name.c_str());
          if (PyErr_Occurred())
            PyErr_Print();
          AssertThrow(pModule, ExcMessage("Failed to load Python module"));

#ifdef ASPECT_USE_FP_EXCEPTIONS
          feenableexcept(FE_DIVBYZERO|FE_INVALID);
#endif

          // Copy the landlab script to the output directory for reproducibility. Only do this on a single rank
          // to avoid having multiple ranks trying to write the same file.
          if (rank == 0)
            {
              const std::string module_filename = (script_path.empty() ? "" : script_path + "/") + script_module_name + ".py";
              std::ifstream python_source(module_filename, std::ios::binary);

              const std::string copied_module_filename = this->get_output_directory() + "original_landlab.py";
              std::ofstream python_copy(copied_module_filename, std::ios::binary);
              AssertThrow(python_copy,
                          ExcMessage("Failed to open destination file for writing: " + copied_module_filename));

              python_copy << python_source.rdbuf();
            }

          // Call Python initialize() function with communicator handle
          PyObject *pArgs;
          if (n_landlab_ranks == 1)
            pArgs = PyTuple_Pack(1, Py_None);
          else
            pArgs = PyTuple_Pack(1, PyLong_FromLong(MPI_Comm_c2f(landlab_communicator)));
          PyObject *pValue = PythonHelper::call_python_function(pModule, "initialize", pArgs);

          Py_DECREF(pArgs);
          Py_DECREF(pValue);
        }
#else
      AssertThrow(false, ExcMessage("To use the 'Landlab' mesh deformation plugin, ASPECT needs to be configured using ASPECT_WITH_LANDLAB=ON."));
#endif
    }



    template <int dim>
    double
    Landlab<dim>::
    boundary_composition (const types::boundary_id boundary_indicator,
                          const Point<dim> &/*position*/,
                          const unsigned int compositional_field) const
    {
#ifdef ASPECT_WITH_LANDLAB
      if (boundary_indicator != this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top"))
        return 0.0;

      if ( this->introspection().compositional_name_exists("sediment_age") &&
           compositional_field == this->introspection().compositional_index_for_name("sediment_age"))
        {
          return this->get_parameters().convert_to_years ? this->get_time()/year_in_seconds : this->get_time();
        }
      return 0.0;
#else
      (void) boundary_indicator;
      (void) compositional_field;
      return 0.0;
#endif
    }



    template <int dim>
    void
    Landlab<dim>::update ()
    {
#ifdef ASPECT_WITH_LANDLAB
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
            PyObject *pValue = PythonHelper::call_python_function(pModule, "set_mesh_information", pArgs);
            Py_DECREF(pArgs);
            Py_DECREF(pValue);
          }

          {
            // get grid:
            PyObject *pArgs = PyTuple_Pack(1, PyLong_FromLong(dim));
            PyObject *pgrid_x = PythonHelper::call_python_function(pModule, "get_grid_x", pArgs);

            // Depending on the ASPECT model geometry and the dimension, we need to
            // include the y coordinates (3D Cartesian), and the z coordinates
            // (3D spherical) of the Landlab grid.
            PyObject *pgrid_y = nullptr;
            PyObject *pgrid_z = nullptr;
            if (dim == 3)
              pgrid_y = PythonHelper::call_python_function(pModule, "get_grid_y", pArgs);
            if (dim == 3 && is_spherical)
              pgrid_z = PythonHelper::call_python_function(pModule, "get_grid_z", pArgs);
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
                // If the geometry is spherical, Landlab's spherical mesh returns the x,y,z coordinates of the surface points
                // in the same Cartesian coordinate system as ASPECT.
                // TODO: Support 2D spherical annulus
                if (is_spherical)
                  {
                    if (dim == 3)
                      surface_points[i] = Point<dim>(data_x[i], data_y[i], data_z[i]);
                    if (dim == 2)
                      AssertThrow(false, ExcMessage("Spherical coordinates in 2D is not yet implemented"));

                  }
                else
                  {
                    const double surface_coordinate = this->get_geometry_model().representative_point(0.0)[dim-1];
                    if (dim == 3)
                      surface_points[i] = Point<dim>(data_x[i], data_y[i], surface_coordinate);
                    if (dim == 2)
                      surface_points[i] = Point<dim>(data_x[i], surface_coordinate);
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
#endif
    }



    template <int dim>
    std::vector<Tensor<1,dim>>
    Landlab<dim>::compute_updated_velocities_at_points (const std::vector<std::vector<double>> &current_solution_at_points) const
    {
#ifdef ASPECT_WITH_LANDLAB
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

          // Add any additional derived outputs to the variable list before sending the data to Landlab.
          evaluate_derived_quantities_at_points(variable_data, variable_names);

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
          PyObject *pValue = PythonHelper::call_python_function(pModule, "update_until", pArgs);

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
              else
                AssertThrow(false, ExcMessage("The gravity vector has a zero norm. To use the Landlab mesh deformation plugin, "
                                              "the gravity vector must have a non-zero norm."));

              velocities[i] = topography_direction * data[i] * one_over_dt;
            }
          // Remove the python object from memory.
          Py_DECREF(pValue);
        }

      // Produce debug output as a vtu file
#if DEAL_II_VERSION_GTE(9,8,0)
#ifdef DEBUG
      {
        static unsigned int output_no = 0;

        DataOutPoints<dim, dim> out;
        std::vector<Point<dim>> real_evaluation_points(this->evaluation_points.size());
        std::vector<std::vector<double>> data(this->evaluation_points.size(), std::vector<double>(dim, 0.0));
        const double one_over_dt = 1.0 / ((this->get_timestep() > 0.0) ? this->get_timestep() : 1.0);

        for (unsigned int i=0; i<this->evaluation_points.size(); ++i)
          {
            // TODO: use mapping to compute real position
            real_evaluation_points[i] = this->evaluation_points[i];
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
#endif

      return velocities;
#else
      (void) current_solution_at_points;
      return {};
#endif
    }



    template <int dim>
    void
    Landlab<dim>::
    evaluate_derived_quantities_at_points (std::vector<std::vector<double>> &variable_data,
                                           std::vector<std::string> &variable_names) const
    {
#ifdef ASPECT_WITH_LANDLAB
      // If the user is not requesting additional quantities, return immediately.
      const unsigned int n_derived_quantities = additional_named_quantities.size();
      if (n_derived_quantities == 0)
        return;

      // Construct a vector for storing the additional derived quantities at each evaluation point.
      const unsigned int n_eval_points = this->evaluation_points.size();
      std::vector<std::vector<double>> derived_quantities_at_points(n_derived_quantities, std::vector<double>(n_eval_points, 0.0));

      std::vector<unsigned int> point_indices(n_eval_points);
      std::iota(point_indices.begin(), point_indices.end(), 0);

      const auto eval_func = [&](const ArrayView<const unsigned int> &values,
                                 const typename Utilities::MPI::RemotePointEvaluation<dim>::CellData &cell_data)
      {
        for (unsigned int derived_quantity_index = 0; derived_quantity_index < n_derived_quantities; ++derived_quantity_index)
          {
            if (additional_named_quantities[derived_quantity_index] == "strain rate")
              for (const auto cell_index : cell_data.cell_indices())
                {
                  const auto cell = cell_data.get_active_cell_iterator(cell_index)->as_dof_handler_iterator(this->get_dof_handler());
                  const ArrayView<const Point<dim>> unit_points = cell_data.get_unit_points(cell_index);

                  // Evaluate the strain rate at the actual requested evaluation points.
                  const Quadrature<dim> quadrature(std::vector<Point<dim>>(unit_points.begin(), unit_points.end()));

                  FEValues<dim> fe_values(this->get_mapping(),
                                          this->get_fe(),
                                          quadrature,
                                          update_gradients);
                  fe_values.reinit(cell);

                  std::vector<SymmetricTensor<2,dim>> strain_rate(unit_points.size());
                  fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients(this->get_solution(), strain_rate);

                  const ArrayView<const unsigned int> local_values(values.data() + cell_data.reference_point_ptrs[cell_index],
                                                                   cell_data.reference_point_ptrs[cell_index + 1] - cell_data.reference_point_ptrs[cell_index]);

                  for (unsigned int i = 0; i < unit_points.size(); ++i)
                    {
                      const unsigned int point_index = local_values[i];
                      derived_quantities_at_points[derived_quantity_index][point_index] =
                        std::sqrt(std::fabs(Utilities::Tensors::consistent_second_invariant_of_deviatoric_tensor(
                                              Utilities::Tensors::consistent_deviator(strain_rate[i]))));
                    }
                }
          }
      };

      this->remote_point_evaluator->template process_and_evaluate<unsigned int, 1>(point_indices, eval_func, /*sort_data*/ true);

      for (unsigned int derived_quantity = 0; derived_quantity < n_derived_quantities; ++derived_quantity)
        {
          variable_data.push_back(derived_quantities_at_points[derived_quantity]);
          variable_names.push_back(additional_named_quantities[derived_quantity]);
        }
#else
      (void) variable_data;
      (void) variable_names;
#endif
    }



    template <int dim>
    void Landlab<dim>::
    compute_initial_deformation_as_constraints(const Mapping<dim> &/*mapping*/,
                                               const DoFHandler<dim> &mesh_deformation_dof_handler,
                                               const types::boundary_id boundary_indicator,
                                               AffineConstraints<double> &constraints) const
    {
#ifdef ASPECT_WITH_LANDLAB
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
          PyObject *pValue = PythonHelper::call_python_function(pModule, "get_initial_topography", pArgs);
          Py_DECREF(pArgs);
          ArrayView<double> data = PythonHelper::numpy_to_array_view(pValue);

          for (size_t i=0; i<data.size(); ++i)
            {
              const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(this->evaluation_points[i]);
              Tensor<1,dim> topography_direction;
              if (gravity.norm() > 0.0)
                topography_direction = -gravity / gravity.norm();
              else
                AssertThrow(false, ExcMessage("The gravity vector has a zero norm. To use the Landlab mesh deformation plugin, "
                                              "the gravity vector must have a non-zero norm."));

              initial_deformation[i] = data[i] * topography_direction;
            }
          Py_DECREF(pValue);
        }

      // 2. Interpolate deformation into a DoF vector:
      LinearAlgebra::Vector initial_deformation_dof_vector = this->interpolate_external_vector_field_to_surface_support_points(initial_deformation);
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
#else
      (void) mesh_deformation_dof_handler;
      (void) boundary_indicator;
      (void) constraints;
#endif
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
                            "Number of ranks to use for the Landlab simulation. Currently, only 1 is supported.");
          prm.declare_entry("Script path", "",
                            Patterns::Anything(),
                            "Path to the Python script to execute. Relative paths and the placeholders "
                            "ASPECT_SOURCE_DIR and ASPECT_BINARY_DIR are allowed.");
          prm.declare_entry("Script name", "",
                            Patterns::Anything(),
                            "Name of the Python module to load (without .py extension).");

          const std::string pattern_of_names = "strain rate";

          prm.declare_entry("List of additional ASPECT quantities", "",
                            Patterns::List(Patterns::Selection(pattern_of_names)),
                            "Comma-separated list of additional ASPECT quantities to send to the Landlab model. Default is 'none', "
                            "and the allowed options are: " + pattern_of_names);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void Landlab<dim>::parse_parameters(ParameterHandler &prm)
    {
#ifdef ASPECT_WITH_LANDLAB
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

          AssertThrow(Utilities::fexists(script_path + script_module_name + ".py"),
                      ExcMessage("The specified script path: " + script_path + script_module_name + ".py" + " for Landlab does not exist."));

          additional_named_quantities = Utilities::split_string_list(prm.get ("List of additional ASPECT quantities"));
          AssertThrow(Utilities::has_unique_entries(additional_named_quantities),
                      ExcMessage("The list of strings for the parameter "
                                 "'MeshDeformation/Landlab/List of additional ASPECT quantities' "
                                 "contains entries more than once. This is not allowed. "
                                 "Please check your parameter file."));
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
#else
      (void) prm;
#endif
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(Landlab,
                                           "landlab",
                                           "A mesh deformation plugin that lets a Python script control the "
                                           "deformation of the surface. It is meant for coupling with the landscape evolution "
                                           "code Landlab, but any other script that provides the necessary functions can be used. "
                                           "It is necessary to have Python and numpy with their C APIs installed and that "
                                           "ASPECT_WITH_PYTHON and ASPECT_WITH_LANDLAB are enabled when ASPECT is configured with "
                                           "CMake. ")
  }
}
