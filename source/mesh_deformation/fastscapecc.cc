#include <iostream>
#include <aspect/global.h>

#include <aspect/mesh_deformation/fastscapecc.h>
#include <aspect/geometry_model/box.h>
#include <deal.II/numerics/vector_tools.h>
#include <aspect/postprocess/visualization.h>
#include <ctime>
#include <aspect/simulator.h>
#include <aspect/simulator/assemblers/interface.h>
#include <aspect/simulator_signals.h>

#include <memory>
#include <algorithm> // For std::copy

#include <aspect/geometry_model/spherical_shell.h>
#include <fastscapelib/grid/healpix_grid.hpp>

// namespace fs = fastscapelib;


namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    void FastScapecc<dim>::initialize ()
    {
        const GeometryModel::Box<dim> *box_geometry
          = dynamic_cast<const GeometryModel::Box<dim>*>(&this->get_geometry_model());

        const GeometryModel::SphericalShell<dim> *spherical_geometry
          = dynamic_cast<const GeometryModel::SphericalShell<dim>*>(&this->get_geometry_model());

        if (geometry_type == GeometryType::Box)
        {
            this->get_pcout() << "Box geometry detected. Initializing FastScape for Box geometry..." << std::endl;

            grid_extent[0].first = box_geometry->get_origin()[0];
            grid_extent[0].second = box_geometry->get_extents()[0];
            grid_extent[1].first = box_geometry->get_origin()[1];
            grid_extent[1].second = box_geometry->get_extents()[1];

            nx = repetitions[0] + 1;
            dx = (grid_extent[0].second) / repetitions[0];

            ny = repetitions[1] + 1;
            dy = (grid_extent[1].second) / repetitions[1];

            x_extent = grid_extent[0].second;
            y_extent = grid_extent[1].second;
            array_size = nx * ny;
        }
        else if (geometry_type == GeometryType::SphericalShell)
        {
            this->get_pcout() << "Spherical Shell geometry detected. Initializing FastScape for Spherical Shell geometry..." << std::endl;

          int nsides =(int) sqrt(48 * std::pow(2, (additional_refinement_levels + surface_refinement_difference) * 2) / 12);
          array_size = nsides;
        }
        else
        {
            AssertThrow(false, ExcMessage("FastScapecc plugin only supports Box or Spherical Shell geometries."));
        }
  }

    template <int dim>
    void
    FastScapecc<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                               AffineConstraints<double> &mesh_velocity_constraints,
                                                               const std::set<types::boundary_id> &boundary_ids) const
    {
      if (this->get_timestep_number() == 0)
        return;

      // std::vector<std::vector<double>> temporary_variables;

      TimerOutput::Scope timer_section(this->get_computing_timer(), "FastScape plugin");

      const unsigned int current_timestep = this->get_timestep_number ();

      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
      std::vector<std::vector<double>> temporary_variables(dim+2, std::vector<double>());

      // Get a quadrature rule that exists only on the corners, and increase the refinement if specified.
      const QIterated<dim-1> face_corners (QTrapezoid<1>(),
                                           static_cast<unsigned int>(std::pow(2,additional_refinement_levels+surface_refinement_difference)));

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        face_corners,
                                        update_values |
                                        update_quadrature_points);

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
            if (cell->face(face_no)->at_boundary())
              {
                if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                  continue;

                std::vector<Tensor<1,dim>> vel(face_corners.size());
                fe_face_values.reinit(cell, face_no);
                fe_face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), vel);

                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                    const Point<dim> vertex = fe_face_values.quadrature_point(corner);

                    // Find what x point we're at. Add 1 or 2 depending on if ghost nodes are used.
                    // Subtract the origin point so that it corresponds to an origin of 0,0 in FastScape.
                    const double indx = 1+(vertex(0) - grid_extent[0].first)/dx;

                    // If our x or y index isn't close to a whole number, then it's likely an artifact
                    // from using an over-resolved quadrature rule, in that case ignore it.
                    if (abs(indx - round(indx)) >= precision)
                      continue;


                    // If we're in 2D, we want to take the values and apply them to every row of X points.
                    if (dim == 2)
                      {
                        for (int ys=0; ys<ny; ++ys)
                          {
                            // FastScape indexes from 1 to n, starting at X and Y = 0, and increases
                            // across the X row. At the end of the row, it jumps back to X = 0
                            // and up to the next X row in increasing Y direction. We track
                            // this to correctly place the variables later on.
                            // Nx*ys effectively tells us what row we are in
                            // and then indx tells us what position in that row.
                            const double index = round(indx)+nx*ys;

                            temporary_variables[0].push_back(vertex(dim-1) - grid_extent[dim-1].second);
                            temporary_variables[1].push_back(index-1);

                            for (unsigned int i=0; i<dim; ++i)
                              {
                                // Always convert to m/yr for FastScape
                                temporary_variables[i+2].push_back(vel[corner][i]*year_in_seconds);
                              }
                          }
                      }
                    // 3D case
                    else
                      {
                        // Because indy only gives us the row we're in, we don't need to add 2 for the ghost node.
                        const double indy = 1+(vertex(1) - grid_extent[1].first)/dx;

                        if (abs(indy - round(indy)) >= precision)
                          continue;

                        const double index = round((indy-1))*nx+round(indx);

                        temporary_variables[0].push_back(vertex(dim-1) - grid_extent[dim-1].second);   //z component
                        temporary_variables[1].push_back(index-1);

                        for (unsigned int i=0; i<dim; ++i)
                          {
                            temporary_variables[i+2].push_back(vel[corner][i]*year_in_seconds);
                          }
                      }
                  }
              }

      // // Vector to hold the velocities that represent the change to the surface.
      std::vector<double> V(array_size);

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
      // Initialize the variables that will be sent to FastScape.
      std::vector<double> h(array_size, std::numeric_limits<double>::max());
      std::vector<double> vx(array_size);
      std::vector<double> vy(array_size);
      std::vector<double> vz(array_size);
      std::vector<double> kf(array_size);
      std::vector<double> kd(array_size);
      std::vector<double> h_old(array_size);

      for (unsigned int i=0; i<temporary_variables[1].size(); ++i)
        {
          int index = static_cast<int>(temporary_variables[1][i]);
          h[index] = temporary_variables[0][i];
          vx[index] = temporary_variables[2][i];
          vz[index] = temporary_variables[dim+1][i];

          if (dim == 2)
            vy[index] = 0;
          else
            vy[index] = temporary_variables[3][i];
        }

      for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
        {

          // First, find out the size of the array a process wants to send.
          MPI_Status status;
          MPI_Probe(p, 42, this->get_mpi_communicator(), &status);
          int incoming_size = 0;
          MPI_Get_count(&status, MPI_DOUBLE, &incoming_size);

          // Resize the array so it fits whatever the process sends.
          for (unsigned int i=0; i<temporary_variables.size(); ++i)
            {
              temporary_variables[i].resize(incoming_size);
            }

          for (unsigned int i=0; i<temporary_variables.size(); ++i)
            MPI_Recv(&temporary_variables[i][0], incoming_size, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

          // Now, place the numbers into the correct place based off the index.
          for (unsigned int i=0; i<temporary_variables[1].size(); ++i)
            {
              int index = static_cast<int>(temporary_variables[1][i]);
              h[index] = temporary_variables[0][i];
              vx[index] = temporary_variables[2][i];
              vz[index] = temporary_variables[dim+1][i];

              // In 2D there are no y velocities, so we set them to zero.
              if (dim == 2 )
                vy[index] = 0;
              else
                vy[index] = temporary_variables[3][i];
            }
        }

      // Initialize kf and kd, and check that there are no empty mesh points due to
      // an improperly set maximum_surface_refinement_level, additional_refinement,
      // and surface_refinement_difference
      // int fastscape_mesh_filled = true;
      for (unsigned int i=0; i<array_size; ++i)
        {
          kf[i] = kff;
          kd[i] = kdd;        
        }

      //Execute Fastscape
      // TimerOutput::Scope timer_section(this->get_computing_timer(), "Execute FastScape");
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
        return;

      this->get_pcout() << "   Initializing FastScape... " << (1 + maximum_surface_refinement_level + additional_refinement_levels) <<
                        " levels, cell size: " << dx << " m." << std::endl;
 

      // Keep initial h values so we can calculate velocity later.
      // In the first timestep, h will be given from other processes.
      // In later timesteps, we copy h directly from FastScape.

      std::mt19937 random_number_generator(fs_seed);
      std::uniform_real_distribution<double> random_distribution(-noise_h, noise_h);
      for (unsigned int i = 0; i < array_size; ++i)
        {
          h_old[i] = h[i];

          // Initialize random noise after h_old is set, so aspect sees this initial topography change.
          // We also only want to apply any of these changes to non-ghost nodes,
          // as the ghost nodes will be changed if necessary when set_ghost_nodes is called.
          if (current_timestep == 1)
            {
              // + or - topography based on the initial noise magnitude.
              const double h_seed = random_distribution(random_number_generator);
              h[i] = h[i] + h_seed;
            }
        }

      const double aspect_timestep_in_years = this->get_timestep() / year_in_seconds;


      // Find a FastScape timestep that is below our maximum timestep.
      unsigned int fastscape_iterations = fastscape_steps_per_aspect_step;
      double fastscape_timestep_in_years = aspect_timestep_in_years / fastscape_iterations;
      while (fastscape_timestep_in_years > maximum_fastscape_timestep)
        {
          fastscape_iterations *= 2;
          fastscape_timestep_in_years *= 0.5;
        }

      // auto triangulate = fastscapelib::triangulation(Points);

      // auto grid = fs.TriMesh(tripoints, tritriangles, outlet_idx);

      // raster grid and boundary conditions
      fastscapelib::raster_boundary_status bs(fastscapelib::node_status::fixed_value);


      // if (geometry_type == GeometryType::Box)
      // {
      //     // grid_box =  std::unique_ptr<fastscapelib::raster_grid<>>(fastscapelib::raster_grid<>::from_length(
      //     //     { static_cast<unsigned long>(nx), static_cast<unsigned long>(ny) }, 
      //     //     { x_extent, y_extent }, 
      //     //     bs
      //     // ));

      //     //grid_box = std::make_unique<fastscapelib::raster_grid<>>({ static_cast<unsigned long>(ny), static_cast<unsigned long>(nx)}, { dy, dx}, 
      //     //    bs);
          

      //     //flow_graph_box = std::make_unique<fastscapelib::flow_graph<fastscapelib::raster_grid<>>>(
      //     //     *grid_box, {
      //     //     fastscapelib::single_flow_router(), 
      //     //     fastscapelib::mst_sink_resolver()
      //     // });
      //     // spl_eroder_box = fastscapelib::make_spl_eroder(*flow_graph_box, 2e-4, 0.4, 1, 1e-5);
      //     // diffusion_eroder_box = fastscapelib::make_diffusion_adi_eroder(*grid_box, 0.01);
      // } 
      // else if (geometry_type == GeometryType::SphericalShell)
      // {

        xt::xarray<fastscapelib::node_status> node_status_array(fastscapelib::node_status::fixed_value);

        // xt::xarray<fastscapelib::node_status> node_status_array(nsides, fastscapelib::node_status::fixed_value);
        auto grid = fastscapelib::healpix_grid<>(nsides, node_status_array,6.371e6);
        auto flow_graph = fastscapelib::flow_graph<fastscapelib::healpix_grid<>>(
            grid, {
            fastscapelib::single_flow_router(), 
            fastscapelib::mst_sink_resolver()}
        );
        auto spl_eroder = fastscapelib::make_spl_eroder(flow_graph, 2e-4, 0.4, 1, 1e-5);
      // }

      xt::xarray<double> uplifted_elevation ;
      xt::xarray<double> drainage_area = xt::zeros<double>(grid.shape());
      xt::xarray<double> sediment_flux = xt::zeros<double>(grid.shape());

            // std::cout << "uplift_rate_in_m_year[1]"<<uplift_rate_in_m_year[100]<<std::endl;

     std::vector<double> uplift_rate_in_m_year(vz.size());
      for (size_t i = 0; i < vz.size(); ++i) {
          uplift_rate_in_m_year[i] = vz[i];
          //  / year_in_seconds;
      }
      std::vector<std::size_t> shape = { static_cast<unsigned long>(nx), static_cast<unsigned long>(ny) };
      auto uplift_rate = xt::adapt(uplift_rate_in_m_year, shape);


      auto elevation = xt::adapt(h, shape);
      auto elevation_old = xt::adapt(h_old, shape);

      // xt::xarray<double> elevation(vec_shape.shape(), h);
      // xt::xarray<double> elevation_old(vec_shape.shape(), h_old);


      // std::vector<double> uplift_rate_in_m_year(vz.size());
      // for (size_t i = 0; i < vz.size(); ++i) {
      //     uplift_rate_in_m_year[i] = vz[i] / year_in_seconds;
      // }

      // run model
      //
      //double dt = 2e4;

      for (unsigned int fastscape_iteration = 0; fastscape_iteration < fastscape_iterations; ++fastscape_iteration)
        {
          std::cout << "Fastscape ite "<<fastscape_iteration<<std::endl;
          // apply uplift
          uplifted_elevation = elevation + fastscape_timestep_in_years * uplift_rate;
          // flow routing
          flow_graph.update_routes(uplifted_elevation);
          // flow accumulation (drainage area)
          flow_graph.accumulate(drainage_area, 1.0);
          // apply channel erosion then hillslope diffusion
          auto spl_erosion = spl_eroder.erode(uplifted_elevation, drainage_area, fastscape_timestep_in_years);

          //calculate the cumulated erosion flux
          sediment_flux = flow_graph.accumulate(spl_erosion);
          // if (geometry_type == GeometryType::Box)
          // {
          // auto diff_erosion = diffusion_eroder_box.erode(uplifted_elevation - spl_erosion, fastscape_timestep_in_years);
          // update topography
          // auto elevation = uplifted_elevation - spl_erosion - diff_erosion;
          // }else if (geometry_type == GeometryType::SphericalShell)
          elevation = uplifted_elevation - spl_erosion;
          }

          std::vector<double> uplifted_erosion_std(uplifted_elevation.begin(), uplifted_elevation.end());
          std::vector<double> elevation_std(elevation.begin(), elevation.end());
          std::vector<double> elevation_old_std(elevation_old.begin(), elevation_old.end());

          std::vector<double> uplift_rate_std(uplift_rate.begin(), uplift_rate.end());       
        //}

      // Find out our velocities from the change in height.
      // Where V is a vector of array size that exists on all processes.
        
        std::cout<<"Timestep : "<<this->get_timestep()<<std::endl;
      for (unsigned int i=0; i<array_size; ++i)
        {
          V[i] = (elevation_std[i] - elevation_old_std[i])/ (this->get_timestep()/ year_in_seconds);
        }
          MPI_Bcast(&V[0], array_size, MPI_DOUBLE, 0, this->get_mpi_communicator());
      }
      else{
          for (unsigned int i=0; i<temporary_variables.size(); ++i)
            MPI_Ssend(&temporary_variables[i][0], temporary_variables[1].size(), MPI_DOUBLE, 0, 42, this->get_mpi_communicator());
          
          MPI_Bcast(&V[0], array_size, MPI_DOUBLE, 0, this->get_mpi_communicator()); 
      }  

      // Get the sizes needed for a data table of the mesh velocities.
      TableIndices<dim> size_idx;
      for (unsigned int d=0; d<dim; ++d)
        {
          size_idx[d] = table_intervals[d]+1;
        }

      // Initialize a table to hold all velocity values that will be interpolated back to ASPECT.
      Table<dim,double> velocity_table = fill_data_table(V, size_idx, nx, ny);

      // As our grid_extent variable end points do not account for the change related to an origin
      // not at 0, we adjust this here into an interpolation extent.
      std::array<std::pair<double,double>,dim> interpolation_extent;
      for (unsigned int i=0; i<dim; ++i)
        {
          interpolation_extent[i].first = grid_extent[i].first;
          interpolation_extent[i].second = (grid_extent[i].second + grid_extent[i].first);
        }
    

      //Functions::InterpolatedUniformGridData<dim> *velocities;
      Functions::InterpolatedUniformGridData<dim> velocities (interpolation_extent,
                                                              table_intervals,
                                                              velocity_table);

      VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
        [&](const Point<dim> &p) -> double
      {
        return velocities.value(p);
      },
      dim-1,
      dim);

      VectorTools::interpolate_boundary_values (mesh_deformation_dof_handler,
                                                *boundary_ids.begin(),
                                                vector_function_object,
                                                mesh_velocity_constraints);
    }


    template <int dim>
    Table<dim,double>
    FastScapecc<dim>::fill_data_table(std::vector<double> &values,
                                      TableIndices<dim> &size_idx,
                                      const int &nx,
                                      const int &ny) const
    {
      // Create data table based off of the given size.
      Table<dim,double> data_table;
      data_table.TableBase<dim,double>::reinit(size_idx);
      TableIndices<dim> idx;
//      Assert(values.size() == data_table.size()[0]*data_table.size()[1]*data_table.size()[2],
//                     ExcMessage("The size of the data table does not match the size of the values.");

      // Loop through the data table and fill it with the velocities from FastScape.

      // Indexes through z, y, and then x.
      for (unsigned int k=0; k<data_table.size()[2]; ++k)
        {
          idx[2] = k;

          for (unsigned int i=0; i<data_table.size()[1]; ++i)
            {
              idx[1] = i;

              for (unsigned int j=0; j<data_table.size()[0]; ++j)
                {
                  idx[0] = j;

                  // Convert back to m/s.
                  data_table(idx) = values[nx*i+j] / year_in_seconds;
                    

                }
            }
        }
      return data_table;
    }



    template <int dim>
    void FastScapecc<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
          // Declare parameters for the Box geometry
          if (prm.get("Model name") == "box")
          {
              prm.enter_subsection("Box");
              {
                  prm.declare_entry("X repetitions", "1", Patterns::Integer(1),
                                    "Number of cells in the X direction.");
                  prm.declare_entry("Y repetitions", "1", Patterns::Integer(1),
                                    "Number of cells in the Y direction.");
                  prm.declare_entry("Z repetitions", "1", Patterns::Integer(1),
                                    "Number of cells in the Z direction.");
              }
              prm.leave_subsection();  // End of Box
          }
          // Declare parameters for the Spherical shell geometry
          else if (prm.get("Model name") == "spherical shell")
          {
              prm.enter_subsection("Spherical shell");
              {
                  prm.declare_entry("Inner radius", "3481000", Patterns::Double(0),
                                    "The inner radius of the spherical shell.");
                  prm.declare_entry("Outer radius", "6336000", Patterns::Double(0),
                                    "The outer radius of the spherical shell.");
                  prm.declare_entry("Opening angle", "360", Patterns::Double(0, 360),
                                    "The opening angle of the spherical shell in degrees.");
              }
              prm.leave_subsection();  
          }
      }
      prm.leave_subsection();  




      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Fastscapecc");
        {
          prm.declare_entry("Number of steps", "10",
                            Patterns::Integer(),
                            "Number of steps per ASPECT timestep");
          prm.declare_entry("Maximum timestep", "10e3",
                            Patterns::Double(0),
                            "Maximum timestep for FastScape. Units: $\\{yrs}$");
          prm.declare_entry("Additional fastscape refinement levels", "0",
                            Patterns::Integer(),
                            "How many levels above ASPECT FastScape should be refined.");
          prm.declare_entry ("Use center slice for 2d", "false",
                             Patterns::Bool (),
                             "If this is set to true, then a 2D model will only consider the "
                             "center slice FastScape gives. If set to false, then aspect will"
                             "average the mesh along Y excluding the ghost nodes.");
          prm.declare_entry("Fastscape seed", "1000",
                            Patterns::Integer(),
                            "Seed used for adding an initial noise to FastScape topography based on the initial noise magnitude.");
          prm.declare_entry("Maximum surface refinement level", "1",
                            Patterns::Integer(),
                            "This should be set to the highest ASPECT refinement level expected at the surface.");
          prm.declare_entry("Surface refinement difference", "0",
                            Patterns::Integer(),
                            "The difference between the lowest and highest refinement level at the surface. E.g., if three resolution "
                            "levels are expected, this would be set to 2.");
          prm.declare_entry("Y extent in 2d", "100000",
                            Patterns::Double(),
                            "FastScape Y extent when using a 2D ASPECT model. Units: $\\{m}$");
          prm.declare_entry ("Use velocities", "true",
                             Patterns::Bool (),
                             "Flag to use FastScape advection and uplift.");
          prm.declare_entry("Precision", "0.001",
                            Patterns::Double(),
                            "Precision value for how close a ASPECT node must be to the FastScape node for the value to be transferred.");
          prm.declare_entry("Initial noise magnitude", "5",
                            Patterns::Double(),
                            "Maximum topography change from the initial noise. Units: $\\{m}$");

          prm.enter_subsection ("Boundary conditions");
          {
            prm.declare_entry ("Front", "1",
                               Patterns::Integer (0, 1),
                               "Front (bottom) boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Right", "1",
                               Patterns::Integer (0, 1),
                               "Right boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Back", "1",
                               Patterns::Integer (0, 1),
                               "Back (top) boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Left", "1",
                               Patterns::Integer (0, 1),
                               "Left boundary condition, where 1 is fixed and 0 is reflective.");
          }
          prm.leave_subsection();

          prm.enter_subsection ("Erosional parameters");
          {
            prm.declare_entry("Drainage area exponent", "0.4",
                              Patterns::Double(),
                              "Exponent for drainage area.");
            prm.declare_entry("Slope exponent", "1",
                              Patterns::Double(),
                              "The  slope  exponent  for  SPL (n).  Generally  m/n  should  equal  approximately 0.4");
            prm.declare_entry("Bedrock river incision rate", "1e-5",
                              Patterns::Double(),
                              "River incision rate for bedrock in the Stream Power Law. Units: $\\{m^(1-2*drainage_area_exponent)/yr}$");
            prm.declare_entry("Sediment river incision rate", "-1",
                              Patterns::Double(),
                              "River incision rate for sediment in the Stream Power Law. -1 sets this to the bedrock river incision rate. Units: $\\{m^(1-2*drainage_area_exponent)/yr}$ ");
            prm.declare_entry("Bedrock diffusivity", "1e-2",
                              Patterns::Double(),
                              "Transport coefficient (diffusivity) for bedrock. Units: $\\{m^2/yr}$ ");
            prm.declare_entry("Sediment diffusivity", "-1",
                              Patterns::Double(),
                              "Transport coefficient (diffusivity) for sediment. -1 sets this to the bedrock diffusivity. Units: $\\{m^2/yr}$");
            prm.declare_entry("Elevation factor", "1",
                              Patterns::Double(),
                              "Amount to multiply kf and kd by past given orographic elevation control.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void FastScapecc<dim>::parse_parameters(ParameterHandler &prm)
    {
      end_time = prm.get_double ("End time");
      if (prm.get_bool ("Use years in output instead of seconds") == true)
        end_time *= year_in_seconds;

      // prm.enter_subsection("Geometry model");
      // {
      //   prm.enter_subsection("Box");
      //   {
      //     repetitions[0] = prm.get_integer ("X repetitions");
      //     if (dim >= 2)
      //       {
      //         repetitions[1] = prm.get_integer ("Y repetitions");
      //       }
      //     if (dim >= 3)
      //       {
      //         repetitions[dim-1] = prm.get_integer ("Z repetitions");
      //       }
      //   }
      //   prm.leave_subsection();
      // }
      // prm.leave_subsection();



      prm.enter_subsection("Geometry model");
      {
          // Parse parameters for the Box geometry
          if (prm.get("Model name") == "box")
          {
              prm.enter_subsection("Box");
              {
                  repetitions[0] = prm.get_integer("X repetitions");
                  repetitions[1] = prm.get_integer("Y repetitions");
                  if (dim == 3)
                      repetitions[2] = prm.get_integer("Z repetitions");
              }
              prm.leave_subsection();  // End of Box
          }
          // Parse parameters for the Spherical shell geometry
          else if (prm.get("Model name") == "spherical shell")
          {
              prm.enter_subsection("Spherical shell");
              {
                  inner_radius = prm.get_double("Inner radius");
                  outer_radius = prm.get_double("Outer radius");
                  opening_angle = prm.get_double("Opening angle");
              }
              prm.leave_subsection();  
          }
      }
      prm.leave_subsection();  

      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("Fastscapecc");
        {
          fastscape_steps_per_aspect_step = prm.get_integer("Number of steps");
          maximum_fastscape_timestep = prm.get_double("Maximum timestep");
          additional_refinement_levels = prm.get_integer("Additional fastscape refinement levels");
          center_slice = prm.get_bool("Use center slice for 2d");
          fs_seed = prm.get_integer("Fastscape seed");
          maximum_surface_refinement_level = prm.get_integer("Maximum surface refinement level");
          surface_refinement_difference = prm.get_integer("Surface refinement difference");
          y_extent_2d = prm.get_double("Y extent in 2d");
          precision = prm.get_double("Precision");
          noise_h = prm.get_double("Initial noise magnitude");

          if (!this->convert_output_to_years())
            {
              maximum_fastscape_timestep /= year_in_seconds;
            }

          prm.enter_subsection("Boundary conditions");
          {
            bottom = prm.get_integer("Front");
            right = prm.get_integer("Right");
            top = prm.get_integer("Back");
            left = prm.get_integer("Left");
          }
          prm.leave_subsection();

          prm.enter_subsection("Erosional parameters");
          {
            m = prm.get_double("Drainage area exponent");
            n = prm.get_double("Slope exponent");
            kfsed = prm.get_double("Sediment river incision rate");
            kff = prm.get_double("Bedrock river incision rate");
            kdsed = prm.get_double("Sediment diffusivity");
            kdd = prm.get_double("Bedrock diffusivity");

            // if (!this->convert_output_to_years())
            //   {
            //     kff *= year_in_seconds;
            //     kdd *= year_in_seconds;
            //     kfsed *= year_in_seconds;
            //     kdsed *= year_in_seconds;
            //   }

          }
          prm.leave_subsection();

        }
        prm.leave_subsection();
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
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(FastScapecc,
                                           "fastscapecc",
                                           "A plugin, which prescribes the surface mesh to "
                                           "deform according to an analytically prescribed "
                                           "function. Note that the function prescribes a "
                                           "deformation velocity, i.e. the return value of "
                                           "this plugin is later multiplied by the time step length "
                                           "to compute the displacement increment in this time step. "
                                           "The format of the "
                                           "functions follows the syntax understood by the "
                                           "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
