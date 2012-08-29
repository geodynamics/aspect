/*
 Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */

#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/numerics/data_out.h>
#include <aspect/global.h>
#include <aspect/postprocess/tracer.h>
#include <aspect/simulator.h>

#ifdef DEAL_II_HAVE_HDF5
#include <hdf5.h>
#endif

namespace aspect
{
  namespace Postprocess
  {
    // Initialize a particle using the data in an MPI_Particle struct
    template <int dim>
    Particle<dim>::Particle(const MPI_Particle<dim> &particle_data) : _is_local(false)
    {
      for (unsigned int d=0; d<dim; ++d)
        {
          _pos(d) = particle_data.pos[d];
          _pos0(d) = particle_data.pos0[d];
          _velocity(d) = particle_data.velocity[d];
        }

      _id = particle_data.id;
    }

    // Write the data for this particle into an MPI structure
    template <int dim>
    MPI_Particle<dim> Particle<dim>::mpi_data(void) const
    {
      MPI_Particle<dim> new_data;

      for (unsigned int d=0; d<dim; ++d)
        {
          new_data.pos[d] = _pos(d);
          new_data.pos0[d] = _pos0(d);
          new_data.velocity[d] = _velocity(d);
        }
      new_data.id = _id;

      return new_data;
    }

    // Write the data for this particle into an HDF5 structure
    template <int dim>
    HDF5_Particle<dim> Particle<dim>::hdf5_data(void) const
    {
      HDF5_Particle<dim> new_data;

      for (unsigned int d=0; d<dim; ++d)
        {
          new_data.pos[d] = _pos(d);
          new_data.velocity[d] = _velocity(d);
        }
      new_data.id = _id;

      return new_data;
    }

    // TODO: add better error checking to MPI calls
    template <int dim>
    void ParticleSet<dim>::setup_mpi(MPI_Comm comm_world)
    {
      int       block_lens[4];
      MPI_Aint    indices[4];
      MPI_Datatype  old_types[4];
      int       res;

      block_lens[0] = dim;    // Position
      block_lens[1] = dim;    // Original position
      block_lens[2] = dim;    // Velocity
      block_lens[3] = 1;      // ID

      indices[0] = 0;
      indices[1] = indices[0]+sizeof(double)*dim;
      indices[2] = indices[1]+sizeof(double)*dim;
      indices[3] = indices[2]+sizeof(double)*dim;

      old_types[0] = MPI_DOUBLE;
      old_types[1] = MPI_DOUBLE;
      old_types[2] = MPI_DOUBLE;
      old_types[3] = MPI_INT;

      res = MPI_Type_struct(4, block_lens, indices, old_types, &particle_type);
      if (res != MPI_SUCCESS) exit(-1);

      res = MPI_Type_commit(&particle_type);
      if (res != MPI_SUCCESS) exit(-1);

      // Determine the size of the MPI comm world
      MPI_Comm_size(comm_world, &world_size);
      MPI_Comm_rank(comm_world, &self_rank);

      // Initialize send/recv structures appropriately
      num_send = new int[world_size];
      num_recv = new int[world_size];
      send_offset = new int[world_size];
      recv_offset = new int[world_size];
      send_reqs = new MPI_Request[world_size];
      recv_reqs = new MPI_Request[world_size];
    }

    template <int dim>
    ParticleSet<dim>::~ParticleSet(void)
    {
      if (world_size) MPI_Type_free(&particle_type);

      if (num_send) delete num_send;
      if (num_recv) delete num_recv;
      if (send_offset) delete send_offset;
      if (recv_offset) delete recv_offset;
      if (send_reqs) delete send_reqs;
      if (recv_reqs) delete recv_reqs;
    }


    template <int dim>
    std::pair<std::string,std::string> ParticleSet<dim>::execute (TableHandler &statistics)
    {
      std::string     result_string = "done.", data_file_name;
      bool            output_data;

      output_data = false;
      if (!initialized)
        {
          next_data_output_time = this->get_time();
          setup_mpi(this->get_mpi_communicator());
          initialized = true;
          generate_particles(this->get_triangulation(), num_initial_tracers);
          this->get_triangulation().signals.post_refinement.connect(std_cxx1x::bind(&ParticleSet<dim>::mesh_changed,
                                                                                    std_cxx1x::ref(*this)));
        }
      advect_particles(this->get_triangulation(),
                       this->get_timestep(),
                       this->get_solution(),
                       this->get_dof_handler(),
                       this->get_mapping());

      if (data_output_format != "none" && this->get_time() >= next_data_output_time)
        {
          set_next_data_output_time (this->get_time());
          data_file_name = output_particle_data(this->get_output_directory(), this->get_triangulation());
          output_data = true;
        }
      if (output_data) result_string += " Wrote particle data: " + data_file_name + ".";

      triangulation_changed = false;

      return std::make_pair("Advecting particles...", result_string);
    }

    // TODO: determine file format, write this function
    /*template <int dim>
    void ParticleSet<dim>::read_particles_from_file(parallel::distributed::Triangulation<dim> &triangulation,
                                                    std::string filename) {
        std::ifstream pos_file(filename.c_str());

        if (!pos_file.is_open()) exit(-1);  // error out better than this

        while (!pos_file.eof()) {

        }

        pos_file.close();
    }*/

    // Generate a set of particles in the specified triangulation
    // TODO: fix the numbering scheme so we have exactly the right number of particles for all processor configurations
    // TODO: fix the particle system so it works even with processors assigned 0 cells
    template <int dim>
    void ParticleSet<dim>::generate_particles(const parallel::distributed::Triangulation<dim> &triangulation,
                                              double total_particles)
    {
      double      total_volume, local_volume, subdomain_fraction, start_fraction, end_fraction;
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator  it;

      global_sum_particles = total_particles;
      // Calculate the number of particles in this domain as a fraction of total volume
      total_volume = local_volume = 0;
      for (it=triangulation.begin_active(); it!=triangulation.end(); ++it)
        {
          double cell_volume = it->measure();
          AssertThrow (cell_volume != 0, ExcMessage ("Found cell with zero volume."));
          if (it->is_locally_owned()) local_volume += cell_volume;
        }

      // Sum the local volumes over all nodes
      MPI_Allreduce(&local_volume, &total_volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      // Assign this subdomain the appropriate fraction
      subdomain_fraction = local_volume/total_volume;

      // Sum the subdomain fractions so we don't miss particles from rounding and to create unique IDs
      MPI_Scan(&subdomain_fraction, &end_fraction, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      start_fraction = end_fraction-subdomain_fraction;

      // Calculate start and end IDs so there are no gaps
      const unsigned int  start_id = static_cast<unsigned int>(std::floor(start_fraction*total_particles));
      const unsigned int  end_id   = static_cast<unsigned int>(std::floor(end_fraction*total_particles));
      const unsigned int  subdomain_particles = end_id - start_id;

      generate_particles_in_subdomain(triangulation, subdomain_particles, start_id);
    }

    // Generate a set of particles uniformly distributed within the specified triangulation.
    // This is done using "roulette wheel" style selection weighted by cell size.
    // We do cell-by-cell assignment of particles because the decomposition of the mesh may
    // result in a highly non-rectangular local mesh which makes uniform particle distribution difficult.
    template <int dim>
    void ParticleSet<dim>::generate_particles_in_subdomain(const parallel::distributed::Triangulation<dim> &triangulation,
                                                           unsigned int num_particles,
                                                           unsigned int start_id)
    {
      unsigned int          i, d, v, num_tries, cur_id;
      double                total_volume, roulette_spin;
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator  it;
      std::map<double, LevelInd>        roulette_wheel;
      const unsigned int n_vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;
      Point<dim>            pt, max_bounds, min_bounds;
      LevelInd              select_cell;

      // Create the roulette wheel based on volumes of local cells
      total_volume = 0;
      for (it=triangulation.begin_active(); it!=triangulation.end(); ++it)
        {
          if (it->is_locally_owned())
            {
              // Assign an index to each active cell for selection purposes
              total_volume += it->measure();
              // Save the cell index and level for later access
              roulette_wheel.insert(std::make_pair(total_volume, std::make_pair(it->level(), it->index())));
            }
        }

      // Pick cells and assign particles at random points inside them
      cur_id = start_id;
      for (i=0; i<num_particles; ++i)
        {
          // Select a cell based on relative volume
          roulette_spin = total_volume*drand48();
          select_cell = roulette_wheel.lower_bound(roulette_spin)->second;
          it = typename parallel::distributed::Triangulation<dim>::active_cell_iterator(&triangulation, select_cell.first, select_cell.second);

          // Get the bounds of the cell defined by the vertices
          for (d=0; d<dim; ++d)
            {
              min_bounds[d] = INFINITY;
              max_bounds[d] = -INFINITY;
            }
          for (v=0; v<n_vertices_per_cell; ++v)
            {
              pt = it->vertex(v);
              for (d=0; d<dim; ++d)
                {
                  min_bounds[d] = fmin(pt[d], min_bounds[d]);
                  max_bounds[d] = fmax(pt[d], max_bounds[d]);
                }
            }

          // Generate random points in these bounds until one is within the cell
          num_tries = 0;
          while (num_tries < 100)
            {
              for (d=0; d<dim; ++d)
                {
                  pt[d] = drand48()*(max_bounds[d]-min_bounds[d]) + min_bounds[d];
                }
              try
                {
                  if (it->point_inside(pt)) break;
                }
              catch (...)
                {
                  // Debugging output, remove when Q4 mapping 3D sphere problem is resolved
                  //std::cerr << "ooo eee ooo aaa aaa " << pt << " " << select_cell.first << " " << select_cell.second << std::endl;
                  //for (int z=0;z<8;++z) std::cerr << "V" << z <<": " << it->vertex(z) << ", ";
                  //std::cerr << std::endl;
                  MPI_Abort(MPI_COMM_WORLD, 1);
                }
              num_tries++;
            }
          AssertThrow (num_tries < 100, ExcMessage ("Couldn't generate particle (unusual cell shape?)."));

          // Add the generated particle to the set
          Particle<dim> new_particle(pt, cur_id);
          particles.insert(std::make_pair(select_cell, new_particle));

          cur_id++;
        }
    }

    // Finds the cell the particle is contained in and returns the corresponding level/index
    template <int dim>
    LevelInd ParticleSet<dim>::find_cell(const Mapping<dim> &mapping,
                                         const parallel::distributed::Triangulation<dim> &triangulation,
                                         Particle<dim> &particle,
                                         LevelInd cur_cell)
    {
      typename parallel::distributed::Triangulation<dim>::cell_iterator         it, found_cell;
      typename parallel::distributed::Triangulation<dim>::active_cell_iterator  ait;
      LevelInd    res;

      // First check the last recorded cell since particles will generally stay in the same area
      if (!triangulation_changed)
        {
          found_cell = typename parallel::distributed::Triangulation<dim>::cell_iterator(&triangulation, cur_cell.first, cur_cell.second);
          if (found_cell != triangulation.end() && found_cell->point_inside(particle.getPosition()) && found_cell->active())
            {
              // If the cell is active, we're at the finest level of refinement and can finish
              particle.setLocal(found_cell->is_locally_owned());
              return cur_cell;
            }
        }

      // Check all the cells on level 0 and recurse down
      for (it=triangulation.begin(0); it!=triangulation.end(0); ++it)
        {
          res = recursive_find_cell(triangulation, particle, std::make_pair(it->level(), it->index()));
          if (res.first != -1 && res.second != -1) return res;
        }

      // If we couldn't find it there, we need to check the active cells
      // since the particle could be on a curved boundary not included in the
      // coarse grid
      for (ait=triangulation.begin_active(); ait!=triangulation.end(); ++ait)
        {
          if (ait->point_inside(particle.getPosition()))
            {
              particle.setLocal(ait->is_locally_owned());
              return std::make_pair(ait->level(), ait->index());
            }
        }

      // If it failed all these tests, the particle is outside the mesh
      particle.setLocal(false);
      return std::make_pair(-1, -1);
    }

    // Recursively determines which cell the given particle belongs to
    // Returns true if the particle is in the specified cell and sets the particle
    // cell information appropriately, false otherwise
    template <int dim>
    LevelInd ParticleSet<dim>::recursive_find_cell(const parallel::distributed::Triangulation<dim> &triangulation,
                                                   Particle<dim> &particle,
                                                   LevelInd cur_cell)
    {
      typename parallel::distributed::Triangulation<dim>::cell_iterator it, found_cell, child_cell;
      unsigned int    child_num;
      LevelInd        res, child_li;

      // If the particle is in the specified cell
      found_cell = typename parallel::distributed::Triangulation<dim>::cell_iterator(&triangulation, cur_cell.first, cur_cell.second);
      if (found_cell != triangulation.end() && found_cell->point_inside(particle.getPosition()))
        {
          // If the cell is active, we're at the finest level of refinement and can finish
          if (found_cell->active())
            {
              particle.setLocal(found_cell->is_locally_owned());
              return cur_cell;
            }
          else
            {
              // Otherwise we need to search deeper
              for (child_num=0; child_num<found_cell->n_children(); ++child_num)
                {
                  child_cell = found_cell->child(child_num);
                  child_li = LevelInd(child_cell->level(), child_cell->index());
                  res = recursive_find_cell(triangulation, particle, child_li);
                  if (res.first != -1 && res.second != -1) return res;
                }
            }
        }

      // If we still can't find it, return false
      return LevelInd(-1, -1);
    }

    template <int dim>
    void ParticleSet<dim>::get_particle_velocities(const parallel::distributed::Triangulation<dim> &triangulation,
                                                   const TrilinosWrappers::MPI::BlockVector &solution,
                                                   const DoFHandler<dim> &dh,
                                                   const Mapping<dim> &mapping,
                                                   std::vector<Point<dim> > &velocities)
    {
      Vector<double>                single_res(dim+2);
      std::vector<Vector<double> >  result;
      Point<dim>                    velocity;
      unsigned int                  i, num_particles;
      LevelInd                      cur_cell;
      typename ParticleMap::iterator  it;
      typename DoFHandler<dim>::active_cell_iterator  found_cell;
      std::vector<Point<dim> >      particle_points;

      // Prepare the field function
      Functions::FEFieldFunction<dim, DoFHandler<dim>, TrilinosWrappers::MPI::BlockVector> fe_value(dh, solution, mapping);

      // Get the velocity for each cell at a time so we can take advantage of knowing the active cell
      for (it=particles.begin(); it!=particles.end();)
        {
          // Get the current cell
          cur_cell = it->first;

          // Resize the vectors to the number of particles in this cell
          num_particles = particles.count(cur_cell);
          particle_points.resize(num_particles);
          result.resize(num_particles, single_res);

          // Get a vector of the point positions in this cell
          i=0;
          while (it != particles.end() && it->first == cur_cell)
            {
              particle_points[i++] = it->second.getPosition();
              it++;
            }

          // Get the cell the particle is in
          found_cell = typename DoFHandler<dim>::active_cell_iterator(&triangulation, cur_cell.first, cur_cell.second, &dh);

          // Interpolate the velocity field for each of the particles
          fe_value.set_active_cell(found_cell);
          fe_value.vector_value_list(particle_points, result);

          // Copy the resulting velocities to the appropriate vector
          for (i=0; i<num_particles; ++i)
            {
              for (int d=0; d<dim; ++d) velocity(d) = result[i](d);
              velocities.push_back(velocity);
            }
        }
    }

    template <int dim>
    unsigned int ParticleSet<dim>::get_global_particle_count(void)
    {
      unsigned int    local_particles = particles.size();
      unsigned int    global_particles;
      int             res;

      res = MPI_Allreduce(&local_particles, &global_particles, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      if (res != MPI_SUCCESS) exit(-1);
      return global_particles;
    }

    // Ensures that particles are not lost in the simulation
    template <int dim>
    void ParticleSet<dim>::check_particle_count(void)
    {
      unsigned int    global_particles = get_global_particle_count();

      AssertThrow (global_particles==global_sum_particles, ExcMessage ("Particle count changed."));
    }

    template <int dim>
    void ParticleSet<dim>::move_particles_back_in_mesh(const Mapping<dim> &mapping,
                                                       const parallel::distributed::Triangulation<dim> &triangulation)
    {
      // TODO: fix this to work with arbitrary meshes
    }

    /*
        Transfer particles that have crossed domain boundaries to other processors
        Because domains can change drastically during mesh refinement, particle transfer occurs as follows:
            - Each domain finds particles that have fallen outside it
            - If the new particle position is in a known domain (e.g. artificial cells), send the particle only to that domain
            - If the new particle position is not in a known domain (e.g. due to mesh refinement), send the particle to all domains
            - The position of each of these particles is broadcast to the specified domains
            - Each domain determines which of the broadcast particles is in itself, keeps these and deletes the others
            - TODO: handle particles outside any domain
            - TODO: if we know the domain of a particle (e.g. bordering domains), send it only to that domain
    */
    template <int dim>
    void ParticleSet<dim>::send_recv_particles(const Mapping<dim> &mapping,
                                               const parallel::distributed::Triangulation<dim> &triangulation)
    {
      typename ParticleMap::iterator  it;
      typename parallel::distributed::Triangulation<dim>::cell_iterator found_cell;
      int                 i, rank;
      std::vector<Particle<dim> >     send_particles;
      typename std::vector<Particle<dim> >::iterator     sit;
      MPI_Particle<dim>         *send_data, *recv_data;

      // Go through the particles and take out those which need to be moved to another processor
      for (it=particles.begin(); it!=particles.end();)
        {
          if (!it->second.isLocal())
            {
              send_particles.push_back(it->second);
              particles.erase(it++);
            }
          else
            {
              ++it;
            }
        }

      // Determine the total number of particles we will send to other processors
      total_send = send_particles.size();
      for (rank=0; rank<world_size; ++rank)
        {
          if (rank != self_rank) num_send[rank] = total_send;
          else num_send[rank] = 0;
          send_offset[rank] = 0;
        }

      // Notify other processors how many particles we will be sending
      MPI_Alltoall(num_send, 1, MPI_INT, num_recv, 1, MPI_INT, MPI_COMM_WORLD);

      total_recv = 0;
      for (rank=0; rank<world_size; ++rank)
        {
          recv_offset[rank] = total_recv;
          total_recv += num_recv[rank];
        }

      // Allocate space for sending and receiving particle data
      send_data = new MPI_Particle<dim>[total_send];
      recv_data = new MPI_Particle<dim>[total_recv];

      // Copy the particle data into the send array
      for (i=0,sit=send_particles.begin(); sit!=send_particles.end(); ++sit,++i)
        {
          send_data[i] = sit->mpi_data();
        }

      // Exchange the particle data between domains
      MPI_Alltoallv(send_data, num_send, send_offset, particle_type,
                    recv_data, num_recv, recv_offset, particle_type,
                    MPI_COMM_WORLD);

      int put_in_domain = 0;
      // Put the received particles into the domain if they are in the triangulation
      for (i=0; i<total_recv; ++i)
        {
          Particle<dim>       recv_particle(recv_data[i]);
          LevelInd            found_cell;
          found_cell = find_cell(mapping, triangulation, recv_particle, std::make_pair(-1,-1));
          if (recv_particle.isLocal())
            {
              put_in_domain++;
              particles.insert(std::make_pair(found_cell, recv_particle));
            }
        }

      delete send_data;
      delete recv_data;
    }

    // Advect particles based on the velocity field using a RK2 method.
    // Currently investigating whether this is the optimal method or not.
    // Will likely be replaced by a hybrid method.
    template <int dim>
    void ParticleSet<dim>::advect_particles(const parallel::distributed::Triangulation<dim> &triangulation,
                                            double timestep,
                                            const TrilinosWrappers::MPI::BlockVector &solution,
                                            const DoFHandler<dim> &dh,
                                            const Mapping<dim> &mapping)
    {
      typename ParticleMap::iterator  it;
      std::vector<Point<dim> >        velocities;
      Point<dim>                      cur_pos, orig_pos, new_pos, old_vel, new_vel;
      unsigned int                    i;
      LevelInd                        found_cell;
      ParticleMap                     tmp_map;

      // Find the cells that the particles moved to
      tmp_map.clear();
      for (it=particles.begin(); it!=particles.end();)
        {
          found_cell = find_cell(mapping, triangulation, it->second, it->first);
          if (found_cell != it->first)
            {
              tmp_map.insert(std::make_pair(found_cell, it->second));
              particles.erase(it++);
            }
          else ++it;
        }
      particles.insert(tmp_map.begin(),tmp_map.end());

      // If the triangulation changed, we may need to move particles between processors
      if (triangulation_changed)
        {
          send_recv_particles(mapping, triangulation);
        }

      // If particles fell out of the mesh, put them back in at the closest point in the mesh
      move_particles_back_in_mesh(mapping, triangulation);

      // Starting out, particles must know which cells they belong to
      // Using this we can quickly interpolate the velocities
      get_particle_velocities(triangulation, solution, dh, mapping, velocities);

      // Move the particles based on the velocity at their current position
      for (it=particles.begin(),i=0; it!=particles.end(); ++it,++i)
        {
          // Save current position and velocity
          cur_pos = it->second.getPosition();
          it->second.setVelocity(velocities[i]);
          it->second.setOriginalPos(cur_pos);

          // Move the particle
          cur_pos += timestep*velocities[i];
          it->second.setPosition(cur_pos);
        }

      // Find the cells that the particles moved to
      tmp_map.clear();
      for (it=particles.begin(); it!=particles.end();)
        {
          found_cell = find_cell(mapping, triangulation, it->second, it->first);
          if (found_cell != it->first)
            {
              tmp_map.insert(std::make_pair(found_cell, it->second));
              particles.erase(it++);
            }
          else ++it;
        }
      particles.insert(tmp_map.begin(),tmp_map.end());

      // If particles fell out of the mesh, put them back in at the closest point in the mesh
      move_particles_back_in_mesh(mapping, triangulation);

      // Swap particles between processors if needed
      send_recv_particles(mapping, triangulation);

      // Get the velocities at the new positions
      velocities.clear();
      get_particle_velocities(triangulation, solution, dh, mapping, velocities);

      // Estimate actual position based on old velocity/position and new velocity/position
      for (it=particles.begin(),i=0; it!=particles.end(); ++it,++i)
        {
          orig_pos = it->second.getOriginalPos();
          cur_pos = it->second.getPosition();
          old_vel = it->second.getVelocity();
          new_vel = velocities[i];
          new_pos = orig_pos + timestep * 0.5 * (old_vel + new_vel);
          it->second.setPosition(new_pos);
        }

      // Find the cells that the particles moved to
      tmp_map.clear();
      for (it=particles.begin(); it!=particles.end();)
        {
          found_cell = find_cell(mapping, triangulation, it->second, it->first);
          if (found_cell != it->first)
            {
              tmp_map.insert(std::make_pair(found_cell, it->second));
              particles.erase(it++);
            }
          else ++it;
        }
      particles.insert(tmp_map.begin(),tmp_map.end());

      // If particles fell out of the mesh, put them back in at the closest point in the mesh
      move_particles_back_in_mesh(mapping, triangulation);

      // Swap particles between processors if needed
      send_recv_particles(mapping, triangulation);

      // Ensure we didn't lose any particles
      check_particle_count();
    }

    // Output the particle locations to a specified data file type
    template <int dim>
    std::string ParticleSet<dim>::output_particle_data(const std::string &output_dir, const parallel::distributed::Triangulation<dim> &triangulation)
    {
      typename ParticleMap::iterator  it;
      MPI_Particle<dim>       particle_data;
      unsigned int            d, i;
      std::string             output_file_prefix, output_path_prefix, full_filename;
      DataOut<dim>            data_out;
      unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      unsigned int nproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

      output_file_prefix = "particle-" + Utilities::int_to_string (data_out_index, 5);
      output_path_prefix = output_dir + output_file_prefix;
      if (data_output_format == "ascii")
        {
          full_filename = output_path_prefix + "." + Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4) + ".txt";
          std::ofstream output (full_filename.c_str());
          if (!output) std::cout << "ERROR: proc " << triangulation.locally_owned_subdomain() << " could not create " << full_filename << std::endl;
          output << "# ";
          for (d=0; d<dim; ++d) output << "pos[" << d << "] ";
          for (d=0; d<dim; ++d) output << "vel[" << d << "] ";
          output << "id\n";
          for (it=particles.begin(); it!=particles.end(); ++it)
            {
              particle_data = it->second.mpi_data();
              for (d=0; d<dim; ++d) output << particle_data.pos[d] << " ";
              for (d=0; d<dim; ++d) output << particle_data.velocity[d] << " ";
              output << particle_data.id;
              output << "\n";
            }

          output.close();

          data_out_index++;
        }
      else if (data_output_format == "vtu")
        {
          const std::string filename = (output_file_prefix +
                                        "." +
                                        Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4) +
                                        ".vtu");
          const std::string full_filename = (output_dir + filename);

          std::ofstream output (full_filename.c_str());
          if (!output) std::cout << "ERROR: proc " << myid << " could not create " << filename << std::endl;

          // Write VTU file XML
          output << "<?xml version=\"1.0\"?>\n";
          output << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
          output << "  <UnstructuredGrid>\n";
          output << "    <Piece NumberOfPoints=\"" << particles.size() << "\" NumberOfCells=\"" << particles.size() << "\">\n";

          // Go through the particles on this domain and print the position of each one
          // TODO: can we change this to dim-dimensions rather than 3?
          output << "      <Points>\n";
          output << "        <DataArray name=\"Position\" type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
          for (it=particles.begin(); it!=particles.end(); ++it)
            {
              output << "          ";
              for (d=0; d<3; ++d)
                {
                  if (d < dim) output << it->second.getPosition()[d] << " ";
                  else output << "0.0 ";
                }
              output << "\n";
            }
          output << "        </DataArray>\n";
          output << "      </Points>\n";

          // Write cell related data (empty)
          output << "      <Cells>\n";
          output << "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
          for (i=1; i<=particles.size(); ++i) output << "          " << i-1 << "\n";
          output << "        </DataArray>\n";
          output << "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
          for (i=1; i<=particles.size(); ++i) output << "          " << i << "\n";
          output << "        </DataArray>\n";
          output << "        <DataArray type=\"UInt8\" Name=\"types\" Format=\"ascii\">\n";
          for (i=1; i<=particles.size(); ++i) output << "          1\n";
          output << "        </DataArray>\n";
          output << "      </Cells>\n";

          // Write point related data (id, velocity, etc)
          output << "      <PointData Scalars=\"scalars\">\n";

          output << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
          for (it=particles.begin(); it!=particles.end(); ++it)
            {
              output << "          ";
              for (d=0; d<3; ++d)
                {
                  if (d < dim) output << it->second.getVelocity()[d] << " ";
                  else output << "0.0 ";
                }
              output << "\n";
            }
          output << "        </DataArray>\n";
          output << "        <DataArray type=\"Float64\" Name=\"id\" Format=\"ascii\">\n";
          for (it=particles.begin(); it!=particles.end(); ++it)
            {
              output << "          " << it->second.getID() << "\n";
            }
          output << "        </DataArray>\n";
          output << "      </PointData>\n";

          output << "    </Piece>\n";
          output << "  </UnstructuredGrid>\n";
          output << "</VTKFile>\n";

          output.close();

          // Write the parallel pvtu and pvd files on the root process
          if (myid == 0)
            {
              const std::string pvtu_filename = (output_path_prefix + ".pvtu");

              std::ofstream pvtu_output (pvtu_filename.c_str());
              if (!pvtu_output) std::cout << "ERROR: Could not create " << filename << std::endl;

              pvtu_output << "<?xml version=\"1.0\"?>\n";
              pvtu_output << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
              pvtu_output << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
              pvtu_output << "    <PPoints>\n";
              pvtu_output << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
              pvtu_output << "    </PPoints>\n";
              pvtu_output << "    <PPointData Scalars=\"scalars\">\n";
              pvtu_output << "      <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
              pvtu_output << "      <PDataArray type=\"Float64\" Name=\"id\" format=\"ascii\"/>\n";
              pvtu_output << "    </PPointData>\n";
              for (i=0; i<nproc; ++i)
                {
                  pvtu_output << "    <Piece Source=\"" << output_file_prefix << "." << Utilities::int_to_string(i, 4) << ".vtu\"/>\n";
                }
              pvtu_output << "  </PUnstructuredGrid>\n";
              pvtu_output << "</VTKFile>\n";
              pvtu_output.close();

              times_and_pvtu_names.push_back(std::pair<double,std::string>(this->get_time(), output_file_prefix+".pvtu"));
              const std::string pvd_master_filename = (this->get_output_directory() + "particle.pvd");
              std::ofstream pvd_master (pvd_master_filename.c_str());
              data_out.write_pvd_record (pvd_master, times_and_pvtu_names);
            }
        }
#ifdef DEAL_II_HAVE_HDF5
      // TODO: error checking for H5 calls
      else if (data_output_format == "hdf5")
        {
          hid_t   h5_file_id;
          hid_t   plist_id, xfer_plist_id, position_dataset_id, velocity_dataset_id, pid_dataset_id;
          hid_t   file_dataspace_id, pos_file_dataspace_id, vel_file_dataspace_id, pid_file_dataspace_id;
          hid_t   dim_dataspace_id, one_dim_ds_id;
          hid_t   dim_mem_ds_id, one_dim_mem_ds_id;
          hid_t   pattr_id;
          herr_t  status;
          int     i;
          hsize_t global_particle_count, dims[2], offset[2], count[2];
          unsigned int  mpi_offset, mpi_count;

          double  *pos_data, *vel_data, *id_data;

          std::string h5_filename = output_path_prefix+".h5";

          // Create parallel file access
          plist_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef H5_HAVE_PARALLEL
          H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif

          // Create the file
          h5_file_id = H5Fcreate(h5_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

          // Create the file dataspace descriptions
          global_particle_count = get_global_particle_count();
          dims[0] = global_particle_count;
          dims[1] = 3;
          one_dim_ds_id = H5Screate_simple(1, dims, NULL);
          dim_dataspace_id = H5Screate_simple(2, dims, NULL);

          // Create the datasets
          position_dataset_id = H5Dcreate(h5_file_id, "nodes", H5T_NATIVE_DOUBLE, dim_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          velocity_dataset_id = H5Dcreate(h5_file_id, "velocity", H5T_NATIVE_DOUBLE, dim_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
          pid_dataset_id = H5Dcreate(h5_file_id, "id", H5T_NATIVE_DOUBLE, one_dim_ds_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

          // Close the file dataspaces
          H5Sclose(dim_dataspace_id);
          H5Sclose(one_dim_ds_id);

          // Just in case we forget what they are
          count[0] = 1;
          file_dataspace_id = H5Screate_simple (1, count, NULL);
          pattr_id = H5Acreate(h5_file_id, "Ermahgerd! Pertecrs!", H5T_NATIVE_INT, file_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
          H5Aclose(pattr_id);
          H5Sclose(file_dataspace_id);

          // Count the number of particles on this process and get the offset among all processes
          mpi_count = particles.size();
          MPI_Scan(&mpi_count, &mpi_offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
          count[0] = mpi_count;
          count[1] = 3;
          offset[0] = mpi_offset-mpi_count;
          offset[1] = 0;

          // Select the appropriate dataspace for this process
          one_dim_mem_ds_id = H5Screate_simple(1, count, NULL);
          dim_mem_ds_id = H5Screate_simple(2, count, NULL);
          pos_file_dataspace_id = H5Dget_space(position_dataset_id);
          vel_file_dataspace_id = H5Dget_space(velocity_dataset_id);
          pid_file_dataspace_id = H5Dget_space(pid_dataset_id);

          // And select the hyperslabs from each dataspace
          status = H5Sselect_hyperslab(pos_file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
          status = H5Sselect_hyperslab(vel_file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
          status = H5Sselect_hyperslab(pid_file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

          // Create property list for collective dataset write
          xfer_plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef H5_HAVE_PARALLEL
          H5Pset_dxpl_mpio(xfer_plist_id, H5FD_MPIO_COLLECTIVE);
#endif

          // Read the local particle data
          pos_data = new double[3*particles.size()];
          vel_data = new double[3*particles.size()];
          id_data = new double[particles.size()];

          for (i=0,it=particles.begin(); it!=particles.end(); ++i,++it)
            {
              for (d=0; d<dim; ++d)
                {
                  pos_data[i*3+d] = it->second.getPosition()(d);
                  vel_data[i*3+d] = it->second.getVelocity()(d);
                }
              if (dim < 3)
                {
                  pos_data[i*3+2] = 0;
                  vel_data[i*3+2] = 0;
                }
              id_data[i] = it->second.getID();
            }

          // Write particle data to the HDF5 file
          H5Dwrite(position_dataset_id, H5T_NATIVE_DOUBLE, dim_mem_ds_id, pos_file_dataspace_id, xfer_plist_id, pos_data);
          H5Dwrite(velocity_dataset_id, H5T_NATIVE_DOUBLE, dim_mem_ds_id, vel_file_dataspace_id, xfer_plist_id, vel_data);
          H5Dwrite(pid_dataset_id, H5T_NATIVE_DOUBLE, one_dim_mem_ds_id, pid_file_dataspace_id, xfer_plist_id, id_data);

          // Clear allocated resources
          delete pos_data;
          delete vel_data;
          delete id_data;
          status = H5Pclose(xfer_plist_id);
          status = H5Sclose(one_dim_mem_ds_id);
          status = H5Sclose(dim_mem_ds_id);
          status = H5Sclose(pos_file_dataspace_id);
          status = H5Sclose(vel_file_dataspace_id);
          status = H5Sclose(pid_file_dataspace_id);
          status = H5Dclose(position_dataset_id);
          status = H5Dclose(velocity_dataset_id);
          status = H5Dclose(pid_dataset_id);
          status = H5Pclose(plist_id);
          status = H5Fclose(h5_file_id);

          // Record and output XDMF info on root process
          if (myid == 0)
            {
              std::string local_h5_filename = output_file_prefix+".h5";
              XDMFEntry   entry(local_h5_filename, this->get_time(), global_particle_count, 0, 3);
              DataOut<dim> data_out;
              const std::string xdmf_filename = (this->get_output_directory() + "particle.xdmf");

              entry.add_attribute("velocity", 3);
              entry.add_attribute("id", 1);
              xdmf_entries.push_back(entry);

              data_out.write_xdmf_file(xdmf_entries, xdmf_filename.c_str(), MPI_COMM_WORLD);
            }
          data_out_index++;
        }
#endif

      return output_path_prefix;
    }

    template <int dim>
    void
    ParticleSet<dim>::set_next_data_output_time (const double current_time)
    {
      // if output_interval is positive, then set the next output interval to
      // a positive multiple; we need to interpret output_interval either
      // as years or as seconds
      if (data_output_interval > 0)
        {
          if (this->convert_output_to_years() == true)
            next_data_output_time = std::ceil(current_time / (data_output_interval * year_in_seconds)) *
                                    (data_output_interval * year_in_seconds);
          else
            next_data_output_time = std::ceil(current_time / (data_output_interval)) *
                                    (data_output_interval);
        }
    }

    template <int dim>
    void
    ParticleSet<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Tracers");
        {
          prm.declare_entry ("Number of tracers", "1e3",
                             Patterns::Double (0),
                             "Total number of tracers to create (not per processor or per element).");
          prm.declare_entry ("Time between data output", "1e8",
                             Patterns::Double (0),
                             "The time interval between each generation of "
                             "output files. A value of zero indicates that "
                             "output should be generated every time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry("Data output format", "none",
                            Patterns::Selection("none|"
                                                "ascii|"
                                                "vtu|"
                                                "hdf5"
                                               ),
                            "File format to output raw particle data in.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    ParticleSet<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Tracers");
        {
          num_initial_tracers = prm.get_double ("Number of tracers");
          data_output_interval = prm.get_double ("Time between data output");
          data_output_format = prm.get("Data output format");
#ifndef DEAL_II_HAVE_HDF5
          AssertThrow (data_output_format == "hdf5",
                       ExcMessage ("deal.ii was not compiled with HDF5 support, "
                                   "so HDF5 output is not possible. Please "
                                   "recompile deal.ii with HDF5 support turned on."));
#endif
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ParticleSet,
                                  "tracers",
                                  "Postprocessor that propagates tracer particles based on the "
                                  "velocity field.")
  }
}
