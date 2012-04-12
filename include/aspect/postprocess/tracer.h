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

#ifndef __aspect__postprocess_tracer_h
#define __aspect__postprocess_tracer_h

#include <aspect/postprocess/interface.h>

namespace aspect
{
  namespace Postprocess
  {
    // Structure used to transfer data over MPI
    template <int dim>
    struct MPI_Particle
    {
      double          pos[dim];
      double          pos0[dim];
      double          velocity[dim];
      unsigned int    id;
    };

    // Typedef of cell level/index pair
    typedef std::pair<int, int> LevelInd;

    template <int dim>
    class Particle
    {
      private:
        // Current particle position
        Point<dim>    _pos;
        // Previous position
        Point<dim>    _pos0;
        // Velocity at previous position
        Point<dim>    _velocity;

        // Globally unique ID of particle
        unsigned int  _id;

        // Flag indicating if this particle has moved off the local subdomain
        bool          _is_local;

      public:
        Particle(const Point<dim> &new_pos, const int &new_id) : _pos(new_pos), _pos0(), _velocity(), _id(new_id), _is_local(true) {};

        Particle(const MPI_Particle<dim> &particle_data);

        void setPosition(Point<dim> new_pos)
        {
          _pos = new_pos;
        };
        Point<dim> getPosition(void) const
        {
          return _pos;
        };

        void setOriginalPos(Point<dim> orig_pos)
        {
          _pos0 = orig_pos;
        };
        Point<dim> getOriginalPos(void) const
        {
          return _pos0;
        };

        void setVelocity(Point<dim> new_velocity)
        {
          _velocity = new_velocity;
        };
        Point<dim> getVelocity(void) const
        {
          return _velocity;
        };

        unsigned int getID(void) const
        {
          return _id;
        };

        void setLocal(bool new_local)
        {
          _is_local = new_local;
        };
        bool isLocal(void) const
        {
          return _is_local;
        };

        MPI_Particle<dim> mpi_data(void) const;
    };

    // MPI tag for particle transfers
    const int           PARTICLE_XFER_TAG = 382;

    template <int dim>
    class ParticleSet : public Interface<dim>, public SimulatorAccess<dim>
    {
      private:
        typedef std::multimap<LevelInd, Particle<dim> > ParticleMap;

        // Whether this set has been initialized yet or not
        bool                            initialized;

        // If the triangulation was changed (e.g. through refinement), in which
        // case we treat all recorded particle level/index values as invalid
        bool                            triangulation_changed;

        // Number of initial particles to create
        // Use a double rather than int since doubles can represent up to 2^52
        double                          num_initial_tracers;

        // Interval between output (in years if appropriate
        // simulation parameter is set, otherwise seconds)
        double                          output_interval;

        // Records time for next output to occur
        double                          next_output_time;

        // Set of particles currently in the local domain, organized by
        // the level/index of the cell they are in
        ParticleMap                     particles;

        // Index of output file
        unsigned int                    out_index;

        // Total number of particles in simulation
        double                          global_sum_particles;

        // MPI related variables
        // MPI registered datatype encapsulating the MPI_Particle struct
        MPI_Datatype                    particle_type;

        // Size and rank in the MPI communication world
        int                             world_size;
        int                             self_rank;

        // Buffers count, offset, request data for sending and receiving particles
        int                             *num_send, *num_recv;
        int                             total_send, total_recv;
        int                             *send_offset, *recv_offset;
        MPI_Request                     *send_reqs, *recv_reqs;

        void generate_particles_in_subdomain(const parallel::distributed::Triangulation<dim> &triangulation,
                                             unsigned int num_particles,
                                             unsigned int start_id);

        LevelInd recursive_find_cell(const parallel::distributed::Triangulation<dim> &triangulation,
                                     Particle<dim> &particle,
                                     LevelInd cur_cell);

        void mesh_changed(void)
        {
          triangulation_changed = true;
        };

      public:
        ParticleSet(void) : initialized(false), triangulation_changed(true), next_output_time(std::numeric_limits<double>::quiet_NaN()), out_index(0), world_size(0), self_rank(0), num_send(NULL), num_recv(NULL), send_offset(NULL), recv_offset(NULL), send_reqs(NULL), recv_reqs(NULL) {};
        ~ParticleSet(void);

        virtual std::pair<std::string,std::string> execute (TableHandler &statistics);

        void setup_mpi(MPI_Comm comm_world);
        void read_particles_from_file(parallel::distributed::Triangulation<dim> &triangulation,
                                      std::string filename);
        void generate_particles(const parallel::distributed::Triangulation<dim> &triangulation,
                                double total_particles);
        std::string output_particles(const std::string &output_dir,
                                     const parallel::distributed::Triangulation<dim> &triangulation);
        void advect_particles(const parallel::distributed::Triangulation<dim> &triangulation,
                              double timestep,
                              const TrilinosWrappers::MPI::BlockVector &solution,
                              const DoFHandler<dim> &dh,
                              const Mapping<dim> &mapping);
        void move_particles_back_in_mesh(const Mapping<dim> &mapping,
                                         const parallel::distributed::Triangulation<dim> &triangulation);
        LevelInd find_cell(const Mapping<dim> &mapping,
                           const parallel::distributed::Triangulation<dim> &triangulation,
                           Particle<dim> &particle,
                           LevelInd cur_cell);
        void send_recv_particles(const Mapping<dim> &mapping,
                                 const parallel::distributed::Triangulation<dim> &triangulation);
        void get_particle_velocities(const parallel::distributed::Triangulation<dim> &triangulation,
                                     const TrilinosWrappers::MPI::BlockVector &solution,
                                     const DoFHandler<dim> &dh,
                                     const Mapping<dim> &mapping,
                                     std::vector<Point<dim> > &velocities);
        void check_particle_count(void);

        void set_next_output_time (const double current_time);


        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };
  }
}

#endif
