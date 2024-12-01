/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.
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

#include <aspect/global.h>
#include <aspect/postprocess/crystal_preferred_orientation.h>
#include <aspect/particle/property/crystal_preferred_orientation.h>
#include <aspect/utilities.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <cstdio>
#include <unistd.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    CrystalPreferredOrientation<dim>::CrystalPreferredOrientation ()
      :
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      last_output_time (std::numeric_limits<double>::quiet_NaN()),
      output_file_number (numbers::invalid_unsigned_int),
      group_files(0),
      write_in_background_thread(false)
    {}



    template <int dim>
    CrystalPreferredOrientation<dim>::~CrystalPreferredOrientation ()
    {
      // make sure a thread that may still be running in the background,
      // writing data, finishes
      if (background_thread_main.joinable())
        background_thread_main.join ();

      if (background_thread_content_raw.joinable())
        background_thread_content_raw.join ();

      if (background_thread_content_draw_volume_weighting.joinable())
        background_thread_content_draw_volume_weighting.join ();
    }



    template <int dim>
    void
    CrystalPreferredOrientation<dim>::initialize ()
    {
      const unsigned int my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      this->random_number_generator.seed(random_number_seed+my_rank);
    }



    template <int dim>
    std::list<std::string>
    CrystalPreferredOrientation<dim>::required_other_postprocessors () const
    {
      return {"particles"};
    }



    template <int dim>
    // We need to pass the arguments by value, as this function can be called on a separate thread:
    void CrystalPreferredOrientation<dim>::writer (const std::string &filename,
                                                   const std::string &temporary_output_location,
                                                   const std::string &file_contents,
                                                   const bool compress_contents)
    {
      std::string tmp_filename = filename;
      if (temporary_output_location != "")
        {
          tmp_filename = temporary_output_location + "/aspect.tmp.XXXXXX";

          // Create the temporary file; get at the actual filename
          // by using a C-style string that mkstemp will then overwrite
          std::vector<char> tmp_filename_x (tmp_filename.size()+1);
          std::strcpy(tmp_filename_x.data(), tmp_filename.c_str());
          const int tmp_file_desc = mkstemp(tmp_filename_x.data());
          tmp_filename = tmp_filename_x.data();

          // If we failed to create the temp file, just write directly to the target file.
          // We also provide a warning about this fact. There are places where
          // this fails *on every node*, so we will get a lot of warning messages
          // into the output; in these cases, just writing multiple pieces to
          // std::cerr will produce an unreadable mass of text; rather, first
          // assemble the error message completely, and then output it atomically
          if (tmp_file_desc == -1)
            {
              const std::string x = ("***** WARNING: could not create temporary file <"
                                     +
                                     tmp_filename
                                     +
                                     ">, will output directly to final location. This may negatively "
                                     "affect performance. (On processor "
                                     + Utilities::int_to_string(Utilities::MPI::this_mpi_process (MPI_COMM_WORLD))
                                     + ".)\n");

              std::cerr << x << std::flush;

              tmp_filename = filename;
            }
          else
            close(tmp_file_desc);
        }

      std::ofstream out(tmp_filename.c_str(), std::ofstream::binary);

      AssertThrow (out, ExcMessage(std::string("Trying to write to file <") +
                                   filename +
                                   ">, but the file can't be opened!"));

      // now write and then move the tmp file to its final destination
      // if necessary
      if (compress_contents)
        {
          namespace bio = boost::iostreams;

          std::stringstream origin(file_contents);

          bio::filtering_streambuf<bio::input> compress_out;
          compress_out.push(bio::zlib_compressor());
          compress_out.push(origin);
          bio::copy(compress_out, out);
        }
      else
        {
          out << file_contents;
        }

      out.close ();

      if (tmp_filename != filename)
        {
          std::string command = std::string("mv ") + tmp_filename + " " + filename;
          int error = system(command.c_str());

          AssertThrow(error == 0,
                      ExcMessage("Could not move " + tmp_filename + " to "
                                 + filename + ". On processor "
                                 + Utilities::int_to_string(Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)) + "."));
        }
    }



    template <int dim>
    std::pair<std::string,std::string>
    CrystalPreferredOrientation<dim>::execute (TableHandler &statistics)
    {
      unsigned int particle_manager_index = numbers::invalid_unsigned_int;
      for (unsigned int i = 0; i < this->n_particle_managers(); ++i)
        {
          if (this->get_particle_manager(i).get_property_manager().template has_matching_active_plugin<Particle::Property::CrystalPreferredOrientation<dim>>())
            {
              Assert(particle_manager_index == numbers::invalid_unsigned_int,
                     ExcMessage("Multiple particle worlds with CPO data found. This is not supported."));

              particle_manager_index = i;
            }
        }

      AssertThrow(particle_manager_index != numbers::invalid_unsigned_int,
                  ExcMessage("No CPO particle properties found. Make sure that the CPO particle property plugin is selected."));

      // get particle data and a reference to the CPO particle property
      const Particle::Property::Manager<dim> &manager = this->get_particle_manager(particle_manager_index).get_property_manager();
      const Particle::Property::ParticleHandler<dim> &particle_handler = this->get_particle_manager(particle_manager_index).get_particle_handler();
      const Particle::Property::ParticlePropertyInformation &property_information = manager.get_data_info();
      const Particle::Property::CrystalPreferredOrientation<dim> &cpo_particle_property =
        manager.template get_matching_active_plugin<Particle::Property::CrystalPreferredOrientation<dim>>();

      const bool cpo_elastic_decomposition_plugin_exists = manager.plugin_name_exists("elastic tensor decomposition");

      const unsigned int n_grains = cpo_particle_property.get_number_of_grains();
      const unsigned int n_minerals = cpo_particle_property.get_number_of_minerals();

      // if this is the first time we get here, set the last output time
      // to the current time - output_interval. this makes sure we
      // always produce data during the first time step
      if (std::isnan(last_output_time))
        last_output_time = this->get_time() - output_interval;

      // If it's not time to generate an output file or we do not write output
      // return early.
      if (this->get_time() < last_output_time + output_interval && this->get_time() != end_time)
        return std::make_pair("","");

      if (output_file_number == numbers::invalid_unsigned_int)
        output_file_number = 0;
      else
        ++output_file_number;

      // Now prepare everything for writing the output
      std::string particle_file_prefix_main = this->get_output_directory() +  "particles_cpo/particles-" + Utilities::int_to_string (output_file_number, 5);
      std::string particle_file_prefix_content_raw = this->get_output_directory() +  "particles_cpo/CPO-" + Utilities::int_to_string (output_file_number, 5);
      std::string particle_file_prefix_content_draw_volume_weighting = this->get_output_directory() +  "particles_cpo/weighted_CPO-" + Utilities::int_to_string (output_file_number, 5);

      std::stringstream string_stream_main;
      std::stringstream string_stream_content_raw;
      std::stringstream string_stream_content_draw_volume_weighting;

      string_stream_main << "id x y" << (dim == 3 ? " z" : "") << " olivine_deformation_type"
                         << (cpo_elastic_decomposition_plugin_exists ? (std::string(" full_norm_square ")
                                                                        + "triclinic_norm_square_p1 triclinic_norm_square_p2 triclinic_norm_square_p3 "
                                                                        + "monoclinic_norm_square_p1 monoclinic_norm_square_p2 monoclinic_norm_square_p3 "
                                                                        + "orthohombic_norm_square_p1 orthohombic_norm_square_p2 orthohombic_norm_square_p3 "
                                                                        + "tetragonal_norm_square_p1 tetragonal_norm_square_p2 tetragonal_norm_square_p3 "
                                                                        + "hexagonal_norm_square_p1 hexagonal_norm_square_p2 hexagonal_norm_square_p3 "
                                                                        + "isotropic_norm_square") : "") << std::endl;

      const unsigned int cpo_data_position = property_information.n_fields() == 0
                                             ?
                                             0
                                             :
                                             property_information.get_position_by_field_name("cpo mineral 0 type");

      std::vector<std::vector<Tensor<2,3>>> rotation_matrices (n_minerals, {n_grains,Tensor<2,3>()});
      std::vector<std::vector<std::array<double,3>>> euler_angles(n_minerals, {n_grains,{{0}}});

      // write unweighted header
      if (write_raw_cpo.size() != 0)
        {
          string_stream_content_raw << "id" << " " << std::setprecision(12);
          for (unsigned int property_i = 0; property_i < write_raw_cpo.size(); ++property_i)
            {
              switch (write_raw_cpo[property_i].second)
                {
                  case Output::VolumeFraction:
                    string_stream_content_raw << "mineral_" << write_raw_cpo[property_i].first << "_volume_fraction" << " ";
                    break;

                  case Output::RotationMatrix:
                    string_stream_content_raw << "mineral_" << write_raw_cpo[property_i].first << "_RM_0" << " "
                                              << "mineral_" << write_raw_cpo[property_i].first << "_RM_1" << " "
                                              << "mineral_" << write_raw_cpo[property_i].first << "_RM_2" << " "
                                              << "mineral_" << write_raw_cpo[property_i].first << "_RM_3" << " "
                                              << "mineral_" << write_raw_cpo[property_i].first << "_RM_4" << " "
                                              << "mineral_" << write_raw_cpo[property_i].first << "_RM_5" << " "
                                              << "mineral_" << write_raw_cpo[property_i].first << "_RM_6" << " "
                                              << "mineral_" << write_raw_cpo[property_i].first << "_RM_7" << " "
                                              << "mineral_" << write_raw_cpo[property_i].first << "_RM_8" << " ";
                    break;

                  case Output::EulerAngles:
                    string_stream_content_raw << "mineral_" << write_raw_cpo[property_i].first << "_EA_phi" << " "
                                              << "mineral_" << write_raw_cpo[property_i].first << "_EA_theta" << " "
                                              << "mineral_" << write_raw_cpo[property_i].first << "_EA_z" << " ";
                    break;
                  default:
                    Assert(false, ExcMessage("Internal error: raw CPO postprocess case not found."));
                    break;
                }
            }
          string_stream_content_raw << std::endl;
        }


      // write weighted header
      if (write_draw_volume_weighted_cpo.size() != 0)
        {
          string_stream_content_draw_volume_weighting << "id" << " ";
          for (unsigned int property_i = 0; property_i < write_draw_volume_weighted_cpo.size(); ++property_i)
            {
              switch (write_draw_volume_weighted_cpo[property_i].second)
                {
                  case Output::VolumeFraction:
                    string_stream_content_draw_volume_weighting << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_volume_fraction" << " ";
                    break;

                  case Output::RotationMatrix:
                    string_stream_content_draw_volume_weighting << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_RM_0" << " "
                                                                << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_RM_1" << " "
                                                                << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_RM_2" << " "
                                                                << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_RM_3" << " "
                                                                << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_RM_4" << " "
                                                                << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_RM_5" << " "
                                                                << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_RM_6" << " "
                                                                << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_RM_7" << " "
                                                                << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_RM_8" << " ";
                    break;

                  case Output::EulerAngles:
                    string_stream_content_draw_volume_weighting << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_EA_phi" << " "
                                                                << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_EA_theta" << " "
                                                                << "mineral_" << write_draw_volume_weighted_cpo[property_i].first << "_EA_z" << " ";
                    break;
                  default:
                    Assert(false, ExcMessage("Internal error: raw CPO postprocess case not found."));
                    break;
                }
            }
          string_stream_content_draw_volume_weighting << std::endl;

        }

      for (const auto &particle: particle_handler)
        {
          AssertThrow(particle.has_properties(),
                      ExcMessage("No particle properties found. Make sure that the CPO particle property plugin is selected."));

          const unsigned int id = particle.get_id();
          const ArrayView<const double> properties = particle.get_properties();

          const Point<dim> position = particle.get_location();

          // always make a vector of rotation matrices
          for (unsigned int mineral = 0; mineral < n_minerals; ++mineral)
            {
              rotation_matrices[mineral].resize(n_grains);
              for (unsigned int i_grain = 0; i_grain < n_grains; ++i_grain)
                {
                  rotation_matrices[mineral][i_grain] = cpo_particle_property.get_rotation_matrix_grains(
                                                          cpo_data_position,
                                                          properties,
                                                          mineral,
                                                          i_grain);
                }
            }

          const unsigned int lpo_hex_data_position = property_information.n_fields() == 0 || cpo_elastic_decomposition_plugin_exists == false
                                                     ?
                                                     0
                                                     :
                                                     property_information.get_position_by_field_name("cpo elastic axis e1");

          // write main file
          string_stream_main << id << " " << position << " " << properties[cpo_data_position];

          if (cpo_elastic_decomposition_plugin_exists == true)
            {
              string_stream_main << " " << properties[lpo_hex_data_position+12] << " " << properties[lpo_hex_data_position+13]
                                 << " " << properties[lpo_hex_data_position+14] << " " << properties[lpo_hex_data_position+15]
                                 << " " << properties[lpo_hex_data_position+16] << " " << properties[lpo_hex_data_position+17]
                                 << " " << properties[lpo_hex_data_position+18] << " " << properties[lpo_hex_data_position+19]
                                 << " " << properties[lpo_hex_data_position+20] << " " << properties[lpo_hex_data_position+21]
                                 << " " << properties[lpo_hex_data_position+22] << " " << properties[lpo_hex_data_position+23]
                                 << " " << properties[lpo_hex_data_position+24] << " " << properties[lpo_hex_data_position+25]
                                 << " " << properties[lpo_hex_data_position+26] << " " << properties[lpo_hex_data_position+27]
                                 << " " << properties[lpo_hex_data_position+28];
            }
          string_stream_main << std::endl;

          // write content file
          if (compute_raw_euler_angles == true)
            {
              euler_angles.resize(n_minerals);
              for (unsigned int mineral = 0; mineral < n_minerals; ++mineral)
                {
                  euler_angles[mineral].resize(n_grains);
                  for (unsigned int i_grain = 0; i_grain < n_grains; ++i_grain)
                    {
                      euler_angles[mineral][i_grain] = Utilities::zxz_euler_angles_from_rotation_matrix(
                                                         rotation_matrices[mineral][i_grain]);
                    }
                }
            }

          if (write_raw_cpo.size() != 0)
            {
              // write unweighted data
              for (unsigned int grain = 0; grain < n_grains; ++grain)
                {
                  string_stream_content_raw << id << " ";
                  for (unsigned int property_i = 0; property_i < write_raw_cpo.size(); ++property_i)
                    {
                      switch (write_raw_cpo[property_i].second)
                        {
                          case Output::VolumeFraction:
                            string_stream_content_raw << cpo_particle_property.get_volume_fractions_grains(
                                                        cpo_data_position,
                                                        properties,
                                                        write_raw_cpo[property_i].first,
                                                        grain) << " ";
                            break;

                          case Output::RotationMatrix:
                            string_stream_content_raw << rotation_matrices[write_raw_cpo[property_i].first][grain]<< " ";
                            break;

                          case Output::EulerAngles:
                            Assert(compute_raw_euler_angles == true,
                                   ExcMessage("Internal error: writing out raw Euler angles, without them being computed."));
                            string_stream_content_raw << euler_angles[write_raw_cpo[property_i].first][grain][0] << " "
                                                      <<  euler_angles[write_raw_cpo[property_i].first][grain][1] << " "
                                                      <<  euler_angles[write_raw_cpo[property_i].first][grain][2] << " ";
                            break;
                          default:
                            Assert(false, ExcMessage("Internal error: raw CPO postprocess case not found."));
                            break;
                        }
                    }
                  string_stream_content_raw << std::endl;
                }

            }

          if (write_draw_volume_weighted_cpo.size() != 0)
            {
              std::vector<std::vector<dealii::Tensor<2,3>>> weighted_rotation_matrices;
              std::vector<std::vector<std::array<double,3>>> weighted_euler_angles;

              weighted_euler_angles.resize(n_minerals);
              weighted_rotation_matrices.resize(n_minerals);
              std::vector<std::vector<double>> volume_fractions_grains(n_minerals,std::vector<double>(n_grains,-1.));
              for (unsigned int mineral = 0; mineral < n_minerals; ++mineral)
                {
                  for (unsigned int i_grain = 0; i_grain < n_grains; ++i_grain)
                    {
                      weighted_euler_angles[mineral].resize(n_grains);
                      weighted_rotation_matrices[mineral].resize(n_grains);
                      volume_fractions_grains[mineral][i_grain] = cpo_particle_property.get_volume_fractions_grains(
                                                                    cpo_data_position,
                                                                    properties,
                                                                    mineral,
                                                                    i_grain);
                    }
                  weighted_rotation_matrices[mineral] = Utilities::rotation_matrices_random_draw_volume_weighting(volume_fractions_grains[mineral], rotation_matrices[mineral], n_grains, this->random_number_generator);

                  Assert(weighted_rotation_matrices[mineral].size() == euler_angles[mineral].size(),
                         ExcMessage("Weighted rotation matrices vector (size = " + std::to_string(weighted_rotation_matrices[mineral].size()) +
                                    ") has different size from input angles (size = " + std::to_string(euler_angles[mineral].size()) + ")."));

                  for (unsigned int i_grain = 0; i_grain < n_grains; ++i_grain)
                    {
                      weighted_euler_angles[mineral][i_grain] = Utilities::zxz_euler_angles_from_rotation_matrix(
                                                                  weighted_rotation_matrices[mineral][i_grain]);
                    }
                }
              string_stream_content_draw_volume_weighting << std::endl;

              // write data
              for (unsigned int grain = 0; grain < n_grains; ++grain)
                {
                  string_stream_content_draw_volume_weighting << id << " ";
                  for (unsigned int property_i = 0; property_i < write_draw_volume_weighted_cpo.size(); ++property_i)
                    {
                      switch (write_draw_volume_weighted_cpo[property_i].second)
                        {
                          case Output::VolumeFraction:
                            string_stream_content_draw_volume_weighting << volume_fractions_grains[write_draw_volume_weighted_cpo[property_i].first][grain] << " ";
                            break;

                          case Output::RotationMatrix:
                            string_stream_content_draw_volume_weighting << weighted_rotation_matrices[write_draw_volume_weighted_cpo[property_i].first][grain] << " ";
                            break;

                          case Output::EulerAngles:
                            Assert(compute_raw_euler_angles == true,
                                   ExcMessage("Internal error: writing out raw Euler angles, without them being computed."));
                            string_stream_content_draw_volume_weighting << weighted_euler_angles[write_draw_volume_weighted_cpo[property_i].first][grain][0] << " "
                                                                        <<  weighted_euler_angles[write_draw_volume_weighted_cpo[property_i].first][grain][1] << " "
                                                                        <<  weighted_euler_angles[write_draw_volume_weighted_cpo[property_i].first][grain][2] << " ";
                            break;

                          default:
                            Assert(false, ExcMessage("Internal error: raw CPO postprocess case not found."));
                            break;
                        }
                    }
                  string_stream_content_draw_volume_weighting << std::endl;
                }
            }
        }

      std::string filename_main = particle_file_prefix_main + "." + Utilities::int_to_string(dealii::Utilities::MPI::this_mpi_process (this->get_mpi_communicator()),4) + ".dat";
      std::string filename_raw = particle_file_prefix_content_raw + "." + Utilities::int_to_string(dealii::Utilities::MPI::this_mpi_process (this->get_mpi_communicator()),4) + ".dat";
      std::string filename_draw_volume_weighting = particle_file_prefix_content_draw_volume_weighting + "." + Utilities::int_to_string(dealii::Utilities::MPI::this_mpi_process (this->get_mpi_communicator()),4) + ".dat";

      std::unique_ptr<std::string> file_contents_main = std::make_unique<std::string>(string_stream_main.str());
      std::unique_ptr<std::string> file_contents_raw = std::make_unique<std::string>(string_stream_content_raw.str());
      std::unique_ptr<std::string> file_contents_draw_volume_weighting = std::make_unique<std::string>(string_stream_content_draw_volume_weighting.str());

      if (write_in_background_thread)
        {
          // Wait for all previous write operations to finish, should
          // any be still active,
          if (background_thread_main.joinable())
            background_thread_main.join ();

          // then continue with writing the main file
          background_thread_main
            = std::thread([ my_filename = std::move(filename_main),
                            my_temporary_output_location = temporary_output_location,
                            my_file_contents = std::move(file_contents_main)]()
          {
            writer (my_filename, my_temporary_output_location, *my_file_contents, false);
          });

          if (write_raw_cpo.size() != 0)
            {
              // Wait for all previous write operations to finish, should
              // any be still active,
              if (background_thread_content_raw.joinable())
                background_thread_content_raw.join ();

              // then continue with writing our own data.
              background_thread_content_raw
                = std::thread([ my_filename = std::move(filename_raw),
                                my_temporary_output_location = temporary_output_location,
                                my_file_contents = std::move(file_contents_raw),
                                my_compress_cpo_data_files = compress_cpo_data_files]()
              {
                writer (my_filename, my_temporary_output_location, *my_file_contents, my_compress_cpo_data_files);
              });
            }

          if (write_draw_volume_weighted_cpo.size() != 0)
            {
              // Wait for all previous write operations to finish, should
              // any be still active,
              if (background_thread_content_draw_volume_weighting.joinable())
                background_thread_content_draw_volume_weighting.join ();

              // then continue with writing our own data.
              background_thread_content_draw_volume_weighting
                = std::thread([ my_filename = std::move(filename_draw_volume_weighting),
                                my_temporary_output_location = temporary_output_location,
                                my_file_contents = std::move(file_contents_draw_volume_weighting),
                                my_compress_cpo_data_files = compress_cpo_data_files]()
              {
                writer (my_filename, my_temporary_output_location, *my_file_contents, my_compress_cpo_data_files);
              });
            }
        }
      else
        {
          writer(filename_main,temporary_output_location,*file_contents_main, false);
          if (write_raw_cpo.size() != 0)
            writer(filename_raw,temporary_output_location,*file_contents_raw, compress_cpo_data_files);
          if (write_draw_volume_weighted_cpo.size() != 0)
            writer(filename_draw_volume_weighting,temporary_output_location,*file_contents_draw_volume_weighting, compress_cpo_data_files);
        }


      // up the next time we need output
      set_last_output_time (this->get_time());

      const std::string &particles_cpo_output = particle_file_prefix_content_raw;

      // record the file base file name in the output file
      statistics.add_value ("Particle CPO file name", particles_cpo_output);
      return std::make_pair("Writing particle cpo output:", particles_cpo_output);
    }



    template <int dim>
    void
    CrystalPreferredOrientation<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then update the last supposed output
      // time
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
        }
    }



    template <int dim>
    typename CrystalPreferredOrientation<dim>::Output
    CrystalPreferredOrientation<dim>::string_to_output_enum(const std::string &string)
    {
      // olivine volume fraction, olivine rotation matrix, olivine Euler angles, enstatite volume fraction, enstatite rotation matrix, enstatite Euler angles
      if (string == "volume fraction")
        return Output::VolumeFraction;
      if (string == "rotation matrix")
        return Output::RotationMatrix;
      if (string == "Euler angles")
        return Output::EulerAngles;
      return Output::not_found;
    }



    template <int dim>
    void
    CrystalPreferredOrientation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Crystal Preferred Orientation");
        {
          prm.declare_entry ("Time between data output", "1e8",
                             Patterns::Double (0),
                             "The time interval between each generation of "
                             "output files. A value of zero indicates that "
                             "output should be generated every time step.\n\n"
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");

          prm.declare_entry ("Random number seed", "1",
                             Patterns::Integer (0),
                             "The seed used to generate random numbers. This will make sure that "
                             "results are reproducible as long as the problem is run with the "
                             "same amount of MPI processes. It is implemented as final seed = "
                             "random number seed + MPI Rank. ");

          prm.declare_entry ("Write in background thread", "false",
                             Patterns::Bool(),
                             "File operations can potentially take a long time, blocking the "
                             "progress of the rest of the model run. Setting this variable to "
                             "`true' moves this process into background threads, while the "
                             "rest of the model continues.");

          prm.declare_entry ("Temporary output location", "",
                             Patterns::Anything(),
                             "On large clusters it can be advantageous to first write the "
                             "output to a temporary file on a local file system and later "
                             "move this file to a network file system. If this variable is "
                             "set to a non-empty string it will be interpreted as a "
                             "temporary storage location.");

          prm.declare_entry ("Write out raw cpo data",
                             "olivine volume fraction,olivine Euler angles,enstatite volume fraction,enstatite Euler angles",
                             Patterns::List(Patterns::Anything()),
                             "A list containing what particle cpo data needs "
                             "to be written out after the particle id. This writes out the raw "
                             "cpo data files for each MPI process. It can write out the following data: "
                             "olivine volume fraction, olivine rotation matrix, olivine Euler angles, "
                             "enstatite volume fraction, enstatite rotation matrix, enstatite Euler angles. \n"
                             "Note that the rotation matrix and Euler angles both contain the same "
                             "information, but in a different format. Euler angles are recommended "
                             "over the rotation matrix since they only require to write 3 values instead "
                             "of 9. If the list is empty, this file will not be written."
                             "Furthermore, the entries will be written out in the order given, "
                             "and if entries are entered multiple times, they will be written "
                             "out multiple times.");

          prm.declare_entry ("Write out draw volume weighted cpo data",
                             "olivine Euler angles,enstatite Euler angles",
                             Patterns::List(Patterns::Anything()),
                             "A list containing the what part of the random draw volume "
                             "weighted particle cpo data needs to be written out after "
                             "the particle id. after using a random draw volume weighting. "
                             "The random draw volume weigthing uses a uniform random distribution "
                             "This writes out the raw cpo data files for "
                             "each MPI process. It can write out the following data: "
                             "olivine volume fraction, olivine rotation matrix, olivine Euler angles, "
                             "enstatite volume fraction, enstatite rotation matrix, enstatite Euler angles. \n"
                             "Note that the rotation matrix and Euler angles both contain the same "
                             "information, but in a different format. Euler angles are recommended "
                             "over the rotation matrix since they only require to write 3 values instead "
                             "of 9. If the list is empty, this file will not be written. "
                             "Furthermore, the entries will be written out in the order given, "
                             "and if entries are entered multiple times, they will be written "
                             "out multiple times.");
          prm.declare_entry ("Compress cpo data files", "true",
                             Patterns::Bool(),
                             "Whether to compress the raw and weighted cpo data output files with zlib.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

    }



    template <int dim>
    void
    CrystalPreferredOrientation<dim>::parse_parameters (ParameterHandler &prm)
    {
      end_time = prm.get_double ("End time");
      if (this->convert_output_to_years())
        end_time *= year_in_seconds;

      unsigned int n_minerals;

      prm.enter_subsection("Particles");
      {
        prm.enter_subsection("Crystal Preferred Orientation");
        {
          prm.enter_subsection("Initial grains");
          {
            // Static variable of CPO has not been initialize yet, so we need to get it directly.
            n_minerals = dealii::Utilities::split_string_list(prm.get("Minerals")).size();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Crystal Preferred Orientation");
        {
          output_interval = prm.get_double ("Time between data output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;

          random_number_seed = prm.get_integer ("Random number seed");

          AssertThrow(this->get_parameters().run_postprocessors_on_nonlinear_iterations == false,
                      ExcMessage("Postprocessing nonlinear iterations in models with "
                                 "particles is currently not supported."));

          aspect::Utilities::create_directory (this->get_output_directory() + "particles_cpo/",
                                               this->get_mpi_communicator(),
                                               true);

          write_in_background_thread = prm.get_bool("Write in background thread");
          temporary_output_location = prm.get("Temporary output location");

          if (temporary_output_location != "")
            {
              // Check if a command-processor is available by calling system() with a
              // null pointer. System is guaranteed to return non-zero if it finds
              // a terminal and zero if there is none (like on the compute nodes of
              // some cluster architectures, e.g. IBM BlueGene/Q)
              AssertThrow(system((char *)nullptr) != 0,
                          ExcMessage("Usage of a temporary storage location is only supported if "
                                     "there is a terminal available to move the files to their final location "
                                     "after writing. The system() command did not succeed in finding such a terminal."));
            }

          std::vector<std::string> write_raw_cpo_list = Utilities::split_string_list(prm.get("Write out raw cpo data"));
          write_raw_cpo.resize(write_raw_cpo_list.size());
          bool found_euler_angles = false;
          for (unsigned int i = 0; i < write_raw_cpo_list.size(); ++i)
            {
              std::vector<std::string> split_raw_cpo_instructions = Utilities::split_string_list(write_raw_cpo_list[i],':');

              AssertThrow(split_raw_cpo_instructions.size() == 2,
                          ExcMessage("Value \""+ write_raw_cpo_list[i] +"\", set in \"Write out raw cpo data\", is not a correct option "
                                     + "because it should contain a mineral identification and a output specifier separated by a colon (:). This entry "
                                     + "does not follow those rules."));

              // get mineral number
              std::vector<std::string> mineral_instructions = Utilities::split_string_list(split_raw_cpo_instructions[0],' ');
              AssertThrow(mineral_instructions.size() == 2,
                          ExcMessage("Value \""+ write_raw_cpo_list[i] +"\", set in \"Write out raw cpo data\", is not a correct option "
                                     + "because the mineral identification part should contain two elements, the word mineral and a number."));

              AssertThrow(mineral_instructions[0] == "mineral",
                          ExcMessage("Value \""+ write_raw_cpo_list[i] +"\", set in \"Write out raw cpo data\", is not a correct option "
                                     + "because the mineral identification part should start with the word mineral and it starts with \""
                                     + mineral_instructions[0] + "\"."));

              int mineral_number = Utilities::string_to_int(mineral_instructions[1]);
              Assert(mineral_number >= 0, ExcMessage("Internal error: mineral_number is negative: " + std::to_string(mineral_number) + "."));

              AssertThrow((unsigned int) mineral_number < n_minerals,
                          ExcMessage("Value \""+ write_raw_cpo_list[i] +"\", set in \"Write out raw cpo data\", is not a correct option "
                                     + "because the mineral number (" + std::to_string(mineral_number) + ") is larger than the number of minerals "
                                     + "provided in the CPO subsection (" + std::to_string(n_minerals) + ")."));

              // get mineral fabric/deformation type
              Output cpo_fabric_instruction = string_to_output_enum(split_raw_cpo_instructions[1]);

              AssertThrow(cpo_fabric_instruction != Output::not_found,
                          ExcMessage("Value \""+ write_raw_cpo_list[i] +"\", set in \"Write out raw cpo data\", is not a correct option."));

              if (cpo_fabric_instruction == Output::EulerAngles)
                found_euler_angles = true;

              write_raw_cpo[i] = std::make_pair(mineral_number,cpo_fabric_instruction);
            }

          std::vector<std::string> write_draw_volume_weighted_cpo_list = Utilities::split_string_list(prm.get("Write out draw volume weighted cpo data"));
          write_draw_volume_weighted_cpo.resize(write_draw_volume_weighted_cpo_list.size());
          bool found_rotation_matrix = false;
          for (unsigned int i = 0; i < write_draw_volume_weighted_cpo_list.size(); ++i)
            {
              std::vector<std::string> split_draw_volume_weighted_cpo_instructions = Utilities::split_string_list(write_draw_volume_weighted_cpo_list[i],':');

              AssertThrow(split_draw_volume_weighted_cpo_instructions.size() == 2,
                          ExcMessage("Value \""+ write_draw_volume_weighted_cpo_list[i] +"\", set in \"Write out draw volume weighted cpo data\", is not a correct option "
                                     + "because it should contain a mineral identification and a output specifier separated by a colon (:). This entry "
                                     + "does not follow those rules."));

              // get mineral number
              std::vector<std::string> mineral_instructions = Utilities::split_string_list(split_draw_volume_weighted_cpo_instructions[0],' ');
              AssertThrow(mineral_instructions.size() == 2,
                          ExcMessage("Value \""+ write_draw_volume_weighted_cpo_list[i] +"\", set in \"Write out draw volume weighted cpo data\", is not a correct option "
                                     + "because the mineral identification part should contain two elements, the word mineral and a number."));

              AssertThrow(mineral_instructions[0] == "mineral",
                          ExcMessage("Value \""+ write_draw_volume_weighted_cpo_list[i] +"\", set in \"Write out draw volume weighted cpo data\", is not a correct option "
                                     + "because the mineral identification part should start with the word mineral and it starts with \""
                                     + mineral_instructions[0] + "\"."));

              int mineral_number = Utilities::string_to_int(mineral_instructions[1]);
              Assert(mineral_number >= 0, ExcMessage("Internal error: mineral_number is negative: " + std::to_string(mineral_number) + "."));

              AssertThrow((unsigned int) mineral_number < n_minerals,
                          ExcMessage("Value \""+ write_raw_cpo_list[i] +"\", set in \"Write out raw cpo data\", is not a correct option "
                                     + "because the mineral number (" + std::to_string(mineral_number) + ") is larger than the number of minerals "
                                     + "provided in the CPO subsection (" + std::to_string(n_minerals) + ")."));

              // get mineral fabric/deformation type
              Output cpo_fabric_instruction = string_to_output_enum(split_draw_volume_weighted_cpo_instructions[1]);

              AssertThrow(cpo_fabric_instruction != Output::not_found,
                          ExcMessage("Value \""+ write_draw_volume_weighted_cpo_list[i] +"\", set in \"Write out draw volume weighted cpo data\", is not a correct option."));

              if (cpo_fabric_instruction == Output::RotationMatrix)
                found_rotation_matrix = true;

              write_draw_volume_weighted_cpo[i] = std::make_pair(mineral_number,cpo_fabric_instruction);
            }

          if (write_draw_volume_weighted_cpo_list.size() != 0 || found_euler_angles == true)
            compute_raw_euler_angles = true;
          else
            compute_raw_euler_angles = false;

          if (write_draw_volume_weighted_cpo_list.size() != 0 && found_rotation_matrix == true)
            compute_weighted_rotation_matrix = true;
          else
            compute_weighted_rotation_matrix = false;

          compress_cpo_data_files = prm.get_bool("Compress cpo data files");
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
    ASPECT_REGISTER_POSTPROCESSOR(CrystalPreferredOrientation,
                                  "crystal preferred orientation",
                                  "A Postprocessor that writes out CPO specific particle data."
                                  "It can write out the CPO data as it is stored (raw) and/or as a"
                                  "random draw volume weighted representation. The latter one"
                                  "is recommended for plotting against real data. For both representations"
                                  "the specific output fields and their order can be set."
                                  "The work of this postprocessor should better be done by the main particles "
                                  "postprocessor, however we need to be able to process the data before outputting it, "
                                  "which does not work with that postprocessor. If this is added to the other "
                                  "postprocessor in the future this one becomes obsolete.")
  }
}
