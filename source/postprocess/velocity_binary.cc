#include <aspect/postprocess/velocity_binary.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/lac/block_vector.h>
#ifdef DEAL_II_WITH_ZLIB
#include <zlib.h>
#endif

namespace aspect {
    namespace Postprocess {
        template <int dim>
        VelocityBinary<dim>::VelocityBinary()
        {}

        template <int dim>
       VelocityBinary<dim>::~VelocityBinary()
        {}

       template <int dim>
       void VelocityBinary<dim>::initialize() {
          // Interface::initialize();
       }

       template <int dim>
       std::pair<std::string, std::string> VelocityBinary<dim>::execute(TableHandler &statistics) {
          std::string fileName = "solution-" + Utilities::int_to_string(this->get_timestep_number(), 5) + ".mesh";
          //std::string fileName = "solution-" + Utilities::int_to_string(this->get_timestep_number(), 5) + "-" + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()), 4);
      
          LinearAlgebra::BlockVector solution;
          solution = this->get_solution(); 
          parallel::distributed::SolutionTransfer<dim, LinearAlgebra::BlockVector> sol_trans(this->get_dof_handler());
          sol_trans.prepare_serialization (solution);
          this->get_triangulation().save(fileName.c_str());
          
//          std::ofstream output (this->get_output_directory() + fileName, std::ios::binary);
          
//          aspect::oarchive oa (output);
 
  //         std::ofstream output (this->get_output_directory() + fileName);

/*          AssertThrow (output,
                        ExcMessage (std::string("Could not open output file <"
                                                +
                                                fileName
                                                +
                                                ">.")));
*/
//          dealii::BlockVector<double> solution(this->get_solution());
//          solution.block_write(output);
 
//            LinearAlgebra::BlockVector solution;
//
//            solution = this->get_solution();
//            LinearAlgebra::Vector s = solution.block(0);
//
            /*for(BlockVectorBase<dealii::TrilinosWrappers::MPI::Vector>::iterator itr = solution.begin(); itr != solution.end(); itr++) {

                output.write(reinterpret_cast<char*>(&(*itr)), sizeof(*itr));

            }*/

 //         output.close();

/*            test.print(output2);
            output2.close();
*/
   //     ofs.write( reinterpret_cast<char*>( &pi ), sizeof pi );
        // Close the file to unlock it
   //     ofs.close();

        // Use a new object so we don't have to worry
        // about error states in the old object
//        ifstream ifs( "atest.txt", ios::binary );
//        double read;
//
//        if ( ifs ) {
//            ifs.read( reinterpret_cast<char*>( &read ), sizeof read );
//            cout << read << '\n';
//        }

//            parallel::distributed::SolutionTransfer<dim, LinearAlgebra::BlockVector>
//                    system_trans (this->get_dof_handler());
//            system_trans.prepare_serialization (solution);
//
//            // save general information This calls the serialization functions on all
//            // processes (so that they can take additional action, if necessary, see
//            // the manual) but only writes to the restart file on process 0
//            {
//                std::ostringstream oss;
//
//                // serialize into a stringstream
//                aspect::oarchive oa (oss);
//                oa << (*this);
//
//                // compress with zlib and write to file on the root processor
//#ifdef DEAL_II_WITH_ZLIB
//                if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
//                {
//                    uLongf compressed_data_length = compressBound (oss.str().length());
//                    std::vector<char *> compressed_data (compressed_data_length);
//                    int err = compress2 ((Bytef *) &compressed_data[0],
//                                         &compressed_data_length,
//                                         (const Bytef *) oss.str().data(),
//                                         oss.str().length(),
//                                         Z_BEST_COMPRESSION);
//                    (void)err;
//                    Assert (err == Z_OK, ExcInternalError());
//
//                    // build compression header
//                    const uint32_t compression_header[4]
//                            = { 1,                                   /* number of blocks */
//                                (uint32_t)oss.str().length(), /* size of block */
//                                (uint32_t)oss.str().length(), /* size of last block */
//                                (uint32_t)compressed_data_length
//                            }; /* list of compressed sizes of blocks */
//
//                    std::ofstream f ((fileName + ".z").c_str());
//                    f.write((const char *)compression_header, 4 * sizeof(compression_header[0]));
//                    f.write((char *)&compressed_data[0], compressed_data_length);
//                }
//#else
//                AssertThrow (false,
//                   ExcMessage ("You need to have deal.II configured with the 'libz' "
//                               "option to support checkpoint/restart, but deal.II "
//                               "did not detect its presence when you called 'cmake'."));
//#endif
//
//            }

           //solution.print(output);
           // aspect::oarchive ar(output);
           // solution.;
           // solution.serialize(ar,0);
           //output.close();

            return std::make_pair("Writing binary output to: ", fileName);
       }

     /* template <int dim>
       void VelocityBinary<dim>::declare_parameters(ParameterHandler &prm) {
           prm.enter_subsection("Postprocess");
           {
               prm.enter_subsection("Raw data");
               {
                   prm.declare_entry ("Time between data output", "0",
                                      Patterns::Double (0),
                                      "The time interval between "
               }
           }
       }

        template <int dim>
       void VelocityBinary<dim>::parse_parameters(ParameterHandler &prm) {

       }*/

/*     template<int dim>
       template<class Archive>
       void VelocityBinary<dim>::serialize(Archive &ar, const unsigned int version) {
          ar << (*this);
       }*/
   }
}


namespace aspect
{
    namespace Postprocess
    {
        ASPECT_REGISTER_POSTPROCESSOR(VelocityBinary,
                                      "velocity binary",
                                      "A Postprocessor that output the velocity solution data per timestep.")
    }
}



