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


#include <aspect/simulator.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>

#if ASPECT_USE_SHARED_LIBS==1
#  include <dlfcn.h>
#  ifdef ASPECT_HAVE_LINK_H
#    include <link.h>
#  endif
#endif



// get the value of a particular parameter from the input
// file. return an empty string if not found
std::string
get_last_value_of_parameter(const std::string &parameter_filename,
                            const std::string &parameter_name)
{
  std::string return_value;

  std::ifstream x_file(parameter_filename.c_str());
  while (x_file)
    {
      // get one line and strip spaces at the front and back
      std::string line;
      std::getline(x_file, line);
      while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
        line.erase(0, 1);
      while ((line.size() > 0)
             && (line[line.size() - 1] == ' ' || line[line.size() - 1] == '\t'))
        line.erase(line.size() - 1, std::string::npos);
      // now see whether the line starts with 'set' followed by multiple spaces
      // if not, try next line
      if (line.size() < 4)
        continue;

      if ((line[0] != 's') || (line[1] != 'e') || (line[2] != 't')
          || !(line[3] == ' ' || line[3] == '\t'))
        continue;

      // delete the "set " and then delete more spaces if present
      line.erase(0, 4);
      while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
        line.erase(0, 1);
      // now see whether the next word is the word we look for
      if (line.find(parameter_name) != 0)
        continue;

      line.erase(0, parameter_name.size());
      while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
        line.erase(0, 1);

      // we'd expect an equals size here
      if ((line.size() < 1) || (line[0] != '='))
        continue;

      // trim the equals sign at the beginning and possibly following spaces
      // as well as spaces at the end
      line.erase(0, 1);
      while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
        line.erase(0, 1);
      while ((line.size() > 0) && (line[line.size()-1] == ' ' || line[line.size()-1] == '\t'))
        line.erase(line.size()-1, std::string::npos);

      // the rest should now be what we were looking for
      return_value = line;
    }

  return return_value;
}


// extract the dimension in which to run ASPECT from the
// parameter file. this is something that we need to do
// before processing the parameter file since we need to
// know whether to use the dim=2 or dim=3 instantiation
// of the main classes
unsigned int
get_dimension(const std::string &parameter_filename)
{
  const std::string dimension = get_last_value_of_parameter(parameter_filename, "Dimension");
  if (dimension.size() > 0)
    return dealii::Utilities::string_to_int (dimension);
  else
    return 2;
}


#if ASPECT_USE_SHARED_LIBS==1

#ifdef ASPECT_HAVE_LINK_H
// collect the names of the shared libraries linked to by this program. this
// function is a callback for the dl_iterate_phdr() function we call below
int get_names_of_shared_libs (struct dl_phdr_info *info,
                              size_t size,
                              void *data)
{
  reinterpret_cast<std::set<std::string>*>(data)->insert (info->dlpi_name);
  return 0;
}
#endif


// make sure the list of shared libraries we currently link with
// has deal.II only once
void validate_shared_lib_list (const bool before_loading_shared_libs)
{
#ifdef ASPECT_HAVE_LINK_H
  // get the list of all shared libs we currently link against
  std::set<std::string> shared_lib_names;
  dl_iterate_phdr(get_names_of_shared_libs, &shared_lib_names);

  // find everything that is interesting
  std::set<std::string> dealii_shared_lib_names;
  for (std::set<std::string>::const_iterator p = shared_lib_names.begin();
       p != shared_lib_names.end(); ++p)
    if (p->find ("libdeal_II") != std::string::npos)
      dealii_shared_lib_names.insert (*p);

  // produce an error if we load deal.II more than once
  if (dealii_shared_lib_names.size() != 1)
    {
      std::ostringstream error;
      error << "........................................................\n"
            << "ASPECT currently links against different versions of the\n"
            << "deal.II library, namely the ones at these locations:\n";
      for (std::set<std::string>::const_iterator p = dealii_shared_lib_names.begin();
           p != dealii_shared_lib_names.end(); ++p)
        error << "  " << *p << '\n';
      error << "This can not work.\n\n";

      if (before_loading_shared_libs)
        error << "Since this is happening already before opening additional\n"
              << "shared libraries, this means that something must have gone\n"
              << "wrong when you configured deal.II and/or ASPECT. Please\n"
              << "contact the mailing lists for help.\n";
      else
        error << "Since this is happening after opening additional shared\n"
              << "library plugins, this likely means that you have compiled\n"
              << "ASPECT in release mode and the plugin in debug mode, or the\n"
              << "other way around. Please re-compile the plugin in the same\n"
              << "mode as ASPECT.\n";

      error << "........................................................\n";

      // if not success, then throw an exception: ExcMessage on processor 0,
      // QuietException on the others
      if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        {
          AssertThrow (false, dealii::ExcMessage (error.str()));
        }
      else
        throw aspect::QuietException();
    }
#else
  // simply mark the argument as read, to avoid compiler warnings
  (void)before_loading_shared_libs;
#endif
}


#endif


// retrieve a list of shared libraries from the parameter file and
// dlopen them so that we can load plugins declared in them
void possibly_load_shared_libs (const std::string &parameter_filename)
{
  using namespace dealii;


  const std::string shared_libs
    = get_last_value_of_parameter(parameter_filename,
                                  "Additional shared libraries");
  if (shared_libs.size() > 0)
    {
#if ASPECT_USE_SHARED_LIBS==1
      // check up front whether the list of shared libraries is internally
      // consistent or whether we link, for whatever reason, with both the
      // debug and release versions of deal.II
      validate_shared_lib_list (true);

      const std::vector<std::string>
      shared_libs_list = Utilities::split_string_list (shared_libs);

      for (unsigned int i=0; i<shared_libs_list.size(); ++i)
        {
          if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
            std::cout << "Loading shared library <"
                      << shared_libs_list[i]
                      << ">" << std::endl;

          void *handle = dlopen (shared_libs_list[i].c_str(), RTLD_LAZY);
          AssertThrow (handle != NULL,
                       ExcMessage (std::string("Could not successfully load shared library <")
                                   + shared_libs_list[i] + ">. The operating system reports "
                                   + "that the error is this: <"
                                   + dlerror() + ">."));

          // check again whether the list of shared libraries is
          // internally consistent or whether we link with both the
          // debug and release versions of deal.II. this may happen if
          // the plugin was compiled against the debug version of
          // deal.II but aspect itself against the release version, or
          // the other way around
          validate_shared_lib_list (false);

          // on systems where we can detect that both libdeal_II.so and
          // libdeal_II.g.so is loaded, the test above function above will
          // throw an exception and we will terminate. on the other hand, on
          // systems where we can't detect this we should at least mitigate
          // some of the ill effects -- in particular, make sure that
          // deallog is set to use the desired output depth since otherwise
          // we get lots of output from the linear solvers
          deallog.depth_console(0);
        }

      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        std::cout << std::endl;
#else
      if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        {
          std::cerr << std::endl << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
          std::cerr << "You can not load plugins through additional shared libraries " << std::endl
                    << "on systems where you link ASPECT as a static executable."
                    << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
        }
      std::exit (1);
#endif
    }
}


/**
 *  Look up break line sign (\\) at the end of a line and merge this line with the next one.
 *  Return the result as a string in which all lines of the input file are
 *  separated by \n characters, unless the corresponding lines ended
 *  in backslashes.
 */
std::string
expand_backslashes (const std::string &filename)
{
  std::string result;

  unsigned int need_empty_lines = 0;

  std::ifstream input (filename.c_str());
  while (input)
    {
      // get one line and strip spaces at the back
      std::string line;
      std::getline(input, line);
      while ((line.size() > 0)
             && (line[line.size() - 1] == ' ' || line[line.size() - 1] == '\t'))
        line.erase(line.size() - 1, std::string::npos);

      // if the line ends in a backslash, add it without the backslash to
      // the buffer and increase the counter for the number of lines we have
      // just concatenated
      if ((line.size() > 0) && (line[line.size()-1] == '\\'))
        {
          result += line.substr(0, line.size()-1);
          ++need_empty_lines;
        }
      else
        // if it doesn't end in a newline, concatenate the current line
        // with what we have in the buffer and add the \n character
        {
          result += line;
          result += '\n';

          // if we have just added a line (not ending in a backslash)
          // to something that was obtained by addressing backslashes,
          // then add some empty lines to make sure that the line
          // counter is still correct at least for all lines that don't
          // end in a backslash (so that we can ensure that errors
          // message propagating out of ParameterHandler)
          for (; need_empty_lines>0; --need_empty_lines)
            result += '\n';
        }
    }

  // finally return whatever we have in the buffer
  return result;
}


/**
 * Let ParameterHandler parse the input file, here given as a string.
 * Since ParameterHandler unconditionally writes to the screen when it
 * finds something it doesn't like, we get massive amounts of output
 * in parallel computations since every processor writes the same
 * stuff to screen. To avoid this, let processor 0 parse the input
 * first and, if necessary, produce its output. Only if this
 * succeeds, also let the other processors read their input.
 *
 * In case of an error, we need to abort all processors without them
 * having read their data. This is done by throwing an exception of the
 * special class aspect::QuietException that we can catch in main() and terminate
 * the program quietly without generating other output.
 */
void
parse_parameters (const std::string &input_as_string,
                  dealii::ParameterHandler  &prm)
{
  // try reading on processor 0
  bool success = true;
  if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    success = prm.read_input_from_string(input_as_string.c_str());

  // broadcast the result. we'd like to do this with a bool
  // data type but MPI_C_BOOL is not part of old MPI standards.
  // so, do the broadcast in integers
  {
    int isuccess = (success ? 1 : 0);
    MPI_Bcast (&isuccess, 1, MPI_INT, 0, MPI_COMM_WORLD);
    success = (isuccess == 1);
  }

  // if not success, then throw an exception: ExcMessage on processor 0,
  // QuietException on the others
  if (success == false)
    {
      if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
        {
          AssertThrow(false, dealii::ExcMessage ("Invalid input parameter file."));
        }
      else
        throw aspect::QuietException();
    }

  // otherwise, processor 0 was ok reading the data, so we can expect the
  // other processors will be ok as well
  if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) != 0)
    {
      success = prm.read_input_from_string(input_as_string.c_str());
      AssertThrow(success, dealii::ExcMessage ("Invalid input parameter file."));
    }
}


int main (int argc, char *argv[])
{
  using namespace dealii;

  // disable the use of thread. if that is not what you want,
  // use numbers::invalid_unsigned_int instead of 1 to use as many threads
  // as deemed useful by TBB
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, /*n_threads =*/ 1);

  try
    {
      deallog.depth_console(0);
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        print_aspect_header(std::cout);

      if (argc < 2)
        {
          // print usage info only on processor 0
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            std::cout << "\tUsage: ./aspect <parameter_file.prm>" << std::endl;
          return 2;
        }

      // see which parameter file to use
      std::string parameter_filename = argv[1];

      // verify that it can be opened
      {
        std::ifstream parameter_file(parameter_filename.c_str());
        if (!parameter_file)
          {
            const std::string message = (std::string("Input parameter file <")
                                         + parameter_filename + "> not found.");
            if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
              AssertThrow(false, ExcMessage (message));
            return 3;
          }
      }
      // now that we know that the file can, at least in principle, be read
      // try to determine the dimension we want to work in. the default
      // is 2, but if we find a line of the kind "set Dimension = ..."
      // then the last such line wins
      const unsigned int dim = get_dimension(parameter_filename);

      // do the same with lines potentially indicating shared libs to
      // be loaded. these shared libs could contain additional module
      // instantiations for geometries, etc, that would then be
      // available as part of the possible parameters of the input
      // file, so they need to be loaded before we even start processing
      // the parameter file
      possibly_load_shared_libs (parameter_filename);

      // now switch between the templates that code for 2d or 3d. it
      // would be nicer if we didn't have to duplicate code, but the
      // following needs to be known at compile time whereas the dimensionality
      // is only read at run-time
      ParameterHandler prm;

      const std::string input_as_string = expand_backslashes (parameter_filename);
      switch (dim)
        {
          case 2:
          {
            aspect::Simulator<2>::declare_parameters(prm);
            parse_parameters (input_as_string, prm);

            aspect::Simulator<2> flow_problem(MPI_COMM_WORLD, prm);
            flow_problem.run();

            break;
          }

          case 3:
          {
            aspect::Simulator<3>::declare_parameters(prm);
            parse_parameters (input_as_string, prm);

            aspect::Simulator<3> flow_problem(MPI_COMM_WORLD, prm);
            flow_problem.run();

            break;
          }

          default:
            AssertThrow((dim >= 2) && (dim <= 3),
                        ExcMessage ("ASPECT can only be run in 2d and 3d but a "
                                    "different space dimension is given in the parameter file."));
        }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (aspect::QuietException &)
    {
      // Quietly treat an exception used on processors other than
      // root when we already know that processor 0 will generate
      // an exception. We do this to avoid creating too much
      // (duplicate) screen output.
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
