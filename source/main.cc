/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/simulator.h>
#include <aspect/utilities.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/revision.h>
#include <csignal>
#include <string>
#include <thread>
#include <regex>

#ifdef DEBUG
#ifdef ASPECT_USE_FP_EXCEPTIONS
#include <cfenv>
#endif
#endif

#if ASPECT_USE_SHARED_LIBS==1
#  include <dlfcn.h>
#  ifdef ASPECT_HAVE_LINK_H
#    include <link.h>
#  endif
#endif

// This define has to be in exactly one translation unit and sets up the catch testing framework
#define CATCH_CONFIG_RUNNER

// work-around for clang 6 error:
// "error: no member named 'uncaught_exceptions' in namespace 'std'"
// see https://github.com/catchorg/Catch2/issues/1201
#define CATCH_INTERNAL_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS
#define CATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS

#include <catch.hpp>


// get the value of a particular parameter from the contents of the input
// file. return an empty string if not found
std::string
get_last_value_of_parameter(const std::string &parameters,
                            const std::string &parameter_name)
{
  std::string return_value;

  std::istringstream x_file(parameters);
  while (x_file)
    {
      // Get one line and then match a regex to it that matches the parameter
      // we are looking for. Before we do that, strip spaces from the front
      // and back of the line:
      std::string line;
      std::getline(x_file, line);

      while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
        line.erase(0, 1);
      while ((line.size() > 0)
             && (line[line.size() - 1] == ' ' || line[line.size() - 1] == '\t'))
        line.erase(line.size() - 1, std::string::npos);

      std::match_results<std::string::const_iterator> matches;
      const std::string regex = "set[ \t]+" + parameter_name + "[ \t]*=[ \t]*(.*)";
      if (std::regex_match(line, matches, std::regex(regex)))
        {
          // Since the line as a whole matched, the 'matches' variable needs to
          // contain two entries: [0] denotes the whole string, and [1] the
          // one that was matched by the '(.*)' expression.
          Assert (matches.size() == 2, dealii::ExcInternalError());
          return_value = std::string(matches[1].first, matches[1].second);
        }
    }

  return return_value;
}


// Extract the dimension in which to run ASPECT from the
// the contents of the parameter file. This is something that
// we need to do before processing the parameter file since we
// need to know whether to use the dim=2 or dim=3 instantiation
// of the main classes.
//
// This function is essentially the first part of ASPECT to look at the input
// file, so if something is wrong with it, this is the place to generate good
// error messages.
unsigned int
get_dimension(const std::string &parameters)
{
  const std::string dimension = get_last_value_of_parameter(parameters, "Dimension");
  if (dimension.size() > 0)
    {
      // A common problem is that people have .prm files that were generated
      // on Windows, but then run this on Linux where the line endings are
      // different. This is pernicious because it means that the conversion
      // of a string such as "2\r" to an integer fails, but if we print
      // this string, it comes out completely garbled because it contains
      // a carriage-return without a newline -- so the error message looks
      // like this:
      //
      //    >.  While reading the dimension from the input file, ASPECT found a string that can not be converted to an integer: <2
      //
      // Note how the end of the error message overwrites the beginning
      // of the line.
      //
      // To avoid this kind of error, specifically test up front that the
      // text in question does not contain '\r' characters. If we are on
      // linux, then this kind of character would means that the line endings
      // are wrong. On the other hand, if we are on windows, then the
      // getline command we have used in finding 'dimension' would have
      // filtered it out. So its presence points to a problem.

      AssertThrow (dimension.find('\r') == std::string::npos,
                   dealii::ExcMessage ("It appears that your input file uses Windows-style "
                                       "line endings ('\\r\\n') but you are running on a system where "
                                       "the C++ run time environment expects input files to have "
                                       "Unix-style line endings ('\\n'). You need to convert your "
                                       "input file to use the correct line endings before running "
                                       "ASPECT with it."));
      try
        {
          return dealii::Utilities::string_to_int (dimension);
        }
      catch (...)
        {
          AssertThrow (false,
                       dealii::ExcMessage("While reading the dimension from the input file, "
                                          "ASPECT found a string that can not be converted to "
                                          "an integer: <" + dimension + ">."));
          return 0; // we should never get here.
        }
    }
  else
    return 2;
}



#if ASPECT_USE_SHARED_LIBS==1

#ifdef ASPECT_HAVE_LINK_H
// collect the names of the shared libraries linked to by this program. this
// function is a callback for the dl_iterate_phdr() function we call below
int get_names_of_shared_libs (struct dl_phdr_info *info,
                              size_t,
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
  for (const auto &p : shared_lib_names)
    if (p.find ("libdeal_II") != std::string::npos ||
        p.find ("libdeal.ii") != std::string::npos)
      dealii_shared_lib_names.insert (p);

  // produce an error if we load deal.II more than once
  if (dealii_shared_lib_names.size() != 1)
    {
      std::ostringstream error;
      error << "........................................................\n"
            << "ASPECT currently links against different versions of the\n"
            << "deal.II library, namely the ones at these locations:\n";
      for (const auto &p : dealii_shared_lib_names)
        error << "  " << p << '\n';
      error << "This can not work.\n\n";

      if (before_loading_shared_libs)
        error << "Since this is happening already before opening additional\n"
              << "shared libraries, this means that something must have gone\n"
              << "wrong when you configured deal.II and/or ASPECT. Please\n"
              << "contact the forum for help.\n";
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
void possibly_load_shared_libs (const std::string &parameters)
{
  using namespace dealii;


  const std::string shared_libs
    = get_last_value_of_parameter(parameters,
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

      for (const auto &shared_lib : shared_libs_list)
        {
          // The user can specify lib{target}.so, lib{target}.debug.so, or lib{target}.release.so but
          // we need to load the correct file depending on our compilation mode. We will try to make
          // it work regardless of what the users specified:
          std::string filename = shared_lib;

          auto delete_if_ends_with = [](std::string &a, const std::string &b)
          {
            if (a.size()<b.size())
              return;
            if (0==a.compare(a.size()-b.size(), b.size(), b))
              a.erase(a.size()-b.size(), std::string::npos);
          };

          delete_if_ends_with(filename, ".debug.so");
          delete_if_ends_with(filename, ".release.so");
          delete_if_ends_with(filename, ".so");

#ifdef DEBUG
          filename.append(".debug.so");
#else
          filename.append(".release.so");
#endif

          if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
            std::cout << "Loading shared library <"
                      << filename
                      << '>' << std::endl;


          void *handle = dlopen (filename.c_str(), RTLD_LAZY);
          AssertThrow (handle != nullptr,
                       ExcMessage (std::string("Could not successfully load shared library <")
                                   + filename + ">. The operating system reports "
                                   + "that the error is this: <"
                                   + dlerror() +
                                   ">. Did you call 'cmake' and then compile "
                                   "the plugin library you are trying to load, and did "
                                   "you check the spelling of the library's name? "
                                   "Are you running ASPECT in a directory so that the path "
                                   "to the library in question is as specified in the "
                                   ".prm file?"));

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
 * Read until we reach the end of the given stream, and return the contents
 * so obtained in a std::string object that contains the individual
 * lines read separated by `\n` characters.
 */
std::string
read_until_end (std::istream &input)
{
  std::string result;
  while (input)
    {
      std::string line;
      std::getline(input, line);

      result += line + '\n';
    }
  return result;
}



/**
 * Take the name of a parameter file and return all parameters in that file
 * as a string. If @p parameter_file_name is "--" read the parameters from
 * std::cin instead.
 */
std::string
read_parameter_file(const std::string &parameter_file_name,
                    MPI_Comm comm)
{
  using namespace dealii;

  std::string input_as_string;
  const bool i_am_proc_0 = (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

  if (parameter_file_name != "--")
    {
      if (i_am_proc_0 == true &&
          aspect::Utilities::fexists(parameter_file_name) == false &&
          (parameter_file_name=="parameter-file.prm"
           || parameter_file_name=="parameter_file.prm"))
        {
          std::cerr << "***          You should not take everything literally!          ***\n"
                    << "*** Please pass the name of an existing parameter file instead. ***" << std::endl;
          exit(1);
        }

      input_as_string = aspect::Utilities::read_and_distribute_file_content(parameter_file_name, comm);

      // Search and replace include directives in the input file.
      std::match_results<std::string::const_iterator> matches;
      const std::string search_regex = "(?:^|\n)[ \t]*include[ \t]+(.*?)[ \t]*(?:#|\n|$)";
      const std::string replace_regex = "(?:^|\n)[ \t]*include[ \t]+.*";

      unsigned int n_included_files = 0;

      while (std::regex_search(input_as_string, matches, std::regex(search_regex)))
        {
          // Make sure we are not circularly including files. This is not easily possible
          // by making sure included files are unique, because we may have multiple
          // files including the same file in a non-circular way. So we just limit
          // the number of included files to a reasonable number.
          AssertThrow(n_included_files < 15,
                      dealii::ExcMessage("Too many included files in parameter file. You likely have a circular include."));

          // Since the line as a whole matched, the 'matches' variable needs to
          // contain two entries: [0] denotes the whole line, and [1] the
          // part that was matched by the '(.*?)' expression.
          Assert (matches.size() == 2, dealii::ExcInternalError());

          const std::string included_filename(matches.str(1));
          const std::string prefix = "\n# Included content from " + included_filename + ":\n\n";

          // Expand ASPECT_SOURCE_DIR in the included file name, but not the content of the file
          // (to keep the formatting of all parameter files intact, which we will copy into the output directory).
          const std::string expanded_filename = aspect::Utilities::expand_ASPECT_SOURCE_DIR(matches.str(1));

          // Prepend a newline character to the included file content, because if the include directive
          // is not in the first line, we will replace one newlince character
          // from the original input file in the regex_replace below.
          const std::string included_file_content = aspect::Utilities::read_and_distribute_file_content(expanded_filename, comm);

          // Replace the include line with the content of the included file. Note that we only replace the first
          // include line we find (there may be several, which we will replace in subsequent iterations).
          input_as_string = std::regex_replace(input_as_string,
                                               std::regex(replace_regex),
                                               prefix + included_file_content,
                                               std::regex_constants::format_first_only);

          ++n_included_files;
        }
    }
  else
    {
      // As stated in the help string, treat "--" as special: as is common
      // on unix, treat it as a way to read input from stdin.
      // Unfortunately, if you do
      //    echo "abc" | mpirun -np 4 ./aspect
      // then only MPI process 0 gets the data. so we have to
      // read it there, then broadcast it to the other processors
      if (i_am_proc_0)
        {
          input_as_string = read_until_end (std::cin);
          int size = input_as_string.size()+1;
          int ierr = MPI_Bcast (&size,
                                1,
                                MPI_INT,
                                /*root=*/0, MPI_COMM_WORLD);
          AssertThrowMPI(ierr);
          ierr = MPI_Bcast (const_cast<char *>(input_as_string.c_str()),
                            size,
                            MPI_CHAR,
                            /*root=*/0, MPI_COMM_WORLD);
          AssertThrowMPI(ierr);
        }
      else
        {
          // on this side, read what processor zero has broadcast about
          // the size of the input file. then create a buffer to put the
          // text in, get it from processor 0, and copy it to
          // input_as_string
          int size;
          int ierr = MPI_Bcast (&size, 1,
                                MPI_INT,
                                /*root=*/0, MPI_COMM_WORLD);
          AssertThrowMPI(ierr);

          std::vector<char> p (size);
          ierr = MPI_Bcast (p.data(), size,
                            MPI_CHAR,
                            /*root=*/0, MPI_COMM_WORLD);
          AssertThrowMPI(ierr);
          input_as_string = p.data();
        }
    }

  return input_as_string;
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
    try
      {
        prm.parse_input_from_string(input_as_string);
      }
    catch (const dealii::ExceptionBase &e)
      {
        success = false;
        e.print_info(std::cerr);
        std::cerr << std::endl;
      }


  // broadcast the result. we'd like to do this with a bool
  // data type but MPI_C_BOOL is not part of old MPI standards.
  // so, do the broadcast in integers
  {
    int isuccess = (success ? 1 : 0);
    const int ierr = MPI_Bcast (&isuccess, 1, MPI_INT, 0, MPI_COMM_WORLD);
    AssertThrowMPI(ierr);
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
      prm.parse_input_from_string(input_as_string);
    }
}



/**
 * Print help text
 */
void print_help()
{
  std::cout << "Usage: ./aspect [args] <parameter_file.prm>   (to read from an input file)\n"
            << "    or ./aspect [args] --                     (to read parameters from stdin)\n"
            << std::endl;
  std::cout << "    optional arguments [args]:\n"
            << "       -h, --help             (for this usage help)\n"
            << "       -v, --version          (for information about library versions)\n"
            << "       -j, --threads          (to use multi-threading)\n"
            << "       --output-json          (print parameters in JSON format to standard output and exit)\n"
            << "       --output-xml           (print parameters in XML format to standard output and exit)\n"
            << "       --output-plugin-graph  (write a representation of all plugins to standard output and exit)\n"
            << "       --validate             (parse parameter file and exit or report errors)\n"
            << "       --test                 (run the unit tests from unit_tests/, run --test -h for more info)\n"
            << std::endl;
}



// hook into SIGABRT/SIGFPE and kill off the program
void signal_handler(int signal)
{
  if (signal == SIGABRT)
    {
      std::cerr << "SIGABRT received\n";
    }
  else if (signal == SIGFPE)
    {
      std::cerr << "SIGFPE received\n";
    }
  else
    {
      std::cerr << "Unexpected signal " << signal << " received\n";
    }

  // Kill the program without performing any other cleanup, which would likely
  // lead to a deadlock.
  std::_Exit(EXIT_FAILURE);
}



template <int dim>
void
run_simulator(const std::string &raw_input_as_string,
              const std::string &input_as_string,
              const bool output_json,
              const bool output_xml,
              const bool output_plugin_graph,
              const bool validate_only)
{
  using namespace dealii;

  ParameterHandler prm;
  const bool i_am_proc_0 = (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);
  aspect::Simulator<dim>::declare_parameters(prm);

  if (validate_only)
    {
      try
        {
          parse_parameters (input_as_string, prm);
        }
      catch (...)
        {
          throw aspect::QuietException();
        }
      if (i_am_proc_0)
        std::cout << "The provided parameter file is valid."
                  "\n\n"
                  "Note: This validation only checks parameter file syntax errors, like typos\n"
                  "in keywords or parameter names, and that each parameter value satisfies a\n"
                  "basic check by itself. However, it may miss more nuanced errors that\n"
                  "are only checked when the model actually begins running. In particular,\n"
                  "checks that involve two or more parameters can not be verified at this\n"
                  "stage of an ASPECT run. Examples for such errors that can not already\n"
                  "be reported here are: (i) That every boundary is assigned exactly one type\n"
                  "of boundary condition; (ii) that parameters that take a list with values for\n"
                  "each composition field receive a list of the correct size."
                  << std::endl;
      return;
    }

  parse_parameters (input_as_string, prm);

  if (output_json)
    {
      if (i_am_proc_0)
        prm.print_parameters(std::cout, ParameterHandler::JSON);
    }
  else if (output_xml)
    {
      if (i_am_proc_0)
        prm.print_parameters(std::cout, ParameterHandler::XML);
    }
  else if (output_plugin_graph)
    {
      aspect::Simulator<dim> simulator(MPI_COMM_WORLD, prm);
      if (i_am_proc_0)
        simulator.write_plugin_graph (std::cout);
    }
  else
    {
      aspect::Simulator<dim> simulator(MPI_COMM_WORLD, prm);
      if (i_am_proc_0)
        {
          // write create output/original.prm containing exactly what we got
          // started with:
          std::string output_directory = prm.get ("Output directory");
          if (output_directory.size() == 0)
            output_directory = "./";
          else if (output_directory[output_directory.size()-1] != '/')
            output_directory += "/";

          std::ofstream file(output_directory + "original.prm");
          file << raw_input_as_string;

          // If using the Geodynamic World Builder, create output/original.wb
          // containing the exact file used to create the world:
          std::string world_builder_file = prm.get("World builder file");
          if (world_builder_file != "")
            {
              // TODO: We just want to make a copy of the file, but to do this
              // platform-independently we need the C++17 filesystem header.
              // Update this when we require C++17
              std::ifstream  wb_source(world_builder_file, std::ios::binary);
              std::ofstream  wb_destination(output_directory + "original.wb",   std::ios::binary);

              wb_destination << wb_source.rdbuf();
            }
        }

      simulator.run();
    }
}



int main (int argc, char *argv[])
{
  using namespace dealii;

#ifdef DEBUG
#ifdef ASPECT_USE_FP_EXCEPTIONS
  // Some implementations seem to not initialize the floating point exception
  // bits to zero. Make sure we start from a clean state.
  feclearexcept(FE_DIVBYZERO|FE_INVALID);

  // enable floating point exceptions
  feenableexcept(FE_DIVBYZERO|FE_INVALID);
#endif
#endif

  std::string prm_name = "";
  bool output_json         = false;
  bool output_xml          = false;
  bool output_plugin_graph = false;
  bool output_version      = false;
  bool output_help         = false;
  bool use_threads         = false;
  bool run_unittests       = false;
  bool validate_only       = false;
  int current_argument = 1;

  // Loop over all command line arguments. Handle a number of special ones
  // starting with a dash, and then take the first non-special one as the
  // name of the input file. We will later check that there are no further
  // arguments left after that.
  while (current_argument<argc)
    {
      const std::string arg = argv[current_argument];
      ++current_argument;
      if (arg == "--output-json")
        {
          output_json = true;
        }
      else if (arg == "--output-xml")
        {
          output_xml = true;
        }
      else if (arg == "--output-plugin-graph")
        {
          output_plugin_graph = true;
        }
      else if (arg=="-h" || arg =="--help")
        {
          output_help = true;
        }
      else if (arg=="-v" || arg =="--version")
        {
          output_version = true;
        }
      else if (arg=="-j" || arg =="--threads")
        {
          use_threads = true;
        }
      else if (arg == "--test")
        {
          run_unittests = true;
          break;
        }
      else if (arg == "--validate")
        {
          validate_only = true;
        }
      else
        {
          // Not a special argument, so we assume that this is the .prm
          // filename (or "--"). We can now break out of this loop because
          // we are not going to parse arguments passed after the filename.
          prm_name = arg;
          break;
        }
    }

  // There might be remaining arguments for PETSc, only hand those over to
  // the MPI initialization, but not the ones we parsed above.
  int n_remaining_arguments = argc - current_argument;
  char **remaining_arguments = (n_remaining_arguments > 0) ? &argv[current_argument] : nullptr;

  try
    {
      // Note: we initialize this class inside the try/catch block and not
      // before, so that the destructor of this instance can react if we are
      // currently unwinding the stack if an unhandled exception is being
      // thrown to avoid MPI deadlocks.
      Utilities::MPI::MPI_InitFinalize mpi_initialization(n_remaining_arguments,
                                                          remaining_arguments,
                                                          (use_threads ?
                                                           numbers::invalid_unsigned_int :
                                                           1));

      if (run_unittests)
        {
          // Construct new_argc, new_argv from argc, argv for catch without
          // the "--test" arg, so we can control catch from the command
          // line. It turns out catch needs argv[0] to be the executable name
          // so we can not use remaining_arguments from above.
          int new_argc = n_remaining_arguments + 1;
          std::vector<char *> args; // use to construct a new argv of type char **
          args.emplace_back(argv[0]);
          for (int i=0; i<n_remaining_arguments; ++i)
            args.emplace_back(argv[i+current_argument]);
          char **new_argv = args.data();

          // Finally run catch
          return Catch::Session().run(new_argc, new_argv);
        }


      deallog.depth_console(0);

      const bool i_am_proc_0 = (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

      if (i_am_proc_0)
        {
          // Output header, except for a clean output for xml or plugin graph
          if (!output_xml && !output_json && !output_plugin_graph && !validate_only)
            print_aspect_header(std::cout);

          if (output_help)
            print_help();

          // If only help or version is requested, we are done.
          if (output_help || output_version)
            return 0;
        }
      else
        {
          // If only help or version is requested, we are done.
          if (output_help || output_version)
            return 0;

          // We hook into the abort handler on ranks != 0 to avoid an MPI
          // deadlock. The deal.II library will call std::abort() when an
          // Assert is triggered, which can lead to a deadlock because it
          // runs the things that are associated with atexit() which may
          // itself trigger MPI communication. The same happens for other
          // signals we may trigger, such as floating point exceptions
          // (SIGFPE).
          //
          // We work around this by immediately calling _Exit in the
          // signal handler and thus aborting the program without running
          // cleanup functions set via atexit(). This is only necessary on
          // rank != 0 for some reason.
          std::signal(SIGABRT, signal_handler);
          std::signal(SIGFPE, signal_handler);
        }

      // If no parameter given or somebody gave additional parameters,
      // show help and exit.
      if ((prm_name == "") || (current_argument < argc))
        {
          if (i_am_proc_0)
            print_help();
          return 0;
        }

      // See where to read input from, then do the reading and
      // put the contents of the input into a string.
      const std::string raw_input_as_string = read_parameter_file(prm_name, MPI_COMM_WORLD);

      // Replace $ASPECT_SOURCE_DIR in the input so that include statements
      // like "include $ASPECT_SOURCE_DIR/tests/bla.prm" work.
      const std::string input_as_string = aspect::Utilities::expand_ASPECT_SOURCE_DIR(raw_input_as_string);

      // Determine the dimension we want to work in. the default
      // is 2, but if we find a line of the kind "set Dimension = ..."
      // then the last such line wins.
      const unsigned int dim = get_dimension(input_as_string);

      // Do the same with lines potentially indicating shared libs to
      // be loaded. These shared libs could contain additional module
      // instantiations for geometries, etc, that would then be
      // available as part of the possible parameters of the input
      // file, so they need to be loaded before we even start processing
      // the parameter file.
      possibly_load_shared_libs (input_as_string);

      // Now switch between the templates that start the model for 2d or 3d.
      switch (dim)
        {
          case 2:
          {
            run_simulator<2>(raw_input_as_string,input_as_string,output_json,output_xml,output_plugin_graph,validate_only);
            break;
          }
          case 3:
          {
            run_simulator<3>(raw_input_as_string,input_as_string,output_json,output_xml,output_plugin_graph,validate_only);
            break;
          }
          default:
            AssertThrow((dim >= 2) && (dim <= 3),
                        ExcMessage ("ASPECT can only be run in 2d and 3d but a "
                                    "different space dimension is given in the parameter file."));
        }
    }
  catch (ExceptionBase &exc)
    {
      // report name of the deal.II exception:
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception '" << exc.get_exc_name() << "'"
                << " on rank " << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                << " on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception"
                << " on rank " << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                << " on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
  catch (aspect::QuietException &)
    {
      // Quietly treat an exception used on processors other than root
      // when we already know that processor 0 will generate an
      // exception. We do this to avoid creating too much (duplicate)
      // screen output. Note that QuietException is not derived from
      // std::exception, so the order of this and the previous 'catch'
      // block does not matter.
      //
      // Sleep a few seconds before aborting. This allows text output from
      // other ranks to be printed before the MPI implementation might kill
      // the computation.
      std::this_thread::sleep_for(std::chrono::seconds(5));

      MPI_Abort(MPI_COMM_WORLD, 1);
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

      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }

  return 0;
}
