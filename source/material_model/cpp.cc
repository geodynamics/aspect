/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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


#include <aspect/material_model/cpp.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    bool
    CPP<dim>::
    is_compressible () const
    {
      return is_compressible_param;
    }

    template <int dim>
    double
    CPP<dim>::
    reference_viscosity () const
    {
      return reference_viscosity_param;
    }

    template <int dim>
    void
    CPP<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      user_eval (in, out, this);
    }

    /*
     * Attempt to locate an indenter command to beautify the generated
     * C++ code.
     *
     * @return Command line tool for indenting c++ source code, in a format
     * such that its only missing argument is the file name, to be inserted
     * at the end of the string.
     */
    std::string
    find_indenter()
    {
      if (execute("which astyle >/dev/null") == 0)
        if (execute("[ -f '" ASPECT_SOURCE_DIR "/doc/astyle.rc' ]") == 0)
          return "astyle --options=" ASPECT_SOURCE_DIR "/doc/astyle.rc";
        else
          return "astyle";
      else if (execute("which astyle2.04 >/dev/null") == 0)
        if (execute("[ -f '" ASPECT_SOURCE_DIR "/doc/astyle.rc' ]") == 0)
          return "astyle2.04 --options=" ASPECT_SOURCE_DIR "/doc/astyle.rc";
        else
          return "astyle2.04";
      else if (execute("which indent >/dev/null") == 0)
        return "indent";
      else if (execute("which clang-format >/dev/null") == 0)
        return "clang-format";

      return "";
    }

    /*
     * Make sure user didn't specify 'none' along with other nonlinear
     * dependencies.
     */
    void
    assert_valid_deps (std::vector<std::string> list, std::string param_name)
    {
      AssertThrow(list.size() == 1
                  || (std::find(list.begin(),list.end(),"none") == list.end()),
                  ExcMessage("The property 'none' in the parameter '"+param_name+"' "
                             "cannot be combined with any other options."));
    }

    /*
     * Write a simple CMakeLists.txt for the generated code, to simplify
     * linking with Aspect and Deal.II. Basically, treat the code as a plugin,
     * since it occupies the same namespace, and should expose most of the same
     * interfaces as a native plugin.
     */
    void
    generate_cmakelists (const std::string cdir)
    {
      std::ofstream cmfile;
      cmfile.open((cdir + "CMakeLists.txt").c_str(), std::ios_base::out);
      cmfile << "CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)\n"
             << "\n"
             << "FIND_PACKAGE(Aspect QUIET HINTS "
             << ASPECT_SOURCE_DIR <<" "
             << ASPECT_CURRENT_BINARY_DIR << " "
             << "  ${Aspect_DIR} $ENV{ASPECT_DIR})\n"
             << "\n"
             << "IF (NOT Aspect_FOUND)\n"
             << "  MESSAGE(FATAL_ERROR \"\\n\"\n"
             << "    \"Could not find a valid ASPECT build/installation directory. \"\n"
             << "    \"Please specify the directory where you are building ASPECT by\\n\"\n"
             << "    \"setting the environment variable ASPECT_DIR in your shell.\")\n"
             << "ENDIF ()\n"
             << "\n"
             << "DEAL_II_INITIALIZE_CACHED_VARIABLES()\n"
             << "\n"
             << "SET(TARGETmaterial \"material\")\n"
             << "PROJECT(${TARGETmaterial})\n"
             << "\n"
             << "ADD_LIBRARY(${TARGETmaterial} SHARED material.cc)\n"
             << "ASPECT_SETUP_PLUGIN(${TARGETmaterial})\n";
      cmfile.close();
    }

    /*
     * Check whether an arbitrary code snippet references a particular
     * variable name.
     *
     * Performs a regular expression search on the code for the
     * requested variable name, making sure it's surrounded by
     * characters that are not allowed in C++ variable names.
     *
     * @return True if `code_snippet' references the variable named in
     * `dep', otherwise returns False.
     */
    inline
    bool
    depends (const std::string code_snippet,
             const std::string dep)
    {
      const std::regex exp ("[^_[:alnum:]]("+dep+")[^_[:alnum:]]");
      return std::regex_search (code_snippet,exp);
    }

    /*
     * To avoid unnecessary computational work, we make each generated
     * function with the minimum number of parameters for it to operate.
     * For example, if the user code never queries the 'strain_rate'
     * variable, we don't need to bother passing it (or a pointer to it)
     * as an argument.
     *
     * This function generates the parameter list for each function declaration
     * and its corresponding function call, by checking whether the function
     * contains a reference to each of the values to which it is entitled access.
     *
     * @return A pair of strings representing the function declaration
     * parameter list, and the function call variable list, respectively.
     */
    inline
    std::pair<std::string, std::string>
    mk_declarations (const unsigned int dim,
                     const std::string code_snippet)
    {
      const std::string dimstr  = Utilities::int_to_string(dim);
      std::vector<std::string> function_declaration_arguments;
      std::vector<std::string> args;
      if (depends(code_snippet,"position"))
        {
          function_declaration_arguments.push_back ("const Point<"+dimstr+"> &position");
          args.push_back ("in.position[_i]");
        }
      if (depends(code_snippet,"temperature"))
        {
          function_declaration_arguments.push_back ("const double &temperature");
          args.push_back ("in.temperature[_i]");
        }
      if (depends(code_snippet,"pressure"))
        {
          function_declaration_arguments.push_back ("const double &pressure");
          args.push_back ("in.pressure[_i]");
        }
      if (depends(code_snippet,"pressure_gradient"))
        {
          function_declaration_arguments.push_back ("const Tensor<1,"+dimstr+"> &pressure_gradient");
          args.push_back ("in.pressure_gradient[_i]");
        }
      if (depends(code_snippet,"composition"))
        {
          function_declaration_arguments.push_back ("const std::vector<double> &composition");
          args.push_back ("in.composition[_i]");
        }
      if (depends(code_snippet,"strain_rate"))
        {
          function_declaration_arguments.push_back ("const SymmetricTensor<2,"+dimstr+"> &strain_rate");
          args.push_back ("in.strain_rate[_i]");
        }
      if (depends(code_snippet,"simulator"))
        {
          function_declaration_arguments.push_back ("const ::aspect::SimulatorAccess<"+dimstr+"> *simulator");
          args.push_back ("simulator");
        }

      std::string function_declaration_argumentlist = "";
      std::string function_call_argumentlist = "";
      for (unsigned int i=0; i<function_declaration_arguments.size(); ++i)
        {
          const bool islast = (i == function_declaration_arguments.size()-1);
          function_declaration_argumentlist += function_declaration_arguments[i] + (!islast ? ",\n" : "");
          function_call_argumentlist += args[i] + (!islast ? ",\n" : "");
        }
      return std::pair<std::string,std::string> (function_declaration_argumentlist, function_call_argumentlist);
    }

    /*
     * Generates source code for a function that evaluates material properties.
     * It takes code snippets for evaluating each material property, as properties
     * of a `UserCode` struct, and writes code to a file named `fname`. If
     * `indenter` is not blank, it will attempt to beautify the resulting source
     * code.
     */
    void
    generate_src (const unsigned int dim,
                  const UserCode &user_code,
                  const std::string &indenter,
                  const std::string &fname)
    {
      // The functions for evaluating each field each have a dynamically chosen
      // set of parameters, based on which parameters are actually used.
      //
      // Build up the complete external source code by generating each function
      // declaration and use commands dynamically.
      std::string function_defs = "";
      std::string function_calls = "";

      // "Update" function
      if (user_code.update_function != "")
        {
          std::pair<std::string,std::string> arglist = mk_declarations(dim, user_code.update_function);
          function_defs += "void update (" + arglist.first + ") {\n" + user_code.update_function + ";\n}\n";
          function_calls += "local.update (" + arglist.second + ");\n";
        }

      // Function for compositional field reaction terms.
      if (user_code.reaction_function != "")
        {
          std::pair<std::string,std::string> arglist = mk_declarations(dim, user_code.reaction_function);
          const std::string sep = (arglist.first != "") ? ",\n" : "";
          function_defs += "double reaction (const unsigned int c" + sep + arglist.first + ") {\n" + user_code.reaction_function + ";\n}\n";
          function_calls += "for (unsigned int c=0; c<in.composition[_i].size(); ++c)\n";
          function_calls += "  out.reaction_terms[_i][c] = local.reaction (c" + sep + arglist.second + ");\n";
        }

      // All other functions are basically the same, so just create them via macro...
#define MKFUNCTIONS(_OUTPUTVAR,_FUNCTION) \
  do { \
      if (_FUNCTION != "") \
        { \
          std::pair<std::string,std::string> arglist = mk_declarations(dim, _FUNCTION); \
          function_defs += "double " #_OUTPUTVAR " (" + arglist.first + ") {\n" + _FUNCTION + ";\n}\n"; \
          function_calls += "out." #_OUTPUTVAR "[_i] = local." #_OUTPUTVAR " (" + arglist.second + ");\n"; \
        } \
    } while (false)

      MKFUNCTIONS(viscosities,user_code.viscosity_function);
      MKFUNCTIONS(densities,user_code.density_function);
      MKFUNCTIONS(thermal_conductivities,user_code.thermal_conductivity_function);
      MKFUNCTIONS(thermal_expansion_coefficients,user_code.thermal_expansivity_function);
      MKFUNCTIONS(specific_heat,user_code.specific_heat_function);
      MKFUNCTIONS(compressibilities,user_code.compressibility_function);
      MKFUNCTIONS(entropy_derivative_pressure,user_code.entropy_derivative_p_function);
      MKFUNCTIONS(entropy_derivative_temperature,user_code.entropy_derivative_t_function);

      // Write the code to an output file. Make liberal use of semicolons
      // after user code, since it's easy to forget them in parameter files,
      // and extra ones don't bother non-pedantic compilers.
      std::ofstream srcfile;
      srcfile.open((fname).c_str(), std::ios_base::out);
      for (unsigned int i = 0; i<user_code.includes.size(); ++i)
        srcfile << "#include <" << user_code.includes[i] << ">\n";
      srcfile << "#include <aspect/material_model/interface.h>\n"
              << "#include <aspect/simulator.h>\n"
              << "#include <aspect/simulator_access.h>\n"
              << "namespace aspect {\n"
              << "  namespace MaterialModel {\n"
              << "    using namespace dealii;\n"
              << "    class Local_\n"
              << "    {\n"
              << "      public:\n"
              << "        const unsigned int dim = "<<dim<<";\n"
              <<          user_code.variable_defs << ";\n"
              <<          function_defs
              << "    } local;\n"
              << "\n"
              << "    extern \"C\" void\n"
              << "    eval (const MaterialModel::MaterialModelInputs<"<<dim<<"> &in,\n"
              << "          MaterialModel::MaterialModelOutputs<"<<dim<<"> &out,\n"
              << "          const ::aspect::SimulatorAccess<"<<dim<<"> *simulator)\n"
              << "    {\n"
              << "      for (unsigned int _i=0; _i<in.position.size(); ++_i)\n"
              << "        {\n"
              <<            function_calls
              << "        }\n"
              << "    }\n"
              << "  }\n"
              << "}\n";
      srcfile.close();

      // Try (not very hard) to properly indent the generated code.
      // Ignore nonzero exit codes.
      if (indenter != "")
        execute(indenter + " '" + fname + "'");
    }


    template <int dim>
    void
    CPP<dim>::declare_parameters (ParameterHandler &prm)
    {
      const std::string pattern_of_deps
        = "none|compositional fields|pressure|strain rate|temperature";
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("C plus plus");
        {
          prm.declare_entry ("Is compressible", "false",
                             Patterns::Bool (),
                             "Whether or not the model is compressible.");
          prm.declare_entry ("Reference viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the reference viscosity $\\eta$. Units: $kg/m/s$.");
          prm.declare_entry ("Includes", "",
                             Patterns::List (Patterns::Anything()),
                             "List additional library headers to include in material model. "
                             "For example: `set Includes =  stdio.h, math.h`. By default "
                             "only the headers required to operate in a MaterialModel are "
                             "included.");
          prm.declare_entry ("Global variable definitions", "",
                             Patterns::Anything (),
                             "Define variables (and whatever else you want) to be accessible "
                             "by all other functions.");
          prm.declare_entry ("Update function", "",
                             Patterns::Anything (),
                             "Contents of a function that can modify global variables "
                             "defined in 'Material model/C plus plus/Global variable "
                             "definitions'. It is invoked before any of the evaluation "
                             "functions are called at each model position."
                             "\n\nThe state variables `position', `temperature', `pressure', "
                             "`pressure_gradient', `velocity', `composition', `strain_rate' are "
                             "available within this function, as well as a SimulatorAccess object "
                             "pointer called `simulator'.");
          prm.declare_entry ("Viscosity function", "",
                             Patterns::Anything (),
                             "Contents of a function for calculating viscosity. It has access to "
                             "all variables defined in 'Material model/C plus plus/Global variable "
                             "definitions'."
                             "\n\nThe state variables `position', `temperature', `pressure', "
                             "`pressure_gradient', `velocity', `composition', `strain_rate' are "
                             "available within this function, as well as a SimulatorAccess object "
                             "pointer called `simulator'.");
          prm.declare_entry ("Density function", "",
                             Patterns::Anything (),
                             "Contents of a function for calculating density. It has access to "
                             "all variables defined in 'Material model/C plus plus/Global variable "
                             "definitions'."
                             "\n\nThe state variables `position', `temperature', `pressure', "
                             "`pressure_gradient', `velocity', `composition', `strain_rate' are "
                             "available within this function, as well as a SimulatorAccess object "
                             "pointer called `simulator'.");
          prm.declare_entry ("Thermal conductivity function", "",
                             Patterns::Anything (),
                             "Contents of a function for calculating thermal conductivity. It has access to "
                             "all variables defined in 'Material model/C plus plus/Global variable "
                             "definitions'."
                             "\n\nThe state variables `position', `temperature', `pressure', "
                             "`pressure_gradient', `velocity', `composition', `strain_rate' are "
                             "available within this function, as well as a SimulatorAccess object "
                             "pointer called `simulator'.");
          prm.declare_entry ("Thermal expansivity function", "",
                             Patterns::Anything (),
                             "Contents of a function for calculating thermal expansivity. It has access to "
                             "all variables defined in 'Material model/C plus plus/Global variable "
                             "definitions'."
                             "\n\nThe state variables `position', `temperature', `pressure', "
                             "`pressure_gradient', `velocity', `composition', `strain_rate' are "
                             "available within this function, as well as a SimulatorAccess object "
                             "pointer called `simulator'.");
          prm.declare_entry ("Specific heat function", "",
                             Patterns::Anything (),
                             "Contents of a function for calculating specific heat. It has access to "
                             "all variables defined in 'Material model/C plus plus/Global variable "
                             "definitions'."
                             "\n\nThe state variables `position', `temperature', `pressure', "
                             "`pressure_gradient', `velocity', `composition', `strain_rate' are "
                             "available within this function, as well as a SimulatorAccess object "
                             "pointer called `simulator'.");
          prm.declare_entry ("Compressibility function", "",
                             Patterns::Anything (),
                             "Contents of a function for calculating compressibility. It has access to "
                             "all variables defined in 'Material model/C plus plus/Global variable "
                             "definitions'."
                             "\n\nThe state variables `position', `temperature', `pressure', "
                             "`pressure_gradient', `velocity', `composition', `strain_rate' are "
                             "available within this function, as well as a SimulatorAccess object "
                             "pointer called `simulator'.");
          prm.declare_entry ("Entropy derivative pressure function", "",
                             Patterns::Anything (),
                             "Contents of a function for calculating entropy derivative with respect to pressure. It has access to "
                             "all variables defined in 'Material model/C plus plus/Global variable "
                             "definitions'."
                             "\n\nThe state variables `position', `temperature', `pressure', "
                             "`pressure_gradient', `velocity', `composition', `strain_rate' are "
                             "available within this function, as well as a SimulatorAccess object "
                             "pointer called `simulator'.");
          prm.declare_entry ("Entropy derivative temperature function", "",
                             Patterns::Anything (),
                             "Contents of a function for calculating entropy derivative with respect to temperature. It has access to "
                             "all variables defined in 'Material model/C plus plus/Global variable "
                             "definitions'."
                             "\n\nThe state variables `position', `temperature', `pressure', "
                             "`pressure_gradient', `velocity', `composition', `strain_rate' are "
                             "available within this function, as well as a SimulatorAccess object "
                             "pointer called `simulator'.");
          prm.declare_entry ("Reaction function", "",
                             Patterns::Anything (),
                             "Contents of a function for calculating reaction terms. It has access to "
                             "all variables defined in 'Material model/C plus plus/Global variable "
                             "definitions'. It receives an integer called `c' that denotes the particular "
                             "compositional field for which to calculate a reaction rate."
                             "\n\nThe state variables `position', `temperature', `pressure', "
                             "`pressure_gradient', `velocity', `composition', `strain_rate' are "
                             "available within this function, as well as a SimulatorAccess object "
                             "pointer called `simulator'.");
          prm.declare_entry ("List of viscosity dependencies", "none",
                             Patterns::MultipleSelection(pattern_of_deps),
                             "A comma separated list of nonlinear dependencies for viscosity. "
                             "The following options are available:\n\n"
                             +
                             pattern_of_deps);
          prm.declare_entry ("List of density dependencies", "none",
                             Patterns::MultipleSelection(pattern_of_deps),
                             "A comma separated list of nonlinear dependencies for density. "
                             "The following options are available:\n\n"
                             +
                             pattern_of_deps);
          prm.declare_entry ("List of compressibility dependencies", "none",
                             Patterns::MultipleSelection(pattern_of_deps),
                             "A comma separated list of nonlinear dependencies for compressibility. "
                             "The following options are available:\n\n"
                             +
                             pattern_of_deps);
          prm.declare_entry ("List of specific heat dependencies", "none",
                             Patterns::MultipleSelection(pattern_of_deps),
                             "A comma separated list of nonlinear dependencies for specific heat. "
                             "The following options are available:\n\n"
                             +
                             pattern_of_deps);
          prm.declare_entry ("List of thermal conductivity dependencies", "none",
                             Patterns::MultipleSelection(pattern_of_deps),
                             "A comma separated list of nonlinear dependencies for thermal conductivity. "
                             "The following options are available:\n\n"
                             +
                             pattern_of_deps);
          prm.declare_entry ("Indenter", find_indenter(),
                             Patterns::Anything (),
                             "Command line tool for indenting generated C++ code. Useful for debugging.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    CPP<dim>::parse_parameters (ParameterHandler &prm)
    {
#if ASPECT_SAFE==0
#if ASPECT_USE_SHARED_LIBS==1
      prm.enter_subsection ("Material model");
      {
        prm.enter_subsection ("C plus plus");
        {

          // Read global plugin settings
          is_compressible_param               = prm.get_bool ("Is compressible");
          reference_viscosity_param           = prm.get_double ("Reference viscosity");
          const std::string indenter          = prm.get ("Indenter");

          UserCode user_code;
          // Read user-defined code snippets
          user_code.includes                       = Utilities::split_string_list(prm.get ("Includes"));
          user_code.variable_defs                  = prm.get ("Global variable definitions");
          user_code.update_function                = prm.get ("Update function");
          user_code.viscosity_function             = prm.get ("Viscosity function");
          user_code.density_function               = prm.get ("Density function");
          user_code.thermal_conductivity_function  = prm.get ("Thermal conductivity function");
          user_code.thermal_expansivity_function   = prm.get ("Thermal expansivity function");
          user_code.specific_heat_function         = prm.get ("Specific heat function");
          user_code.compressibility_function       = prm.get ("Compressibility function");
          user_code.entropy_derivative_p_function  = prm.get ("Entropy derivative pressure function");
          user_code.entropy_derivative_t_function  = prm.get ("Entropy derivative temperature function");
          user_code.reaction_function              = prm.get ("Reaction function");

          // Read nonlinear dependencies
          const std::vector<std::string> viscosity_deps =
            Utilities::split_string_list (prm.get ("List of viscosity dependencies"));
          const std::vector<std::string> density_deps =
            Utilities::split_string_list (prm.get ("List of density dependencies"));
          const std::vector<std::string> compressibility_deps =
            Utilities::split_string_list (prm.get ("List of compressibility dependencies"));
          const std::vector<std::string> specific_heat_deps =
            Utilities::split_string_list (prm.get ("List of specific heat dependencies"));
          const std::vector<std::string> thermal_conductivity_deps =
            Utilities::split_string_list (prm.get ("List of thermal conductivity dependencies"));

          // Assert that user didn't mix 'none' with other nonlinear dependencies
          assert_valid_deps (viscosity_deps,
                             "Material model/C plus plus/List of viscosity dependencies");
          assert_valid_deps (density_deps,
                             "Material model/C plus plus/List of density dependencies");
          assert_valid_deps (compressibility_deps,
                             "Material model/C plus plus/List of compressibility dependencies");
          assert_valid_deps (specific_heat_deps,
                             "Material model/C plus plus/List of specific heat dependencies");
          assert_valid_deps (thermal_conductivity_deps,
                             "Material model/C plus plus/List of thermal conductivity dependencies");

          // Declare dependencies
          this->model_dependence.viscosity              = NonlinearDependence::uninitialized;
          this->model_dependence.density                = NonlinearDependence::uninitialized;
          this->model_dependence.compressibility        = NonlinearDependence::uninitialized;
          this->model_dependence.specific_heat          = NonlinearDependence::uninitialized;
          this->model_dependence.thermal_conductivity   = NonlinearDependence::uninitialized;

          for (unsigned int i=0; i<viscosity_deps.size(); ++i)
            this->model_dependence.viscosity |= str2dep(viscosity_deps[i]);
          for (unsigned int i=0; i<density_deps.size(); ++i)
            this->model_dependence.density |= str2dep(density_deps[i]);
          for (unsigned int i=0; i<compressibility_deps.size(); ++i)
            this->model_dependence.compressibility |= str2dep(compressibility_deps[i]);
          for (unsigned int i=0; i<specific_heat_deps.size(); ++i)
            this->model_dependence.specific_heat |= str2dep(specific_heat_deps[i]);
          for (unsigned int i=0; i<thermal_conductivity_deps.size(); ++i)
            this->model_dependence.thermal_conductivity |= str2dep(thermal_conductivity_deps[i]);

          // Choose filenames and create any necessary directories
          std::string cdir = this->get_output_directory() + "cpp_material/";
          std::string sfname = cdir + "material.cc";
          std::string ofname = cdir + "libmaterial.so";
          Utilities::create_directory(cdir, MPI_COMM_WORLD, true);

          // Generate and compile code
          int error;
          if (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
            {
              generate_src (dim, user_code, indenter, sfname);

              generate_cmakelists (cdir);

              std::cout << "Compiling material model <" + sfname + ">.\n" << std::endl;
              error = execute ("cd '" + cdir + "' && cmake . > build.log 2>&1 && make >> build.log 2>&1");
              MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
              AssertThrow (error == 0,
                           ExcMessage (std::string("Compiling material model failed. See build log at <")
                                       + cdir + "build.log>."));
            }
          else
            {
              MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
              if (error != 0)
                throw aspect::QuietException();
            }

          // Load function via dlopen
          void *handle = dlopen (ofname.c_str(), RTLD_LAZY);
          AssertThrow (handle != nullptr,
                       ExcMessage (std::string("Could not load compiled material model ")
                                   + "from shared library <" + ofname + ">, based on generated "
                                   + "source code <" + sfname + ">. The operating system "
                                   + "returned error: <" + dlerror() + ">."));

          user_eval = (eval_t) dlsym(handle, "eval");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
#else
      // No shared libraries
      AssertThrow(false, ExcMessage("The 'cpp' material model requires "
                                    "that Aspect can load shared libraries. "
                                    "Please check your build configuration."));
#endif
#else
      // ASPECT_SAFE == 1
      AssertThrow(false, ExcMessage("The 'cpp' material model requires "
                                    "the ability to execute arbitrary code, "
                                    "which is disabled by the ASPECT_SAFE "
                                    "variable. Please check your build "
                                    "configuration."));
#endif
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(CPP,
                                   "cpp",
                                   "A material model whose evaluate function is dynamically "
                                   "compiled from code defined in section 'Material model/C "
                                   "plus plus' of the parameter file.")
  }
}
