/*
  Copyright (C) 2018-2024 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef WORLD_BUILDER_PARAMETERS_H
#define WORLD_BUILDER_PARAMETERS_H

#include <map>
#include <memory>
#include <vector>

#include "rapidjson/schema.h"
#include "world_builder/point.h"
#include "world_builder/types/unsigned_int.h"

namespace WorldBuilder
{
  namespace Types
  {
    class Interface;
    template<unsigned int dim>
    class Point;
    class Double;
    class String;
    class Segment;
    class Array;
    class Bool;
    class UnsignedInt;
    class Int;
  } // namespace Types

  namespace Features
  {
    class Interface;
  } // namespace Features

  namespace CoordinateSystems
  {
    class Interface;
  } // namespace CoordinateSystems

  namespace GravityModel
  {
    class Interface;
  } // namespace GravityModel

  class World;

  /**
   * A class to hold all the parameters needed by the world builder. Internally
   * it holds all values in the form of vectors of class Types. Values can be
   * entered in two ways into this class. The first way is through the
   * load_entry function which load the value from the provided world builder
   * file. The second way is through the set_entry function, through which
   * values can be directly entered into the parameter class. Values can be
   * retrieved through the get functions which take the name of the value
   * with which it was set. It is also required for bot loading, setting and
   * getting values to do it in the correct subsection. Subsections can be
   * entered with the enter_subsection function and left with the
   * leave_subsection function. The current path can be retrieved through the
   * function get_current_path() and get_current_path_without_arrays().
   */
  class Parameters
  {
    public:
      /**
       * Constructor
       * \param world A reference to the World class
       */
      Parameters(World &world);

      /**
       * Destructor
       */
      ~Parameters();

      /**
       * Initializes the parameter file
       * \param filename A string with the path to the world builder file
       * \param has_output_dir A bool indicating whether the world builder may write out information.
       * \param output_dir A string with the path to the directory where it can output information if allowed by has_output_dir
       */
      void initialize(std::string &filename, bool has_output_dir = false, const std::string &output_dir = "");

      /**
       * A generic get function to retrieve setting from the parameter file.
       * Note that this is dependent on the current path/subsection which you are in.
       * \param name The name of the entry to retrieved
       * @see path
       * @see enter_subsection()
       * @see leave_subsection()
       */
      template<class T>
      T get(const std::string &name);

      /**
       * A specialized version of get which can return vectors/arrays.
       * \param name The name of the entry to retrieved
       */
      template<class T>
      std::vector<T> get_vector(const std::string &name);

      std::vector<std::vector<double>> get_vector_or_double(const std::string &name);

      /**
       * A specialized version of get which can return a value at points type.
       * \param name The name of the entry to retrieved
       * \param name additional points to be added to the list at either the default value or at the value of a single value array in the list
       */
      std::pair<std::vector<double>,std::vector<double>>
                                                      get(const std::string &name,
                                                          const std::vector<Point<2> > &addition_points = {});

      /**
       * A specialized version of get which can return a values at times type.
       * \param name The name of the entry to retrieved
       */
      std::pair<std::vector<double>,std::vector<double>> get_value_at_array(const std::string &name);

      /**
       * A specialized version of get which can return vectors/arrays.
       * This version is designed for the plugin system.
       * \param name The name of the entry to retrieved
       */
      template<class T, class A, class B, class C>
      std::vector<T> get_vector(const std::string &name, std::vector<std::shared_ptr<A> > &, std::vector<std::shared_ptr<B> > &, std::vector<std::shared_ptr<C> > &);

      /**
       * A specialized version of get which can return unique pointers.
       * \param name The name of the entry to retrieved
       */
      template<class T>
      std::unique_ptr<T> get_unique_pointer(const std::string &name);

      /**
       * A specialized version of get which can return unique pointers as an argument
       * and returns a bool to indicate whether it was successful or not.
       * Note that this function will erase all information in the vector.
       * \param name The name of the entry to retrieved
       * \param vector A vector of unique pointers.
       */
      template<class T>
      bool
      get_unique_pointers(const std::string &name, std::vector<std::unique_ptr<T> > &vector);

      /**
       * A specialized version of get which can return shared pointers as an argument
       * and returns a bool to indicate whether it was successful or not.
       * Note that this function will erase all information in the vector.
       * \param name The name of the entry to retrieved
       * \param vector A vector of shared pointers.
       */
      template<class T>
      bool
      get_shared_pointers(const std::string &name, std::vector<std::shared_ptr<T> > & /*vector*/);

      /**
       * Checks for the existence of an entry in the parameter file.
       * Return true when an entry is specified and false when it is not.
       * This is independent of whether an entry has been declared or not.
       * The main intended usage is to check whether the user has provided
       * the specified entry in the user supplied parameters file, since
       * the get functions may use default values.
       * \param name The name of the entry to be checked.
       */
      bool
      check_entry(const std::string &name) const;

      /**
       * Declares the existence an entry in the parameters class.
       * Default values are supplied by the type.
       * \param name The name of the entry to be declared
       * \param type The type of entry (e.g. Double, Array, etc.)
       * \param documentation A string containing information about this parameter.
       */
      void declare_entry(const std::string &name,
                         const Types::Interface &type,
                         const std::string &documentation);


      /**
       * This function is used to enter a subsection. It appends to the path
       * variable. This action is revesed by the leave subsection function.
       * \param name The name of the subsection to be entered.
       * @see path
       * @see leave_subsection()
       */
      void enter_subsection(const std::string &name);

      /**
       * This function is used to leave a subsection by removing the last
       * element of the path variable. It reverses the action of the enter
       * subsection function.
       * @see path
       * @see enter_subsection()
       */
      void leave_subsection();

      /**
       * A utilities function for declaring plugin model entries. This always contains a model declaration entry with the plugin name.
       * @param model_group_name The name of the model group which is declared.
       * @param parent_name The name of the parent declaration group.
       * @param declaration_map A map containing plugin names and plugin declaration functions
       * @param required_entries A vector containing what entries should be required from the user. Default value is empty.
       * @param extra_declarations A vector containing extra declarations common to all plugins in this group. Default value is empty.
       */
      void
      declare_model_entries(const std::string &model_group_name,
                            const std::string &parent_name,
                            const std::map<std::string, void ( *)(Parameters &,const std::string &)> &declare_map,
                            const std::vector<std::string> &required_entries = {},
                            const std::vector<std::tuple<std::string,const WorldBuilder::Types::Interface &, std::string> > &extra_declarations = {});


      /**
       * A reference to the World class. This is needed to create the features.
       */
      World &world;

      /**
       * This variable stores what path separator is used in the property tree
       * and in this class.
       */
      const std::string path_separator = ".";

      /**
       * This variable stores the path in a vector of strings.
       * @see enter_subsection()
       * @see leave_subsection()
       */
      std::vector<std::string> path;

      rapidjson::Document declarations;
      rapidjson::Document parameters;




      /**
       * A vector containing all the pointers to the features. This vector is
       * responsible for the features and has ownership over them. Therefore
       * unique pointers are used.
       * @see Features
       */
      std::vector<std::unique_ptr<WorldBuilder::Features::Interface> > features;

      /**
       * A pointers to the corodinate system. This variable is responsible for
       * the coordinate system and has ownership over it. Therefore a unique
       * pointer are used.
       * @see CoordinateSystem
       */
      std::unique_ptr<WorldBuilder::CoordinateSystems::Interface> coordinate_system;

      /**
       * A pointers to the gravity model. This variable is responsible for
       * the gravity model and has ownership over it. Therefore a unique
       * pointer are used.
       * @see CoordinateSystem
       */
      std::unique_ptr<WorldBuilder::GravityModel::Interface> gravity_model;

      /**
       * This function return the current path as stored in the path variable
       * as a string in json pointer format.
       * \return std::string
       */
      std::string get_full_json_path(size_t max_size = std::numeric_limits<size_t>::max()) const;

      /**
       * todo: Warning: do not use before declarations is filled.
       * This function return the current path as stored in the path variable
       * as a string in json pointer format.
       * \return std::string
       */
      std::string get_full_json_schema_path() const;

      /**
       * This function return the current path as stored in the path variable
       * as a string.
       * \return std::string
       */
      //std::string get_full_path() const;

      /**
       * This function return the current path as stored in the path variable
       * as a string, but the arrays are striped. This is useful for working
       * with the boost property tree.
       * \return std::string
       */
      //std::string get_full_path_without_arrays() const;

    private:


      /**
       * This is used for the get relative path functions. It stores how many
       * top entries of the path should be ignored.
       */
      size_t path_level;

      /**
       * A function which returns the relative path, which is the full path
       * minus the path_level top entries..
       * @see get_current_path()
       */
      std::string get_relative_path() const;

      /**
       * A function which returns the relative path, which is the full path
       * minus the path_level top entries., without names for the arrays.
       * @see get_current_path_without_arrays()
       */
      std::string get_relative_path_without_arrays() const;
  };
} // namespace WorldBuilder
#endif
