/*
  Copyright (C) 2018 by the authors of the World Builder code.

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

#ifndef _world_builder_parameters_h
#define _world_builder_parameters_h

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>


#include <rapidjson/document.h>
#include "rapidjson/schema.h"

#include <world_builder/point.h>

namespace WorldBuilder
{
  namespace Types
  {
    class Interface;
    template<int dim>
    class Point;
    class Double;
    class String;
    class Segment;
    class Array;
    class Bool;
    class UnsignedInt;
  }

  namespace Features
  {
    class Interface;
  }

  namespace CoordinateSystems
  {
    class Interface;
  }

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
       * \param filename A string with the path to the world builder file
       * \param world A reference to the World class
       */
      Parameters(World &world);

      /**
       * Destructor
       */
      ~Parameters();

      /**
       * Todo
       */
      void initialize(std::string &filename, bool has_output_dir = false, std::string output_dir = "");

      /**
       * Todo
       */
      template<class T>
      T get(const std::string &name);

      /**
       * Todo
       */
      template<class T>
      std::vector<T> get_vector(const std::string &name);

      /**
       * Todo
       */
      template<class T, class A, class B>
      std::vector<T> get_vector(const std::string &name, std::vector<std::shared_ptr<A> > &, std::vector<std::shared_ptr<B> > &);

      /**
       * Todo
       */
      template<class T>
      std::unique_ptr<T> get_unique_pointer(const std::string &name);

      /**
       * Todo
       */
      template<class T>
      bool
      get_unique_pointers(const std::string &name, std::vector<std::unique_ptr<T> > &);

      /**
       * Todo
       */
      template<class T>
      bool
      get_shared_pointers(const std::string &name, std::vector<std::shared_ptr<T> > &);

      /**
       * Todo
       */
      bool
      check_entry(const std::string &name) const;

      /**
       * Todo
       */
      void declare_entry(const std::string name,
                         const Types::Interface &type,
                         const std::string documentation);


      /**
       * This function is used to enter a subsection. It appends to the path
       * variable. This action is revesed by the leave subsection function.
       * \param name The name of the subsection to be entered.
       * @see path
       * @see leave_subsection()
       */
      void enter_subsection(const std::string name);

      /**
       * This function is used to leave a subsection by removing the last
       * element of the path variable. It reverses the action of the enter
       * subsection function.
       * @see path
       * @see enter_subsection()
       */
      void leave_subsection();


      /**
       * A reference to the World class. This is needed to create the features.
       */
      World &world;

      /**
       * This variable stores what path separtor is used in the property tree
       * and in this class.
       */
      const std::string path_seperator = ".";

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
       * This function return the current path as stored in the path variable
       * as a string in json pointer format.
       * \return std::string
       */
      std::string get_full_json_path(unsigned int max_size = std::numeric_limits<unsigned int>::max()) const;

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
      unsigned int path_level;

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
}
#endif
