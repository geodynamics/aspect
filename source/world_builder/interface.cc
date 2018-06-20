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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
 */


#include <aspect/global.h>
#include <aspect/world_builder/interface.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx11/tuple.h>
#include <deal.II/base/utilities.h>
#include <list>
#include <sstream>
#include <regex>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>


namespace aspect
{
  namespace WorldBuilder
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    void
    Interface<dim>::initialize ()
    {}

    template <int dim>
    void
    Interface<dim>::
    declare_parameters (WorldBuilderParameterHandler &)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &)
    {}




    // -------------------------------- Deal with registering world generator models and automating
    // -------------------------------- their setup and selection at run time

    namespace
    {
      std_cxx11::tuple
      <void *,
      void *,
      internal::Plugins::PluginList<Interface<2> >,
      internal::Plugins::PluginList<Interface<3> > > registered_plugins;
    }



    template <int dim>
    void
    Manager<dim>::register_world_builder (const std::string &name,
                                          const std::string &description,
                                          void (*declare_parameters_function) (ParameterHandler &),
                                          Interface<dim> *(*factory_function) ())
    {
      std_cxx11::get<dim>(registered_plugins).register_plugin (name,
                                                               description,
                                                               declare_parameters_function,
                                                               factory_function);
    }

    namespace
    {
      std::string
      mangle (const std::string &s)
      {
        std::string u;

        // reserve the minimum number of characters we will need. it may
        // be more but this is the least we can do
        u.reserve (s.size());

        // see if the name is special and if so mangle the whole thing
        const bool mangle_whole_string = (s == "value");

        // for all parts of the string, see
        // if it is an allowed character or
        // not
        for (unsigned int i=0; i<s.size(); ++i)
          {
            static const std::string allowed_characters
            ("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");

            if ((! mangle_whole_string)
                &&
                (allowed_characters.find (s[i]) != std::string::npos))
              u.push_back (s[i]);
            else
              {
                u.push_back ('_');
                static const char hex[16]
                  = { '0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};
                u.push_back (hex[static_cast<unsigned char>(s[i])/16]);
                u.push_back (hex[static_cast<unsigned char>(s[i])%16]);
              }
          }

        return u;
      }

      std::string
      demangle(const std::string &s)
      {
        std::string u;
        u.reserve(s.size());

        for (unsigned int i = 0; i < s.size(); ++i)
          if (s[i] != '_')
            u.push_back(s[i]);
          else
            {
              Assert(i + 2 < s.size(),
                     ExcMessage("Trying to demangle an invalid string."));

              unsigned char c = 0;
              switch (s[i + 1])
                {
                  case '0':
                    c = 0 * 16;
                    break;
                  case '1':
                    c = 1 * 16;
                    break;
                  case '2':
                    c = 2 * 16;
                    break;
                  case '3':
                    c = 3 * 16;
                    break;
                  case '4':
                    c = 4 * 16;
                    break;
                  case '5':
                    c = 5 * 16;
                    break;
                  case '6':
                    c = 6 * 16;
                    break;
                  case '7':
                    c = 7 * 16;
                    break;
                  case '8':
                    c = 8 * 16;
                    break;
                  case '9':
                    c = 9 * 16;
                    break;
                  case 'a':
                    c = 10 * 16;
                    break;
                  case 'b':
                    c = 11 * 16;
                    break;
                  case 'c':
                    c = 12 * 16;
                    break;
                  case 'd':
                    c = 13 * 16;
                    break;
                  case 'e':
                    c = 14 * 16;
                    break;
                  case 'f':
                    c = 15 * 16;
                    break;
                  default:
                    Assert(false, ExcInternalError());
                }
              switch (s[i + 2])
                {
                  case '0':
                    c += 0;
                    break;
                  case '1':
                    c += 1;
                    break;
                  case '2':
                    c += 2;
                    break;
                  case '3':
                    c += 3;
                    break;
                  case '4':
                    c += 4;
                    break;
                  case '5':
                    c += 5;
                    break;
                  case '6':
                    c += 6;
                    break;
                  case '7':
                    c += 7;
                    break;
                  case '8':
                    c += 8;
                    break;
                  case '9':
                    c += 9;
                    break;
                  case 'a':
                    c += 10;
                    break;
                  case 'b':
                    c += 11;
                    break;
                  case 'c':
                    c += 12;
                    break;
                  case 'd':
                    c += 13;
                    break;
                  case 'e':
                    c += 14;
                    break;
                  case 'f':
                    c += 15;
                    break;
                  default:
                    Assert(false, ExcInternalError());
                }

              u.push_back(static_cast<char>(c));

              // skip the two characters
              i += 2;
            }

        return u;
      }

      std::string indent(int level)
      {
        std::string s;
        for (int i = 0; i<level; i++) s += "  ";
        return s;
      }

      void printTree(boost::property_tree::ptree &pt, int level)
      {
        if (pt.empty())
          {
            std::cout << "\"" << pt.data() << "\"";
          }

        else
          {
            if (level) std::cout << std::endl;

            std::cout << indent(level) << "{" << std::endl;

            for (boost::property_tree::ptree::iterator pos = pt.begin(); pos != pt.end();)
              {
                std::cout << indent(level + 1) << "\"" << pos->first << "\": ";

                printTree(pos->second, level + 1);
                ++pos;
                if (pos != pt.end())
                  {
                    std::cout << ",";
                  }
                std::cout << std::endl;
              }

            std::cout << indent(level) << " }";
          }
        std::cout << std::endl;
        return;
      }

      /**
       * This function loops through the source tree which should contain all the
       * declared variables. It turns the empty subsections from source, which represent
       * arrays, into subsections called 'array[0]' where the zero is the number of the
       * arry entry and adds it to destination. It furthermore adds a section size in that
       * which contains the size of the array and it's pattern.
       */
      void
      make_arrays_recursively(const boost::property_tree::ptree &source,
                              const std::string           &current_path,
                              const char                  path_separator,
                              const unsigned int          indent,
                              std::vector<std::unique_ptr<const Patterns::PatternBase> > &patterns,
                              boost::property_tree::ptree       &destination)
      {

        typedef typename boost::property_tree::ptree::key_type::value_type Ch;
        typedef typename std::basic_string<Ch> Str;

        // Value or object or array
        if (indent > 0 && source.empty())
          {
            // This is a value, write it into destination
            destination.add(current_path + path_separator + "value",source.template get_value<Str>());

          }
        else if (indent > 0 && source.count(Str()) == source.size())
          {
            // This is an array
            // Add a subsection called size which contains the size of the array
            // and the pattern interger.
            unsigned int array_number = 0;
            for (boost::property_tree::ptree::const_iterator it = source.begin(); it != source.end(); ++it)
              {
                make_arrays_recursively(it->second, current_path + path_separator + mangle("array[") + std::to_string(array_number) + mangle("]"), path_separator, indent + 1, patterns, destination);
                array_number++;
              }
            destination.add(current_path + path_separator + "size" + path_separator + "value",std::to_string(array_number));
            patterns.reserve (patterns.size() + 1);
            patterns.emplace_back (Patterns::Integer().clone());
            destination.add (current_path + path_separator + "size" + path_separator + "pattern",
                             static_cast<unsigned int>(patterns.size()-1));

          }
        else
          {
            // This is an object
            typename boost::property_tree::ptree::const_iterator it = source.begin();
            for (; it != source.end(); ++it)
              {
                make_arrays_recursively(it->second,
                                        (current_path == "" ?
                                         mangle(it->first) :
                                         current_path + path_separator + mangle(it->first)),
                                        path_separator,
                                        indent + 1,
                                        patterns,
                                        destination);
              }
          }
      }

      /**
       * This function loops through the source and default arrays and adds the
       * source values to the destination tree when available, otherwise it uses
       * the default values.
       */
      void
      set_default_arrays_recursively(const boost::property_tree::ptree &source,
                                     const boost::property_tree::ptree &defaults,
                                     const std::string           &current_path,
                                     const std::string           &current_defaults_path,
                                     const char                  path_separator,
                                     const unsigned int          indent,
                                     boost::property_tree::ptree       &destination)
      {
        typedef typename boost::property_tree::ptree::key_type::value_type Ch;
        typedef typename std::basic_string<Ch> Str;

        // Value or object or array
        if (indent > 0 && defaults.empty())
          {
            // This is a value
            destination.add(current_path,defaults.template get_value<Str>(current_defaults_path));

          }
        else if (indent > 0
                 && defaults.get_optional<std::string>("array")
                 && source.get_optional<unsigned int>(current_path + path_separator + "size" + path_separator + "value")
                )
          {
            // This is an array
            // Add a subsection called size which contains the size of the array
            // with the source values.
            const std::string full_path_size = current_path + path_separator + "size";
            unsigned int array_size = source.get<unsigned int>(full_path_size + path_separator + "value");
            unsigned int array_pattern = source.get<unsigned int>(full_path_size + path_separator + "pattern");

            destination.add(full_path_size + path_separator + "value",array_size);
            destination.add(full_path_size + path_separator + "pattern",array_pattern);

            // Now continue for the array itself
            typename boost::property_tree::ptree::const_iterator it = defaults.begin();
            for (unsigned int array_number = 0; array_number < array_size; array_number++)
              {
                set_default_arrays_recursively(source,
                                               it->second,
                                               current_path + path_separator + mangle("array[") + std::to_string(array_number)+mangle("]"),
                                               current_path + path_separator + "array",
                                               path_separator,
                                               indent + 1,
                                               destination);
              }
          }
        else
          {
            // This is an object
            typename boost::property_tree::ptree::const_iterator it = defaults.begin();
            for (; it != defaults.end(); ++it)
              {
                std::string new_path = (current_path == "" ?
                                        it->first :
                                        current_path + path_separator + it->first);

                set_default_arrays_recursively(source,
                                               it->second,
                                               new_path,
                                               new_path,
                                               path_separator,
                                               indent + 1,
                                               destination);
              }
          }
      }


      /**
       * This function read a processed source file and if the values are set
       * this function overwrites the value in the destination tree. Furthermore,
       * this function checks the if the right patterns is used for the given value.
       * This function therefore expects destination to contain all the default values.
       */
      void
      read_recursively(const boost::property_tree::ptree &source,
                       const std::string                 &current_path,
                       const char                         path_separator,
                       const std::vector<std::unique_ptr<const Patterns::PatternBase> > &patterns,
                       boost::property_tree::ptree       &destination)
      {
        for (boost::property_tree::ptree::const_iterator p = source.begin();
             p != source.end(); ++p)
          {
            // a sub-tree must either be a parameter node, array, comment or a subsection
            if (p->second.get_optional<std::string>("value") && boost::algorithm::to_lower_copy(p->first) != "comment")
              {
                // make sure we have a corresponding entry in the destination
                // object as well
                const std::string full_path
                  = (current_path == ""
                     ?
                     p->first
                     :
                     current_path + path_separator + p->first);

                const std::string new_value
                  = p->second.get<std::string>("value");

                AssertThrow (destination.get_optional<std::string> (full_path)
                             &&
                             destination.get_optional<std::string>
                             (full_path + path_separator + "value"),
                             ParameterHandler::ExcEntryUndeclared (demangle(full_path)));

                // first make sure that the new entry actually satisfies its
                // constraints
                const unsigned int pattern_index
                  = destination.get<unsigned int> (full_path +
                                                   path_separator +
                                                   "pattern");

                AssertThrow (patterns[pattern_index]->match(new_value),
                             ParameterHandler::ExcInvalidEntryForPatternXML
                             // XML entries sometimes have extra surrounding
                             // newlines
                             (Utilities::trim(new_value),
                              p->first,
                              patterns[pattern_index]->description()));

                // set the found parameter in the destination argument
                destination.put (full_path + path_separator + "value",
                                 new_value);

                // this node might have sub-nodes in addition to "value", such as
                // "default_value", "documentation", etc. we might at some point
                // in the future want to make sure that if they exist that they
                // match the ones in the 'destination' tree
              }
            else if (p->second.get_optional<std::string>("alias"))
              {
                // it is an alias node. alias nodes are static and there is
                // nothing to do here (but the same applies as mentioned in the
                // comment above about the static nodes inside parameter nodes
              }
            else
              {
                // it must be a subsection
                read_recursively (p->second,
                                  (current_path == "" ?
                                   p->first :
                                   current_path + path_separator + p->first),
                                  path_separator,
                                  patterns,
                                  destination);
              }
          }
      }
    }


    WorldBuilderParameterHandler::WorldBuilderParameterHandler(const char path_separator_)
      : path_separator(path_separator_)
    { }

    void
    WorldBuilderParameterHandler::declare_entry (const std::string           &entry,
                                                 const std::string           &default_value,
                                                 const Patterns::PatternBase &pattern,
                                                 const std::string           &documentation)
    {
      std::string full_path = get_current_full_path(entry);
      wb_tree.put (full_path + path_separator + "value",
                   default_value);
      wb_tree.put (full_path + path_separator + "default_value",
                   default_value);
      wb_tree.put (full_path + path_separator + "documentation",
                   documentation);

      patterns.reserve (patterns.size() + 1);
      patterns.emplace_back (pattern.clone());
      wb_tree.put (full_path + path_separator + "pattern",
                   static_cast<unsigned int>(patterns.size()-1));
      // also store the description of
      // the pattern. we do so because we
      // may wish to export the whole
      // thing as XML or any other format
      // so that external tools can work
      // on the parameter file; in that
      // case, they will have to be able
      // to re-create the patterns as far
      // as possible
      wb_tree.put (full_path + path_separator +
                   "pattern_description",
                   patterns.back()->description());

      // as documented, do the default value checking at the very end
      AssertThrow (pattern.match (default_value),
                   ExcValueDoesNotMatchPattern (default_value, pattern.description()));
    }

    void
    WorldBuilderParameterHandler::set_tree(boost::property_tree::ptree &tree)
    {
      wb_tree = tree;
      patterns.clear();
    }

    boost::property_tree::ptree &
    WorldBuilderParameterHandler::get_tree()
    {
      return wb_tree;
    }

    void
    WorldBuilderParameterHandler::enter_subsection (const std::string &subsection)
    {
      subsection_path.push_back (subsection);
    }


    void
    WorldBuilderParameterHandler::leave_subsection ()
    {
      // assert there is a subsection that
      // we may leave
      Assert (subsection_path.size() != 0, ExcMessage("Already at top level"));

      if (subsection_path.size() > 0)
        subsection_path.pop_back ();
    }

    std::string
    WorldBuilderParameterHandler::get(const std::string &entry_string) const
    {
      // assert that the entry is indeed declared
      if (boost::optional<std::string> value
          = wb_tree.get_optional<std::string> ((get_current_path() == "" ? "" : get_current_path() + path_separator) + mangle(entry_string) + path_separator + "value"))
        {
          return value.get();
        }
      else
        {
          Assert (false, ExcMessage("Entry undeclared: " + entry_string + ", path = " + (get_current_path() == "" ? "" : get_current_path() + path_separator)  + mangle(entry_string) + path_separator + "value"));
          return "";
        }
    }


    std::vector<std::string>
    WorldBuilderParameterHandler::get_array(const std::string &entry_string)
    {
      std::vector<std::string> tmp;

      const unsigned int number_of_objects = wb_tree.get<unsigned int>(get_current_path() + path_separator + "size");

      tmp.resize(number_of_objects);
      for (unsigned int array_number = 0; array_number < number_of_objects; array_number++)
        {
          if (boost::optional<std::string> value = wb_tree.get<std::string>(get_current_path() + path_separator + mangle("array[") + std::to_string(array_number)+mangle("].value")))
            {
              Assert(value.get() != "", ExcMessage("value is empty!"));
              tmp[array_number] = value.get();
            }
          else
            {
              Assert (false, ExcMessage("Entry undeclared: " + entry_string + ", path = " + (get_current_path() == "" ? "" : get_current_path() + path_separator)  + mangle(entry_string) + path_separator + "value"));
              return tmp;
            }
        }
      return tmp;
    }

    std::vector<std::vector<std::string> >
    WorldBuilderParameterHandler::get_double_array(const std::string &entry_string)
    {
      std::vector<std::vector<std::string> > tmp;
      const unsigned int number_of_objects_i = wb_tree.get<unsigned int>(get_current_path()
                                                                         + path_separator + entry_string
                                                                         + path_separator + "size"
                                                                         + path_separator + "value");
      tmp.resize(number_of_objects_i);
      for (unsigned int counter_i = 0; counter_i < number_of_objects_i; counter_i++)
        {
          const unsigned int number_of_objects_j = wb_tree.get<unsigned int>(get_current_path()
                                                                             + path_separator + entry_string
                                                                             + path_separator + mangle("array[") + std::to_string(counter_i)+mangle("]")
                                                                             + path_separator + "size"
                                                                             + path_separator + "value");

          tmp[counter_i].resize(number_of_objects_j);
          for (unsigned int counter_j = 0; counter_j < number_of_objects_j; counter_j++)
            {
              if (boost::optional<std::string> value = wb_tree.get<std::string>(get_current_path()
                                                                                + path_separator + entry_string
                                                                                + path_separator + mangle("array[") + std::to_string(counter_i) + mangle("]")
                                                                                + path_separator + mangle("array[") + std::to_string(counter_j)+mangle("]")
                                                                                + path_separator + "value"))
                {
                  Assert(value.get() != "", ExcMessage("value is empty in " + get_current_path()
                                                       + path_separator + entry_string
                                                       + path_separator + "array[" + std::to_string(counter_i)+"]"
                                                       + path_separator + mangle("array[") + std::to_string(counter_j)+mangle("]")
                                                       + path_separator + "value"));
                  tmp[counter_i][counter_j] = value.get();
                }
              else
                {
                  Assert (false, ExcMessage("Entry undeclared: " + entry_string + ", path = " + (get_current_path() == "" ? "" : get_current_path() + path_separator)  + mangle(entry_string) + path_separator + "value"));
                  return tmp;
                }
            }
        }
      return tmp;
    }

    /*std::vector<std::vector<std::string> >
    WorldBuilderParameterHandler::get_array(const std::string &entry_string)
    {
      std::vector<std::vector<double> > tmp(s.size());
      for (unsigned int i = 0; i < s.size(); ++i)
        tmp[i] = get_array(s[i]);
      return tmp;
    }*/

    std::string
    WorldBuilderParameterHandler::get_current_path () const
    {
      if (subsection_path.size() > 0)
        {
          std::string p = mangle(subsection_path[0]);
          for (unsigned int i=1; i<subsection_path.size(); ++i)
            {
              p += path_separator;
              p += mangle(subsection_path[i]);
            }
          return p;
        }
      else
        return "";
    }

    std::string
    WorldBuilderParameterHandler::get_current_full_path(const std::string &name) const
    {
      std::string path = get_current_path();
      if (path.empty() == false)
        path += path_separator;

      path += mangle(name);

      return path;
    }

    template <int dim>
    WorldBuilderParameterHandler Manager<dim>::ph;

    template <int dim>
    Manager<dim>::Manager ()
    {}

    template <int dim>
    void
    Manager<dim>::initialize ()
    {

    }

    template <>
    double
    Manager<2>::initial_temperature(const dealii::Point<2> &/*original_position*/) const
    {
      AssertThrow(false,ExcMessage("This has not been implmented in 2d."))
      return 0;
    }

    template <int dim>
    double
    Manager<dim>::initial_temperature(const dealii::Point<dim> &original_position) const
    {
      const GeometryModel::Interface<dim> *geometry_model= &this->get_geometry_model();

      Point<dim> external_position;
      Point<3> internal_position;
      double temperature;
      if (dim==2)
        {
          //Todo, convert from a 2d point of a cross section to a 3d point for the world generator
          return 2.0;
        }
      else if (dim==3)
        {
          external_position = Utilities::convert_array_to_point<dim>(geometry_model->cartesian_to_natural_coordinates(original_position));//dynamic_cast<const GeometryModel::EllipsoidalChunk<dim>&>(this->get_geometry_model()).get_manifold().pull_back(position); // longitude, lattitude, height


          internal_position(0) = external_position(0);
          internal_position(1) = external_position(1);
          internal_position(2) = external_position(2);

          Point<3> position_rotated(cos(rotation_angle) * (internal_position[0] - rotation_point[0]) - sin(rotation_angle) * (internal_position[1] - rotation_point[1]) + rotation_point[0],
                                    sin(rotation_angle) * (internal_position[0] - rotation_point[0]) + cos(rotation_angle) * (internal_position[1] - rotation_point[1]) + rotation_point[1],
                                    internal_position[2]);
          Point<3> position = geometry_model->natural_to_cartesian_coordinates(Utilities::convert_point_to_array<dim>(reinterpret_cast<dealii::Point<3>&>(position_rotated)));
          Utilities::NaturalCoordinate<3> natural_position = typename Utilities::NaturalCoordinate<3>::NaturalCoordinate(position,this->get_geometry_model());
          const double gravity_norm=dynamic_cast<const GravityModel::Interface<dim>&>(this->get_gravity_model()).gravity_vector(position).norm();

          const double depth = geometry_model->depth(position);

          // mantle background temperature (adiabatic)
          const double mantle_temperature = potential_mantle_temperature + (((potential_mantle_temperature * thermal_expansion_coefficient_alfa * gravity_norm) / specific_heat_Cp) * 1000.0) * ((depth) / 1000.0);

          temperature = mantle_temperature;

          if (std::fabs(depth) < 1e-8)
            {
              /*
               * We force the top layer to have a uniform
               * surface temperature.
               */
              temperature = surface_temperature;
            }
          else
            {
              // Load objects in the right order
              for (unsigned int i_object = 0; i_object < ordered_modules.size(); ++i_object)
                {
                  temperature = ordered_modules[i_object]->temperature(natural_position,i_object,depth,temperature);
                }
            }
        }
      return temperature;
    }

    template <>
    double
    Manager<2>::initial_composition(const dealii::Point<2> &/*original_position*/, const unsigned int /*n_comp*/) const
    {
      AssertThrow(false,ExcMessage("This has not been implmented in 2d."))
      return 0;
    }

    template <int dim>
    double
    Manager<dim>::initial_composition(const dealii::Point<dim> &original_position, const unsigned int n_comp) const
    {

      const GeometryModel::Interface<dim> *geometry_model= &this->get_geometry_model();

      Point<dim> external_position;
      Point<3> internal_position;

      double composition = 0;
      if (dim==2)
        {
          //Todo, convert from a 2d point of a cross section to a 3d point for the world generator
          return 2.0;
        }
      else if (dim==3)
        {
          external_position = Utilities::convert_array_to_point<dim>(geometry_model->cartesian_to_natural_coordinates(original_position));//dynamic_cast<const GeometryModel::EllipsoidalChunk<dim>&>(this->get_geometry_model()).get_manifold().pull_back(position); // longitude, lattitude, height

          internal_position(0) = external_position(0);
          internal_position(1) = external_position(1);
          internal_position(2) = external_position(2);

          Point<3> position_rotated(cos(rotation_angle) * (internal_position[0] - rotation_point[0]) - sin(rotation_angle) * (internal_position[1] - rotation_point[1]) + rotation_point[0],
                                    sin(rotation_angle) * (internal_position[0] - rotation_point[0]) + cos(rotation_angle) * (internal_position[1] - rotation_point[1]) + rotation_point[1],
                                    internal_position[2]);

          Point<3> position = geometry_model->natural_to_cartesian_coordinates(Utilities::convert_point_to_array<3>(reinterpret_cast<dealii::Point<3>&>(position_rotated)));
          Utilities::NaturalCoordinate<3> natural_position = typename Utilities::NaturalCoordinate<3>::NaturalCoordinate(position,this->get_geometry_model());

          double depth = geometry_model->depth(position);

          // Load objects
          for (unsigned int i_object = 0; i_object < ordered_modules.size(); ++i_object)
            {
              composition = ordered_modules[i_object]->composition(natural_position,n_comp,i_object,depth,composition);
            }
        }
      return composition;
    }

    template <int dim>
    void
    Manager<dim>::declare_parameters (WorldBuilderParameterHandler &ph)
    {
      /**
       * This works as much as possible as the normal parameter handler, with the difference that
       * this parameter handler also supports arrays. To support arrays, the word 'array' and the
       * word 'array[]' with any number in between the square brackets are reserved. If you declare
       * a section or a value here with the name `array` you can ask for it with the name 'array[]'
       * with a number between the square brackets between zero and the value stored in size.
       */
      ph.enter_subsection("Surface rotation point");
      //ph.enter_subsection("array");
      {
        ph.declare_entry("array","2", Patterns::Double(),"TODO");
      }
      ph.leave_subsection();

      ph.declare_entry("Surface rotation angle","3", Patterns::Double(),"TODO");
      ph.declare_entry("Minimum parts per distance unit","3", Patterns::Integer(),"TODO");
      ph.declare_entry("Minimum distance points","1e-5", Patterns::Double(),"TODO");

      ph.declare_entry("Surface temperature", "273.15", Patterns::Double(), "Set the surface temperature constant in degree Kelvin.");
      ph.declare_entry("Potential mantle temperature", "1600.0", Patterns::Double(), "Set the temperature constant in degree Kelvin.");
      ph.declare_entry("Reference specific heat", "1250", Patterns::Double(), "Set the reference specific heat.");
      ph.declare_entry("Thermal expansion coefficient", "3.5e-5", Patterns::Double(), "Set the thermal expansion coefficient.");

      ph.enter_subsection("Surface objects");
      {
        ph.enter_subsection("array");
        {
          ph.declare_entry("name","default_name", Patterns::Anything(),"TODO");
          ph.enter_subsection("coordinates");
          {
            ph.enter_subsection("array");
            {
              ph.declare_entry("array","0", Patterns::Double(),"TODO");
            }
            ph.leave_subsection();
          }
          ph.leave_subsection();
        }
        ph.leave_subsection();
      }
      ph.leave_subsection();
      ParameterHandler prm;
      std_cxx11::get<dim>(registered_plugins).declare_parameters (prm);
    }

    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      // Get direct references to manager parameters
      boost::property_tree::ptree &wb_tree = ph.get_tree();

      // Get world builder file and check wether it exists
      const std::string world_builder_file = prm.get("World builder file");
      AssertThrow(access( world_builder_file.c_str(), F_OK ) != -1,
                  ExcMessage("Could not find the world builder file at the specified location: " + world_builder_file));

      // make sure the three is empty before filling it with the
      // declared parameters.
      wb_tree.clear();
      declare_parameters (ph);


      // The declared entries are stored in wb_tree, put them into the
      // default values tree, so we can start processing it.
      boost::property_tree::ptree default_value_tree = wb_tree;
      wb_tree.clear();

      // Now read in the world builder file into a file stream and
      // put it into a boost property tree.
      std::ifstream json_input_stream(world_builder_file.c_str());
      boost::property_tree::ptree original_json_tree;
      boost::property_tree::json_parser::read_json (json_input_stream, original_json_tree);

      /**
       * This function loops through the source tree which should contain all the
       * declared variables. It turns the empty subsections from source, which represent
       * arrays, into subsections called 'array[0]' where the zero is the number of the
       * arry entry and adds it to destination. It furthermore adds a section size in that
       * which contains the size of the array and it's pattern.
       */
      boost::property_tree::ptree processed_tree;
      make_arrays_recursively (original_json_tree,
                               "",
                               '.',
                               0,
                               ph.patterns,
                               processed_tree);

      /**
       * This function loops through the source and default arrays and adds the
       * source values to the destination tree when available, otherwise it uses
       * the default values.
       */
      set_default_arrays_recursively(processed_tree,
                                     default_value_tree,
                                     "",
                                     "",
                                     '.',
                                     0,
                                     wb_tree);


      /**
       * This function read a processed source file and if the values are set
       * this function overwrites the value in the destination tree. Furthermore,
       * this function checks the if the right patterns is used for the given value.
       * This function therefore expects destination to contain all the default values.
       */
      read_recursively(processed_tree, "", '.',ph.patterns, wb_tree);

      /**
       * This is just for the review to show how to get to values.
       */
      rotation_angle = Utilities::string_to_double(ph.get("Surface rotation angle"));

      ph.enter_subsection("Surface rotation point");
      {
        rotation_point.resize(Utilities::string_to_double(ph.get("size")));
        rotation_point[0] = Utilities::string_to_double(ph.get("array[0]"));
        rotation_point[1] = Utilities::string_to_double(ph.get("array[1]"));
      }
      ph.leave_subsection();

      potential_mantle_temperature = Utilities::string_to_double(ph.get("Potential mantle temperature"));
      surface_temperature = Utilities::string_to_double(ph.get("Surface temperature"));
      specific_heat_Cp = Utilities::string_to_double(ph.get("Reference specific heat"));
      thermal_expansion_coefficient_alfa = Utilities::string_to_double(ph.get("Thermal expansion coefficient"));

      // Find out what modules/plugins are being used.
      std::vector<std::string> module_names;
      ph.enter_subsection("Surface objects");
      {
        for (unsigned int i = 0; i < (unsigned int) Utilities::string_to_int(ph.get("size")); i++)
          {
            ph.enter_subsection("array[" + std::to_string(i) + "]");
            {
              ph.enter_subsection("module");
              {
                const std::string module_name = boost::algorithm::trim_copy(boost::algorithm::to_lower_copy(ph.get("name")));
                AssertThrow (module_name != "" && module_name != "unspecified",
                             ExcMessage ("Please specify the module you want to use (It is not allowed to be empty)."));

                module_names.push_back(module_name);
              }
              ph.leave_subsection();
            }
            ph.leave_subsection();
          }
      }
      ph.leave_subsection();

      std::vector<std::string> tmp_module_names = module_names;
      std::sort(tmp_module_names.begin(), tmp_module_names.end());
      std::vector<std::string>::iterator it;
      it = std::unique(tmp_module_names.begin(), tmp_module_names.end());
      tmp_module_names.resize(distance(tmp_module_names.begin(),it));
      modules.resize(tmp_module_names.size());

      for (unsigned int i=0; i<tmp_module_names.size(); ++i)
        {
          modules[i] = (std_cxx11::shared_ptr<Interface<dim> >(
                          std_cxx11::get<dim>(registered_plugins).create_plugin (module_names[i],
                                                                                 "World builder model::model name")));
        }


      for (auto &module: modules)
        module->parse_parameters (prm);

      for (unsigned int i=0; i<module_names.size(); ++i)
        {
          // find where that module is stored
          for (unsigned int j=0; j<tmp_module_names.size(); j++)
            {
              if (module_names[i] == tmp_module_names[j])
                {
                  ordered_modules.push_back(modules[j]);
                }
            }
        }

    }


    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.declare_entry ("World builder file", "",
                         Patterns::FileName (),
                         "The location of the world builder JSON file."
                         "For modules, select one of the following:\n\n"
                         +
                         std_cxx11::get<dim>(registered_plugins).get_description_string());


    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<WorldBuilder::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<WorldBuilder::Interface<2> >::plugins = 0;

      template <>
      std::list<internal::Plugins::PluginList<WorldBuilder::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<WorldBuilder::Interface<3> >::plugins = 0;
    }
  }

  namespace WorldBuilder
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>;

    ASPECT_REGISTER_WORLD_BUILDER_AS_INITIAL_TEMPERATURE_MODEL(Manager,
                                                               "world builder",
                                                               "Temperature is prescribed as an adiabatic "
                                                               "profile with upper and lower thermal boundary layers, "
                                                               "whose ages are given as input parameters.")


    ASPECT_REGISTER_WORLD_BUILDER_AS_INITIAL_COMPOSITION_MODEL(Manager,
                                                               "world builder",
                                                               "Temperature is prescribed as an adiabatic "
                                                               "profile with upper and lower thermal boundary layers, "
                                                               "whose ages are given as input parameters.")

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
