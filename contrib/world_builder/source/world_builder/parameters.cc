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


#include "world_builder/features/continental_plate_models/composition/interface.h"
#include "world_builder/features/continental_plate_models/grains/interface.h"
#include "world_builder/features/continental_plate_models/temperature/interface.h"
#include "world_builder/features/fault.h"
#include "world_builder/features/mantle_layer_models/composition/interface.h"
#include "world_builder/features/mantle_layer_models/grains/interface.h"
#include "world_builder/features/mantle_layer_models/temperature/interface.h"
#include "world_builder/features/oceanic_plate_models/composition/interface.h"
#include "world_builder/features/oceanic_plate_models/grains/interface.h"
#include "world_builder/features/oceanic_plate_models/temperature/interface.h"
#include "world_builder/features/plume_models/composition/interface.h"
#include "world_builder/features/plume_models/grains/interface.h"
#include "world_builder/features/plume_models/temperature/interface.h"
#include "world_builder/features/subducting_plate.h"
#include "world_builder/gravity_model/interface.h"
#include "world_builder/types/object.h"
#include "world_builder/utilities.h"

#include "rapidjson/error/en.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/latexwriter.h"
#include "rapidjson/mystwriter.h"
#include "rapidjson/prettywriter.h"

#include <fstream>
#include <memory>

using namespace rapidjson;

namespace
{
  void remove_key(rapidjson::Value &value, const char *key)
  {
    if (value.IsObject())
      {
        auto it = value.FindMember(key);
        if (it != value.MemberEnd())
          value.RemoveMember(it);

        for (auto &member : value.GetObject())
          remove_key(member.value, key);

      }
    else if (value.IsArray())
      {
        for (auto &element : value.GetArray())
          remove_key(element, key);
      }
  }
}

namespace WorldBuilder
{
  Parameters::Parameters(World &world_)
    :
    world(world_)
  {
  }

  Parameters::~Parameters()
    = default;

  void Parameters::initialize(std::string &filename, bool has_output_dir, const std::string &output_dir)
  {

    if (has_output_dir)
      {
        StringBuffer buffer;
        std::ofstream file;

        // write out json schema
        file.open (output_dir + "world_builder_declarations.schema.json");
        WBAssertThrow(file.is_open(), "Error: Could not open file '" + output_dir + "world_builder_declarations.schema.json' for string the json declarations.");
        PrettyWriter<StringBuffer, UTF8<>, UTF8<>, CrtAllocator, kWriteNanAndInfFlag> json_writer(buffer);
        declarations.Accept(json_writer);
        file << buffer.GetString();
        file.close();
        buffer.Clear();

        // remove Snippets so they don't appear in the documentation:
        remove_key(declarations, "defaultSnippets");

        // write out declarations
        file.open (output_dir + "world_builder_declarations.tex");
        WBAssertThrow(file.is_open(), "Error: Could not open file '" + output_dir + "world_builder_declarations.tex' for string the tex declarations.");

        LatexWriter<StringBuffer, UTF8<>, UTF8<>, CrtAllocator, kWriteNanAndInfFlag> tex_writer(buffer);
        declarations.Accept(tex_writer);
        file << buffer.GetString();
        file.close();
        buffer.Clear();

        // write out declarations (open)
        file.open (output_dir + "world_builder_declarations_open.md");
        WBAssertThrow(file.is_open(), "Error: Could not open file '" + output_dir + "world_builder_declarations_open.md' for string the tex declarations.");

        MySTWriter<StringBuffer, UTF8<>, UTF8<>, CrtAllocator, kWriteNanAndInfFlag> myst_writer_open(buffer, true);
        declarations.Accept(myst_writer_open);
        file << buffer.GetString();
        file.close();
        buffer.Clear();

        // write out declarations (closed)
        file.open (output_dir + "world_builder_declarations_closed.md");
        WBAssertThrow(file.is_open(), "Error: Could not open file '" + output_dir + "world_builder_declarations_closed.md' for string the tex declarations.");

        MySTWriter<StringBuffer, UTF8<>, UTF8<>, CrtAllocator, kWriteNanAndInfFlag> myst_writer_closed(buffer, false);
        declarations.Accept(myst_writer_closed);
        file << buffer.GetString();
        file.close();
        buffer.Clear();
      }

    path_level =0;
    // Now read in the world builder file into a stringstream and
    // put it into a the rapidjson document
    std::stringstream json_input_stream(WorldBuilder::Utilities::read_and_distribute_file_content(filename));
    rapidjson::IStreamWrapper isw(json_input_stream);

    // relaxing syntax by allowing comments () for now, maybe also allow trailing commas and (kParseTrailingCommasFlag) and nan's, inf etc (kParseNanAndInfFlag)?
    //WBAssertThrow(!parameters.ParseStream<kParseCommentsFlag>(isw).HasParseError(), "Parsing errors world builder file");

    WBAssertThrowExc(!(parameters.ParseStream<kParseCommentsFlag | kParseNanAndInfFlag>(isw).HasParseError()), std::ifstream json_input_stream_error(filename.c_str()); ,
                     "Parsing errors world builder file: Error(offset " << static_cast<unsigned>(parameters.GetErrorOffset())
                     << "): " << GetParseError_En(parameters.GetParseError()) << std::endl << std::endl
                     << " Showing 50 chars before and after: "
                     << std::string((std::istreambuf_iterator<char>(json_input_stream_error.seekg(0, json_input_stream_error.beg))),
                                    std::istreambuf_iterator<char>()).substr(static_cast<unsigned>(parameters.GetErrorOffset()) <= 50
                                                                             ?
                                                                             0
                                                                             :
                                                                             static_cast<unsigned>(parameters.GetErrorOffset()) - 50, 100
                                                                            ) << std::endl << std::endl
                     << " Showing 5 chars before and after: "
                     << std::string((std::istreambuf_iterator<char>(json_input_stream_error.seekg(0, json_input_stream_error.beg))),
                                    std::istreambuf_iterator<char>()).substr(static_cast<unsigned>(parameters.GetErrorOffset()) <= 5
                                                                             ? 0
                                                                             :
                                                                             static_cast<unsigned>(parameters.GetErrorOffset())-5,
                                                                             (static_cast<unsigned>(parameters.GetErrorOffset()) + 10 > json_input_stream_error.seekg(0,std::ios::end).tellg()
                                                                              ?
                                                                              static_cast<unsigned>(json_input_stream.tellg())-static_cast<unsigned>(parameters.GetErrorOffset())
                                                                              :
                                                                              10)
                                                                            ));

    WBAssertThrow(parameters.IsObject(), "World builder file is is not an object.");


    const SchemaDocument schema(declarations);
    SchemaValidator validator(schema);

    if (!parameters.Accept(validator))
      {
        // Input JSON is invalid according to the schema
        // Output diagnostic information
        StringBuffer buffer;
        std::stringstream string;
        validator.GetInvalidSchemaPointer().StringifyUriFragment(buffer);
        string << "Invalid schema: " << buffer.GetString() << std::endl;
        string << "Invalid keyword: " << validator.GetInvalidSchemaKeyword();
        buffer.Clear();
        validator.GetInvalidDocumentPointer().StringifyUriFragment(buffer);
        string << "Invalid schema: " << buffer.GetString() << std::endl;
        PrettyWriter<StringBuffer> writer(buffer);
        validator.GetError().Accept(writer);
        WBAssertThrow(false, string.str() << "Error document: " << std::endl << buffer.GetString());
      }
  }

  void
  Parameters::declare_entry(const std::string &name,
                            const Types::Interface &type,
                            const std::string &documentation)
  {
    type.write_schema(*this,name,documentation);
  }

  bool
  Parameters::check_entry(const std::string &name) const
  {
    return Pointer((this->get_full_json_path() + "/" + name).c_str()).Get(parameters) != nullptr;
  }


  template<>
  std::string
  Parameters::get(const std::string &name)
  {
    const std::string base = this->get_full_json_path();
    const Value *value = Pointer((base + "/" + name).c_str()).Get(parameters);

#ifdef debug
    bool required = false;
    if (Pointer((base + "/required").c_str()).Get(declarations) != NULL)
      {
        for (auto &v : Pointer((base + "/required").c_str()).Get(declarations)->GetArray())
          {
            if (v.GetString() == name)
              {
                required = true;
              }
          }
      }

    WBAssert(value != NULL || required == false,
             "Internal error: Value \"" << base << '/' << name << "/type\" not found in the input file, while it was set as required.");
#endif
    if (value == nullptr)
      {
        value = Pointer((get_full_json_schema_path() + "/" + name + "/default value").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << base + "/" + name + "/default value");
      }

    return value->GetString();
  }

  template<>
  double
  Parameters::get(const std::string &name)
  {
    const std::string base = this->get_full_json_path();
    const Value *value = Pointer((base + "/" + name).c_str()).Get(parameters);

#ifdef debug
    bool required = false;
    if (Pointer((base + "/required").c_str()).Get(declarations) != NULL)
      {
        for (auto &v : Pointer((base + "/required").c_str()).Get(declarations)->GetArray())
          {
            if (v.GetString() == name)
              {
                required = true;
              }
          }
      }

    WBAssert(value != NULL || required == false,
             "Internal error: Value \"" << base << '/' << name << "/type\" not found in the input file, while it was set as required.");
#endif
    if (value == nullptr)
      {
        value = Pointer((get_full_json_schema_path() + "/" + name + "/default value").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << get_full_json_schema_path() + "/" + name + "/default value, for value: " << base + "/" + name);
      }

    double return_value;
    try
      {
        return_value = value->GetDouble();
      }
    catch (...)
      {
        WBAssertThrow(false, "Could not convert values of " << base << " into doubles.");
      }
    return return_value;
  }

  template<>
  size_t
  Parameters::get(const std::string &name)
  {
    const std::string base = this->get_full_json_path();
    const Value *value = Pointer((base + "/" + name).c_str()).Get(parameters);

#ifdef debug
    bool required = false;
    if (Pointer((base + "/required").c_str()).Get(declarations) != NULL)
      {
        for (auto &v : Pointer((base + "/required").c_str()).Get(declarations)->GetArray())
          {
            if (v.GetString() == name)
              {
                required = true;
              }
          }
      }

    WBAssert(value != NULL || required == false,
             "Internal error: Value \"" << base << '/' << name << "/type\" not found in the input file, while it was set as required.");
#endif
    if (value == nullptr)
      {
        value = Pointer((get_full_json_schema_path() + "/" + name + "/default value").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << base + "/" + name + "/default value");
      }

    return value->GetUint();
  }

  template<>
  unsigned int
  Parameters::get(const std::string &name)
  {
    const std::string base = this->get_full_json_path();
    const Value *value = Pointer((base + "/" + name).c_str()).Get(parameters);

#ifdef debug
    bool required = false;
    if (Pointer((base + "/required").c_str()).Get(declarations) != NULL)
      {
        for (auto &v : Pointer((base + "/required").c_str()).Get(declarations)->GetArray())
          {
            if (v.GetString() == name)
              {
                required = true;
              }
          }
      }

    WBAssert(value != NULL || required == false,
             "Internal error: Value \"" << base << '/' << name << "/type\" not found in the input file, while it was set as required.");
#endif
    if (value == nullptr)
      {
        value = Pointer((get_full_json_schema_path() + "/" + name + "/default value").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << base + "/" + name + "/default value");
      }

    return value->GetUint();
  }

  template<>
  int
  Parameters::get(const std::string &name)
  {
    const std::string base = this->get_full_json_path();
    const Value *value = Pointer((base + "/" + name).c_str()).Get(parameters);

#ifdef debug
    bool required = false;
    if (Pointer((base + "/required").c_str()).Get(declarations) != NULL)
      {
        for (auto &v : Pointer((base + "/required").c_str()).Get(declarations)->GetArray())
          {
            if (v.GetString() == name)
              {
                required = true;
              }
          }
      }

    WBAssert(value != NULL || required == false,
             "Internal error: Value \"" << base << '/' << name << "/type\" not found in the input file, while it was set as required.");
#endif
    if (value == nullptr)
      {
        value = Pointer((get_full_json_schema_path() + "/" + name + "/default value").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << base + "/" + name + "/default value");
      }

    return value->GetInt();
  }

  template<>
  bool
  Parameters::get(const std::string &name)
  {
    const std::string base = this->get_full_json_path();
    const Value *value = Pointer((base + "/" + name).c_str()).Get(parameters);

#ifdef debug
    bool required = false;
    if (Pointer((base + "/required").c_str()).Get(declarations) != NULL)
      {
        for (auto &v : Pointer((base + "/required").c_str()).Get(declarations)->GetArray())
          {
            if (v.GetString() == name)
              {
                required = true;
              }
          }
      }

    WBAssert(value != NULL || required == false,
             "Internal error: Value \"" << base << '/' << name << "/type\" not found in the input file, while it was set as required.");
#endif
    if (value == nullptr)
      {
        value = Pointer((get_full_json_schema_path() + "/" + name + "/default value").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << base + "/" + name + "/default value");
      }

    return value->GetBool();
  }


  template<>
  Point<2>
  Parameters::get(const std::string &name)
  {

    const std::string strict_base = this->get_full_json_path();
    const Value *array = Pointer((strict_base + "/" + name).c_str()).Get(parameters);

#ifdef debug
    bool required = false;
    if (Pointer((strict_base + "/required").c_str()).Get(declarations) != NULL)
      {
        for (auto &v : Pointer((strict_base + "/required").c_str()).Get(declarations)->GetArray())
          {
            if (v.GetString() == name)
              {
                required = true;
              }
          }
      }

    WBAssert(array != NULL || required == false,
             "Internal error: Value \"" << strict_base << '/' << name << "/type\" not found in the input file, while it was set as required.");
#endif
    if (array != nullptr)
      {
        const std::string base = strict_base + "/" + name;
        //let's assume that the file is correct, because it has been checked with the json schema.
        // So there are exactly two values.
        double value1;
        double value2;

        try
          {
            value1 = Pointer((base + "/0").c_str()).Get(parameters)->GetDouble();
            value2 = Pointer((base + "/1").c_str()).Get(parameters)->GetDouble();
          }
        catch (...)
          {
            WBAssertThrow(false, "Could not convert values of " << base << " into Point<2>, because it could not convert the sub-elements into doubles.");
          }
        return Point<2>(value1,value2,this->coordinate_system->natural_coordinate_system());
      }
    WBAssertThrow(false, "default values not implemented in get<Point<2> >. Looked in: " + strict_base + "/" << name);

    return {invalid};;
  }

  template<>
  std::vector<bool>
  Parameters::get_vector(const std::string &name)
  {
    std::vector<bool> vector;
    const std::string strict_base = this->get_full_json_path();
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr)
      {
        Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

        for (size_t i = 0; i < array->Size(); ++i )
          {
            const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));

            vector.push_back(Pointer(base.c_str()).Get(parameters)->GetBool());
          }
      }
    else
      {
        const Value *value = Pointer((this->get_full_json_schema_path()  + "/" + name + "/minItems").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the minItems value at: "
                      << this->get_full_json_schema_path() + "/" + name + "/minItems value");

        const size_t min_size = value->GetUint();

        const bool default_value = Pointer((this->get_full_json_schema_path()  + "/" + name + "/items/default value").c_str()).Get(declarations)->GetBool();

        // set to min size
        for (size_t i = 0; i < min_size; ++i)
          {
            vector.push_back(default_value);
          }
      }
    return vector;
  }


  std::pair<std::vector<double>,std::vector<double>>
                                                  Parameters::get(const std::string &name,
                                                                  const std::vector<Point<2> > &addition_points)
  {
    // There are four cases:
    // 1. No value provided: use the default value everywhere. Return first with one value and second with size 0.
    // 2. One double provided: use the default value everywhere. Return first with one value and second with size 0.
    // 3. One value in a double array and no points provided: use that value everywhere. Return first with one value and second with size 0.
    // 4. Other: fill the vectors with the default value and addition points and then add new point.
    //    If a value without points is encountered, the additional points are used.
    std::pair<std::vector<double>,std::vector<double>> result;

    const std::string strict_base = this->get_full_json_path();

    // start with adding the additional points with the default value
    // to do this we need the default value
    double default_value = 0;
    bool is_array = true;
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr && Pointer((strict_base + "/" + name).c_str()).Get(parameters)->IsArray())
      {
        const std::string value_def_path = get_full_json_schema_path() + "/" + name + "/oneOf/1/items/items/anyOf/0/default value";
        Value *value_def = Pointer(value_def_path.c_str()).Get(declarations);
        WBAssertThrow(value_def != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << value_def_path);

        // Since the default value is set in the code, if it fails it is an internal error, not a user error.
        // So no try/catch needed.
        default_value = value_def->GetDouble();
      }
    else
      {
        is_array = false;
        Value *value_def = Pointer((get_full_json_schema_path() + "/" + name + "/oneOf/0/default value").c_str()).Get(declarations);
        WBAssertThrow(value_def != nullptr,
                      "internal error: could not retrieve the default value at: "
                      <<get_full_json_schema_path() + "/" + name + "/oneOf/0/default value");


        // Since the default value is set in the code, if it fails it is an internal error, not a user error.
        // So no try/catch needed.
        default_value = value_def->GetDouble();
      }


    // check if there is a user defined value
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr)
      {
        // there is a user defined value, so either case 2, 3 or 4.
        if (is_array)
          {
            Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

            if (array->Size() == 1
                && Pointer((strict_base + "/" + name + "/0/1").c_str()).Get(parameters) == nullptr)
              {
                // case 2: Return first with one value and second with size 0.
                double value = 0;
                try
                  {
                    value = Pointer((strict_base + "/" + name + "/0/0").c_str()).Get(parameters)->GetDouble();
                  }
                catch (...)
                  {
                    WBAssertThrow(false, "Could not convert values of " << strict_base << "/" << name << "/0/0 into a double. "
                                  << "The provided value was \"" <<  Pointer((strict_base + "/" + name + "/0/0").c_str()).Get(parameters)->GetString() << "\".");
                  }
                result.first.emplace_back(value);
              }
            else
              {
                // case 3: fill the vectors with the default value and addition points and then add new point.
                //         If a value without points is encountered, the additional points are used.

                // first fill with additional points at default value
                for (const auto &addition_point : addition_points)
                  {
                    result.first.emplace_back(default_value);
                    result.second.emplace_back(addition_point[0]);
                    result.second.emplace_back(addition_point[1]);
                  }

                // second, go through all the points in order
                for (size_t i = 0; i < array->Size(); ++i )
                  {
                    // now parse a single value_at_point.
                    const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));
                    // Let's assume that the file is correct, because it has been checked with the json schema.
                    // So there are exactly two values, the first value is a double, the second an array of 2d arrays (points).

                    // Get the double
                    double value;
                    Value *value_pointer = Pointer((base + "/0").c_str()).Get(parameters);

                    WBAssertThrow(value_pointer != nullptr, "internal error: this should not happen.");

                    try
                      {
                        value = value_pointer->GetDouble();
                      }
                    catch (...)
                      {
                        WBAssertThrow(false, "Could not convert values of " << base << "/0 into doubles. "
                                      << "The provided value was \"" <<  Pointer((base + "/0").c_str()).Get(parameters)->GetString() << "\".");
                      }

                    // now get the array of points.
                    Value *coordinates_array = Pointer((base + "/1").c_str()).Get(parameters);
                    if (coordinates_array != nullptr)
                      {
                        for (size_t coordinate_i = 0; coordinate_i < coordinates_array->Size(); ++coordinate_i )
                          {
                            // Let's assume that the file is correct, because it has been checked with the json schema.
                            // That means that there are exactly two values per item
                            double coordinate_0;
                            double coordinate_1;
                            try
                              {
                                coordinate_0 = Pointer((base + "/1/" + std::to_string(coordinate_i) + "/0").c_str()).Get(parameters)->GetDouble();
                              }
                            catch (...)
                              {
                                WBAssertThrow(false, "Could not convert values of " << base + "/1/" + std::to_string(coordinate_i) + "/0"
                                              << " into a Point<2> array, because it could not convert the 1st sub-elements into doubles. "
                                              << "The provided value was \""
                                              <<  Pointer((base + "/1/" + std::to_string(coordinate_i) + "/0").c_str()).Get(parameters)->GetString()
                                              << "\".");
                              }
                            try
                              {
                                coordinate_1 = Pointer((base + "/1/" + std::to_string(coordinate_i) + "/1").c_str()).Get(parameters)->GetDouble();
                              }
                            catch (...)
                              {
                                WBAssertThrow(false, "Could not convert values of " << base + "/1/" + std::to_string(coordinate_i) + "/1"
                                              << " into a Point<2> array, because it could not convert the 2nd sub-elements into doubles. "
                                              << "The provided value was \""
                                              <<  Pointer((base + "/1/" + std::to_string(coordinate_i) + "/1").c_str()).Get(parameters)->GetString()
                                              << "\".");
                              }

                            coordinate_0 *= (coordinate_system->natural_coordinate_system() == CoordinateSystem::spherical ? Consts::PI / 180.0 : 1.);
                            coordinate_1 *= (coordinate_system->natural_coordinate_system() == CoordinateSystem::spherical ? Consts::PI / 180.0 : 1.);

                            bool found_same_point = false;
                            unsigned int coordinate_pair_i = 0;
                            for (; coordinate_pair_i < result.second.size(); coordinate_pair_i+=2)
                              {
                                if (Utilities::approx(result.second[coordinate_pair_i],coordinate_0) && Utilities::approx(result.second[coordinate_pair_i+1],coordinate_1))
                                  {
                                    found_same_point = true;
                                    break;
                                  }
                              }
                            if (found_same_point)
                              {
                                // set the value to the new value
                                result.first[coordinate_pair_i/2] = value;
                              }
                            else
                              {
                                // add a new point
                                result.first.emplace_back(value);
                                result.second.emplace_back(coordinate_0);
                                result.second.emplace_back(coordinate_1);
                              }
                          }

                      }
                    else
                      {
                        // no points are provided, so use the value to fill put in the additional points
                        // at that value. Since we know that all the additional points are at the start
                        // we can easily overwrite the values.
                        for (unsigned int addition_point_i = 0; addition_point_i < addition_points.size(); ++addition_point_i)
                          {
                            result.first[addition_point_i] = value;
                          }
                      }
                  }
              }
          }
        else
          {
            // case 2: there one value, not an array
            double value = 0;
            try
              {
                value = Pointer((strict_base + "/" + name).c_str()).Get(parameters)->GetDouble();
              }
            catch (...)
              {
                WBAssertThrow(false, "Could not convert values of " << strict_base << "/" << name << " into a double. "
                              << "The provided value was \"" <<  Pointer((strict_base + "/" + name).c_str()).Get(parameters)->GetString() << "\".");
              }
            result.first.emplace_back(value);
          }
      }
    else
      {
        // there is no user defined value. Case one: return the default value and no points
        result.first.emplace_back(default_value);
      }

    return result;
  }


  std::pair<std::vector<double>,std::vector<double>> Parameters::get_value_at_array(const std::string &name)
  {
    // There are two cases:
    // 1. One double provided: use the default value everywhere. Return first with one value and second with size 0.
    // 2. Other: fill the vectors with the default value and addition points and then add new point.
    //    If a value without points is encountered, the additional points are used.
    std::pair<std::vector<double>,std::vector<double>> result;

    const std::string strict_base = this->get_full_json_path();

    // start with adding the additional points with the default value
    // to do this we need the default value
    double default_value = 0;
    bool is_array = true;
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr && Pointer((strict_base + "/" + name).c_str()).Get(parameters)->IsArray())
      {
        const std::string value_def_path = get_full_json_schema_path() + "/" + name + "/oneOf/1/items/items/anyOf/0/default value";
        Value *value_def = Pointer(value_def_path.c_str()).Get(declarations);
        WBAssertThrow(value_def != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << value_def_path);

        // Since the default value is set in the code, if it fails it is an internal error, not a user error.
        // So no try/catch needed.
        default_value = value_def->GetDouble();
      }
    else
      {
        is_array = false;
        Value *value_def = Pointer((get_full_json_schema_path() + "/" + name + "/oneOf/0/default value").c_str()).Get(declarations);
        WBAssertThrow(value_def != nullptr,
                      "internal error: could not retrieve the default value at: "
                      <<get_full_json_schema_path() + "/" + name + "/oneOf/0/default value");


        // Since the default value is set in the code, if it fails it is an internal error, not a user error.
        // So no try/catch needed.
        default_value = value_def->GetDouble();
      }


    // check if there is a user defined value
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr)
      {
        // there is a user defined value, so either case 2, 3 or 4.
        if (is_array)
          {
            Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

            // go through all the points in order
            for (size_t i = 0; i < array->Size(); ++i )
              {
                // now parse a single value_at_point.
                const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));
                // Let's assume that the file is correct, because it has been checked with the json schema.
                // So there are exactly two values, the first value is a double, the second an array of 2d arrays (points).

                // Get the double
                double value;
                Value *value_pointer = Pointer((base + "/0").c_str()).Get(parameters);

                WBAssertThrow(value_pointer != nullptr, "internal error: this should not happen.");

                try
                  {
                    value = value_pointer->GetDouble();
                  }
                catch (...)
                  {
                    WBAssertThrow(false, "Could not convert values of " << base << "/0 into doubles. "
                                  << "The provided value was \"" <<  Pointer((base + "/0").c_str()).Get(parameters)->GetString() << "\".");
                  }

                // now get the array of points.
                Value *array_of_doubles = Pointer((base + "/1").c_str()).Get(parameters);
                double value_in_array;
                if (array_of_doubles != nullptr)
                  {
                    for (size_t coordinate_i = 0; coordinate_i < array_of_doubles->Size(); ++coordinate_i )
                      {
                        Value *coordinate_j_array = Pointer((base + "/1/" + std::to_string(coordinate_i)).c_str()).Get(parameters);
                        for (size_t coordinate_j = 0; coordinate_j < coordinate_j_array->Size(); ++coordinate_j)
                          {
                            // Let's assume that the file is correct, because it has been checked with the json schema.
                            // That means that there are exactly two values per item
                            try
                              {
                                value_in_array = Pointer((base + "/1/" + std::to_string(coordinate_i) + "/" + std::to_string(coordinate_j)).c_str()).Get(parameters)->GetDouble();
                              }
                            catch (...)
                              {
                                WBAssertThrow(false, "Could not convert values of " << base + "/1/" + std::to_string(coordinate_i) + "/" + std::to_string(coordinate_j)
                                              << " into a Point<2> array, because it could not convert the 1st sub-elements into doubles. "
                                              << "The provided value was \""
                                              <<  Pointer((base + "/1/" + std::to_string(coordinate_i) + "/" + std::to_string(coordinate_j)).c_str()).Get(parameters)->GetString()
                                              << "\".");
                              }

                            result.second.emplace_back(value_in_array);
                          }
                      }
                    result.first.emplace_back(value);
                  }
              }
          }
        else
          {
            // case 1: there one value, not an array
            double value = 0;
            try
              {
                value = Pointer((strict_base + "/" + name).c_str()).Get(parameters)->GetDouble();
              }
            catch (...)
              {
                WBAssertThrow(false, "Could not convert values of " << strict_base << "/" << name << " into a double. "
                              << "The provided value was \"" <<  Pointer((strict_base + "/" + name).c_str()).Get(parameters)->GetString() << "\".");
              }
            result.first.emplace_back(0.0);
            result.second.emplace_back(value);
          }
      }
    else
      {
        // there is no user defined value. Case one: return the default value and no points
        result.first.emplace_back(0.0);
        result.second.emplace_back(default_value);
      }

    return result;
  }


  template<>
  std::vector<Point<2> >
  Parameters::get_vector(const std::string &name)
  {
    std::vector<Point<2> > vector;
    const std::string strict_base = this->get_full_json_path();
    WBAssertThrow(Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr,"Error: " << name
                  << " was not defined in " << strict_base << ", schema path: " << this->get_full_json_schema_path() << ".");

    Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

    for (size_t i = 0; i < array->Size(); ++i )
      {
        const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));
        //let's assume that the file is correct, because it has been checked with the json schema.
        // So there are exactly two values.
        double value1;
        double value2;

        try
          {
            value1 = Pointer((base + "/0").c_str()).Get(parameters)->GetDouble();
            value2 = Pointer((base + "/1").c_str()).Get(parameters)->GetDouble();
          }
        catch (...)
          {
            WBAssertThrow(false, "Could not convert values of " << base << " into a Point<2> array, because it could not convert the sub-elements into doubles.");
          }
        vector.emplace_back(value1,value2,this->coordinate_system->natural_coordinate_system());
      }

    return vector;
  }

  template<>
  std::vector<std::array<double,3> >
  Parameters::get_vector(const std::string &name)
  {
    std::vector<std::array<double,3> > vector;
    const std::string strict_base = this->get_full_json_path();
    WBAssertThrow(Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr,"Error: " << name
                  << " was not defined in " << strict_base << ", schema path: " << this->get_full_json_schema_path() << ".");

    Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

    for (size_t i = 0; i < array->Size(); ++i )
      {
        const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));
        //let's assume that the file is correct, because it has been checked with the json schema.
        // So there are exactly three values.
        try
          {
            const double value1 = Pointer((base + "/0").c_str()).Get(parameters)->GetDouble();
            const double value2 = Pointer((base + "/1").c_str()).Get(parameters)->GetDouble();
            const double value3 = Pointer((base + "/2").c_str()).Get(parameters)->GetDouble();
            vector.push_back({{value1,value2,value3}});
          }
        catch (...)
          {
            WBAssertThrow(false, "Could not convert values of " << base << " into doubles.");
          }
      }


    return vector;
  }

  std::vector<std::vector<double>>
                                Parameters::get_vector_or_double(const std::string &name)
  {
    // There are two cases:
    // 1. One double provided: use the default value everywhere. Return first with one value and second with size 0.
    // 2. std::vector<std::vector<double>> provided: output this array.
    std::vector<std::vector<double>> result;

    const std::string strict_base = this->get_full_json_path();

    // start with adding the additional points with the default value
    // to do this we need the default value
    double default_value = 0;
    bool is_array = true;
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr && Pointer((strict_base + "/" + name).c_str()).Get(parameters)->IsArray())
      {
        const std::string value_def_path = get_full_json_schema_path() + "/" + name + "/oneOf/1/items/items/default value";
        Value *value_def = Pointer(value_def_path.c_str()).Get(declarations);
        WBAssertThrow(value_def != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << value_def_path);

        // Since the default value is set in the code, if it fails it is an internal error, not a user error.
        // So no try/catch needed.
        default_value = value_def->GetDouble();
      }
    else
      {
        is_array = false;
        Value *value_def = Pointer((get_full_json_schema_path() + "/" + name + "/oneOf/0/default value").c_str()).Get(declarations);
        WBAssertThrow(value_def != nullptr,
                      "internal error: could not retrieve the default value at: "
                      <<get_full_json_schema_path() + "/" + name + "/oneOf/0/default value");


        // Since the default value is set in the code, if it fails it is an internal error, not a user error.
        // So no try/catch needed.
        default_value = value_def->GetDouble();
      }


    // check if there is a user defined value
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr)
      {
        // there is a user defined value, so either case 2, 3 or 4.
        if (is_array)
          {
            Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

            // go through all the 1d-arrays within the 2d-array
            for (size_t i = 0; i < array->Size(); ++i )
              {
                // now parse a single value_at_point.
                const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));

                // now get the array of points.
                // Value *coordinates_array = Pointer((base + "/0").c_str()).Get(parameters);
                Value *coordinates_array = Pointer((base).c_str()).Get(parameters);
                double value_in_array;
                std::vector<double> array_of_values;
                if (coordinates_array != nullptr)
                  {
                    for (size_t coordinate_i = 0; coordinate_i < coordinates_array->Size(); ++coordinate_i )
                      {
                        // Let's assume that the file is correct, because it has been checked with the json schema.
                        // That means that there are exactly two values per item
                        try
                          {
                            value_in_array = Pointer(((base + "/") + std::to_string(coordinate_i)).c_str()).Get(parameters)->GetDouble();
                          }
                        catch (...)
                          {
                            WBAssertThrow(false, "Could not convert values of " << base + "/" + std::to_string(i) + "/" + std::to_string(coordinate_i) + "/" + std::to_string(coordinate_i)
                                          << " into a Point<2> array, because it could not convert the 1st sub-elements into doubles. "
                                          << "The provided value was \""
                                          <<  Pointer((base + "/" + std::to_string(i) + "/" + std::to_string(coordinate_i)).c_str()).Get(parameters)->GetString()
                                          << "\".");
                          }
                        array_of_values.emplace_back(value_in_array);
                      }
                    result.emplace_back(array_of_values);

                  }
              }
          }
        else
          {
            // case 1: there's one value, not an array
            double value = 0;
            try
              {
                value = Pointer((strict_base + "/" + name).c_str()).Get(parameters)->GetDouble();
              }
            catch (...)
              {
                WBAssertThrow(false, "Could not convert values of " << strict_base << "/" << name << " into a double. "
                              << "The provided value was \"" <<  Pointer((strict_base + "/" + name).c_str()).Get(parameters)->GetString() << "\".");
              }
            result.emplace_back(std::vector<double> {value});
          }
      }
    else
      {
        // there is no user defined value. Case one: return the default value and no points
        result.emplace_back(std::vector<double> {default_value});
      }

    return result;
  }



  template<>
  std::vector<std::array<std::array<double,3>,3> >
  Parameters::get_vector(const std::string &name)
  {
    std::vector<std::array<std::array<double,3>,3>  > vector;
    const std::string strict_base = this->get_full_json_path();
    WBAssertThrow(Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr,"Error: " << name
                  << " was not defined in " << strict_base << ", schema path: " << this->get_full_json_schema_path() << ".");

    Value *array1 = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

    for (size_t i = 0; i < array1->Size(); ++i )
      {
        const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));
        Value *array2 = Pointer((base).c_str()).Get(parameters);

        // Not sure why cppcheck it is generating the warning
        // Filed a question at: https://sourceforge.net/p/cppcheck/discussion/general/thread/429759f85e/
        // cppcheck-suppress constStatement
        std::array<std::array<double,3>,3> array;
        WBAssertThrow(array2->Size() == 3, "Array " << i << " is supposed to be a 3x3 array, but the outer array dimensions is "
                      << array2->Size() << '.');
        for (size_t j = 0; j < array2->Size(); ++j )
          {
            const std::string base_extended = base + "/" + std::to_string(j);

            WBAssertThrow(Pointer((base_extended).c_str()).Get(parameters)->Size() == 3,
                          "Array " << i << " is supposed to be a 3x3 array, but the inner array dimensions of "
                          << j << " is " << Pointer((base_extended).c_str()).Get(parameters)->Size() << '.');
            double value1;
            double value2;
            double value3;

            try
              {
                value1 = Pointer((base_extended + "/0").c_str()).Get(parameters)->GetDouble();
                value2 = Pointer((base_extended + "/1").c_str()).Get(parameters)->GetDouble();
                value3 = Pointer((base_extended + "/2").c_str()).Get(parameters)->GetDouble();
              }
            catch (...)
              {
                WBAssertThrow(false, "Could not convert values of " << base << " into doubles.");
              }
            array[j][0] = value1;
            array[j][1] = value2;
            array[j][2] = value3;
          }
        vector.push_back(array);
      }

    return vector;
  }


  template<>
  std::vector<std::vector<Point<2> > >
  Parameters::get_vector(const std::string &name)
  {
    std::vector<std::vector<Point<2> >  > vector;
    const std::string strict_base = this->get_full_json_path();
    WBAssertThrow(Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr,"Error: " << name
                  << " was not defined in " << strict_base << ", schema path: " << this->get_full_json_schema_path() << ".");

    Value *array1 = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

    for (size_t i = 0; i < array1->Size(); ++i )
      {
        const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));
        Value *array2 = Pointer((base).c_str()).Get(parameters);

        // Not sure why cppcheck it is generating the warning
        // Filed a question at: https://sourceforge.net/p/cppcheck/discussion/general/thread/429759f85e/
        // cppcheck-suppress constStatement
        std::vector<Point<2> > sub_vector(array2->Size(),Point<2>(NaN::DSNAN,NaN::DSNAN,coordinate_system->natural_coordinate_system()));
        for (size_t j = 0; j < array2->Size(); ++j )
          {
            const std::string base_extended = base + "/" + std::to_string(j);

            WBAssertThrow(Pointer((base_extended).c_str()).Get(parameters)->Size() == 2,
                          "Array " << i << " is supposed to be a 2d point, but the inner array dimensions of "
                          << j << " is " << Pointer((base_extended).c_str()).Get(parameters)->Size() << '.');
            double value1;
            double value2;

            try
              {
                value1 = Pointer((base_extended + "/0").c_str()).Get(parameters)->GetDouble();
                value2 = Pointer((base_extended + "/1").c_str()).Get(parameters)->GetDouble();
              }
            catch (...)
              {
                WBAssertThrow(false, "Could not convert values of " << base << " into doubles, because it could not covert the sub-elements into doubles.");
              }
            sub_vector[j][0] = value1;
            sub_vector[j][1] = value2;
          }
        vector.push_back(sub_vector);
      }

    return vector;
  }


  template<>
  std::vector<Objects::Segment<Features::SubductingPlateModels::Temperature::Interface,
      Features::SubductingPlateModels::Composition::Interface,
      Features::SubductingPlateModels::Grains::Interface> >
      Parameters::get_vector(const std::string &name,
                             std::vector<std::shared_ptr<Features::SubductingPlateModels::Temperature::Interface> > &default_temperature_models,
                             std::vector<std::shared_ptr<Features::SubductingPlateModels::Composition::Interface> > &default_composition_models,
                             std::vector<std::shared_ptr<Features::SubductingPlateModels::Grains::Interface> > &default_grains_models)
  {
    using namespace Features::SubductingPlateModels;
    std::vector<Objects::Segment<Temperature::Interface,Composition::Interface,Grains::Interface> > vector;
    this->enter_subsection(name);
    const std::string strict_base = this->get_full_json_path();
    WBAssertThrow(Pointer((strict_base).c_str()).Get(parameters) != nullptr,"Error: " << name
                  << " was not defined in " << strict_base << ", schema path: " << this->get_full_json_schema_path() << ".");

    // get the array of segments
    Value *array = Pointer((strict_base).c_str()).Get(parameters);

    for (size_t i = 0; i < array->Size(); ++i )
      {
        this->enter_subsection(std::to_string(i));
        const std::string base = this->get_full_json_path();
        // get one segment
        // length
        const double length = Pointer((base + "/length").c_str()).Get(parameters)->GetDouble();

        // get thickness
        Value *point_array = Pointer((base  + "/thickness").c_str()).Get(parameters);
        Point<2> thickness(invalid);
        WBAssertThrow(point_array != nullptr,"error: " << name
                      << " was not defined in " << strict_base << "/thickness, schema path: " << this->get_full_json_schema_path() << "/thickness.");

        if (point_array->Size() == 1)
          {
            // There is only one value, set it for both elements
            const double local0 = Pointer((base + "/thickness/0").c_str()).Get(parameters)->GetDouble();
            thickness = Point<2>(local0,local0,invalid);
          }
        else
          {
            const double local0 = Pointer((base + "/thickness/0").c_str()).Get(parameters)->GetDouble();
            const double local1 = Pointer((base + "/thickness/1").c_str()).Get(parameters)->GetDouble();
            thickness = Point<2>(local0,local1,invalid);
          }


        // get top truncation (default is 0,0)
        point_array = Pointer((base  + "/top truncation").c_str()).Get(parameters);
        Point<2> top_truncation(invalid);
        if (point_array != nullptr)
          {
            if (point_array->Size() == 1)
              {
                // There is only one value, set it for both elements
                const double local0 = Pointer((base + "/top truncation/0").c_str()).Get(parameters)->GetDouble();
                top_truncation = Point<2>(local0,local0,invalid);
              }
            else
              {
                const double local0 = Pointer((base + "/top truncation/0").c_str()).Get(parameters)->GetDouble();
                const double local1 = Pointer((base + "/top truncation/1").c_str()).Get(parameters)->GetDouble();
                top_truncation = Point<2>(local0,local1,invalid);
              }
          }
        // get thickness
        point_array = Pointer((base  + "/angle").c_str()).Get(parameters);
        Point<2> angle(invalid);
        WBAssertThrow(point_array != nullptr,"error: " << name
                      << " was not defined in " << strict_base << "/angle, schema path: " << this->get_full_json_schema_path() << "/angle.");

        if (point_array->Size() == 1)
          {
            // There is only one value, set it for both elements
            const double local0 = Pointer((base + "/angle/0").c_str()).Get(parameters)->GetDouble();
            angle = Point<2>(local0,local0,invalid);
          }
        else
          {
            const double local0 = Pointer((base + "/angle/0").c_str()).Get(parameters)->GetDouble();
            const double local1 = Pointer((base + "/angle/1").c_str()).Get(parameters)->GetDouble();
            angle = Point<2>(local0,local1,invalid);
          }


        // Get temperature models
        std::vector<std::shared_ptr<Temperature::Interface> > temperature_models;

        //This is a value to look back in the path elements.
        size_t searchback = 0;
        if (!this->get_shared_pointers<Temperature::Interface>("temperature models", temperature_models) ||
            Pointer((base + "/temperature model default entry").c_str()).Get(parameters) != nullptr)
          {
            temperature_models = default_temperature_models;

            // find the default value, which is the closest to the current path
            for (searchback = 0; searchback < path.size(); ++searchback)
              {
                if (Pointer((this->get_full_json_path(path.size()-searchback) + "/temperature models").c_str()).Get(parameters) != nullptr)
                  {
                    break;
                  }
              }

            // if we can not find default value for the temperature model, skip it
            if (searchback < path.size())
              {

                // copy the value, this unfortunately removes it.
                Value value1 = Value(Pointer((this->get_full_json_path(path.size()-searchback) + "/temperature models").c_str()).Get(parameters)->GetArray());

                // now copy it
                Value value2;
                value2.CopyFrom(value1, parameters.GetAllocator());

                // now we should have 2x the same value, so put it back and place it in the correct location.
                Pointer((this->get_full_json_path(path.size()-searchback) + "/temperature models").c_str()).Set(parameters, value1);//.Get(parameters)->Set("temperature models", value1, parameters.GetAllocator());

                Pointer((base).c_str()).Get(parameters)->AddMember("temperature models", value2, parameters.GetAllocator());
                Pointer((base + "/temperature model default entry").c_str()).Set(parameters,true);
              }
          }

        // now do the same for compositions
        std::vector<std::shared_ptr<Composition::Interface> > composition_models;
        if (!this->get_shared_pointers<Composition::Interface>("composition models", composition_models) ||
            Pointer((base + "/composition model default entry").c_str()).Get(parameters) != nullptr)
          {
            composition_models = default_composition_models;


            // find the default value, which is the closest to the current path
            for (searchback = 0; searchback < path.size(); ++searchback)
              {
                if (Pointer((this->get_full_json_path(path.size()-searchback) + "/composition models").c_str()).Get(parameters) != nullptr)
                  {
                    break;
                  }
              }

            // if we can not find default value for the temperature model, skip it
            if (searchback < path.size())
              {

                // copy the value, this unfortunately removes it.
                Value value1 = Value(Pointer((this->get_full_json_path(path.size()-searchback) + "/composition models").c_str()).Get(parameters)->GetArray());

                // now copy it
                Value value2;
                value2.CopyFrom(value1, parameters.GetAllocator());

                // now we should have 2x the same value, so put it back and place it in the correct location.
                Pointer((this->get_full_json_path(path.size()-searchback) + "/composition models").c_str()).Set(parameters, value1);//.Get(parameters)->Set("temperature models", value1, parameters.GetAllocator());

                Pointer((base).c_str()).Get(parameters)->AddMember("composition models", value2, parameters.GetAllocator());
                Pointer((base + "/composition model default entry").c_str()).Set(parameters,true);
              }
          }

        // now do the same for grains
        std::vector<std::shared_ptr<Grains::Interface> > grains_models;
        if (!this->get_shared_pointers<Grains::Interface>("grains models", grains_models) ||
            Pointer((base + "/grains model default entry").c_str()).Get(parameters) != nullptr)
          {
            grains_models = default_grains_models;


            // find the default value, which is the closest to the current path
            for (searchback = 0; searchback < path.size(); ++searchback)
              {
                if (Pointer((this->get_full_json_path(path.size()-searchback) + "/grains models").c_str()).Get(parameters) != nullptr)
                  {
                    break;
                  }
              }

            // if we can not find default value for the temperature model, skip it
            if (searchback < path.size())
              {

                // copy the value, this unfortunately removes it.
                Value value1 = Value(Pointer((this->get_full_json_path(path.size()-searchback) + "/grains models").c_str()).Get(parameters)->GetArray());

                // now copy it
                Value value2;
                value2.CopyFrom(value1, parameters.GetAllocator());

                // now we should have 2x the same value, so put it back and place it in the correct location.
                Pointer((this->get_full_json_path(path.size()-searchback) + "/grains models").c_str()).Set(parameters, value1);//.Get(parameters)->Set("temperature models", value1, parameters.GetAllocator());

                Pointer((base).c_str()).Get(parameters)->AddMember("grains models", value2, parameters.GetAllocator());
                Pointer((base + "/grains model default entry").c_str()).Set(parameters,true);
              }
          }
        vector.emplace_back(length, thickness, top_truncation, angle, temperature_models, composition_models, grains_models);

        this->leave_subsection();
      }

    this->leave_subsection();
    return vector;
  }


  template<>
  std::vector<Objects::Segment<Features::FaultModels::Temperature::Interface,Features::FaultModels::Composition::Interface, Features::FaultModels::Grains::Interface> >
  Parameters::get_vector(const std::string &name,
                         std::vector<std::shared_ptr<Features::FaultModels::Temperature::Interface> > &default_temperature_models,
                         std::vector<std::shared_ptr<Features::FaultModels::Composition::Interface> > &default_composition_models,
                         std::vector<std::shared_ptr<Features::FaultModels::Grains::Interface> > &default_grains_models)
  {
    using namespace Features::FaultModels;
    std::vector<Objects::Segment<Temperature::Interface,Composition::Interface,Grains::Interface> > vector;
    this->enter_subsection(name);
    const std::string strict_base = this->get_full_json_path();
    WBAssertThrow(Pointer((strict_base).c_str()).Get(parameters) != nullptr,"error: " << name
                  << " was not defined in " << strict_base << ", schema path: " << this->get_full_json_schema_path() << ".");

    // get the array of segments
    Value *array = Pointer((strict_base).c_str()).Get(parameters);

    for (size_t i = 0; i < array->Size(); ++i )
      {
        this->enter_subsection(std::to_string(i));
        const std::string base = this->get_full_json_path();
        // get one segment
        // length
        const double length = Pointer((base + "/length").c_str()).Get(parameters)->GetDouble();

        // get thickness
        Value *point_array = Pointer((base  + "/thickness").c_str()).Get(parameters);
        Point<2> thickness(invalid);
        WBAssertThrow(point_array != nullptr,"error: " << name
                      << " was not defined in " << strict_base << "/thickness, schema path: " << this->get_full_json_schema_path() << "/thickness.");

        if (point_array->Size() == 1)
          {
            // There is only one value, set it for both elements
            const double local0 = Pointer((base + "/thickness/0").c_str()).Get(parameters)->GetDouble();
            thickness = Point<2>(local0,local0,invalid);
          }
        else
          {
            const double local0 = Pointer((base + "/thickness/0").c_str()).Get(parameters)->GetDouble();
            const double local1 = Pointer((base + "/thickness/1").c_str()).Get(parameters)->GetDouble();
            thickness = Point<2>(local0,local1,invalid);
          }


        // get top truncation (default is 0,0)
        point_array = Pointer((base  + "/top truncation").c_str()).Get(parameters);
        Point<2> top_truncation(invalid);
        if (point_array != nullptr)
          {
            if (point_array->Size() == 1)
              {
                // There is only one value, set it for both elements
                const double local0 = Pointer((base + "/top truncation/0").c_str()).Get(parameters)->GetDouble();
                top_truncation = Point<2>(local0,local0,invalid);
              }
            else
              {
                const double local0 = Pointer((base + "/top truncation/0").c_str()).Get(parameters)->GetDouble();
                const double local1 = Pointer((base + "/top truncation/1").c_str()).Get(parameters)->GetDouble();
                top_truncation = Point<2>(local0,local1,invalid);
              }
          }
        // get thickness
        point_array = Pointer((base  + "/angle").c_str()).Get(parameters);
        Point<2> angle(invalid);
        WBAssertThrow(point_array != nullptr,"error: " << name
                      << " was not defined in " << strict_base << "/angle, schema path: " << this->get_full_json_schema_path() << "/angle.");

        if (point_array->Size() == 1)
          {
            // There is only one value, set it for both elements
            const double local0 = Pointer((base + "/angle/0").c_str()).Get(parameters)->GetDouble();
            angle = Point<2>(local0,local0,invalid);
          }
        else
          {
            const double local0 = Pointer((base + "/angle/0").c_str()).Get(parameters)->GetDouble();
            const double local1 = Pointer((base + "/angle/1").c_str()).Get(parameters)->GetDouble();
            angle = Point<2>(local0,local1,invalid);
          }


        // Get temperature models
        std::vector<std::shared_ptr<Temperature::Interface> > temperature_models;

        //This is a value to look back in the path elements.
        size_t searchback = 0;
        if (!this->get_shared_pointers<Temperature::Interface>("temperature models", temperature_models) ||
            Pointer((base + "/temperature model default entry").c_str()).Get(parameters) != nullptr)
          {
            temperature_models = default_temperature_models;

            // find the default value, which is the closest to the current path
            for (searchback = 0; searchback < path.size(); ++searchback)
              {
                if (Pointer((this->get_full_json_path(path.size()-searchback) + "/temperature models").c_str()).Get(parameters) != nullptr)
                  {
                    break;
                  }
              }

            // if we can not find default value for the temperature model, skip it
            if (searchback < path.size())
              {

                // copy the value, this unfortunately removes it.
                Value value1 = Value(Pointer((this->get_full_json_path(path.size()-searchback) + "/temperature models").c_str()).Get(parameters)->GetArray());

                // now copy it
                Value value2;
                value2.CopyFrom(value1, parameters.GetAllocator());

                // now we should have 2x the same value, so put it back and place it in the correct location.
                Pointer((this->get_full_json_path(path.size()-searchback) + "/temperature models").c_str()).Set(parameters, value1);//.Get(parameters)->Set("temperature models", value1, parameters.GetAllocator());

                Pointer((base).c_str()).Get(parameters)->AddMember("temperature models", value2, parameters.GetAllocator());
                Pointer((base + "/temperature model default entry").c_str()).Set(parameters,true);
              }
          }

        // now do the same for compositions
        std::vector<std::shared_ptr<Composition::Interface> > composition_models;
        if (!this->get_shared_pointers<Composition::Interface>("composition models", composition_models) ||
            Pointer((base + "/composition model default entry").c_str()).Get(parameters) != nullptr)
          {
            composition_models = default_composition_models;


            // find the default value, which is the closest to the current path
            for (searchback = 0; searchback < path.size(); ++searchback)
              {
                if (Pointer((this->get_full_json_path(path.size()-searchback) + "/composition models").c_str()).Get(parameters) != nullptr)
                  {
                    break;
                  }
              }

            // if we can not find default value for the temperature model, skip it
            if (searchback < path.size())
              {

                // copy the value, this unfortunately removes it.
                Value value1 = Value(Pointer((this->get_full_json_path(path.size()-searchback) + "/composition models").c_str()).Get(parameters)->GetArray());

                // now copy it
                Value value2;
                value2.CopyFrom(value1, parameters.GetAllocator());

                // now we should have 2x the same value, so put it back and place it in the correct location.
                Pointer((this->get_full_json_path(path.size()-searchback) + "/composition models").c_str()).Set(parameters, value1);//.Get(parameters)->Set("temperature models", value1, parameters.GetAllocator());

                Pointer((base).c_str()).Get(parameters)->AddMember("composition models", value2, parameters.GetAllocator());
                Pointer((base + "/composition model default entry").c_str()).Set(parameters,true);
              }
          }

        // now do the same for grains
        std::vector<std::shared_ptr<Grains::Interface> > grains_models;
        if (!this->get_shared_pointers<Grains::Interface>("grains models", grains_models) ||
            Pointer((base + "/grains model default entry").c_str()).Get(parameters) != nullptr)
          {
            grains_models = default_grains_models;


            // find the default value, which is the closest to the current path
            for (searchback = 0; searchback < path.size(); ++searchback)
              {
                if (Pointer((this->get_full_json_path(path.size()-searchback) + "/grains models").c_str()).Get(parameters) != nullptr)
                  {
                    break;
                  }
              }

            // if we can not find default value for the temperature model, skip it
            if (searchback < path.size())
              {

                // copy the value, this unfortunately removes it.
                Value value1 = Value(Pointer((this->get_full_json_path(path.size()-searchback) + "/grains models").c_str()).Get(parameters)->GetArray());

                // now copy it
                Value value2;
                value2.CopyFrom(value1, parameters.GetAllocator());

                // now we should have 2x the same value, so put it back and place it in the correct location.
                Pointer((this->get_full_json_path(path.size()-searchback) + "/grains models").c_str()).Set(parameters, value1);//.Get(parameters)->Set("temperature models", value1, parameters.GetAllocator());

                Pointer((base).c_str()).Get(parameters)->AddMember("grains models", value2, parameters.GetAllocator());
                Pointer((base + "/grains model default entry").c_str()).Set(parameters,true);
              }
          }
        vector.emplace_back(length, thickness, top_truncation, angle, temperature_models, composition_models, grains_models);

        this->leave_subsection();
      }

    this->leave_subsection();
    return vector;
  }

  template<>
  std::vector<double>
  Parameters::get_vector(const std::string &name)
  {
    std::vector<double> vector;
    const std::string strict_base = this->get_full_json_path();
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr)
      {
        Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

        for (size_t i = 0; i < array->Size(); ++i )
          {
            const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));

            vector.push_back(Pointer(base.c_str()).Get(parameters)->GetDouble());
          }
      }
    else
      {
        const Value *value = Pointer((this->get_full_json_schema_path()  + "/" + name + "/minItems").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the minItems value at: "
                      << this->get_full_json_schema_path() + "/" + name + "/minItems");

        const size_t min_size = value->GetUint();

        value = Pointer((this->get_full_json_schema_path()  + "/" + name + "/items/default value").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << this->get_full_json_schema_path() + "/" + name + "/default value");

        const double default_value = value->GetDouble();

        // set to min size
        for (size_t i = 0; i < min_size; ++i)
          {
            vector.push_back(default_value);
          }
      }
    return vector;
  }

  template<>
  std::vector<size_t>
  Parameters::get_vector(const std::string &name)
  {
    std::vector<size_t> vector;
    const std::string strict_base = this->get_full_json_path();
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr)
      {
        Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

        for (size_t i = 0; i < array->Size(); ++i )
          {
            const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));

            vector.push_back(Pointer(base.c_str()).Get(parameters)->GetUint());
          }
      }
    else
      {
        const Value *value = Pointer((this->get_full_json_schema_path()  + "/" + name + "/minItems").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the minItems value at: "
                      << this->get_full_json_schema_path() + "/" + name + "/minItems");

        const size_t min_size = value->GetUint();

        value = Pointer((this->get_full_json_schema_path()  + "/" + name + "/items/default value").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << this->get_full_json_schema_path() + "/" + name + "/default value");

        const size_t default_value = value->GetUint();

        // set to min size
        for (size_t i = 0; i < min_size; ++i)
          {
            vector.push_back(default_value);
          }
      }
    return vector;
  }

  template<>
  std::vector<unsigned int>
  Parameters::get_vector(const std::string &name)
  {
    std::vector<unsigned int> vector;
    const std::string strict_base = this->get_full_json_path();
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr)
      {
        Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

        for (size_t i = 0; i < array->Size(); ++i )
          {
            const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));

            vector.push_back(Pointer(base.c_str()).Get(parameters)->GetUint());
          }
      }
    else
      {
        const Value *value = Pointer((this->get_full_json_schema_path()  + "/" + name + "/minItems").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the minItems value at: "
                      << this->get_full_json_schema_path() + "/" + name + "/minItems value");

        const size_t min_size = value->GetUint();

        const unsigned int default_value = Pointer((this->get_full_json_schema_path()  + "/" + name + "/items/default value").c_str()).Get(declarations)->GetUint();

        // set to min size
        for (size_t i = 0; i < min_size; ++i)
          {
            vector.push_back(default_value);
          }
      }
    return vector;
  }

  template<class T>
  std::unique_ptr<T>
  Parameters::get_unique_pointer(const std::string &name)
  {
    const std::string base = this->get_full_json_path();
    Value *value = Pointer((base + "/" + name + "/model").c_str()).Get(parameters);

#ifdef debug
    bool required = false;
    if (Pointer((base + "/required").c_str()).Get(declarations) != NULL)
      for (auto &v : Pointer((base + "/required").c_str()).Get(declarations)->GetArray())
        {
          if (v.GetString() == name)
            {
              required = true;
            }
        }

    WBAssert(value != NULL || required == false,
             "Internal error: Value \"" << base << '/' << name << "/model\" not found in the input file, while it was set as required.");
#endif
    if (value == nullptr)
      {
        value = Pointer((get_full_json_schema_path() + "/" + name + "/default value").c_str()).Get(declarations);
        WBAssertThrow(value != nullptr,
                      "internal error: could not retrieve the default value at: "
                      << base + "/" + name + "/default value. Make sure the value has been declared.");
      }

    return T::create(value->GetString(),&world);
  }

  template<class T>
  bool
  Parameters::get_unique_pointers(const std::string &name, std::vector<std::unique_ptr<T> > &vector)
  {
    vector.resize(0);
    const std::string strict_base = this->get_full_json_path();
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr)
      {
        Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

        for (size_t i = 0; i < array->Size(); ++i )
          {
            const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));

            const std::string value = Pointer((base + "/model").c_str()).Get(parameters)->GetString();

            vector.push_back(std::move(T::create(value, &world)));
          }
      }
    else
      {
        return false;
      }

    return true;
  }

  template<>
  bool
  Parameters::get_unique_pointers(const std::string &name, std::vector<std::unique_ptr<Features::SubductingPlate> > &vector)
  {
    vector.resize(0);
    const std::string strict_base = this->get_full_json_path();
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr)
      {
        Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

        for (size_t i = 0; i < array->Size(); ++i )
          {
            vector.push_back(std::make_unique<Features::SubductingPlate>(&world));
          }
      }
    else
      {
        return false;
      }

    return true;
  }

  template<>
  bool
  Parameters::get_unique_pointers(const std::string &name, std::vector<std::unique_ptr<Features::Fault> > &vector)
  {
    vector.resize(0);
    const std::string strict_base = this->get_full_json_path();
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr)
      {
        Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

        for (size_t i = 0; i < array->Size(); ++i )
          {
            vector.push_back(std::make_unique<Features::Fault>(&world));
          }
      }
    else
      {
        return false;
      }

    return true;
  }



  template<class T>
  bool
  Parameters::get_shared_pointers(const std::string &name, std::vector<std::shared_ptr<T> > &vector)
  {
    vector.resize(0);
    const std::string strict_base = this->get_full_json_path();
    if (Pointer((strict_base + "/" + name).c_str()).Get(parameters) != nullptr)
      {
        Value *array = Pointer((strict_base  + "/" + name).c_str()).Get(parameters);

        for (size_t i = 0; i < array->Size(); ++i )
          {
            const std::string base = (strict_base + "/").append(name).append("/").append(std::to_string(i));

            const std::string value = Pointer((base + "/model").c_str()).Get(parameters)->GetString();

            vector.push_back(std::move(T::create(value, &world)));
          }
      }
    else
      {
        return false;
      }

    return true;
  }

  void
  Parameters::enter_subsection(const std::string &name)
  {
    path.push_back(name);
    //TODO: WBAssert(is path valid?)
  }

  void
  Parameters::leave_subsection()
  {
    path.pop_back();
  }



  void
  Parameters::declare_model_entries(const std::string &model_group_name,
                                    const std::string &parent_name,
                                    const std::map<std::string, void ( *)(Parameters &,const std::string &)> &declare_map,
                                    const std::vector<std::string> &required_entries,
                                    const std::vector<std::tuple<std::string,const WorldBuilder::Types::Interface &, std::string> > &extra_declarations)
  {
    unsigned int counter = 0;
    for (const auto &it : declare_map)
      {
        typedef std::tuple<std::string,const WorldBuilder::Types::Interface &, std::string> DeclareEntry;
        // prevent infinite recursion
        if (it.first != parent_name)
          {
            enter_subsection("oneOf");
            {
              enter_subsection(std::to_string(counter));
              {
                enter_subsection("properties");
                {
                  declare_entry("", Types::Object(required_entries), model_group_name + " object");

                  declare_entry("model", Types::String("",it.first),
                                "The name of the " + model_group_name + " model.");

                  for (DeclareEntry extra_declaration : extra_declarations)
                    {
                      declare_entry(std::get<0>(extra_declaration),std::get<1>(extra_declaration),std::get<2>(extra_declaration));
                    }


                  it.second(*this, parent_name);
                }
                leave_subsection();
              }
              leave_subsection();
            }
            leave_subsection();

            counter++;

          }
      }
  }



  std::string
  Parameters::get_full_json_path(size_t max_size) const
  {
    std::string collapse;
    for (size_t i = 0; i < path.size() && i < max_size; i++)
      {
        collapse +=  "/" + path[i];
      }
    return collapse;
  }

  std::string
  Parameters::get_full_json_schema_path() const
  {
    std::string collapse = "/properties";
    for (size_t i = 0; i < path.size(); i++)
      {
        // first get the type
        //WBAssert(Pointer((collapse + "/" + path[i] + "/type").c_str()).Get(declarations) != NULL, "Internal error: could not find " << collapse + "/" + path[i] + "/type");

        const std::string base_path = Pointer((collapse + "/" + path[i] + "/type").c_str()).Get(declarations) != nullptr
                                      ?
                                      collapse + "/" + path[i]
                                      :
                                      collapse;

        const std::string type = Pointer((base_path + "/type").c_str()).Get(declarations)->GetString();

        if (type == "array")
          {
            // the type is an array. Arrays always have an items, but can also
            // have a oneOf (todo: or anyOf ...). Find out whether this is the case
            //collapse += path[i] + "/items";
            if (Pointer((base_path + "/items/oneOf").c_str()).Get(declarations) != nullptr)
              {
                // it has a structure with oneOf. Find out which of the entries is needed.
                // This means we have to take a sneak peak to figure out how to get to the
                // next value.
                const size_t size = Pointer((base_path + "/items/oneOf").c_str()).Get(declarations)->Size();
#ifdef debug
                bool found = false;
#endif
                size_t index = 0;
                for (; index < size; ++index)
                  {
                    const std::string declarations_string = Pointer((base_path + "/items/oneOf/" + std::to_string(index)
                                                                     + "/properties/model/enum/0").c_str()).Get(declarations)->GetString();

                    // we need to get the json path relevant for the current declaration string
                    // we are interested in, which requires an offset of 2.
                    WBAssert(Pointer((get_full_json_path(i+2) + "/model").c_str()).Get(parameters) != nullptr, "Could not find model in: " << get_full_json_path(i+2) + "/model");
                    const std::string parameters_string = Pointer((get_full_json_path(i+2) + "/model").c_str()).Get(parameters)->GetString();

                    // currently in our case these are always objects, so go directly to find the option we need.
                    if (declarations_string == parameters_string)
                      {
                        // found it for index i;
#ifdef debug
                        found = true;
#endif
                        break;
                      }
                  }
#ifdef debug
                WBAssert(found == true,
                         "Internal error: This is an array with several possible values, "
                         "but could not find the correct value " << collapse + "/" + path[i] + "/items/oneOf");
#endif
                collapse += "/" + path[i] + "/items/oneOf/" + std::to_string(index) + "/properties";
                // add one to i, to skip the array
                ++i;
              }
            else
              {
                collapse = base_path + "/items";
              }
          }
        else if (type == "object")
          {
            // normal objects just have properties, but some have oneOf
            if (Pointer((base_path + "/oneOf").c_str()).Get(declarations) != nullptr)
              {
                // it has a structure with oneOf. Find out which of the entries is needed.
                // This means we have to take a sneak peak to figure out how to get to the
                // next value.
                const size_t size = Pointer((base_path + "/oneOf").c_str()).Get(declarations)->Size();
#ifdef debug
                bool found = false;
#endif
                size_t index = 0;
                for (; index < size; ++index)
                  {
                    const std::string declarations_string = Pointer((base_path + "/oneOf/" + std::to_string(index)
                                                                     + "/properties/model/enum/0").c_str()).Get(declarations)->GetString();

                    // we need to get the json path relevant for the current declaration string
                    // we are interested in, which requires an offset of 2.
                    const Value *value = Pointer((get_full_json_path(i+2) + "/model").c_str()).Get(parameters);


                    if (value == nullptr)
                      {
                        value = Pointer((collapse + "/" + path[i]  + "/default value").c_str()).Get(declarations);
                        WBAssertThrow(value != nullptr,
                                      "internal error: could not retrieve the default value at: "
                                      << collapse + "/" + path[i]  + "/default value, for value: ");
                      }

                    const std::string parameters_string = value->GetString();

                    // currently in our case these are always objects, so go directly to find the option we need.
                    if (declarations_string == parameters_string)
                      {
                        // found it for index i;
#ifdef debug
                        found = true;
#endif
                        break;
                      }
                  }
#ifdef debug
                WBAssert(found == true,
                         "Internal error: This is an array with several possible values, "
                         "but could not find the correct value " << collapse + "/" + path[i] + "/items/oneOf");
#endif
                collapse += "/" + path[i] + "/oneOf/" + std::to_string(index) + "/properties";
                // add one to i, to skip the array
                ++i;

              }
            else
              {
                collapse += "/properties";
              }
          }
        else
          {
            collapse += "/" + path[i];
          }
      }
    return collapse;
  }





  /**
   * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
   * Note that the variable with this name has to be loaded before this function is called.
   */
  template std::unique_ptr<CoordinateSystems::Interface> Parameters::get_unique_pointer<CoordinateSystems::Interface>(const std::string &name);


  /**
   * Returns a vector of pointers to the Gravity Model based on the provided name.
   * Note that the variable with this name has to be loaded before this function is called.
   */
  template std::unique_ptr<GravityModel::Interface> Parameters::get_unique_pointer<GravityModel::Interface>(const std::string &name);



  /**
   * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
   * Note that the variable with this name has to be loaded before this function is called.
   */
  template bool
  Parameters::get_unique_pointers<Features::Interface>(const std::string &name,
                                                       std::vector<std::unique_ptr<Features::Interface> > &vector);


  /**
   * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
   * Note that the variable with this name has to be loaded before this function is called.
   */
  template bool
  Parameters::get_unique_pointers<Features::ContinentalPlateModels::Temperature::Interface>(const std::string &name,
      std::vector<std::unique_ptr<Features::ContinentalPlateModels::Temperature::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::ContinentalPlateModels::Composition::Interface>(const std::string &name,
      std::vector<std::unique_ptr<Features::ContinentalPlateModels::Composition::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::ContinentalPlateModels::Grains::Interface>(const std::string &name,
      std::vector<std::unique_ptr<Features::ContinentalPlateModels::Grains::Interface> > &vector);


  /**
   * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
   * Note that the variable with this name has to be loaded before this function is called.
   */
  template bool
  Parameters::get_unique_pointers<Features::PlumeModels::Temperature::Interface>(const std::string &name,
                                                                                 std::vector<std::unique_ptr<Features::PlumeModels::Temperature::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::PlumeModels::Composition::Interface>(const std::string &name,
                                                                                 std::vector<std::unique_ptr<Features::PlumeModels::Composition::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::PlumeModels::Grains::Interface>(const std::string &name,
                                                                            std::vector<std::unique_ptr<Features::PlumeModels::Grains::Interface> > &vector);



  /**
   * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
   * Note that the variable with this name has to be loaded before this function is called.
   */
  template bool
  Parameters::get_unique_pointers<Features::OceanicPlateModels::Temperature::Interface>(const std::string &name,
      std::vector<std::unique_ptr<Features::OceanicPlateModels::Temperature::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::OceanicPlateModels::Composition::Interface>(const std::string &name,
      std::vector<std::unique_ptr<Features::OceanicPlateModels::Composition::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::OceanicPlateModels::Grains::Interface>(const std::string &name,
      std::vector<std::unique_ptr<Features::OceanicPlateModels::Grains::Interface> > &vector);



  /**
   * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
   * Note that the variable with this name has to be loaded before this function is called.
   */
  template bool
  Parameters::get_unique_pointers<Features::MantleLayerModels::Temperature::Interface>(const std::string &name,
      std::vector<std::unique_ptr<Features::MantleLayerModels::Temperature::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::MantleLayerModels::Composition::Interface>(const std::string &name,
      std::vector<std::unique_ptr<Features::MantleLayerModels::Composition::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::MantleLayerModels::Grains::Interface>(const std::string &name,
                                                                                  std::vector<std::unique_ptr<Features::MantleLayerModels::Grains::Interface> > &vector);



  /**
   * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
   * Note that the variable with this name has to be loaded before this function is called.
   */
  template bool
  Parameters::get_unique_pointers<Features::SubductingPlateModels::Temperature::Interface>(const std::string &name,
      std::vector<std::unique_ptr<Features::SubductingPlateModels::Temperature::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::SubductingPlateModels::Composition::Interface>(const std::string &name,
      std::vector<std::unique_ptr<Features::SubductingPlateModels::Composition::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::SubductingPlateModels::Grains::Interface>(const std::string &name,
      std::vector<std::unique_ptr<Features::SubductingPlateModels::Grains::Interface> > &vector);



  /**
   * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
   * Note that the variable with this name has to be loaded before this function is called.
   */
  template bool
  Parameters::get_unique_pointers<Features::FaultModels::Temperature::Interface>(const std::string &name,
                                                                                 std::vector<std::unique_ptr<Features::FaultModels::Temperature::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::FaultModels::Composition::Interface>(const std::string &name,
                                                                                 std::vector<std::unique_ptr<Features::FaultModels::Composition::Interface> > &vector);
  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_unique_pointers<Features::FaultModels::Grains::Interface>(const std::string &name,
                                                                            std::vector<std::unique_ptr<Features::FaultModels::Grains::Interface> > &vector);



  /**
   * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
   * Note that the variable with this name has to be loaded before this function is called.
   */
  template bool
  Parameters::get_shared_pointers<Features::SubductingPlateModels::Temperature::Interface>(const std::string &name,
      std::vector<std::shared_ptr<Features::SubductingPlateModels::Temperature::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_shared_pointers<Features::SubductingPlateModels::Composition::Interface>(const std::string &name,
      std::vector<std::shared_ptr<Features::SubductingPlateModels::Composition::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_shared_pointers<Features::SubductingPlateModels::Grains::Interface>(const std::string &name,
      std::vector<std::shared_ptr<Features::SubductingPlateModels::Grains::Interface> > &vector);


  /**
   * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
   * Note that the variable with this name has to be loaded before this function is called.
   */
  template bool
  Parameters::get_shared_pointers<Features::FaultModels::Temperature::Interface>(const std::string &name,
                                                                                 std::vector<std::shared_ptr<Features::FaultModels::Temperature::Interface> > &vector);

  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_shared_pointers<Features::FaultModels::Composition::Interface>(const std::string &name,
                                                                                 std::vector<std::shared_ptr<Features::FaultModels::Composition::Interface> > &vector);
  /**
  * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
  * Note that the variable with this name has to be loaded before this function is called.
  */
  template bool
  Parameters::get_shared_pointers<Features::FaultModels::Grains::Interface>(const std::string &name,
                                                                            std::vector<std::shared_ptr<Features::FaultModels::Grains::Interface> > &vector);





} // namespace WorldBuilder

