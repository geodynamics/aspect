// Tencent is pleased to support the open source community by making RapidJSON available.
//
// Copyright (C) 2015 THL A29 Limited, a Tencent company, and Milo Yip. All rights reserved.
//
// Licensed under the MIT License (the "License"); you may not use this file except
// in compliance with the License. You may obtain a copy of the License at
//
// http://opensource.org/licenses/MIT
//
// Unless required by applicable law or agreed to in writing, software distributed
// under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
// CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.

#ifndef RAPIDJSON_MYSTWRITER_H_
#define RAPIDJSON_MYSTWRITER_H_

#include "writer.h"

#include "assert.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#ifdef __GNUC__
RAPIDJSON_DIAG_PUSH
RAPIDJSON_DIAG_OFF(effc++)
#endif

#if defined(__clang__)
RAPIDJSON_DIAG_PUSH
RAPIDJSON_DIAG_OFF(c++98-compat)
#endif

RAPIDJSON_NAMESPACE_BEGIN

#define MAX_PATH_LEVEL 25

//! Writer with indentation and spacing.
/*!
    \tparam OutputStream Type of output os.
    \tparam SourceEncoding Encoding of source string.
    \tparam TargetEncoding Encoding of output stream.
    \tparam StackAllocator Type of allocator for allocating memory of stack.
*/
template<typename OutputStream, typename SourceEncoding = UTF8<>, typename TargetEncoding = UTF8<>, typename StackAllocator = CrtAllocator, unsigned writeFlags = kWriteDefaultFlags>
class MySTWriter : public Writer<OutputStream, SourceEncoding, TargetEncoding, StackAllocator, writeFlags>
{
  public:
    typedef Writer<OutputStream, SourceEncoding, TargetEncoding, StackAllocator, writeFlags> Base;
    typedef typename Base::Ch Ch;

    //! Constructor
    /*! \param os Output stream.
        \param allocator User supplied allocator. If it is null, it will create a private one.
        \param levelDepth Initial capacity of stack.
    */
    explicit MySTWriter(OutputStream &os, bool _make_open = false, StackAllocator *allocator = 0, size_t levelDepth = Base::kDefaultLevelDepth) :
      Base(os, allocator, levelDepth)
    {
      make_open = _make_open;
      level_type.resize(0);
      array_number.resize(0);
      path.resize(0);
    }


    explicit MySTWriter(StackAllocator *allocator = 0, size_t levelDepth = Base::kDefaultLevelDepth) :
      Base(allocator, levelDepth)
    {
      level_type.resize(0);
      array_number.resize(0);
      path.resize(0);
    }

#if RAPIDJSON_HAS_CXX11_RVALUE_REFS
    MySTWriter(MySTWriter &&rhs) :
      Base(std::forward<MySTWriter>(rhs))
    {
      level_type.resize(0);
      array_number.resize(0);
      path.resize(0);
    }
#endif

    /*! @name Implementation of Handler
        \see Handler
    */
    //@{

    bool Null()
    {
      return Base::EndValue(Base::WriteNull());
    }
    bool Bool(bool b)
    {
      Base::EndValue(Base::WriteBool(b));
      Base::os_->Put('\n');
      return true;
    }
    bool Int(int i)
    {
      Base::EndValue(Base::WriteInt(i));
      Base::os_->Put('\n');
      return true;
    }
    bool Uint(unsigned u)
    {
      Base::EndValue(Base::WriteUint(u));
      Base::os_->Put('\n');
      return true;
    }
    bool Int64(int64_t i64)
    {
      Base::EndValue(Base::WriteInt64(i64));
      Base::os_->Put('\n');
      return true;
    }
    bool Uint64(uint64_t u64)
    {
      Base::EndValue(Base::WriteUint64(u64));
      Base::os_->Put('\n');
      return true;
    }
    bool Double(double d)
    {
      Base::EndValue(Base::WriteDouble(d));
      Base::os_->Put('\n');
      return true;
    }

    bool RawNumber(const Ch *str, SizeType length, bool /*copy = false*/)
    {
      RAPIDJSON_ASSERT(str != 0);
      Base::EndValue(Base::WriteString(str, length, false));
      Base::os_->Put('\n');
      return true;
    }

    bool String(const Ch *str, SizeType length, bool /*copy = false*/)
    {
      RAPIDJSON_ASSERT(str != 0);
      if (small_array == true)
        {
          if (first_small_array == true)
            first_small_array = false;
          else
            {
              Base::os_->Put(',');
              Base::os_->Put(' ');
            }
        }
      Base::EndValue(Base::WriteString(str, length, false));

      if (small_array == false)
        Base::os_->Put('\n');
      return true;
    }

#if RAPIDJSON_HAS_STDSTRING
    bool String(const std::basic_string<Ch> &str)
    {
      return String(str.data(), SizeType(str.size()));
    }
#endif

    bool StartObject()
    {
      bool wrote_dropdown = false;
      std::string begin = "";
      if (level_type.size() != 0 && level_type.back() == 1)
        {
          // this is a properties starting, so first clear all the itemize
          // the lvel_type.push_back(1) has already been done by the key: properties

          if (skip_next_push_back == false)
            {
              for (size_t level =MAX_PATH_LEVEL-path.size(); level != 0; level--)
                {
                  begin += ":";
                }
              begin += "{dropdown} " + get_path();

              wrote_dropdown =  true;
              level_type.push_back(0);
            }
          else
            {
              skip_next_push_back = false;
            }
        }
      else if (level_type.size() != 0 && level_type.back() == 2)
        {
          unsigned int number = array_number.back() + 1;
          array_number.back() = number;
          path.push_back(std::to_string(array_number.back()));
          level_type.push_back(0);

          for (size_t level =MAX_PATH_LEVEL-path.size(); level != 0; level--)
            {
              begin += ":";
            }
          begin += "{dropdown} " + get_path();


          wrote_dropdown =  true;
        }
      else if (level_type.size() != 0 && level_type.back() == 5)
        {
          unsigned int number = array_number.back() + 1;
          array_number.back() = number;
          path.push_back(std::to_string(array_number.back()));
          level_type.push_back(0);

          for (size_t level =MAX_PATH_LEVEL-path.size(); level != 0; level--)
            {
              begin += ":";
            }
          begin += "{dropdown} " + get_path();


          wrote_dropdown =  true;
        }
      else
        {
          for (size_t level =MAX_PATH_LEVEL-path.size(); level != 0; level--)
            {
              begin += ":";
            }
          begin += "{dropdown} " + get_path();


          wrote_dropdown =  true;
          if (level_type.size() > 0 && (level_type.back() != 3 || level_type.size() == 0))
            level_type.push_back(0);
        }

      Base::WriteString(begin.c_str(), static_cast<rapidjson::SizeType>(begin.size()), false);

      Base::os_->Put('\n');
      if (level_type.size() == 0 && !make_open)
        {
          std::string string_open = ":open:";
          Base::WriteString(string_open.c_str(), static_cast<rapidjson::SizeType>(string_open.size()), false);
          Base::os_->Put('\n');
        }

      if (wrote_dropdown)
        {
          if (make_open)
            {
              std::string open = ":open:";
              Base::WriteString(open.c_str(), static_cast<rapidjson::SizeType>(open.size()), false);
              Base::os_->Put('\n');
              std::string name = ":name: open" + get_path_underscore() + "";
              Base::WriteString(name.c_str(), static_cast<rapidjson::SizeType>(name.size()), false);
              Base::os_->Put('\n');
              Base::os_->Put('\n');
            }
          else
            {
              std::string name = ":name: closed" + get_path_underscore() + "";
              Base::WriteString(name.c_str(), static_cast<rapidjson::SizeType>(name.size()), false);
              Base::os_->Put('\n');
              Base::os_->Put('\n');
            }
        }

      return true;
    }

    bool Key(const Ch *str, SizeType /*length*/, bool /*copy = false*/)
    {
      std::string item = "";


      std::string key(str);
      if (key == "properties")
        {
          level_type.push_back(1);
          skip_next_push_back = true;
        }
      else if (key == "oneOf")
        {
          level_type.push_back(2);
          path.push_back(key);
        }
      else if (key == "items")
        {
          level_type.push_back(3);
          path.push_back(key);
        }
      else if (key == "enum" || key == "required" )
        {
          level_type.push_back(4);
          item += "- **" + std::string(str) + "**:";
        }
      else if (key == "anyOf")
        {
          level_type.push_back(5);
          path.push_back(key);
        }
      else if (level_type.size() > 0 && (level_type.back() == 1 || level_type.back() == 2 || level_type.back() == 5))
        {
          // the level just below properties, these are the sections
          // add to the path
          path.push_back(key);
        }
      else
        {
          item += "- **" + std::string(str) + "**:";
        }

      Base::WriteString(item.c_str(), static_cast<rapidjson::SizeType>(item.length()), false);

      return true;
    }

#if RAPIDJSON_HAS_STDSTRING
    bool Key(const std::basic_string<Ch> &str)
    {
      return Key(str.data(), SizeType(str.size()));
    }
#endif

    bool EndObject(SizeType /*memberCount = 0*/)
    {
      std::string end = "";

      // If this assert fails it could be an indicator that
      // the world builder library is not fuly linked to the
      // appllication.

      if (level_type.size() > 0 && path.size() > 0 &&
          ( level_type.back() == 3 || level_type.back() == 0))
        {

          for (size_t level =MAX_PATH_LEVEL-path.size(); level != 0; level--)
            {
              WBAssert(path.size() < MAX_PATH_LEVEL, "EA: path size is larger than 90: " << path.size() << ", level = " << level);
              end += ':';
            }

          path.pop_back();
        }

      if (level_type.size() > 0)
        {
          level_type.pop_back();
        }

      Base::WriteString(end.c_str(), static_cast<rapidjson::SizeType>(end.length()), false);

      Base::os_->Put('\n');
      Base::os_->Put('\n');
      return true;
    }

    bool StartArray()
    {
      bool wrote_dropdown = false;
      std::string begin = "";
      if (level_type.size() != 0 && (level_type.back() == 2|| level_type.back() == 3|| level_type.back() == 5 ))
        {
          //begin += "level type" + std::to_string(level_type.back());
          // this is a properties starting, so first clear all the itemize
          // the lvel_type.push_back(1) has already been done by the key: properties

          for (size_t level = MAX_PATH_LEVEL-path.size(); level != 0; level--)
            {
              WBAssert(path.size() < MAX_PATH_LEVEL, "SA: path size is larger than 90: " << path.size() << ", level = " << level);
              begin += ":";
            }
          begin += "{dropdown} " + get_path();

          wrote_dropdown = true;
          array_number.push_back(0);
        }
      else
        {
          begin += "[";
          small_array = true;
          first_small_array = true;
        }

      Base::WriteString(begin.c_str(), static_cast<rapidjson::SizeType>(begin.size()), false);

      if (small_array == false)
        Base::os_->Put('\n');


      if (wrote_dropdown)
        {
          if (make_open)
            {
              std::string open = ":open:";
              Base::WriteString(open.c_str(), static_cast<rapidjson::SizeType>(open.size()), false);
              Base::os_->Put('\n');
              std::string name = ":name: open" + get_path_underscore() + "";
              Base::WriteString(name.c_str(), static_cast<rapidjson::SizeType>(name.size()), false);
              Base::os_->Put('\n');
              Base::os_->Put('\n');
            }
          else
            {
              std::string name = ":name: closed" + get_path_underscore() + "";
              Base::WriteString(name.c_str(), static_cast<rapidjson::SizeType>(name.size()), false);
              Base::os_->Put('\n');
              Base::os_->Put('\n');
            }
        }
      return true;
    }

    bool EndArray(SizeType /*memberCount = 0*/)
    {
      std::string end = "";
      if (level_type.size() != 0 && (level_type.back() == 1 || level_type.back() == 2 || level_type.back() == 3 || level_type.back() == 5))
        {
          array_number.pop_back();
          level_type.pop_back();

          if (path.size() > 0)
            path.pop_back();
        }
      else
        {

          end = "]";

          small_array = false;

          if (level_type.size() > 0 && (level_type.back() == 4))// || level_type.back() == 0))
            {
              level_type.pop_back();
            }
        }

      Base::WriteString(end.c_str(), static_cast<rapidjson::SizeType>(end.size()), false);

      Base::os_->Put('\n');
      return true;
    }

    //@}

    /*! @name Convenience extensions */
    //@{

    //! Simpler but slower overload.
    bool String(const Ch *str)
    {
      return String(str, internal::StrLen(str));
    }
    bool Key(const Ch *str)
    {
      return Key(str, internal::StrLen(str));
    }

    //@}

    //! Write a raw JSON value.
    /*!
        For user to write a stringified JSON as a value.

        \param json A well-formed JSON value. It should not contain null character within [0, length - 1] range.
        \param length Length of the json.
        \param type Type of the root of json.
        \note When using MySTWriter::RawValue(), the result json may not be indented correctly.
    */
    bool RawValue(const Ch *json, size_t length, Type /*type*/)
    {
      RAPIDJSON_ASSERT(json != 0);
      return Base::EndValue(Base::WriteRawValue(json, length));
    }

  protected:

    std::string get_path()
    {
      std::string return_path;
      for (auto string : path)
        {
          return_path += "/" + string;
        }
      return return_path == "" ? "/" : return_path;
    }

    std::string get_path_underscore()
    {
      std::string return_path;
      for (auto string : path)
        {
          return_path += "_" + string;
        }
      const char space  = ' ';
      const char dash = '-';
      std::replace(return_path.begin(), return_path.end(), space, dash);
      return return_path == "" ? "_" : return_path;
    }

    std::vector<unsigned int> level_type;
    std::vector<unsigned int> array_number;
    bool make_open = false;
    bool skip_next_push_back = false;
    bool small_array = false;
    bool first_small_array = false;

  private:
    std::vector<std::string> path;
    // Prohibit copy constructor & assignment operator.
    MySTWriter(const MySTWriter &);
    MySTWriter &operator=(const MySTWriter &);
};

RAPIDJSON_NAMESPACE_END

#if defined(__clang__)
RAPIDJSON_DIAG_POP
#endif

#ifdef __GNUC__
RAPIDJSON_DIAG_POP
#endif

#endif // RAPIDJSON_RAPIDJSON_H_
