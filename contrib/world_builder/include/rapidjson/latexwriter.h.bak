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

#ifndef RAPIDJSON_LATEXWRITER_H_
#define RAPIDJSON_LATEXWRITER_H_

#include "writer.h"

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


//! Writer with indentation and spacing.
/*!
    \tparam OutputStream Type of output os.
    \tparam SourceEncoding Encoding of source string.
    \tparam TargetEncoding Encoding of output stream.
    \tparam StackAllocator Type of allocator for allocating memory of stack.
*/
template<typename OutputStream, typename SourceEncoding = UTF8<>, typename TargetEncoding = UTF8<>, typename StackAllocator = CrtAllocator, unsigned writeFlags = kWriteDefaultFlags>
class LatexWriter : public Writer<OutputStream, SourceEncoding, TargetEncoding, StackAllocator, writeFlags>
{
  public:
    typedef Writer<OutputStream, SourceEncoding, TargetEncoding, StackAllocator, writeFlags> Base;
    typedef typename Base::Ch Ch;

    //! Constructor
    /*! \param os Output stream.
        \param allocator User supplied allocator. If it is null, it will create a private one.
        \param levelDepth Initial capacity of stack.
    */
    explicit LatexWriter(OutputStream &os, StackAllocator *allocator = 0, size_t levelDepth = Base::kDefaultLevelDepth) :
      Base(os, allocator, levelDepth)
    {
      level_type.resize(0);
      array_number.resize(0);
      path.resize(0);
    }


    explicit LatexWriter(StackAllocator *allocator = 0, size_t levelDepth = Base::kDefaultLevelDepth) :
      Base(allocator, levelDepth)
    {
      level_type.resize(0);
      array_number.resize(0);
      path.resize(0);
    }

#if RAPIDJSON_HAS_CXX11_RVALUE_REFS
    LatexWriter(LatexWriter &&rhs) :
      Base(std::forward<LatexWriter>(rhs))
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
      std::string begin = "";
      if (level_type.size() != 0 && level_type.back() == 1)
        {
          // this is a properties starting, so first clear all the itemize
          // the lvel_type.push_back(1) has already been done by the key: properties
          for (unsigned int i = 0; i < itemize_open; ++i)
            {
              begin += "\\end{itemize}";
            }
          itemize_open = 0;
          if (skip_next_push_back == false)
            {
              if (section_level <= 2)
                begin += "\\section{(" + std::to_string(section_level) + ") " + get_path() + "}";
              else if (section_level <= 5)
                begin += "\\subsection{(" + std::to_string(section_level) + ") " + get_path() + "}";
              else if (section_level <= 8)
                begin += "\\subsubsection{(" + std::to_string(section_level) + ") " + get_path() + "}";
              else
                begin += "\\paragraph{(" + std::to_string(section_level) + ") " + get_path() + "}";

              open_new_itemize = true;
              level_type.push_back(0);
            }
          else
            {
              skip_next_push_back = false;
            }
        }
      else if (level_type.size() != 0 && level_type.back() == 2)
        {
          // this is a array starting, so first clear all the itemize
          // the lvel_type.push_back(1) has already been done by the key: properties
          for (unsigned int i = 0; i < itemize_open; ++i)
            {
              begin += "\\end{itemize}";
            }
          itemize_open = 0;
          unsigned int number = array_number.back() + 1;
          array_number.back() = number;
          path.push_back(std::to_string(array_number.back()));
          level_type.push_back(0);
          // have arrrays (/oneOf/1) as new section level
          section_level++;

          if (section_level <= 2)
            begin += "\\section{(" + std::to_string(section_level) + ") " + get_path() + "}";
          else if (section_level <= 5)
            begin += "\\subsection{(" + std::to_string(section_level) + ") " + get_path() + "}";
          else if (section_level <= 8)
            begin += "\\subsubsection{(" + std::to_string(section_level) + ") " + get_path() + "}";
          else
            begin += "\\paragraph{(" + std::to_string(section_level) + ") " + get_path() + "}";

          open_new_itemize = true;
        }
      else
        {
          for (unsigned int i = 0; i < itemize_open; ++i)
            {
              begin += "\\end{itemize}";
            }
          itemize_open = 0;

          //begin += "C";
          if (section_level <= 2)
            begin += "\\section{(" + std::to_string(section_level) + ") " + get_path() + "}";
          else if (section_level <= 5)
            begin += "\\subsection{(" + std::to_string(section_level) + ") " + get_path() + "}";
          else if (section_level <= 8)
            begin += "\\subsubsection{(" + std::to_string(section_level) + ") " + get_path() + "}";
          else
            begin += "\\paragraph{(" + std::to_string(section_level) + ") " + get_path() + "}";

          open_new_itemize = true;
          if (level_type.size() > 0 && (level_type.back() != 3 || level_type.size() == 0))
            level_type.push_back(0);
        }

      Base::WriteString(begin.c_str(), static_cast<rapidjson::SizeType>(begin.size()), false);
      Base::os_->Put('\n');

      return true;
    }

    bool Key(const Ch *str, SizeType /*length*/, bool /*copy = false*/)
    {
      std::string item = "";


      std::string key(str);
      if (key == "properties")
        {
          if (open_new_itemize == true)
            {
              item += "\\begin{itemize}[leftmargin=" + std::to_string(section_level) + "em]";
              itemize_open++;
            }
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
          if (open_new_itemize == true)
            {
              item += "\\begin{itemize}[leftmargin=" + std::to_string(section_level) + "em]";
              itemize_open++;
            }
          section_level++;
          level_type.push_back(3);
          path.push_back(key);
        }
      else if (level_type.size() > 0 && level_type.back() == 1)
        {
          // the level just below properties, these are the sections
          // add to the path

          section_level++;
          path.push_back(key);
        }
      else
        {
          if (open_new_itemize == true)
            {
              item += "\\begin{itemize}[leftmargin=" + std::to_string(section_level) + "em]";
              itemize_open++;

            }
          item += "\\item {\\bf " + std::string(str) + "}: ";
        }

      open_new_itemize = false;
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
      if (itemize_open > 0)
        {
          end = "\\end{itemize}";
          itemize_open--;
        }
      itemize_open = 0;
      // If this assert fails it could be an indicator that
      // the world builder library is not fuly linked to the
      // appllication.
      if (level_type.size() > 0 && (level_type.back() == 0 || level_type.back() == 3))
        section_level = section_level > 0 ? section_level-1 : 0;


      if (level_type.size() > 0 && path.size() > 0 &&
          ( level_type.back() == 3 || level_type.back() == 0))
        {
          path.pop_back();
        }

      if (level_type.size() > 0)
        level_type.pop_back();

      Base::WriteString(end.c_str(), static_cast<rapidjson::SizeType>(end.length()), false);

      return true;
    }

    bool StartArray()
    {
      std::string begin = "";
      if (level_type.size() != 0 && level_type.back() == 2)
        {
          // this is a properties starting, so first clear all the itemize
          // the lvel_type.push_back(1) has already been done by the key: properties
          for (unsigned int i = 0; i < itemize_open; ++i)
            {
              begin += "\\end{itemize}";
            }
          itemize_open = 0;

          open_new_itemize = true;
          array_number.push_back(0);
        }
      else
        {
          begin = "[";
          small_array = true;
          first_small_array = true;
        }

      Base::WriteString(begin.c_str(), static_cast<rapidjson::SizeType>(begin.size()), false);

      if (small_array == false)
        Base::os_->Put('\n');

      return true;
    }

    bool EndArray(SizeType /*memberCount = 0*/)
    {
      std::string end = "";
      if (level_type.size() != 0 && level_type.back() == 2)
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
        }

      Base::WriteString(end.c_str(), static_cast<rapidjson::SizeType>(end.size()), false);

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
        \note When using LatexWriter::RawValue(), the result json may not be indented correctly.
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

    unsigned int itemize_open = 0;
    std::vector<unsigned int> level_type;
    std::vector<unsigned int> array_number;
    bool skip_next_push_back = false;
    bool small_array = false;
    bool first_small_array = false;
    bool open_new_itemize = false;
    unsigned int section_level = 0;

  private:
    std::vector<std::string> path;
    // Prohibit copy constructor & assignment operator.
    LatexWriter(const LatexWriter &);
    LatexWriter &operator=(const LatexWriter &);
};

RAPIDJSON_NAMESPACE_END

#if defined(__clang__)
RAPIDJSON_DIAG_POP
#endif

#ifdef __GNUC__
RAPIDJSON_DIAG_POP
#endif

#endif // RAPIDJSON_RAPIDJSON_H_
