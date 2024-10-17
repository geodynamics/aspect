#ifndef FASTSCAPELIB_CONFIG_HPP
#define FASTSCAPELIB_CONFIG_HPP

#define FASTSCAPELIB_VERSION_MAJOR 0
#define FASTSCAPELIB_VERSION_MINOR 1
#define FASTSCAPELIB_VERSION_PATCH 3

#include <string>

namespace fastscapelib
{
    namespace version
    {
        constexpr int version_major = FASTSCAPELIB_VERSION_MAJOR;
        constexpr int version_minor = FASTSCAPELIB_VERSION_MINOR;
        constexpr int version_patch = FASTSCAPELIB_VERSION_PATCH;
        static const std::string version_str = std::to_string(version_major) + "."
                                               + std::to_string(version_minor) + "."
                                               + std::to_string(version_patch);
    }
}
#endif
