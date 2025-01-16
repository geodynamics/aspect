/**
 * @file
 * @brief Defines a number of constants used in this library.
 */

#pragma once

namespace fastscapelib
{

    template <class T = double>
    struct numeric_constants
    {
        static constexpr T EARTH_RADIUS_METERS = 6.371e6;
    };

}  // namespace fastscapelib
