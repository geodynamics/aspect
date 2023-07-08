/**
 * @file
 * @brief Defines a number of constants used in this library.
 */

#pragma once

namespace fastscapelib
{

    namespace consts
    {

        /*
         * @brief Row index offsets (from a central grid point) of D8
         *  neighbour directions.
         *
         *  The first element correspond to the central grid point itself.
         */
        const short d8_row_offsets[9] = { 0, -1, -1, 0, 1, 1, 1, 0, -1 };


        /**
         * @brief Column index offsets (from a central grid point) of D8
         * neighbour directions.
         *
         * The first element correspond to the central grid point itself.
         */
        const short d8_col_offsets[9] = { 0, 0, -1, -1, -1, 0, 1, 1, 1 };

    }  // namespace consts

}  // namespace fastscapelib
