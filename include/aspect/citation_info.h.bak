/*
  Copyright (C) 2017 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_citation_info_h
#define _aspect_citation_info_h

#include <string>
#include <iostream>

/**
 * A namespace to provide information to the user about stuff to cite.
 */
namespace CitationInfo
{
  /**
   * Get the URL in the format "citing.html?(parameters)" that describes how to
   * cite ASPECT based on the current model you are running.
   */
  const std::string get_url_part ();

  /**
   * Add the paper identified by the given id to the currently used list of
   * papers. See citing.html for the list of ids. For specific features inside
   * ASPECT that have associated publications, call this function if the
   * feature is used in the current computation. For example, if the
   * computation requires melt migration, call <tt>add("melt")</tt>.
   */
  void add (const std::string &id);

  /**
   * Print the info text containing the citation info into the given
   * stream.
   */
  template <class Stream>
  void print_info_block (Stream &stream)
  {
    stream << "-----------------------------------------------------------------------------\n"
           << "-- For information on how to cite ASPECT, see:\n"
           << "--   https://aspect.geodynamics.org/" << get_url_part() << "\n"
           << "-----------------------------------------------------------------------------"
           << std::endl;
  }
}

#endif
