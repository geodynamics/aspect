// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2017-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_citation_info_h
#define _aspect_citation_info_h

#include <string>
#include <iostream>


namespace aspect
{
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

}


#endif
