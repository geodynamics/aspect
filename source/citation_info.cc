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

#include <aspect/citation_info.h>
#include <aspect/revision.h>

#include <iostream>
#include <set>

namespace aspect
{
  namespace CitationInfo
  {
    std::set<std::string> citation_ids;

    const std::string get_url_part ()
    {
      // version:
      std::string url = "citing.html?ver=";
      url += ASPECT_PACKAGE_VERSION;

      // all ids:
      for (const auto &id : citation_ids)
        url += "&" + id + "=1";

      // sha1:
      url += "&sha=";
      url += ASPECT_GIT_SHORTREV;

      // src:
      url += "&src=code";

      return url;
    }

    void add (const std::string &id)
    {
      citation_ids.insert(id);
    }
  }
}
