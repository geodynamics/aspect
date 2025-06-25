/*
  Copyright (C) 2017 - 2018 by the authors of the ASPECT code.

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
