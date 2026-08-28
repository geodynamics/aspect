// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2023 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/postprocess/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    class PointIsInDomain : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &);
    };
  }
}


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    PointIsInDomain<dim>::execute(TableHandler &)
    {
      // The airy isostasy box domain is 1x1 with max 0.05 topography.
      const Point<dim> point_in_domain(0.1, 0.1);
      const Point<dim> point_not_in_domain(0.1, 2.0);
      const Point<dim> point_not_in_domain_2(2.0, 0.01);

      std::ostringstream screen_text;
      screen_text << std::boolalpha << this->get_geometry_model().point_is_in_domain(point_in_domain) << " ";
      screen_text << std::boolalpha << this->get_geometry_model().point_is_in_domain(point_not_in_domain) << " ";
      screen_text << std::boolalpha << this->get_geometry_model().point_is_in_domain(point_not_in_domain_2);

      return std::pair<std::string, std::string> ("Points lie in domain: ",
                                                  screen_text.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(PointIsInDomain,
                                  "point is in domain",
                                  "A postprocessor that tests whether certain points "
                                  "lie within the domain.")
  }
}
