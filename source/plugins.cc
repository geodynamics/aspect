// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2024-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/plugins.h>

namespace aspect
{
  namespace Plugins
  {
    void
    InterfaceBase::initialize ()
    {}



    void
    InterfaceBase::update ()
    {}



    void
    InterfaceBase::
    declare_parameters (dealii::ParameterHandler &)
    {}



    void
    InterfaceBase::parse_parameters (dealii::ParameterHandler &)
    {}


    void
    InterfaceBase::save (std::map<std::string,std::string> &) const
    {}



    void
    InterfaceBase::load (const std::map<std::string,std::string> &)
    {}
  }
}
