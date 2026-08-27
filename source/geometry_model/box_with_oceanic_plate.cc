#include <aspect/geometry_model/box_with_oceanic_plate.h>

#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>

namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    BoxWithOceanicPlate<dim>::
    set_oceanic_plate_boundary_indicators (parallel::distributed::Triangulation<dim> &triangulation) const
    {
      // These are the same two extra indicators returned by
      // get_used_boundary_indicators()/get_symbolic_boundary_names_map():
      // the first two indicators not already used by TwoMergedBoxes. They are
      // exposed to the user as "left lithosphere"/"right lithosphere" (the
      // thin plate strip), while the base class' own indicators 2*dim/2*dim+1
      // (originally called "left/right lithosphere" by TwoMergedBoxes) are
      // exposed as "left upper mantle"/"right upper mantle" instead -- see
      // get_symbolic_boundary_names_map().
      const types::boundary_id left_plate_id  = 2*dim + 2*(dim-1);
      const types::boundary_id right_plate_id = 2*dim + 2*(dim-1) + 1;

      for (const auto &cell : triangulation.active_cell_iterators())
        {
          // face(0) is the left boundary (x = minimum), face(1) is the
          // right boundary (x = maximum), for both 2D and 3D box geometries.
          if (cell->face(0)->at_boundary())
            if (cell->face(0)->vertex(cell->face(0)->n_vertices()-1)[dim-1] > height_oceanic)
              cell->face(0)->set_boundary_id (left_plate_id);

          if (cell->face(1)->at_boundary())
            if (cell->face(1)->vertex(cell->face(1)->n_vertices()-1)[dim-1] > height_oceanic)
              cell->face(1)->set_boundary_id (right_plate_id);
        }
    }



    template <int dim>
    void
    BoxWithOceanicPlate<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &total_coarse_grid) const
    {
      // Build the standard mantle/lithosphere mesh and boundary indicators
      // exactly the way the base class (TwoMergedBoxes) does. This also
      // connects the base class' own post-refinement signal, which keeps
      // reapplying its "left/right lithosphere" (now renamed "left/right
      // upper mantle") split after every refinement.
      TwoMergedBoxes<dim>::create_coarse_mesh (total_coarse_grid);

      // Now carve the thin plate strip -- exposed to the user as "left/right
      // lithosphere" -- out of the top of the base class' indicators. Our
      // check is based purely on absolute height, so it is safe regardless
      // of whether it runs before or after the base class' own signal
      // handler on subsequent refinements.
      set_oceanic_plate_boundary_indicators (total_coarse_grid);

      total_coarse_grid.signals.post_refinement.connect (
        [&]()
      {
        this->set_oceanic_plate_boundary_indicators (total_coarse_grid);
      });
    }



    template <int dim>
    std::set<types::boundary_id>
    BoxWithOceanicPlate<dim>::
    get_used_boundary_indicators () const
    {
      std::set<types::boundary_id> indicators = TwoMergedBoxes<dim>::get_used_boundary_indicators();
      indicators.insert (2*dim + 2*(dim-1));
      indicators.insert (2*dim + 2*(dim-1) + 1);
      return indicators;
    }



    template <int dim>
    std::map<std::string,types::boundary_id>
    BoxWithOceanicPlate<dim>::
    get_symbolic_boundary_names_map () const
    {
      std::map<std::string,types::boundary_id> names = TwoMergedBoxes<dim>::get_symbolic_boundary_names_map();

      // The base class calls the whole depth range from the surface down to
      // "Lithospheric thickness" (indicators 2*dim / 2*dim+1) "left/right
      // lithosphere". We rename that zone -- which is now only the part
      // *below* the thin plate strip -- to "left/right upper mantle", and
      // use "left/right lithosphere" instead for the new, thinner strip at
      // the very top (indicators 2*dim+2*(dim-1) / +1), where the rigid
      // plate itself is meant to sit. "left"/"right" (the mantle below
      // "Lithospheric thickness") are unchanged.
      names.erase ("left lithosphere");
      names.erase ("right lithosphere");
      names["left upper mantle"]  = 2*dim;
      names["right upper mantle"] = 2*dim + 1;
      names["left lithosphere"]   = 2*dim + 2*(dim-1);
      names["right lithosphere"]  = 2*dim + 2*(dim-1) + 1;

      return names;
    }



    template <int dim>
    void
    BoxWithOceanicPlate<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      // Declare all of the parameters of the base class first (X/Y/Z
      // extent, Lithospheric thickness, repetitions, periodicity, ...).
      TwoMergedBoxes<dim>::declare_parameters (prm);

      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box with lithosphere boundary indicators");
        {
          prm.declare_entry ("Oceanic plate thickness", "0.05",
                             Patterns::Double (0.),
                             "The thickness of the strip carved out of the top of the "
                             "'left upper mantle'/'right upper mantle' boundary indicators "
                             "(the part of the domain within 'Lithospheric thickness' of "
                             "the surface) and given the separate 'left lithosphere'/'right "
                             "lithosphere' boundary indicators instead. This does not need "
                             "to coincide exactly with a mesh cell face: the actual location "
                             "of the new boundary indicator follows the nearest cell edge "
                             "on the initial (coarse) mesh to this depth. Must be smaller "
                             "than 'Lithospheric thickness'. Units: \\si{\\meter}.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    BoxWithOceanicPlate<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      // Let the base class parse everything it needs. This fully sets up
      // TwoMergedBoxes' internal state (extents, mesh repetitions,
      // the mantle/lithosphere split height, ...), so every function we
      // do *not* override (depth(), get_extents(), maximal_depth(),
      // representative_point(), natural coordinates, ...) keeps working
      // completely unchanged.
      TwoMergedBoxes<dim>::parse_parameters (prm);

      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box with lithosphere boundary indicators");
        {
          const double y_extent          = prm.get_double ("Y extent");
          const double y_origin          = prm.get_double ("Box origin Y coordinate");
          const double thickness_lith    = prm.get_double ("Lithospheric thickness");
          const double thickness_oceanic = prm.get_double ("Oceanic plate thickness");

          AssertThrow (thickness_oceanic < thickness_lith,
                       ExcMessage ("The 'Oceanic plate thickness' must be smaller than "
                                   "the 'Lithospheric thickness': the 'left/right "
                                   "lithosphere' boundary indicators are carved out of the "
                                   "top part of the 'left/right upper mantle' boundary "
                                   "indicators."));

          // Top of the domain minus the requested oceanic plate thickness.
          // (This mirrors how the base class computes height_lith as the
          // top of the mantle box, i.e. origin + extent - lithospheric
          // thickness, ignoring initial topography for the purposes of
          // assigning boundary indicators.)
          height_oceanic = y_origin + y_extent - thickness_oceanic;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace GeometryModel
  {
    ASPECT_REGISTER_GEOMETRY_MODEL(BoxWithOceanicPlate,
                                   "box with lithosphere and oceanic plate boundary indicators",
                                   "A geometry model that behaves exactly like "
                                   "``box with lithosphere boundary indicators'', but "
                                   "additionally splits the top part of its "
                                   "``left lithosphere''/``right lithosphere'' boundary "
                                   "indicators into a separate pair of indicators, and "
                                   "renames boundaries so that, from the surface down: "
                                   "``left/right lithosphere'' now covers only the thin "
                                   "plate strip given by the ``Oceanic plate thickness'' "
                                   "parameter, ``left/right upper mantle'' covers the "
                                   "remainder of ``Lithospheric thickness'' below that "
                                   "strip (this is what the base model calls ``left/right "
                                   "lithosphere''), and ``left''/``right'' cover the mantle "
                                   "below that, unchanged from the base model. All "
                                   "parameters are declared in the same ``Box with "
                                   "lithosphere boundary indicators'' subsection as the "
                                   "base model, with one addition, ``Oceanic plate "
                                   "thickness''. This allows, for example, prescribing an "
                                   "inflow traction on a thin oceanic plate at the top of "
                                   "the lithospheric boundary (``left/right lithosphere''), "
                                   "a different velocity boundary condition (e.g. tangential "
                                   "or zero velocity) on the upper mantle below it "
                                   "(``left/right upper mantle''), and an outflow traction "
                                   "on the underlying mantle boundary below that "
                                   "(``left''/``right'').")
  }
}
