# This python module contains utility and helper functions 
# intended to make it easier to use pyvista to plot ASPECT solutions.
# For examples of how this module is used, refer to the plot_figure.py 
# scripts in some cookbook and benchmark directories.

from typing import Any, Sequence

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from cmcrameri import cm

__all__ = ["camera_translate","render_mesh_2d"]

def camera_translate(plotter: pv.Plotter, translation_delta: list[float])-> None:
    """Translates the position and focal point of a pyvista camera. 
    
    Moves the camera position and the focal point by the same 3D offset,
    which shifts the rendered image without any rotating or zooming. 

    This is useful for managing whitespace when plotting 2D images with scalar
    bars and titles.

    Args:
        plotter: pyvista.Plotter
            The plotter whose camera to move. The camera's position is modified in place.
        translation_delta: list of floats, ``[dx,dy,dz]``

    Returns:
        None
    """


    pos, focal, view_up = plotter.camera_position
    delta = np.asarray(translation_delta, dtype=float)
    if delta.shape != (3,):
        raise ValueError(f"translation_delta must have 3 components, got shape {delta.shape}")
    plotter.camera_position = [tuple(np.asarray(pos) + delta),
                            tuple(np.asarray(focal) + delta),
                            view_up]

# This dictionary can be updated with more entries for fields as needed
SCALARS_TO_FULLNAME = {
    "T":"Temperature (K)",
}

def render_mesh_2d(
              mesh: pv.pv.DataSet, 
              scalars: str, 
              clim:list[float] | None = None, 
              cmap:plt.colors.ListedColormap | None = cm.vik, 
              arrows_enabled:bool | None = True,
              arrows_field:str | None ="velocity",
              arrows_cmap:plt.colors.ListedColormap | None = cm.imola,
              arrows_scale:float | None = 0.08,
              arrows_normalize:bool | None = None,
              scalarbar_enable:bool | None = True, 
              axisFormat:str | None = "%.0f",
              axisPadding:float | None = 0.015,
              grid_enabled:bool | None = False,
              bounds_enable:bool | None = True,
              sargs: dict[str, Any] | None = None,
              window_size:Sequence[int] = (1024,1024),
              x_title:str | None = "X Axis(m)",
              y_title:str | None = "Y Axis(m)",
              n_x_labels:int | None = 6,
              n_y_labels:int | None = 6,
              background_color="white")-> np.ndarray:
    """Returns a 2D render of a pyvista mesh. Rendering is done offscreen.
    A scalar bar and direction arrow glyphs (indicating velocity, for example) 
    can also be rendered. This function is especially useful when making figures
    composed of separate panels, since you can quickly generate a single panel using
    this function.

    Args:
        mesh: PyVista mesh (e.g. a loaded ``.pvtu`` solution). Must contain the
            ``scalars`` field, and ``arrows_field`` if arrows are enabled.
        scalars: Name of the mesh field to color the surface with.
        clim: ``[min, max]`` color range for the scalar field. By default this value is
            computed from the data's min and max.
        cmap: Colormap for the scalar field. In general 
            this should be a colormap from the cmcrameri package.
        arrows_enabled: Whether or not to draw arrows.
        arrows_field: Name of the vector field the arrows represent. By default,
            arrows are set to render the velocity of the solution.
        arrows_cmap: Colormap for the arrows.
        arrows_scale: Maximum arrow length as a fraction of the domain size. 
            The default is 0.08. 
        arrows_median_normalize: Optional vector array (e.g. the velocity of
            another run) whose median magnitude normalizes arrow sizes, so
            arrows from two panels with different magnitudes compare visually.
        scalarbar_enable: Whether to show the scalar bar.
        sargs: Optional dict of scalar-bar arguments, passed through to pyvista's
            ``scalar_bar_args`` and overrides the pyvista defaults. The title of 
            the scalar bar is set to the full name of the field. The mapping between
            field name and full name is found in the SCALARS_TO_FULLNAME dict. For 
            example, the field `T` in the output data maps to the string `Temperature (K)`.
        window_size: ``(width, height)`` of the rendered image in pixels.
        x_title: Label for the x axis.
        y_title: Label for the y axis.
        n_x_labels: Number of x-axis tick labels.
        n_y_labels: Number of y-axis tick labels.
        background_color: Background color of the rendered image.
        bounds_enable: Draw the bounds around the mesh or not.

    Returns:
        np.ndarray:
            The image is represented as a numpy array.

    """
   
    # These variables are used to scale the arrow glyphs and
    # scalar bar ticks properly. In some examples,
    # like convection-box example, the mesh has a size of 1 so this isn't needed.
    bounds_original = mesh.bounds
    size = max(mesh.bounds[1] - mesh.bounds[0], mesh.bounds[3] - mesh.bounds[2])

    plotter = pv.Plotter(off_screen=True)
    plotter.window_size = window_size
    plotter.set_background(background_color)

    if (clim==None):
        clim = [mesh[scalars].min(), mesh[scalars].max()]

    # Make room for the scalar bar if there is one
    if (scalarbar_enable == True):
        mesh = mesh.scale(0.80)

    # If drawing grid, extract the grid.
    if (grid_enabled):
        mesh_grid = mesh.extract_all_edges()
        if (scalarbar_enable == True):
            mesh_grid.scale(0.80)

    if (arrows_enabled == True):
        # median velocity for arrow scaling
        if (arrows_normalize is None):
            if arrows_field not in mesh.point_data and arrows_field not in mesh.cell_data:
                available = ", ".join(list(mesh.point_data) + list(mesh.cell_data)) or "(none)"
                raise ValueError(
                    f"field to plot as arrows {arrows_field!r} not found in mesh data; "
                    f"available fields: {available}"
                )
            arrows_variable = mesh[arrows_field]
        else:
            arrows_variable = arrows_normalize
        arrows_mag = np.linalg.norm(arrows_variable, axis=1)
        arrows_median = float(np.median(arrows_mag[arrows_mag > 0])) if np.any(arrows_mag > 0) else 1.0
        arrows = mesh.glyph(scale=arrows_field, factor=(arrows_scale*size) / arrows_median, orient=arrows_field, tolerance=0.03)
        glyph_actor = plotter.add_mesh(arrows, cmap=arrows_cmap, show_scalar_bar=False,clim=[0,np.max(arrows_mag)])
        # Setting use_bounds to false here allows the camera to cleanly crop out the arrows
        # when using plotter.camera.tight(). This is noticeable when not drawing
        # scalar bars.
        glyph_actor.use_bounds = False

    

    # This enables you to pass in a different value 
    # for scalar bar arguments.
    if (sargs is None):
        sargs = dict(
            title=SCALARS_TO_FULLNAME.get(scalars),
            height=0.13,
            width=0.75,
            position_x=0.125,
            position_y=0.01,
            title_font_size=40,
            label_font_size=32,
        )

    mesh_actor = plotter.add_mesh(mesh, scalars=scalars, cmap=cmap, clim=clim,
                    lighting=False, show_scalar_bar=scalarbar_enable,scalar_bar_args=sargs)
    if (grid_enabled):
        grid_actor = plotter.add_mesh(mesh_grid, color="ivory_black", style="wireframe",line_width=0.3, opacity=1)
    
    # straight-down 2d camera; parallel_scale fits the render exactly onto the domain bounds
    plotter.camera_position = "xy"
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = size / 2

    if (bounds_enable==True):
        # shows the axes labels on the actual mesh (show_axes shows a widget in 3d space)
        plotter.show_bounds(
        use_3d_text=False,use_2d=True, # These parameters make it easier to plot 2D data
        axes_ranges=bounds_original, # We scaled the mesh earlier, so we need to set the axes ranges to their true values
        location='outer',
        ticks='outside',
        xtitle=x_title,ytitle=y_title,
        grid=None,
        fmt=axisFormat,
        n_xlabels=n_x_labels,
        n_ylabels=n_y_labels,
        padding=axisPadding
        )

    # If drawing the scalar bar, move the camera upwards to
    # make room for the scalar bar under it. Need to move the
    # camera after setting camera_position to "xy", not before.
    if (scalarbar_enable):
        if (sargs.get("vertical")==False):
            camera_translate(plotter,[0,-0.09*size,0])
        else:
            camera_translate(plotter,[-0.09*size,0,0])
    plotter.camera.tight(padding=0.0)
    img = plotter.screenshot(return_img=True)
    plotter.close()
    return img
