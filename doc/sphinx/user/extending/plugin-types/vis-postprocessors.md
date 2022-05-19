# Visualization postprocessors

As mentioned in the previous section, one of the postprocessors that are
already implemented in ASPECT is the
[aspect::Postprocess::Visualization][] class that takes the solution and
outputs it as a collection of files that can then be visualized graphically,
see {ref}`sec:viz`21]. The question is which variables to output:
the solution of the basic equations we solve here is characterized by the
velocity, pressure and temperature; on the other hand, we are frequently
interested in derived, spatially and temporally variable quantities such as
the viscosity for the actual pressure, temperature and strain rate at a given
location, or seismic wave speeds.

ASPECT already implements a good number of such
derived quantities that one may want to visualize. On the other hand, always
outputting *all* of them would yield very large output files, and would
furthermore not scale very well as the list continues to grow. Consequently,
as with the postprocessors described in the previous section, what *can* be
computed is implemented in a number of plugins and what *is* computed is
selected in the input parameter file (see {ref}`sec:3.167][]).

Defining visualization postprocessors works in much the same way as for the
other plugins discussed in this section. Specifically, an implementation of
such a plugin needs to be a class that derives from interface classes, should
by convention be in namespace
`aspect::Postprocess::VisualizationPostprocessors`, and is registered using a
macro, here called `ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR`. Like the
postprocessor plugins, visualization postprocessors can derive from class
[aspect::Postprocess::SimulatorAccess][] if they need to know specifics of the
simulation such as access to the material models and to get access to the
introspection facility outlined in {ref}`sec:1.1][]. A typical example is
the plugin that produces the viscosity as a spatially variable field by
evaluating the viscosity function of the material model using the pressure,
temperature and location of each visualization point (implemented in the
`aspect::Postprocess::VisualizationPostprocessors::Viscosity` class). On the
other hand, a hypothetical plugin that simply outputs the norm of the strain
rate $\sqrt{\varepsilon(\mathbf
  u):\varepsilon(\mathbf u)}$ would not need access to anything but the
solution vector (which the plugin&rsquo;s main function is given as an
argument) and consequently is not derived from the
[aspect::Postprocess::SimulatorAccess][] class.[7]

Visualization plugins can come in two flavors:

-   *Plugins that compute things from the solution in a point-wise way:* The
    classes in this group are derived not only from the respective interface
    class (and possibly the [SimulatorAccess][aspect::SimulatorAccess class]
    class) but also from the deal.II class `DataPostprocessor` or any of the
    classes like `DataPostprocessorScalar` or `DataPostprocessorVector`. These
    classes can be thought of as filters: DataOut will call a function in them
    for every cell and this function will transform the values or gradients of
    the solution and other information such as the location of quadrature
    points into the desired quantity to output. A typical case would be if the
    quantity $g(x)$ you want to output can be written as a function
    $g(x) = G(u(x),\nabla u(x), x, ...)$ in a point-wise sense where $u(x)$ is
    the value of the solution vector (i.e., the velocities, pressure,
    temperature, etc) at an evaluation point. In the context of this program
    an example would be to output the density of the medium as a spatially
    variable function since this is a quantity that for realistic media
    depends point-wise on the values of the solution.

    To sum this, slightly confusing multiple inheritance up, visualization
    postprocessors do the following:

    -   If necessary, they derive from
        [aspect::Postprocess::SimulatorAccess][].

    -   They derive from
        [aspect::Postprocess::VisualizationPostprocessors::Interface][]. The
        functions of this interface class are all already implemented as doing
        nothing in the base class but can be overridden in a plugin.
        Specifically, the following functions exist:

        ```{code-block} c++
        class Interface
            {
              public:
                static
                void
                declare_parameters (ParameterHandler &prm);

                virtual
                void
                parse_parameters (ParameterHandler &prm);

                virtual
                void save (std::map<std::string, std::string> &status_strings) const;

                virtual
                void load (const std::map<std::string, std::string> &status_strings);
            };
        ```

    -   They derive from either the `dealii::DataPostprocessor` class, or the
        simpler to use `dealii::DataPostprocessorScalar` or
        `dealii::DataPostprocessorVector` classes. For example, to derive from
        the second of these classes, the following interface functions has to
        be implemented:

        ```{code-block} c++
        class dealii::DataPostprocessorScalar
            {
              public:
                virtual
                void
                compute_derived_quantities_vector
                  (const std::vector<Vector<double> >              &uh,
                   const std::vector<std::vector<Tensor<1,dim> > > &duh,
                   const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                   const std::vector<Point<dim> >                  &normals,
                   const std::vector<Point<dim> >                  &evaluation_points,
                   std::vector<Vector<double> >                    &computed_quantities) const;
            };
        ```

        What this function does is described in detail in the deal.II
        documentation. In addition, one has to write a suitable constructor to
        call `dealii::DataPostprocessorScalar::DataPostprocessorScalar`.

-   *Plugins that compute things from the solution in a cell-wise way:* The
    second possibility is for a class to not derive from
    `dealii::DataPostprocessor` but instead from the
    [aspect::Postprocess::VisualizationPostprocessors::CellDataVectorCreator][]
    class. In this case, a visualization postprocessor would generate and
    return a vector that consists of one element per cell. The intent of this
    option is to output quantities that are not point-wise functions of the
    solution but instead can only be computed as integrals or other
    functionals on a per-cell basis. A typical case would be error estimators
    that do depend on the solution but not in a point-wise sense; rather, they
    yield one value per cell of the mesh. See the documentation of the
    `CellDataVectorCreator` class for more information.

If all of this sounds confusing, we recommend consulting the implementation of
the various visualization plugins that already exist in the
ASPECT sources, and using them as a template.
