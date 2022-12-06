# Categories of parameters

The parameters that can be provided in the input file can roughly be
categorized into the following groups:

-   Global parameters (see {ref}`parameters:global`): These
    parameters determine the overall behavior of the program. Primarily they
    describe things like the output directory, the end time of the simulation,
    or whether the computation should be resumed from a previously saved
    state.

-   Parameters for certain aspects of the numerical algorithm: These describe,
    for example, the specifics of the spatial discretization. In particular,
    this is the case for parameters concerning the polynomial degree of the
    finite element approximation
    ({ref}`parameters:Discretization`), some details about the
    stabilization
    ({ref}`parameters:Discretization/Stabilization_20parameters`),
    and how adaptive mesh refinement is supposed to work
    ({ref}`parameters:Mesh_20refinement`).

-   Parameters that describe certain global aspects of the equations to be
    solved: This includes, for example, a description if certain terms in the
    model should be omitted or not. See
    {ref}`parameters:Formulation` for the list of parameters
    in this category.

-   Parameters that characterize plugins: Certain behaviors of
    ASPECT are described by what we call *plugins* - self-contained parts of
    the code that describe one particular
    aspect of the simulation. An example would be which of the implemented
    material models to use, and the specifics of this material model. The
    sample parameter file above gives an indication of how this works: within
    a subsection of the file that pertains to the material models, one can
    select one out of several plugins (or, in the case of the postprocessors,
    any number, including none, of the available plugins), and one can then
    specify the specifics of this model in a sub-subsection dedicated to this
    particular model.

    A number of components of ASPECT are
    implemented via plugins. Some of these, together with the sections in
    which their parameters are declared, are the following:

    -   The material model:
        Sections&nbsp;{ref}`parameters:Material_20model` and following.

    -   The geometry: Sections&nbsp;{ref}`parameters:Geometry_20model` and
        following.

    -   The gravity description:
        Sections&nbsp;{ref}`parameters:Gravity_20model` and following.

    -   Initial conditions for the temperature:
        Sections&nbsp;{ref}`parameters:Initial_20temperature_20model` and
        following.

    -   Temperature boundary conditions:
        Sections&nbsp;{ref}`parameters:Boundary_20temperature_20model` and
        following.

    -   Postprocessors: Sections&nbsp;{ref}`parameters:Postprocess` and
        following for most postprocessors, section
        {ref}`parameters:Postprocess/Visualization` and following for
        postprocessors related to visualization.

The details of parameters in each of these categories can be found in the
sections linked to above. Some of them will also be used in the cookbooks in
{ref}`cha:cookbooks`.
