# Limiting postprocessing

ASPECT has a lot of postprocessing
capabilities, from generating graphical output to computing average
temperatures or temperature fluxes. To see what all is possible, take a look
at the `List of postprocessors` parameter that can be set in the input file,
see {ref}`parameters:Postprocess`.

Many of these postprocessors take a non-negligible amount of time. How much
they collectively use can be inferred from the timing report
ASPECT prints periodically among its output, see for
example the output shown in {ref}`sec:cookbooks:convection-box`. So, if your computations
take too long, consider limiting which postprocessors you run to those you
really need. Some postprocessors - for example those that generate
graphical output, see
{ref}`parameters:Postprocess/Visualization` - also allow
you to run them only once every once in a while, rather than at every time
step.
