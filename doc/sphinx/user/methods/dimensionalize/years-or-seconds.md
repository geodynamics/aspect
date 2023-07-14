(sec:methods:dimensionalize:years-or-seconds)=
# Years or seconds?

All internal calculations in ASPECT are performed using time units of seconds.
Input quantities with units of time or velocity are assumed to be in seconds or meters per second, and output quantities with units of time or velocity will also be in seconds or meters per second, unless the input parameter `Use years in output instead of seconds` is `true` (see {ref}`parameters:global`).

This parameter is somewhat deceptively named, as it influences how ASPECT treats inputs as well as outputs.
For example, if `Use years in output instead of seconds` is `true`, input values for `Start time`, `End time`, and `Maximum time step` are assumed to be in years instead of seconds.
When the flag is set, ASPECT converts input time and velocity units to MKS internally, computes solutions, and converts time and velocity outputs back to years and meters per year during postprocessing.

By default, `Use years in output instead of seconds` is `true`, since ASPECT is designed primarily for models described in physical units rather than in non-dimensionalized form, and years are often more intuitive time units for mantle convection problems (see {ref}`sec:methods:dimensionalize`).
For non- dimensional models the flag should be set to `false` since conversions between years and seconds do not make sense for non-dimensional quantities.
