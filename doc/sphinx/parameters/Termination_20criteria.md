(parameters:Termination_20criteria)=
# Termination criteria


## **Subsection:** Termination criteria


(parameters:Termination_20criteria/Checkpoint_20on_20termination)=
### __Parameter name:__ Checkpoint on termination
**Default value:** false

**Pattern:** [Bool]

**Documentation:** Whether to checkpoint the simulation right before termination.

(parameters:Termination_20criteria/End_20step)=
### __Parameter name:__ End step
**Default value:** 100

**Pattern:** [Integer range 0...2147483647 (inclusive)]

**Documentation:** Terminate the simulation once the specified timestep has been reached.

(parameters:Termination_20criteria/Termination_20criteria)=
### __Parameter name:__ Termination criteria
**Default value:** end time

**Pattern:** [MultipleSelection end step|end time|steady state heat flux|steady state temperature|steady state velocity|user request|wall time ]

**Documentation:** A comma separated list of termination criteria that will determine when the simulation should end. Whether explicitly stated or not, the &ldquo;end time&rdquo; termination criterion will always be used.The following termination criteria are available:

&lsquo;end step&rsquo;: Terminate the simulation once the specified timestep has been reached.

&lsquo;end time&rsquo;: Terminate the simulation once the end time specified in the input file has been reached. Unlike all other termination criteria, this criterion is *always* active, whether it has been explicitly selected or not in the input file (this is done to preserve historical behavior of ASPECT, but it also likely does not inconvenience anyone since it is what would be selected in most cases anyway).

&lsquo;steady state heat flux&rsquo;: A criterion that terminates the simulation when the integrated heat flux over a given list of boundaries stays within a certain range for a specified period of time.

The criterion considers the total heat flux over all boundaries listed by their boundary indicators, rather than each boundary separately. As a consequence, if the *sum* of heat fluxes over individual parts of the boundary no longer changes, then this criterion recommends termination, even if the heat flux over individual parts of the boundary continues to change.

&lsquo;steady state temperature&rsquo;: A criterion that terminates the simulation when the global integral of the temperature field stays within a certain range for a specified period of time.

&lsquo;steady state velocity&rsquo;: A criterion that terminates the simulation when the RMS of the velocity field stays within a certain range for a specified period of time.

&lsquo;user request&rsquo;: Terminate the simulation gracefully when a file with a specified name appears in the output directory. This allows the user to gracefully exit the simulation at any time by simply creating such a file using, for example, `touch output/terminate`. The file&rsquo;s location is chosen to be in the output directory, rather than in a generic location such as the ASPECT directory, so that one can run multiple simulations at the same time (which presumably write to different output directories) and can selectively terminate a particular one.

&lsquo;wall time&rsquo;: Terminate the simulation once the wall time limit has reached.

(parameters:Termination_20criteria/Wall_20time)=
### __Parameter name:__ Wall time
**Default value:** 24.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The wall time of the simulation. Unit: hours.

(parameters:Termination_20criteria/Steady_20state_20heat_20flux)=
## **Subsection:** Termination criteria / Steady state heat flux
(parameters:Termination_20criteria/Steady_20state_20heat_20flux/Boundary_20indicators)=
### __Parameter name:__ Boundary indicators
**Default value:**

**Pattern:** [List of <[Anything]> of length 0...4294967295 (inclusive)]

**Documentation:** A comma separated list of names denoting those boundaries that should be taken into account for integrating the heat flux. Note that the plugin will compute the integrated heat flux over these boundaries (instead of taking them into account individually).

The names of the boundaries listed here can either be numbers (in which case they correspond to the numerical boundary indicators assigned by the geometry object), or they can correspond to any of the symbolic names the geometry object may have provided for each part of the boundary. You may want to compare this with the documentation of the geometry model you use in your model.

(parameters:Termination_20criteria/Steady_20state_20heat_20flux/Maximum_20relative_20deviation)=
### __Parameter name:__ Maximum relative deviation
**Default value:** 0.05

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The maximum relative deviation of the heat flux in recent simulation time for the system to be considered in steady state. If the actual deviation is smaller than this number, then the simulation will be terminated.

(parameters:Termination_20criteria/Steady_20state_20heat_20flux/Time_20in_20steady_20state)=
### __Parameter name:__ Time in steady state
**Default value:** 1e7

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The minimum length of simulation time that the system should be in steady state before termination. Note that if the time step size is similar to or larger than this value, the termination criterion will only have very few (in the most extreme case, just two) heat flux values to check. To ensure that a larger number of time steps are included in the check for steady state, this value should be much larger than the time step size. Units: years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Termination_20criteria/Steady_20state_20temperature)=
## **Subsection:** Termination criteria / Steady state temperature
(parameters:Termination_20criteria/Steady_20state_20temperature/Maximum_20relative_20deviation)=
### __Parameter name:__ Maximum relative deviation
**Default value:** 0.05

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The maximum relative deviation of the temperature in recent simulation time for the system to be considered in steady state. If the actual deviation is smaller than this number, then the simulation will be terminated.

(parameters:Termination_20criteria/Steady_20state_20temperature/Time_20in_20steady_20state)=
### __Parameter name:__ Time in steady state
**Default value:** 1e7

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The minimum length of simulation time that the system should be in steady state before termination.Units: years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Termination_20criteria/Steady_20state_20velocity)=
## **Subsection:** Termination criteria / Steady state velocity
(parameters:Termination_20criteria/Steady_20state_20velocity/Maximum_20relative_20deviation)=
### __Parameter name:__ Maximum relative deviation
**Default value:** 0.05

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The maximum relative deviation of the RMS in recent simulation time for the system to be considered in steady state. If the actual deviation is smaller than this number, then the simulation will be terminated.

(parameters:Termination_20criteria/Steady_20state_20velocity/Time_20in_20steady_20state)=
### __Parameter name:__ Time in steady state
**Default value:** 1e7

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** The minimum length of simulation time that the system should be in steady state before termination.Units: years if the &rsquo;Use years in output instead of seconds&rsquo; parameter is set; seconds otherwise.

(parameters:Termination_20criteria/User_20request)=
## **Subsection:** Termination criteria / User request
(parameters:Termination_20criteria/User_20request/File_20name)=
### __Parameter name:__ File name
**Default value:** terminate-aspect

**Pattern:** [FileName (Type: input)]

**Documentation:** The name of a file that, if it exists in the output directory (whose name is also specified in the input file) will lead to termination of the simulation. The file&rsquo;s location is chosen to be in the output directory, rather than in a generic location such as the ASPECT directory, so that one can run multiple simulations at the same time (which presumably write to different output directories) and can selectively terminate a particular one.
