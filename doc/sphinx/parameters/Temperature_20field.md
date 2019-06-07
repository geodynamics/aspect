(parameters:Temperature_20field)=
# Temperature field


## **Parameters in section** Temperature field


(parameters:Temperature_20field/Temperature_20method)=
### __Parameter name:__ Temperature method
**Default value:** field

**Pattern:** [Selection field|prescribed field|prescribed field with diffusion|static ]

**Documentation:** A comma separated list denoting the solution method of the temperature field. Each entry of the list must be one of the currently implemented field types.

These choices correspond to the following methods by which the temperature field gains its values:\begin{itemize}\item ``field'': If the temperature is marked with this method, then its values are computed in each time step by solving the temperature advection-diffusion equation. In other words, this corresponds to the usual notion of a temperature.
\item ``prescribed field'': The value of the temperature is determined in each time step from the material model. If a compositional field is marked with this method, then the value of a specific additional material model output, called the `PrescribedTemperatureOutputs' is interpolated onto the temperature. This field does not change otherwise, it is not advected with the flow.
\item ``prescribed field with diffusion'': If the temperature field is marked this way, the value of a specific additional material model output, called the `PrescribedTemperatureOutputs' is interpolated onto the field, as in the ``prescribed field'' method. Afterwards, the field is diffused based on a solver parameter, the diffusion length scale, smoothing the field. Specifically, the field is updated by solving the equation $(I-l^2 \Delta) T_\text{smoothed} = T_\text{prescribed}$, where $l$ is the diffusion length scale. Note that this means that the amount of diffusion is independent of the time step size, and that the field is not advected with the flow.
\item ``static'': If a temperature field is marked this way, then it does not evolve at all. Its values are simply set to the initial conditions, and will then never change.\end{itemize}