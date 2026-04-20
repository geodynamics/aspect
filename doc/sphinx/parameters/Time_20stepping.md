(parameters:Time_20stepping)=
# Time stepping


## **Subsection:** Time stepping


::::{dropdown} __Parameter:__ {ref}`List of model names<parameters:Time_20stepping/List_20of_20model_20names>`
:name: parameters:Time_20stepping/List_20of_20model_20names
**Default value:**

**Pattern:** [MultipleSelection conduction time step|convection time step|function|repeat on cutback|repeat on nonlinear solver failure ]

**Documentation:** A comma separated list of time stepping plugins that will be used to calculate the time step size. The minimum of the  result of each plugin will be used.

The following plugins are available:

&lsquo;conduction time step&rsquo;: This model computes the conduction time step as the minimum over all cells of $ CFL h^2 \cdot \rho C_p / k$, where k is the thermal conductivity. This plugin will always request advancing to the next time step.

&lsquo;convection time step&rsquo;: This model computes the convection time step as $ CFL / \max \| u \| / h$ over all cells, where $u$ is the velocity and $h$ is the product of mesh size and temperature polynomial degree.

&lsquo;function&rsquo;: This model uses a time step specified in the parameter file specified as a function of time. This plugin will always request advancing to the next time step.

&lsquo;repeat on cutback&rsquo;: This time stepping plugin will detect a situation where the computed time step shrinks by more than a user-controlled factor. In that situation, the previous time step will be repeated with a smaller step size.
A large reduction in time step size typically happens when velocities change abruptly. Repeating the time step ensure properly resolving this event. It is useful to consider setting the "Maximum relative increase in time step" option to avoid repeatedly repeating every other time step.

&lsquo;repeat on nonlinear solver failure&rsquo;: This time stepping plugin will react when the nonlinear solver does not converge in the specified maximum number of iterations and repeats the current timestep with a smaller step size. This plugin is enabled automatically if "Nonlinear solver failure strategy" is set to "cut timestep size".
::::

::::{dropdown} __Parameter:__ {ref}`Minimum time step size<parameters:Time_20stepping/Minimum_20time_20step_20size>`
:name: parameters:Time_20stepping/Minimum_20time_20step_20size
**Default value:** 0.

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** Specify a minimum time step size (or 0 to disable).
::::

(parameters:Time_20stepping/Function)=
## **Subsection:** Time stepping / Function
::::{dropdown} __Parameter:__ {ref}`Function constants<parameters:Time_20stepping/Function/Function_20constants>`
:name: parameters:Time_20stepping/Function/Function_20constants
**Default value:**

**Pattern:** [Anything]

**Documentation:** Sometimes it is convenient to use symbolic constants in the expression that describes the function, rather than having to use its numeric value everywhere the constant appears. These values can be defined using this parameter, in the form &lsquo;var1=value1, var2=value2, ...&rsquo;.

A typical example would be to set this runtime parameter to &lsquo;pi=3.1415926536&rsquo; and then use &lsquo;pi&rsquo; in the expression of the actual formula. (That said, for convenience this class actually defines both &lsquo;pi&rsquo; and &lsquo;Pi&rsquo; by default, but you get the idea.)
::::

::::{dropdown} __Parameter:__ {ref}`Function expression<parameters:Time_20stepping/Function/Function_20expression>`
:name: parameters:Time_20stepping/Function/Function_20expression
**Default value:** 1.0

**Pattern:** [Anything]

**Documentation:** Expression for the time step size as a function of &rsquo;time&rsquo;.Units: years if the &rsquo;Use years instead of seconds&rsquo; parameter is set; seconds otherwise.
::::

::::{dropdown} __Parameter:__ {ref}`Variable names<parameters:Time_20stepping/Function/Variable_20names>`
:name: parameters:Time_20stepping/Function/Variable_20names
**Default value:** time

**Pattern:** [Anything]

**Documentation:** Name for the variable representing the current time.
::::

(parameters:Time_20stepping/Repeat_20on_20cutback)=
## **Subsection:** Time stepping / Repeat on cutback
::::{dropdown} __Parameter:__ {ref}`Cut back amount<parameters:Time_20stepping/Repeat_20on_20cutback/Cut_20back_20amount>`
:name: parameters:Time_20stepping/Repeat_20on_20cutback/Cut_20back_20amount
**Default value:** 0.5

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** A factor that controls the size of the time step when repeating. The default of 0.5 corresponds to 50\% of the original step taken.
::::

::::{dropdown} __Parameter:__ {ref}`Relative repeat threshold<parameters:Time_20stepping/Repeat_20on_20cutback/Relative_20repeat_20threshold>`
:name: parameters:Time_20stepping/Repeat_20on_20cutback/Relative_20repeat_20threshold
**Default value:** 0.2

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** A factor that controls when a step is going to be repeated. If the newly computed step size is smaller than the last step size multiplied by this factor, the step is repeated.
::::

(parameters:Time_20stepping/Repeat_20on_20nonlinear_20solver_20failure)=
## **Subsection:** Time stepping / Repeat on nonlinear solver failure
::::{dropdown} __Parameter:__ {ref}`Cut back factor<parameters:Time_20stepping/Repeat_20on_20nonlinear_20solver_20failure/Cut_20back_20factor>`
:name: parameters:Time_20stepping/Repeat_20on_20nonlinear_20solver_20failure/Cut_20back_20factor
**Default value:** 0.5

**Pattern:** [Double 0...MAX_DOUBLE (inclusive)]

**Documentation:** A factor that controls the size of the time step when repeating. The default of 0.5 corresponds to 50\% of the original step taken.
::::
