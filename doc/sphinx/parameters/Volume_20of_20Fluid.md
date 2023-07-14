(parameters:Volume_20of_20Fluid)=
# Volume of Fluid


## **Subsection:** Volume of Fluid


(parameters:Volume_20of_20Fluid/Number_20initialization_20samples)=
### __Parameter name:__ Number initialization samples
**Default value:** 3

**Pattern:** [Integer range 1...2147483647 (inclusive)]

**Documentation:** Number of divisions per dimension when computing the initial volume fractions.If set to the default of 3 for a 2d model, then initialization will be based on the initialization criterion at $3^2=9$ points within each cell. If the initialization based on a composition style initial condition, a larger value may be desired for better approximation of the initial fluid fractions. Smaller values will suffice in the case of level set initializations due to the presence of more information to better approximate the initial fluid fractions.

(parameters:Volume_20of_20Fluid/Volume_20fraction_20threshold)=
### __Parameter name:__ Volume fraction threshold
**Default value:** 1e-6

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** Minimum significant volume. Fluid fractions below this value are considered to be zero.

(parameters:Volume_20of_20Fluid/Volume_20of_20Fluid_20solver_20tolerance)=
### __Parameter name:__ Volume of Fluid solver tolerance
**Default value:** 1e-12

**Pattern:** [Double 0...1 (inclusive)]

**Documentation:** The relative tolerance up to which the linear system for the Volume of Fluid system gets solved. See &rsquo;Solver parameters/Composition solver tolerance&rsquo; for more details.
