(sec:cookbooks:using-particles)=
# Using particles

Using compositional fields to trace where material has come from or is going to
has many advantages from a computational point of view. For example, the
numerical methods to advect along fields are well developed and we can do so at
a cost that is equivalent to one temperature solve for each of the compositional
fields. Unless you have many compositional fields, this cost is therefore
relatively small compared to the overall cost of a time step. Another advantage
is that the value of a compositional field is well defined at every point within
the domain. On the other hand, compositional fields over time diffuse initially
sharp interfaces, as we have seen in the images of the previous section.

At the same time, the geodynamics community has a history of using particles for
this purpose. Historically, this may have been because it is conceptually
simpler to advect along individual particles rather than whole fields, since it
only requires an ODE integrator rather than the stabilization techniques
necessary to advect fields. They also provide the appearance of no diffusion,
though this is arguable. Leaving aside the debate whether fields or particles are the
way to go, ASPECT supports both: using fields and using particles.

:::{toctree}
cookbooks/composition_passive_particles/doc/composition_passive_particles.md
cookbooks/composition_active_particles/doc/composition_active_particles.md
:::
