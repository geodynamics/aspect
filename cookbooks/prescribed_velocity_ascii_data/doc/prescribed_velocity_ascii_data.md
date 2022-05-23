#### Prescribing internal velocity constraints with ASCII files

*This section was contributed by Bob Myhill*

Building on [\[sec:prescribed-velocities\]][1], the
[cookbooks/prescribed_velocity_ascii_data] directory contains a plugin which
uses an ASCII data file to specify where to prescribe internal velocities.
Velocities are prescribed wherever the field value indicated by the ASCII data
file is greater than 0.5. As before, the plugin is loaded in parameter files
as an additional shared library:

``` prmfile
```

An example parameter file using this plugin can be found at
[cookbooks/prescribed_velocity_ascii_data/prescribed_velocity_ascii_data.prm].
In this file, the velocities are constrained to be zero within the letters
&ldquo;ASPECT&rdquo; (Figure&nbsp;[1]). The part of this file which provides
the location of the ASCII file and the prescribed velocity field function is:

``` prmfile
```

A temperature gradient is applied within the letters, while the temperature
field outside the letters is set to be constant. This initial temperature
field is specified by another ASCII data file:

``` prmfile
```

These two ASCII data files are generated from `aspect_name.png` by the python
file `make_ascii_files_from_png.py`, both of which can be found in the same
directory as the parameter file.

<figure>
<img src="cookbooks/prescribed_velocity_ascii_data/doc/prescribed_velocities_ascii_data_initial_conditions.png" id="fig:prescribed-velocity-ascii-data-init" alt="Initial composition and temperature conditions for the prescribed velocity ascii data cookbook, as described in Section 0.0.1." /><figcaption aria-hidden="true"><em>Initial composition and temperature conditions for the prescribed velocity ascii data cookbook, as described in Section <a href="#sec:prescribed-velocities-ascii-data" data-reference-type="ref" data-reference="sec:prescribed-velocities-ascii-data">0.0.1</a>.</em></figcaption>
</figure>

During the simulation, excess heat diffuses out from the tops of the letters,
and into the bases of the letters. The temperature gradients in the
unconstrained part of the domain then generate convective flow.
Figure&nbsp;[2] illustrates the resulting flow field.

<figure>
<img src="cookbooks/prescribed_velocity_ascii_data/doc/prescribed_velocity_ascii_data.png" id="fig:prescribed-velocity-ascii-data" style="width:80.0%" alt="Convective flow around the letters ASPECT, within which velocities are prescribed to be zero, as described in Section 0.0.1." /><figcaption aria-hidden="true"><em>Convective flow around the letters ASPECT, within which velocities are prescribed to be zero, as described in Section <a href="#sec:prescribed-velocities-ascii-data" data-reference-type="ref" data-reference="sec:prescribed-velocities-ascii-data">0.0.1</a>.</em></figcaption>
</figure>

  [1]: #sec:prescribed-velocities
  [cookbooks/prescribed_velocity_ascii_data]: cookbooks/prescribed_velocity_ascii_data
  [cookbooks/prescribed_velocity_ascii_data/prescribed_velocity_ascii_data.prm]:
    cookbooks/prescribed_velocity_ascii_data/prescribed_velocity_ascii_data.prm
  [1]: #fig:prescribed-velocity-ascii-data-init
  [2]: #fig:prescribed-velocity-ascii-data
