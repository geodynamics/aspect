(sec:cookbooks:prescribed_velocity)=
# Prescribed internal velocity constraints

*This section was contributed by Jonathan Perry-Houts*

In cases where it is desirable to investigate the behavior of one part of the
model domain but the controlling physics of another part is difficult to
capture, such as corner flow in subduction zones, it may be useful to force
the desired behavior in some parts of the model domain and solve for the
resulting flow everywhere else. This is possible through the use of &rsquo;s
&ldquo;signal&rdquo; mechanism, as documented in
Section&nbsp;[\[sec:extending-signals\]][1].

Internally, adds &ldquo;constraints&rdquo; to the finite element system for
boundary conditions and hanging nodes. These are places in the finite element
system where certain solution variables are required to match some prescribed
value. Although it is somewhat mathematically inadmissible to prescribe
constraints on nodes inside the model domain, $\Omega$, it is nevertheless
possible so long as the prescribed velocity field fits in to the finite
element&rsquo;s solution space, and satisfies the other constraints (i.e., is
divergence free).

Using &rsquo;s signals mechanism, we write a shared library which provides a
&ldquo;slot&rdquo; that listens for the signal which is triggered after the
regular model constraints are set, but before they are
&ldquo;distributed.&rdquo;

As an example of this functionality, below is a plugin which allows the user
to prescribe internal velocities with functions in a parameter file:

```{literalinclude} prescribed_velocity.cc
```

The above plugin can be compiled with `cmake . && make` in the
[cookbooks/prescribed_velocity][] directory. It can be loaded in a parameter
file as an &ldquo;Additional shared library.&rdquo; By setting parameters like
those shown below, it is possible to produce many interesting flow fields such
as the ones visualized in (Figure&nbsp;[\[fig:prescribed-velocity\]][2]).

```{literalinclude} corner_flow.prm
```

```{figure-md} fig:quickref
<img src="_static/images/aspect_logo.*" alt="Screenshot"  width="100%"/>

This is the figure caption.
```

&nbsp;

  [1]: #sec:extending-signals
  [cookbooks/prescribed_velocity]: cookbooks/prescribed_velocity
  [2]: #fig:prescribed-velocity
