(sec:run-aspect:parameters-overview:muparser-format)=
# A note on the syntax of formulas in input files

#### A note on the syntax of formulas in input files

Input files have different ways of describing certain things to
ASPECT. For example, you could select a plugin for
the temperature initial values that prescribes a constant temperature, or a
plugin that implements a particular formula for these initial conditions in
C++ in the code of the plugin, or a plugin that allows you to describe this
formula in a symbolic way in the input file (see
{ref}`parameters:Initial_20temperature_20model`19]). An example
of this latter case is this snippet of code discussed in
{ref}`5.2.2][]:

``` prmfile
```

The formulas you can enter here need to use a syntax that is understood by the
functions and classes that interpret what you write. Internally, this is done
using the muparser library, see <http://muparser.beltoforion.de/>. The syntax
is mostly self-explanatory in that it allows to use the usual symbols `x`, `y`
and `z` to reference coordinates (unless a particular plugin uses different
variables, such as the depth), the symbol `t` for time in many situations, and
allows you to use all of the typical mathematical functions such as sine and
cosine. Based on the muparser library, deal.II supports additional functions,
including `|` (the logical OR), `&` (the logical AND), `int()`, `ceil()`,
`floor()`, `cot()`, `csc()`, `sec()`, `pow()`, `log()`, `erfc()`, `rand()`,
and `rand_seed()`. For more detailed information, see
<http://www.dealii.org/developer/doxygen/deal.II/classFunctionParser.html#details>.

A common need for function expression is an if-else-statement, for example
"if $1<x<4$ then output 1, else output 0." The muparser uses
lazy-expression syntax `(if-condition ? true-expression : false-expression)`
for if-else statements. This lazy-expression only evaluates the expression
that meets the if-condition, rather than evaluating both expressions, which
can be useful if one of the expressions is not defined (e.g., has a divide by
zero) when the if-condition is not met. Note it is also possible to use the
syntax `if(condition, true-expression, false-expression)`, but in this case
both expressions are always evaluated. This is inefficient, but in addition
may abort the program with a floating point exception if the expression that
will be discarded has invalid floating point operations (such as a division by
zero, or taking the square root of a negative number) that would ordinarily
not be visible because, after all, the expression should be discarded.
Therefore, the lazy-expression syntax is recommended.

As a simple example using the lazy-expression syntax, the statement "if
$1<x<4$ then output 1, else output 0" can be expressed as
`(1<x && x<4 ? 1 : 0)`. Multiple, nested if-else expressions can also be used.
To extend the simple example, the statement "if $1<x<4$ then, if
2\<y\<3, then output 2, else output 1, else output 0" can be expressed
as `((1<x && x<4) ? ((2<y && y<3) ? 2 : 1) : (0))`.

An example for how to translate nested if-else statements into the
lazy-expression syntax is given in the cookbook example found in
{ref}`5.2.14][]. This cookbook includes a python script that defines
the initial temperature structure using nested if-else statements and shows
how this is then rewritten using the lazy-expression. The cookbook runs a
single time-step to show the outcome of using the function option for the
initial temperature. Quite complex initial conditions can be defined in this
way, however, using something like python to debug these expressions before
defining them in the parameter file is recommended. For more examples of
functions used in parameter files, go to the `cookbooks` directory and use
grep to search for "Function expression" in the parameters files.
You can also search "Function expression" on the
ASPECT github page. For more examples of the syntax
understood, reference the documentation of the muparser library linked to
above.
