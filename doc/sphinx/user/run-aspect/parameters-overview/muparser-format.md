(sec:run-aspect:parameters-overview:muparser-format)=
# How to write Function expressions in parameter files

## Syntax
Input files have different ways of prescribing functional forms used by
ASPECT. In many cases, the user can select parameters that feed into a
fixed functional form defined inside an ASPECT plugin. In many others,
ASPECT allows you to define your own function expression in a symbolic way
in the input file (see {ref}`parameters:Initial_20temperature_20model`).
An example of this latter case is this snippet of code discussed in
{ref}`sec:cookbooks:convection_box_3d`:

```{literalinclude} ../../cookbooks/cookbooks/convection_box_3d/doc/initial.part.prm
```

These function expressions need to use a syntax that can be parsed by ASPECT.
This parsing is done by the muparser library (<http://muparser.beltoforion.de/>).
The syntax is mostly self-explanatory in that it allows you to use the usual
symbols `x`, `y` and `z` for spatial coordinates
(unless a particular plugin uses different variables, such as the depth),
the symbol `t` for time in many situations, and allows you to use many common
mathematical functions such as `||` (the logical OR), `&&` (the logical AND),
`^` (raise x to the power of y), `sin()` and `cos()`
(see <https://beltoforion.de/en/muparser/features.php> for a list of functions).
Based on the muparser library, deal.II supports additional functions,
including `int()`, `ceil()`, `floor()`, `cot()`, `csc()`, `sec()`,
`pow()`, `erfc()`, `rand()`, and `rand_seed()`. For more detailed information, see
<http://www.dealii.org/developer/doxygen/deal.II/classFunctionParser.html#details>.

## Help with MuParser syntax
If you would like to check that muparser is parsing your functions correctly,
and to plot the values, you may consider running your expressions through
`pymuparser` (<https://github.com/bobmyhill/pymuparser>),
which can be installed on your machine using pip
(`python -m pip install pymuparser`). This module also allows users to
define functions not included in MuParser, such as the extended library
provided by deal.II. Python's `eval` function may also be useful,
(<https://docs.python.org/3/library/functions.html#eval>)
but be aware that Python syntax may not be the same as that of deal.II.

## Lazy expressions
If-else-statements are particularly useful in function expressions,
for example "if $1<x<4$ then output 1, else output 0". MuParser uses
lazy-expression syntax `(if-condition ? true-expression : false-expression)`
for if-else statements. This lazy-expression only evaluates the expression
that meets the if-condition, rather than evaluating both expressions, which
can be useful if one of the expressions is not defined (e.g., has a divide by
zero) when the if-condition is not met. Note it is also possible to use the
syntax `if(condition, true-expression, false-expression)`, but in this case
both expressions are always evaluated. This is inefficient, but in addition
may yield avoidable floating point exceptions if the expression that
will be discarded has invalid floating point operations (such as a division by
zero, or taking the square root of a negative number).
Therefore, the lazy-expression syntax is recommended.

Using the lazy-expression syntax, the statement "if
$1<x<4$ then output 1, else output 0" can be expressed as
`(1<x && x<4 ? 1 : 0)`. Multiple, nested if-else expressions can also be used.
For example, the statement
"if $1<x<4$ then, if 2\<y\<3, then output 2, else output 1, else output 0"
can be expressed as `((1<x && x<4) ? ((2<y && y<3) ? 2 : 1) : (0))`.

Nested lazy-expression if-else statements can be found in the cookbook
example at {ref}`sec:cookbooks:muparser_temperature_example`.
This cookbook includes a Python script that shows
how an initial temperature structure can be written using the lazy-expression.
The cookbook runs a single time-step to illustrate the functionality.
We recommend using tools like pymuparser to debug complex expressions before
defining them in the parameter file. More examples of
functions can be found in the `cookbooks` and `tests` directories (use
grep to search for "Function expression" in the prm files).
You can also search "Function expression" on the ASPECT github page.
