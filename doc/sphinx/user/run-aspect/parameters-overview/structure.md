# The structure of parameter files

Most of the run-time behavior of ASPECT is
driven by a parameter file that looks in essence like this:

``` prmfile
```

Some parameters live at the top level, but most parameters are grouped into
subsections. An input parameter file is therefore much like a file system: a
few files live in the root directory; others are in a nested hierarchy of
sub-directories. And just as with files, parameters have both a name (the
thing to the left of the equals sign) and a content (what's to the
right).

All parameters you can list in this input file have been *declared* in
ASPECT. What this means is that you can't just
list anything in the input file, and expect that entries that are unknown are
simply ignored. Rather, if your input file contains a line setting a parameter
that is unknown, you will get an error message. Likewise, all declared
parameters have a description of possible values associated with them -
for example, some parameters must be non-negative integers (the number of
initial refinement steps), can either be true or false (whether the
computation should be resumed from a saved state), or can only be a single
element from a selection (the name of the material model). If an entry in your
input file doesn't satisfy these constraints, it will be rejected at the
time of reading the file (and not when a part of the program actually accesses
the value and the programmer has taken the time to also implement some error
checking at this location). Finally, because parameters have been declared,
you do not *need* to specify a parameter in the input file: if a parameter
isn't listed, then the program will simply use the default provided when
declaring the parameter.

:::{note}
In cases where a parameter requires a significant amount of text, you can end a line in the
input file with a backslash. This indicates that the following line will simply continue to be part
of the text of the current line, in the same way as the C/C++ preprocessor expands lines that
end in backslashes. The underlying implementation always eats whitespace at the beginning of
each continuing line, but not before the backslash. This means that the parameter file
```
set Some parameter = abc\
def
```
is equivalent to
```
set Some parameter = abcdef
```
that is, with no space between `abc` and `def` despite the leading whitespace at the beginning of
the second line. If you do want space between these two parts, you need to add it before the
backslash in the first of the two lines.
