---
orphan: true
---
# MyST Quick reference

(sec-quickref)=
# Heading 1
## Heading 2
### Heading 3
#### Heading 4

Refer to Section {ref}`sec-quickref`.

## Admonitions

:::{admonition} General admonition as warning
:class: warning

Text goes here.
:::

:::{attention}
This is an `attention` admonition.
:::

:::{danger}
This is a `danger` admonition.
:::

:::{error}
This is an `error` admonition.
:::

:::{important}
This is an `important` admonition.
:::

:::{note}
This is a `note` admonition.
:::

:::{tip}
This is a `tip` admonition.
:::

:::{warning}
This is a `warning` admonition.
:::

:::{seealso}
This is a `seealso` admonition.
:::

:::{admonition} TODO
:class: error

This is a custom `TODO` admonition.
:::

## Lists

### Itemized lists

* Level 1
  * Level 2
    * Level 3

### Definition lists

Term 1
: Definition of term 1

Term 2
: Definition of term 2

## Code blocks

```{code-block} c++
---
caption: C++ code block.
emphasize-lines: 3-4
---
int
main(int argc, char* argv[]) {
    // Emphasized lines corresponding to body of main().
    return 0;
}
```

```{code-block} python
---
caption: Python code block.
---
def square(x):
    return x**2
```

```{code-block} console
---
caption: Interactive shell.
---
$ ls
a b c
```

```{code-block} bash
---
caption: Bash code block.
---
# Comment
for i in "a b c"; do
  print $i
done
```

```{code-block} cfg
---
caption: Config code block.
---
# Comment
[pylithapp]
journal.info.problem = 1

[pylithapp.petsc]
ksp_rtol = 1.0e-3
```
Captions are optional.

## Codeblocks from external files

```{literalinclude} ../manual/cookbooks/burnman/doc/material_model.part.prm
```

When you're in the working on documentation for cookbooks, here's the path to the cookbooks folder: ../../../../manual/cookbooks/

And if you're benchmarks documentation: ../../../manual/cookbooks/

## Tables

Please see {numref}`tab:quickref`.
Figure labels cannot contain more than one dash, so we use colons instead.
This is likely a bug.

```{table} Table caption
:name: tab:quickref
|             Header 1 |   Header 2    | Header 3       |
| -------------------: | :-----------: | :------------- |
|        right aligned |   centered    | left aligned   |
| more data, more data | yet more data | even more data |
```

:::{note}
The numref command doesn't work in this particular file because its not a part of the main file structure. It does work in the manual, as long as the file is in a toctree.

## Figures

Please see {numref}`fig:quickref`.
Figure labels cannot contain more than one dash, so we use colons instead.
This is likely a bug.


:::{figure-md} fig:quickref
<img src="_static/images/aspect_logo.*" alt="Screenshot"  width="100px"/>

This is the figure caption.
:::

## Math

Look at {math:numref}`eq:one:two`.

```{math}
:label: eqn:one:two
F = m a
```

(The label is not necessary if the equation doesn't have one, like so):

```{math}
F = ma
```

Environments are cool too:
```{math}
:label: eq:aligned
\begin{align}
  -\nabla \cdot \left[2\eta \left(\varepsilon(\mathbf u) - \frac{1}{3}(\nabla \cdot \mathbf u)\mathbf 1\right)\right] + \nabla p' &=  -\bar \alpha \bar\rho T' \mathbf g & \qquad  & \textrm{in $\Omega$},  \\
  \nabla \cdot (\bar\rho \mathbf u) &= 0  & \qquad  & \textrm{in $\Omega$}.
\end{align}
```

##Links

[Put a link in text](https://en.wikipedia.org/wiki/Gibbs_free_energy) or just leave it like this:
<https://en.wikipedia.org/wiki/Gibbs_free_energy>


## Citations

Traditional citation {cite}`kronbichler:etal:2012`.
Citation as noun {cite:t}`kronbichler:etal:2012` or for multiple {cite:t}`kronbichler:etal:2012,heister_aspect_methods2`.

## Footnotes

Here's how you put a footnote[^footnote1] in a sentence.

[^footnote1]: And here's what the footnote says (put this at the very bottom of the document)

# See here for more

<https://myst-parser.readthedocs.io/en/latest/syntax/syntax.html>
or
<https://jupyterbook.org/en/stable/reference/cheatsheet.html>
