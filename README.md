# FEYNSON

*Feynson* finds symmetries between families of Feynman integrals.
One can use it to reduce the number of families before handing
the integrals over to an IBP relation solver.

In principle IBP solvers perform all the functions of *Feynson*
too. The need for a separate tool is that the construction of
integral families needs to be done before the final expression
is ready for an IBP solver. Additionally, not all IBP solvers
are optimized to receive thousands of (mostly redundant) integral
families; *Feynson* is.

# BUILDING

To build *Feynson*, first install GiNaC [1] and Nauty [2]
libraries, and then run:

    make

A static binary of *Feynson* can also be built. For this you will
need UPX [3], and the static versions of GiNaC and Nauty libraries
(as well as their dependencies). Once those are ready, run:

    make feynson.static

[1] https://www.ginac.de/

[2] https://pallini.di.uniroma1.it/

[3] https://upx.github.io/

# MANUAL

## NAME

`feynson` -- a tool for Feynman integral symmetries.

## SYNOPSYS

`feynson` [`options`] **command** *args* ...

## DESCRIPTION

## COMMANDS

* `feynson` **symmetrize** [`-d`] *spec-file*

  Print a list of momenta substitutions that make symmetries
  between a list of integral families explicit.

  These momenta substitutions will make the set of
  denominators of two integral families exactly match (up
  to a reordering) if the families are isomorphic, and
  will make one a subset of the other if one family is
  isomorphic to a subsector of another family.

  The input specification file should be a list of three
  elements:
  1) a list of all integral families, with each family
     being a list of propagators (e.g. `(l1+l2)^2`);
  2) a list of all loop momenta;
  3) a list of external invariant substitution rules, each
     rule being a list of two elements: a scalar product
     and its substitution (e.g. `{q^2, 1}` or `{p1*p2, s12}`).

  For example: `{ {{(q-l)^2, l^2}, {(q+l)^2, l^2}}, {l}, {} }`.

  Each family that can be mapped to (a subsector of) another
  is guaranteed to be mapped to the first possible family,
  prefering families that are larger or listed earlier.
  If `-d` flag is given, earlier families are prefered
  irrespective of their size.

  Note that if non-trivial invariant substitution rules
  are supplied, it becomes possible that two families are
  identical, but no loop momenta substitution exists to map
  them onto each other. For example, a 1-loop propagator
  with momenta `p1` is equal to a 1-loop propagator with
  `p2`, but only if `p1^2 = p2^2`, in which case no loop
  momenta substitution can make the integrands identical.

  For this reason, it is best to use **symmetrize** with the
  invariant substitution rules set to `{}`, and to fall back
  to **mapping-rules** otherwise.

* `feynson` **mapping-rules** [`-d`] *spec-file*

  Same as **symmetrize**, but instead of printing the loop
  momenta substitutions, produce explicit rules of mapping
  between families: for each family that is symmetric to
  another, print `{fam, {n1, n2, ...}}`, meaning that any
  integral in this family with indices `{i_1, i_2, ...}`
  is equal to an integral in the family number `fam` with
  indices `{i_n1, i_n2, ...}`. For unique families, print `{}`.
  The families are numbered starting at 1. If a given family
  is symmetric to a subfamily, some of the `n` indices will
  be `0`: the convention is that `i_0 = 0`.

* `feynson` **zero-sectors** [`-s`] *spec-file*

  Print a list of all zero sectors of a given integral
  family.

  The input specification file should be a list of four
  elements:
  1) a list of all propagator momenta (e.g. `(l1-q)^2`);
  2) a list of cut flags, `0` for normal propagators, `1`
     for cut propagators;
  3) a list of all loop momenta (e.g. `l1`);
  4) and a list of external invariant substitutions (e.g.
     `{q^2, 1}`).

  For example: `{ {(q-l)^2, l^2}, {0, 0}, {l}, {{q^2,1}} }`.

  The output will be a list of zero sectors, each denoted
  by an integer s=2^{i_1-1} + ... + 2^{i_n-1}, where i_k
  are the indices of denominators that belong to this
  sector (counting from 1).

  If the `-s` ("short") flag is given, the output will
  be shortened by only listing the topmost zero sectors:
  all the remaining zero sectors are their subsectors.

  Every sector that is missing a cut propagator of its
  supersectors will be reported as zero.

* `feynson` **ufx** *spec-file*

  Print Feynman parametrization (U, F, X) of an integral
  defined by a set of propagators.

  The input specification file should be a list of three
  elements:
  1) a list of all propagators, e.g. `(l1-q)^2`;
  2) a list of all loop momenta, e.g. `l1`;
  3) and a list of external invariant substitutions, e.g.
     `{q^2, 1}`.

  For example: `{ {(q-l)^2, l^2}, {l}, {{q^2,1}} }`.

  The output will be a list of three items: the U polynomial,
  the F polynomial, and the list of Feynman parameter
  variables.

## OPTIONS

* `-j` *jobs*

  Parallelize calculations using at most this many workers.

* `-d`

  Prioritize families in the definition order, irrespective of size.

* `-s`

  Shorten the output (in **zero-sectors** and **minimize-family**).

* `-q`

  Print a more quiet log.

* `-h`

  Show this help message.

* `-V`

  Print version information.

## ARGUMENTS

* *spec-file*

  Filename of the input file, with `-` meaning the standard input.

## ENVIRONMENT

`TMPDIR`     Temporary files will be created here.

## AUTHORS

Vitaly Magerya
