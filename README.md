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

`feynson` [options] **command** *args* ...

## DESCRIPTION

## COMMANDS

* `feynson` **symmetrize** *spec-file*

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

  Each family that can be mapped to (a subsector of) another
  is guaranteed to be mapped to the first possible family,
  prefering families that are larger or listed earlier.

* `feynson` **zero-sectors** [-s] *spec-file*

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

  If the -s flag is given, the output will be shortened
  by only listing the topmost zero sectors: all the remaining
  zero sectors are their subsectors.

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

* -j *jobs*

  Parallelize calculations using at most this many workers.

* -s

  Shorten the output (depending on the command).

* -q

  Print a more quiet log.

* -h

  Show this help message.

* -V

  Print version information.

## ARGUMENTS

* *spec-file*

  Filename of the input file, with `-` meaning the standard input.

## ENVIRONMENT

`TMPDIR`     Temporary files will be created here.

## AUTHORS

Vitaly Magerya
