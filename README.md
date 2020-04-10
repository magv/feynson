# FEYNSON

*Feynson* finds symmetries between families of Feynman integrals.
One can use it to reduce the number of families before handing
the integrals over to an IBP relation solver.

# MANUAL

## NAME

`feynson` -- a tool for Feynman integral symmetries.

## SYNOPSYS

`feynson` [options] **command** *args* ...

## DESCRIPTION

## COMMANDS

* **symmetrize** *spec-file*

  Print a list of momenta substitutions that make symmetries
  between a list of integral families explicit.

  These momenta substitutions will make the set of
  denominators of two integral families exactly match (up
  to a reordering) if the families are isomorphic, and
  will make one a subset of the other if one family is
  isomorphic to a subsector of another family.

* **zero-sectors** [-s] *spec-file*

  Print a list of all zero sectors of a given integral
  family.

  The input specification file should be a list of four
  elements:
  1) a list of all propagator momenta (e.g. "(l1-q)^2");
  2) a list of cut flags, "0" for normal propagators, "1"
     for cut propagators;
  3) a list of all loop momenta (e.g. "l1");
  4) and a list of external invariant substitutions (e.g.
     "{q^2, 1}").

  The output will be a list of zero sectors, each denoted
  by an integer s=2^{i_1-1} + ... + 2^{i_n-1}, where i_k
  are the indices of denominators that belong to this
  sector (counting from 1).

  If the -s flag is given, the output will be shortened
  by only listing the topmost zero sectors: all the remaining
  zero sectors are their subsectors.

  Every sector that is missing a cut propagator of its
  supersectors will be reported as zero.

* **ufx** *spec-file*

  Print Feynman parametrization (U, F, X) of an integral
  defined by a set of propagators.

  The input specification file should be a list of three
  elements:
  1) a list of all propagators, e.g. "(l1-q)^2";
  2) a list of all loop momenta, e.g. "l1";
  3) and a list of external invariant substitutions, e.g.
     "{q^2, 1}".

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

  Filename of the input file, with "-" meaning the standard input.

## ENVIRONMENT

`TMPDIR`     Temporary files will be created here.

## AUTHORS

Vitaly Magerya
