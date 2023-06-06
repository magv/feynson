(* Print the error message and stop the computation. Exit with
 * an error code if running in a script; raise an exception when
 * in GUI. *)
Error[msg__] := If[Length[Cases[$CommandLine, "-script"]] > 0,
    Print["ERROR: ", msg]; Exit[1];
    ,
    Print[Style["ERROR: ", Red, Bold], msg]; Throw[$Failed];
]

(* Fail the computation unless a condition is met. Useful for
 * assetions and unit tests. *)
FailUnless[tests___] := Module[{test, idx, result},
  Do[
    test = Extract[Hold[tests], {idx}, HoldForm];
    If[test === HoldForm[Null], Continue[]];
    result = ReleaseHold[test];
    If[result =!= True,
      If[MatchQ[Extract[test, {1,0}, HoldForm], HoldForm[_Symbol]],
        Print["!!! Test: ", Extract[test, {1,0}, HoldForm], " => ", result];
        Print["!!! 1: ", Extract[test, {1,1}, HoldForm]];
        Print["!!! == ", Extract[test, {1,1}]];
        Print["!!! 2: ", Extract[test, {1,2}, HoldForm]];
        Print["!!! == ", Extract[test, {1,2}]];
        ,
        Print["!!! Test: ", test];
        Print["!!!    => ", result];
      ];
      Error["Test failed!"];
    ];
    ,
    {idx, Length[Hold[tests]]}];
];
SetAttributes[FailUnless, {HoldAll}]

(* Call feynson with the given command-line arguments, pass the
 * last parameter on the standard input (Mathematica format),
 * and return the parsed standard output (Mathematica format).
 *)
Feynson[cmd_, lastarg_] := RunThrough[{"./feynson ", cmd}//Flatten//Map[ToString]//Apply[StringJoin], lastarg]

FailUnless[
  Feynson["ufx -q -", { {(q-l)^2, l^2}, {l}, {{q^2, 1}} }] ===
  {x1 + x2, x1 x2, {x1, x2}}]

FailUnless[
  Feynson["ufx -q -", {
      {(l)^2-mm, (l-p1)^2-mm, (l-p1-p2)^2-mm},
      {l},
      {{p1^2, p1sq}, {p2^2, p2sq}, {p1*p2, p3sq/2 - p1sq/2 - p2sq/2}}
  }] ===
  {
    x1 + x2 + x3,
    -mm*x1^2 -mm*x2^2 -mm*x3^2
    -2*mm*x1*x2 -2*mm*x1*x3 -2*mm*x2*x3
    +p1sq*x1*x2 +p3sq*x1*x3 +p2sq*x2*x3,
    {x1, x2, x3}
  }]

FailUnless[
  Feynson["ufx -q -", {
      {(l)^2, (l-p1)^2, (l-p1-p2)^2, (l-p1-p2-p3)^2},
      {l},
      {{p1^2, 0}, {p2^2, 0}, {p3^2, 0}, {p1*p2, s/2}, {p1*p3, t/2}, {p2*p3, u/2}}
    }] ===
  {x1 + x2 + x3 + x4, s x1 x3 + s x1 x4 + t x1 x4 + u x1 x4 + u x2 x4, {x1, x2, x3, x4}}]

FailUnless[
  Feynson["zero-sectors -q -", { {(q-l)^2, l^2}, {0, 0}, {l}, {} }] ===
  {0, 1, 2}]

FailUnless[
  Feynson["zero-sectors -q -s -", { {(q-l)^2, l^2}, {0, 0}, {l}, {} }] ===
  {1, 2}]

FailUnless[
  Feynson["symmetrize -q -", { {{(q-l)^2, l^2}, {(q+l)^2, l^2}}, {l}, {} }] ===
  {{}, {{l, -l}}}]

FailUnless[
  Feynson["mapping-rules -q -", { {{(q-l)^2, l^2}, {(q+l)^2, l^2}}, {l}, {} }] ===
  {{}, {1, {1, 2}}}]

FailUnless[
  Feynson["symmetrize -q -", {{
      {(l-p1)^2, (l+p2)^2},
      {(l)^2, (l+p1+p2)^2},
      {(l)^2, (l-p1)^2, (l-p1-p2)^2},
      {(l)^2, (l-p2)^2, (l-p1-p2)^2}
    }, {l}, {} }] ===
  {{{l, -l+p1}}, {{l, -l}}, {}, {{l, -l + p1 + p2}}}]

FailUnless[
  Feynson["symmetrize -q -d -", {{
      {(l-p1)^2, (l+p2)^2},
      {(l)^2, (l+p1+p2)^2},
      {(l)^2, (l-p1)^2, (l-p1-p2)^2},
      {(l)^2, (l-p2)^2, (l-p1-p2)^2}
    }, {l}, {} }] ===
  {{}, {{l, l-p1}}, {}, {{l, -l + p1 + p2}}}]

FailUnless[
  Feynson["mapping-rules -q -", {{
      {(l-p1)^2, (l+p2)^2},
      {(l)^2, (l+p1+p2)^2},
      {(l)^2, (l-p1)^2, (l-p1-p2)^2},
      {(l)^2, (l-p2)^2, (l-p1-p2)^2}
    }, {l}, {} }] ===
  {{3, {1, 0, 2}}, {3, {1, 0, 2}}, {}, {3, {3, 2, 1}}}]

FailUnless[
  Feynson["mapping-rules -q -d -", {{
      {(l-p1)^2, (l+p2)^2},
      {(l)^2, (l+p1+p2)^2},
      {(l)^2, (l-p1)^2, (l-p1-p2)^2},
      {(l)^2, (l-p2)^2, (l-p1-p2)^2}
    }, {l}, {} }] ===
  {{}, {1, {1, 2}}, {}, {3, {3, 2, 1}}}]
