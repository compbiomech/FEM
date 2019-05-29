
(* Solving the heat equation with the Wolfram Engine *)
(* Written by Dr. Kewei Li, TU Graz*)

Needs["NDSolve`FEM`"];

(*mesh size*)
h = 0.01;

mesh = ToElementMesh[Rectangle[{0, 0}, {2, 2}], MeshOrder -> 1,
   MaxCellMeasure -> h];

sol = NDSolveValue[{2 D[u[x, y], x, x] + 2 D[u[x, y], y, y] + 2 x +
      y^2 == NeumannValue[5.0, x == 2 && 0 <= y <= 2],
    DirichletCondition[u[x, y] == 20, x == 0 && 0 <= y <= 2]},
   u, {x, y} \[Element] mesh];

fig2 = DensityPlot[sol[x, y], {x, y} \[Element] mesh, Mesh -> 11,
   MeshFunctions -> {#3 &}, ColorFunction -> "Rainbow",
   BaseStyle -> {FontWeight -> "Normal", FontSize -> 18},
   PlotRange -> All,
   PlotLegends ->
    BarLegend[Automatic,
     LabelStyle -> {Black, FontSize -> 18}],
   LabelStyle -> Directive[Black]];

Export["heateq2d.pdf", fig2];
