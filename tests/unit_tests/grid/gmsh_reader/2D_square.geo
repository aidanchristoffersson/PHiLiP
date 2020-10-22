//+
SetFactory("OpenCASCADE");
Circle(1) = {-0, 0, 0, 0.5, 0, 2*Pi};
//+
Circle(2) = {0.1, 0.1, 0.0, 0.25, 0, 2*Pi};
//+
Curve Loop(3) = {1};
//+
Curve Loop(4) = {2};
//+
Plane Surface(1) = {3, 4};

// Frontal-Delaunay for quads
Mesh.Algorithm = 8;
Mesh 2;

// 2 : simple full quad
// 3 : blossom full quad
Mesh.RecombinationAlgorithm = 2; // or 3
RecombineMesh;

// 1 : all quads
// 2 : all hex
Mesh.SubdivisionAlgorithm = 1;
RefineMesh;

// (0: none, 1: optimization, 2: elastic+optimization, 3: elastic, 4: fast curving)
Mesh.ElementOrder = 4;
SetOrder 4;
Mesh.HighOrderOptimize = 2;
OptimizeMesh "HighOrder";

Save "2D_square.msh";
