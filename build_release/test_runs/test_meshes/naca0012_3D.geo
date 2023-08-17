Include "naca0012.geo";

// parameters
xmax = 150;
ymax = 100;

spanlength = 4; // width of extrusion

n_spanwise = 1;
n_TE = 32;

// Mesh 000
n_inlet = 3;
n_vertical = 15;
r_vertical = 1/0.8;
n_airfoil = 5;
n_wake = 7;
r_wake = 1/0.7;

// // Mesh 001
// n_inlet = 10;
// n_vertical = 20;
// r_vertical = 1/0.8;
// n_airfoil = 15;
// n_wake = 20;
// r_wake = 1/0.7;

// // Mesh 002
// n_inlet = 12;
// n_vertical = 40;
// r_vertical = 1/0.8;
// n_airfoil = 20;
// n_wake = 30;
// r_wake = 1/0.8;

// // Mesh 003
// n_inlet = 10;
// n_vertical = 60;
// r_vertical = 1/0.85;
// n_airfoil = 30;
// n_wake = 60;
// r_wake = 1/0.9;

// borders
Point(1)={-5, ymax, 0, lc1};
Point(2)={1, ymax, 0, lc1};
Point(3)={xmax, ymax, 0, lc1};
Point(4)={xmax, 0, 0, lc1};
Point(5)={xmax, -ymax, 0, lc1};
Point(6)={1, -ymax, 0, lc1};
Point(7)={-5, -ymax, 0, lc1};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,6};
Line(6)={6,7};

//+
Circle(1004) = {7, 197, 1};
//+
Line(1005) = {1, 214};
//+
Line(1006) = {180, 7};
//+
Line(1007) = {2, 101};
//+
Line(1008) = {101, 6};
//+
Line(1009) = {101, 4};
//+
Split Curve {1003} Point {214};
//+
Split Curve {1001} Point {180};
//+
Transfinite Curve {1013, 1010, 1002} = n_inlet Using Progression 1;
//+
Transfinite Curve {1004} = 3 * n_inlet - 2 Using Progression 1;
//+
Transfinite Curve {1005, 1007, 3} = n_vertical Using Progression 1 / r_vertical;
//+
Transfinite Curve {1006, 1008, 4} = n_vertical Using Progression r_vertical;
//+
Transfinite Curve {6, 1012} = n_airfoil Using Bump 1/5;
//+
Transfinite Curve {1011, 1} = n_airfoil Using Bump 1/5;
//+
Transfinite Curve {1009} = n_wake Using Progression r_wake;
//+
Transfinite Curve {1092} = n_wake Using Progression r_wake;
//+
Transfinite Curve {2} = n_wake Using Progression r_wake;
//+
Transfinite Curve {5} = n_wake Using Progression 1 / r_wake;
//+
Curve Loop(1) = {1004, 1005, -1010, -1002, -1013, 1006};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1005, 1011, -1007, -1};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1012, 1006, -6, -1008};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {1007, 1009, -3, -2};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {1008, -5, -4, -1009};
//+
Plane Surface(5) = {5};
//+
Transfinite Surface {5} = {6, 5, 4, 101};
//+
Transfinite Surface {3} = {101, 180, 7, 6};
//+
Transfinite Surface {1} = {7, 1, 214, 180};
//+
Transfinite Surface {2} = {214, 101, 2, 1};
//+
Transfinite Surface {4} = {101, 4, 3, 2};
//+
Recombine Surface {1, 2, 3, 5, 4};
//+
Mesh.RecombinationAlgorithm = 2;
//+
Extrude {0, 0, spanlength} {
	Surface{1, 2, 3, 4, 5}; 
	Layers {n_spanwise}; 
	Recombine;
}
//+
Transfinite Curve {1057} = n_TE Using Progression 1;
//+
Mesh 3;
//+
RecombineMesh;
//+
Mesh.SubdivisionAlgorithm = 1;
//+
Mesh.SecondOrderLinear = 1;
//+
RefineMesh;
//+
Physical Volume("MeshInterior", 32) = {2, 4, 5, 3, 1};
//+
Physical Surface("Airfoil", 1001) = {1076, 1058, 1040, 1036, 1032};
//+
Physical Surface("Farfield", 1004) = {1024, 1084, 1124, 1110, 1066, 1128, 1106};
//+
Physical Surface("SideWall_z0", 2005) = {1, 2, 3, 4, 5};
//+
Physical Surface("SideWall_z1", 2006) = {1045, 1067, 1089, 1111, 1133};
