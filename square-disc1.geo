//square-disc1

mesh_size  = 1.0;
refined_mesh_size = 0.1;
radio = 1;

Point(1) = {0, 0, 0, mesh_size};
//+
Point(2) = {10, 0, 0, mesh_size};
//+
Point(2) = {10, 0, 0, mesh_size};
//+
Point(3) = {10, 10, 0, mesh_size};
//+
Point(4) = {0, 10, 0, mesh_size};
//+
Point(5) = {5, 2.5, 0};
//+
Point(6) = {5 + radio, 2.5, 0, refined_mesh_size};
//+
Point(7) = {5 - radio, 2.5, 0, refined_mesh_size};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {6, 5, 7};
//+
Circle(6) = {7, 5, 6};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5, 6};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("down", 1) = {1};
//+
Physical Curve("up", 2) = {3};
//+
Physical Curve("left", 3) = {4};
//+
Physical Curve("right", 4) = {1};
Physical Surface("mi superficie", 5) = {1};
