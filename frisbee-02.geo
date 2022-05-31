mesh_size  = 1.0;
refined_mesh_size = 1.0;
radius = 2.0;
Side = 10.0;

Point(1) = {-Side, -Side, 0, mesh_size};
//+
Point(2) = {Side, -Side, 0, mesh_size};
//+
Point(3) = {Side, Side, 0, mesh_size};
//+
Point(4) = {-Side, Side, 0, mesh_size};
//+
Point(5) = {3, 0.2, 0};
//+
Point(6) = {2, 0.2, 0};
//+
Point(7) = {1, 0.25, 0};
//+
Point(8) = {0, 0.255, 0};
//+
Point(9) = {-1, 0.25, 0};
//+
Point(10) = {-2, 0.2, 0};
//+
Point(11) = {-3, 0.2, 0};
//+
Point(12) = {-2.5, -0.2, 0};
//+
Point(13) = {-1, -0.2, 0};
//+
Point(14) = {0, -0.2, 0};
//+
Point(15) = {1, -0.2, 0};
//+
Point(16) = {2.5, -0.2, 0};
//+
Point(17) = {3, -0.9, 0};
//+
Point(18) = {-3, -0.9, 0};
//+
Point(19) = {2.5, -0.9, 0};
//+
Point(20) = {-2.5, -0.9, 0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
BSpline(5) = {14, 15, 16, 19, 17, 5, 6, 7, 8, 9, 10, 11, 18, 20, 12, 13, 14};
//+
Transfinite Curve(5) = 301;
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("down", 1) = {1};
//+
Physical Curve("up", 2) = {3};
//+
Physical Curve("left", 3) = {4};
//+
Physical Curve("right", 4) = {2};
//+
Physical Curve("Ring", 5) = {5};
//+
Physical Surface("Wind",6) = {1};