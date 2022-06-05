//****************
//** Parameters **
//****************
//+
R = 0.8; // Radius of border of the RBC
//+
Lcent = 20; // Length of the central part of the RBC
//+
depr = R/4; //magnitude of the depression of central part of the RBC
//+
th = 2*R/5; //Thickness of the cell
//+
lcm = R/4; //Characteristic length of the mesh
//+
mesh_size = 5; //characteristic length near the bounders
//+
theta = Pi/6;
//**************
//** Geometry **
//**************
//+
x1 = -Lcent/2;
y1 = -R;		 
Point(1) = {x1*Cos(theta) + y1*Sin(theta), -x1*Sin(theta) + y1*Cos(theta), 0};
//+
x2 = -Lcent/2;
y2 = 0;
Point(2) = {x2*Cos(theta) + y2*Sin(theta), -x2*Sin(theta) + y2*Cos(theta), 0};
//+
x3 = -Lcent/2;
y3 = R;
Point(3) = {x3*Cos(theta) + y3*Sin(theta), -x3*Sin(theta) + y3*Cos(theta), 0};
//+
x4 = -R-Lcent/2;
y4 = 0;
Point(4) = {x4*Cos(theta) + y4*Sin(theta), -x4*Sin(theta) + y4*Cos(theta), 0};
//+
Circle(1) = {3, 2, 4};
//+
Circle(2) = {4, 2, 1};
//+
//Top central line (with Spline)
//+
x5 = -Lcent/4;
y5 = R + 1.5*depr;
Point(5) = {x5*Cos(theta) + y5*Sin(theta), -x5*Sin(theta) + y5*Cos(theta), 0, lcm};
//+
x6 = 0;
y6 = R + 2*depr;
Point(6) = {x6*Cos(theta) + y6*Sin(theta), -x6*Sin(theta) + y6*Cos(theta), 0, lcm};
//+
Spline(3) = {3, 5, 6};
//+
//Bottom central line
//+
//inner frisbee
//+
x9 = -Lcent/2;
y9 = -R + th;
Point(9) = {x9*Cos(theta) + y9*Sin(theta), -x9*Sin(theta) + y9*Cos(theta), 0, lcm};
//+
x10 = -Lcent/2;
y10 = 0;
Point(10) = {x10*Cos(theta) + y10*Sin(theta), -x10*Sin(theta) + y10*Cos(theta), 0, lcm};
//+
x11 = -Lcent/2;
y11 = R - th;
Point(11) = {x11*Cos(theta) + y11*Sin(theta), -x11*Sin(theta) + y11*Cos(theta), 0, lcm};
//+
x12 = -R + th-Lcent/2;
y12 = 0;
Point(12) = {x12*Cos(theta) + y12*Sin(theta), -x12*Sin(theta) + y12*Cos(theta), 0, lcm};
//+
Circle(4) = {11, 10, 12};
//+
Circle(5) = {12, 10, 9};
//+
//Top central line (with Spline)
//+
x13 = -Lcent/4;
y13 = R - th*0.7 +  1.5*depr;
Point(13) = {x13*Cos(theta) + y13*Sin(theta), -x13*Sin(theta) + y13*Cos(theta), 0, lcm};
//+
x14 = 0;
y14 = R - th + 2*depr;
Point(14) = {x14*Cos(theta) + y14*Sin(theta), -x14*Sin(theta) + y14*Cos(theta), 0, lcm};
//+
Spline(6) = {11, 13, 14};
//+
x15 = -Lcent/2;
y15 = -R + th/2;
Point(15) = {x15*Cos(theta) + y15*Sin(theta), -x15*Sin(theta) + y15*Cos(theta), 0, lcm};
//+
x16 = th/2-Lcent/2;
y16 = -R + th/2;
Point(16) = {x16*Cos(theta) + y16*Sin(theta), -x16*Sin(theta) + y16*Cos(theta), 0, lcm};
//+
Circle(7) = {1, 15, 16};
//+
Circle(8) = {16, 15, 9};
//+
//sección derecha
//+
x17 = Lcent/2;
y17 = -R;		 
Point(17) = {x17*Cos(theta) + y17*Sin(theta), -x17*Sin(theta) + y17*Cos(theta), 0};
//+
x18 = Lcent/2;
y18 = 0;
Point(18) = {x18*Cos(theta) + y18*Sin(theta), -x18*Sin(theta) + y18*Cos(theta), 0};
//+
x19 = Lcent/2;
y19 = R;
Point(19) = {x19*Cos(theta) + y19*Sin(theta), -x19*Sin(theta) + y19*Cos(theta), 0};
//+
x20 = R+Lcent/2;
y20 = 0;
Point(20) = {x20*Cos(theta) + y20*Sin(theta), -x20*Sin(theta) + y20*Cos(theta), 0};
//+
Circle(9) = {19, 18, 20};
//+
Circle(10) = {20, 18, 17};
//+
//Top central line (with Spline)
//+
x21 = Lcent/4;
y21 = R + 1.5*depr;
Point(21) = {x21*Cos(theta) + y21*Sin(theta), -x21*Sin(theta) + y21*Cos(theta), 0, lcm};
//+
Spline(11) = {19, 21, 6};
//+
//Bottom central line
//+
//inner frisbee
//+
x23 = Lcent/2;
y23 = -R + th;
Point(23) = {x23*Cos(theta) + y23*Sin(theta), -x23*Sin(theta) + y23*Cos(theta), 0, lcm};
//+
x24 = Lcent/2;
y24 = 0;
Point(24) = {x24*Cos(theta) + y24*Sin(theta), -x24*Sin(theta) + y24*Cos(theta), 0, lcm};
//+
x25 = Lcent/2;
y25 = R - th;
Point(25) = {x25*Cos(theta) + y25*Sin(theta), -x25*Sin(theta) + y25*Cos(theta), 0, lcm};
//+
x26 = R - th+Lcent/2;
y26 = 0;
Point(26) = {x26*Cos(theta) + y26*Sin(theta), -x26*Sin(theta) + y26*Cos(theta), 0, lcm};
//+
Circle(12) = {25, 24, 26};
//+
Circle(13) = {26, 24, 23};
//+
//Top central line (with Spline)
//+
x27 = Lcent/4;
y27 = R - th*0.7 +  1.5*depr;
Point(27) = {x27*Cos(theta) + y27*Sin(theta), -x27*Sin(theta) + y27*Cos(theta), 0, lcm};
//+
Spline(14) = {25, 27, 14};
//+
x29 = Lcent/2;
y29 = -R + th/2;
Point(29) = {x29*Cos(theta) + y29*Sin(theta), -x29*Sin(theta) + y29*Cos(theta), 0, lcm};
//+
x30 = -th/2+Lcent/2;
y30 = -R + th/2;
Point(30) = {x30*Cos(theta) + y30*Sin(theta), -x30*Sin(theta) + y30*Cos(theta), 0, lcm};
//+
Circle(15) = {17, 29, 30};
//+
Circle(16) = {30, 29, 23};
//+
//surroundings
//+
Point(31) = {-50, -20, 0, mesh_size};
//+
Point(32) = {50, -20, 0, mesh_size};
//+
Point(33) = {50, 20, 0, mesh_size};
//+
Point(34) = {-50, 20, 0, mesh_size};
//+
Line(17) = {31, 32};
//+
Line(18) = {32, 33};
//+
Line(19) = {33, 34};
//+
Line(20) = {34, 31};
//+
//**********
//** Mesh **
//**********
//Number of divisions on the circular lines
//+
Transfinite Curve {1, 2, 4, 8} = 10*R;
//+
//define surface
//+
Curve Loop(1) = {17, 18, 19, 20};
//+
Curve Loop(2) = {1, 2, 7, 8, -5, -4, 6, -14, 12, 13, -16, -15, -10, -9, 11, -3};
//+
Plane Surface(1) = {1, 2};
//+
//+
Physical Curve("down", 1) = {17};
//+
Physical Curve("up", 2) = {19};
//+
Physical Curve("left", 3) = {20};
//+
Physical Curve("right", 4) = {18};
//+
Physical Curve("Frisbee", 5) = {1, 2, 7, 8, -5, -4, 6, -14, 12, 13, -16, -15, -10, -9, 11, -3};
//+
Physical Surface("Wind",6) = {1};
//+
Mesh 2;
//+
Mesh.SurfaceFaces = 1;
//+
Mesh.Points = 1;