L = 1;
Point(1) = {1, 0, 0};
Point(2) = {0.8, 0.027, 0};
Point(3) = {0.55, 0.058, 0};
Point(4) = {0.3, 0.08, 0};
Point(5) = {0.15, 0.078, 0};
Point(6) = {0, 0.052, 0};
Point(7) = {0, 0, 0};
Point(8) = {0, -0.052, 0};
Point(9) = {0.15, -0.078, 0};
Point(10) = {0.3, -0.08, 0};
Point(11) = {0.55, -0.058, 0};
Point(12) = {0.8, -0.027, 0};

af = newl; BSpline(af) = {1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 1};

Transfinite Curve(af) = 301;
cint = newl; Curve Loop(cint) = {af};

p1 = newp; Point(p1) = {0.5-L, -L, 0};
p2 = newp; Point(p2) = {0.5+L, -L, 0};
p3 = newp; Point(p3) = {0.5+L, L, 0};
p4 = newp; Point(p4) = {0.5-L, L, 0};
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};

cext =  newl; Curve Loop(cext) = {l1, l2, l3, l4};

s1 = news; Plane Surface(s1) = {cext, cint};

Physical Surface("Aire") = {s1};

Mesh 2;


Save "airfoil1.msh";

Mesh.SurfaceFaces = 1;
Mesh.Points = 1;