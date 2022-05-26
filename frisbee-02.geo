L = 10;
Point(1) = {3, 0.2, 0};
Point(2) = {2, 0.2, 0};
Point(3) = {1, 0.25, 0};
Point(4) = {0, 0.255, 0};
Point(5) = {-1, 0.25, 0};
Point(6) = {-2, 0.2, 0};
Point(7) = {-3, 0.2, 0};
Point(8) = {-2.5, -0.2, 0};
Point(9) = {-1, -0.2, 0};
Point(10) = {0, -0.2, 0};
Point(11) = {1, -0.2, 0};
Point(12) = {2.5, -0.2, 0};
Point(13) = {3, -0.9, 0};
Point(14) = {-3, -0.9, 0};
Point(15) = {2.5, -0.9, 0};
Point(16) = {-2.5, -0.9, 0};

af = newl; BSpline(af) = {10, 11, 12, 15, 13, 1, 2, 3, 4, 5, 6, 7, 14, 16, 8, 9, 10};

Transfinite Curve(af) = 301;
cint = newl; Curve Loop(cint) = {af};

p1 = newp; Point(p1) = {-L, -L, 0};
p2 = newp; Point(p2) = {L, -L, 0};
p3 = newp; Point(p3) = {L, L, 0};
p4 = newp; Point(p4) = {-L, L, 0};
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};

cext =  newl; Curve Loop(cext) = {l1, l2, l3, l4};

s1 = news; Plane Surface(s1) = {cext, cint};

Physical Surface("Aire") = {s1};

Mesh 2;


Save "frisbee-02.msh";

Mesh.SurfaceFaces = 1;
Mesh.Points = 1;