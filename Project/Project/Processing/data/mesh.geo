w  = 1;
h  = 1;
r  = w/4;
lc = 0.1;

Point(1) = {-w/2, -h/2, 0, lc};
Point(2) = { w/2, -h/2, 0, lc};
Point(3) = { w/2,  h/2, 0, lc};
Point(4) = {-w/2,  h/2, 0, lc};
Point(5) = {-w/2,  r  , 0, lc};
Point(6) = {   0,  r  , 0, lc};
Point(7) = {   0, -r  , 0, lc};
Point(8) = {-w/2, -r  , 0, lc};
Point(9) = {   0,  0  , 0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Circle(6) = {7,9,6};
Line(7) = {7,8};
Line(8) = {8,1};

Curve Loop(1) = {1:5, -6, 7:8};
Plane Surface(1) = {1};



