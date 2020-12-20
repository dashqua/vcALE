Point(1) = {-.5, 0, -.5, 1.0};
Point(2) = {.5, 0, -.5, 1.0};
Point(3) = {.5, 0, .5, 1.0};
Point(4) = {-.5, 0, .5, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Point(5) = {-10, 0, -2, 1.0};
Point(6) = {10, 0, -2, 1.0};
Point(7) = {10, 0, 2, 1.0};
Point(8) = {-10, 0, 2, 1.0};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line Loop(1) = {5, 6, 7, 8};
Line Loop(2) = {4, 1, 2, 3};
Plane Surface(1) = {1, 2};

Extrude {0, 6, 0} {
  Surface{1}; 
}

Line Loop(3) = {16, 17, 14, 15};
Plane Surface(51) = {3};

Extrude {0, 6, 0} {
  Line{10}; Line{13}; Line{12}; Line{11}; 
}

Line Loop(4) = {52, 64, 60, 56};
Plane Surface(68) = {4};
Surface Loop(1) = {68, 55, 67, 63, 59, 50, 51};
Volume(2) = {1};

Physical Surface("inlet") = {59, 33};
Physical Surface("outlet") = {67, 25};
Physical Surface("top") = {68};
Physical Surface("bottom") = {1};
Physical Surface("sides") = {21, 55, 29, 63};
Physical Surface("interface") = {49, 51, 41, 37, 45};
Physical Volume("channel") = {2, 1};
