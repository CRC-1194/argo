// Parameters to Determine center and radius
x = 0.5;
y = 0.5;
z = 0.5;

r = 0.15;

// Control maximum element size
elemsize = 0.05;

// Define support points
// Center point
Point(1) = {x, y, z};

// Support points in plane parallel to x-y-plane through center
Point(2) = {x-3*r, y, z};
Point(3) = {x, y-2*r, z};
Point(4) = {x+3*r, y, z};
Point(5) = {x, y+2*r, z};

// Support points in plane parallel to x-z-plane through center
Point(6) = {x, y, z-r};
Point(7) = {x, y, z+r};

// Define Ellipse  arcs in plane parallel to x-y-plane
Ellipse (1) = {2, 1, 2, 3};
Ellipse (2) = {3, 1, 4, 4};
Ellipse (3) = {4, 1, 4, 5};
Ellipse (4) = {5, 1, 2, 2};

// Define Ellipse  arcs in plane parallel to x-z-plane
Ellipse (5) = {6, 1, 4, 4};
Ellipse (6) = {4, 1, 4, 7};
Ellipse (7) = {7, 1, 2, 2};
Ellipse (8) = {2, 1, 2, 6};

// Define Ellipse  arcs in plane parallel to y-z-plane
Ellipse (9) = {6, 1, 3, 3};
Ellipse (10) = {3, 1, 3, 7};
Ellipse (11) = {7, 1, 5, 5};
Ellipse (12) = {5, 1, 5, 6};

// Define ellipsoid segments' boundaries:
// Begin in positive sector, iterate counter clockwise, segments with x >= 0
Line Loop(1) = {3, -11, -6};
Line Loop(2) = {3, 12, 5};
Line Loop(3) = {-2, -9, 5};
Line Loop(4) = {-2, 10, -6};

// Begin in sector (x < 0, y > 0, z > 0) and iterate counter clockwise
Line Loop(5) = {-4, -11, 7};
Line Loop(6) = {-4, 12, -8};
Line Loop(7) = {1, -9, -8};
Line Loop(8) = {1, 10, 7};

// Define ellipsoid surface patches
Surface(1) = {1} ;
Surface(2) = {2};
Surface(3) = {3};
Surface(4) = {4};
Surface(5) = {5};
Surface(6) = {6};
Surface(7) = {7};
Surface(8) = {8};

// Test to orient normal vectors consistently
// They should now point inwards
Reverse Surface {1};
Reverse Surface {3};
Reverse Surface {6};
Reverse Surface {8};

// Combine to ellipsoid surface
Physical Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8};

Mesh.CharacteristicLengthMax = elemsize;
Mesh.Algorithm = 6; // Frontal