// -----------------------------------------------------------------------------
//
//  Gmsh GEO TEAM 3 Coil mesh file
//  Author: Karthik Chockalingam
//
// -----------------------------------------------------------------------------
SetFactory("OpenCASCADE");
//Mesh parameters
Mesh.Algorithm = 6;
Mesh.CharacteristicLengthMin = 0.002;
Mesh.CharacteristicLengthMax = 0.020;
Mesh.PartitionOldStyleMsh2 = 1;
Mesh.PreserveNumberingMsh2 = 1;

// Parametric sizes of Geometry
h = 0.002;

// plate parameters
Offset = 0.01;
Hollow = 0.04;

// -----------------------------------------------------------------------------
// Plate Geometry
// -----------------------------------------------------------------------------

// Outer Points on coil
Point(1) = {0.0, 0.0,  0.0,  h};
Point(2) = {0.06, 0.0, 0.0, h};
Point(3) = {0.01, 0.02, 0.0, h};
Point(4) = {0.05, 0.02, 0.0, h};
Point(5) = {0.01, 0.05, 0.0, h};
Point(6) = {0.05, 0.05, 0.0, h};
Point(7) = {0.01, 0.06, 0.0, h};
Point(8) = {0.05, 0.06, 0.0, h};
Point(9) = {0.01, 0.09, 0.0, h};
Point(10) = {0.05, 0.09, 0.0, h};
Point(11) = {0.0, 0.11, 0.0, h};
Point(12) = {0.06, 0.11, 0.0, h};

// Connect lines
Line(1) = {1, 2};  
Line(2) = {2, 12};  
Line(3) = {12, 11};  
Line(4) = {11, 1};  
Line(5) = {3, 4};  
Line(6) = {4, 6};  
Line(7) = {6, 5};  
Line(8) = {5, 3}; 
Line(9) = {7, 8};  
Line(10) = {8, 10};  
Line(11) = {10, 9};  
Line(12) = {9, 7}; 

//Line Loops
Curve Loop(13) = {1, 2, 3, 4};
Curve Loop(14) = {5, 6, 7, 8};
Curve Loop(15) = {9, 10, 11, 12};

Plane Surface(16) = {13};
Plane Surface(17) = {14};
Plane Surface(18) = {15};

// -----------------------------------------------------------------------------
// Coil Geometry
// (Defined before the BooleanDifference so the explicit tags 13-21 are claimed
//  first; otherwise the OpenCASCADE boolean grabs tags 13-20 for its generated
//  boundary entities and collides with these on newer gmsh versions.)
// -----------------------------------------------------------------------------


// Center
Point(13) = {0.030, 0.035, 0.02135, h};

// Points on the circle
Point(14) = {-0.010, 0.035, 0.02135, h};
Point(15) = {0.030, -0.005, 0.02135, h};
Point(16) = {0.070, 0.035, 0.02135, h};
Point(17) = {0.030, 0.075, 0.02135, h};

// Outer Circle arcs
Circle(13) = {14, 13, 15};
Circle(14) = {15, 13, 16};
Circle(15) = {16, 13, 17};
Circle(16) = {17, 13, 14};

// Points on the circle
Point(18) = {0.010, 0.035, 0.02135, h};
Point(19) = {0.030, 0.015, 0.02135, h};
Point(20) = {0.050, 0.035, 0.02135, h};
Point(21) = {0.030, 0.055, 0.02135, h};

// inner Circle arcs
Circle(17) = {18, 13, 19};
Circle(18) = {19, 13, 20};
Circle(19) = {20, 13, 21};
Circle(20) = {21, 13, 18};


Line(101) = {17, 21};
Line(102) = {19, 15};

Curve Loop(16) = {101, -19, -18, 102, 14, 15};
Curve Loop(17) = {101, 20, 17, 102, -13, -16};

Plane Surface(19) = {16};
Plane Surface(20) = {17};

BooleanDifference{ Surface{16}; Delete; }{ Surface{17,18}; Delete; }

// Extrude first-half

Extrude {0, 0, 0.020} {
  Surface{19};
}

// Extrude second-half

Extrude {0, 0, 0.020} {
  Surface{20};
}

// Extrude plate

Extrude {0, 0, 0.00635} {
  Surface{16};
}

Coherence;

Box(4) = {-0.0415, -0.0230, -0.028135, 0.138, 0.156, 0.09762};
BooleanFragments{ Volume{4}; Delete; }{ Volume{1}; Volume{2}; Volume{3}; Delete; }

Physical Volume("First_half")  = {1}; //first half
Physical Volume("Second_half")  = {2}; //second half
Physical Volume("Plate")  = {3}; //plate
Physical Volume("Free_space")  = {4}; //free space

Physical Surface("Cut") = {48};
Physical Surface("Boundary") = {62, 63, 64, 65, 66, 67}; //free space boundary
Physical Surface("Coil_Surface") = {49, 50, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61}; //Coil surface

Mesh.MshFileVersion = 2.2;