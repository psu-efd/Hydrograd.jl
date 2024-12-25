// Gmsh .geo file example for a rectangular domain

// Parameters
Lx = 25;    // Length of the rectangle along the x-axis
Ly = 1;    // Length of the rectangle along the y-axis
nx = 20;   // Number of divisions along x-axis for structured mesh
ny = 1;   // Number of divisions along y-axis for structured mesh

// Define points for the rectangle
Point(1) = {0, 0, 0, 1.0};
Point(2) = {Lx, 0, 0, 1.0};
Point(3) = {Lx, Ly, 0, 1.0};
Point(4) = {0, Ly, 0, 1.0};

// Define lines for the rectangle
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Define a closed loop and surface for the rectangle
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Define the structured mesh on the surface
Transfinite Line{1, 3} = nx + 1;    // x-direction divisions
Transfinite Line{2, 4} = ny + 1;    // y-direction divisions
Transfinite Surface{1};              // Structured mesh for the surface
Recombine Surface{1};                // Recombine to get quadrilateral elements

// Define physical groups (for boundary conditions and regions)
Physical Line("Inlet") = {4};   
Physical Line("Outlet") = {2};  
Physical Line("Side1") = {1};  
Physical Line("Side2") = {3};  
Physical Surface("Channel") = {1};  // Entire domain

// Mesh generation settings (optional)
Mesh.Algorithm = 5;  // Structured mesh algorithm
Mesh.Format = 1;     // Export mesh in MSH format version 2.2

//+
Show "*";
