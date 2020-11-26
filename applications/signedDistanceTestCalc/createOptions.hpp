    argList::addNote
    (
        "Test signed distance calculation"
    );

    argList::addOption
    (
        "npoints", 
        "label",
        "Number of points in each dimension."
    ); 

    argList::addOption
    (
        "rsearch", 
        "scalar",
        "Search radius for octree search."
    ); 

    argList::addOption
    (
        "bbfactor", 
        "scalar",
        "Scale bounding box of the surface by this factor. Trial points are generated in the resulting domain."
    ); 

    argList::addOption
    (
        "surfaceFile",
        "Surface mesh file."
        "STL file containing the interface description. Requires a consistent, inward normal orientation. Use OpenFOAM's 'surfaceOrient' tool for this purpose."
    );
