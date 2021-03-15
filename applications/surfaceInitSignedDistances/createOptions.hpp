    argList::addNote
    (
        "Initialize a signed distance field from an oriented, triangulated surface."
    );

    argList::addOption
    (
        "fieldName", 
        "signedDistance",
        "Name of the signed distance field. Default: signedDistance." 
    ); 

    argList::addOption
    (
        "surfaceFile",
        "Surface mesh file."
        "STL or VTK file containing the interface description. Requires a consistent normal orientation. Use OpenFOAM's 'surfaceOrient' tool for this purpose. Default: surface.stl."
    );

    argList::addOption
    (
        "searchDistanceFactor",
        "scalar",
        "Factor for the search radius computation. Set to zero or negative to compute signed distance in the entire domain rather than in a narrowband. Default: 4.0."
    );

    argList::addBoolOption
    (
        "propagateInsideOutside",
        "Propagate inside/outside information in the entire domain."
    );

    argList::addBoolOption
    (
        "invert",
        "Invert inside / outside as given by the surface normals."
    );
