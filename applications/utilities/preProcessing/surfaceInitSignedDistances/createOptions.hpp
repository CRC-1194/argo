    argList::addNote
    (
        "Initialize a signed distance field from an oriented, triangulated surface or a level set."
    );

    argList::addOption
    (
        "fieldName", 
        "signedDistance",
        "Name of the signed distance field. Default: signedDistance" 
    ); 

    argList::addOption
    (
        "surfaceType",
        "triSurface / levelSet"
        "Surface type, meaning either a triSurface, e.g. from an STL, or a level set. Default: triSurface"
    );

    argList::addOption
    (
        "surfaceFile",
        "Surface mesh file."
        "STL or VTK file containing the interface description. Requires a consistent normal orientation. Use OpenFOAM's 'surfaceOrient' tool for this purpose. Default: surface.stl"
    );

    argList::addOption
    (
        "narrowBandWidth",
        "scalar",
        "Enable narrow band around the interface and set its thickness in number of cells on each side. Set to zero or negative to compute signed distance in the entire domain rather than in a narrowband. Default: -1."
    );

    argList::addOption
    (
        "bulkValue",
        "scalar",
        "Set the distance of cells outside of the narrow band to this value. Default: 0.0."
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

    argList::addBoolOption
    (
        "writeAllFields",
        "Write all fields used for signed distance calculation."
    );
