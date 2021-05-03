    argList::addNote
    (
        "Initialize a volume fraction field from a triangulated surface or level set using the SMCI/A algorithm"
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
        "File containing the interface description. Requires a consistent, inward normal orientation. Use OpenFOAM's 'surfaceOrient' tool for this purpose. Default: surface.stl"
    );

    argList::addOption
    (
        "narrowBandWidth", 
        "scalar",
        "Number of cells in interface normal direction constituting the narrow band. Default: 4.0"
    ); 

    argList::addOption
    (
        "fieldName", 
        "alpha.water",
        "Name of the volume fraction field. Default: alpha.water" 
    ); 

    argList::addOption
    (
        "algorithm", 
        "SMCI/SMCA",
        "Name of the volume fraction calculator. Default: SMCI" 
    ); 

    argList::addOption
    (
        "refinementLevel", 
        "label",
        "Maximum refinement level to be used (SMCA only). Default: -1 (auto mode)"
    ); 

    argList::addBoolOption
    (
        "writeGeometry",
        "Write the intersected geometry (SMCI) or the testrahedral decomposition (SMCA)used to compute the volume fractions"
    );

    argList::addBoolOption
    (
        "invert",
        "Invert the computed volume fraction field with respect to the surface orientation."
    );

    argList::addBoolOption
    (
        "writeAllFields",
        "Write out all fields used by the initialization method."
    );

    argList::addBoolOption
    (
        "checkVolume",
        "Check the volume given by volume fraction compared with the volume of the surface mesh. Works only with closed surfaces."
    );
