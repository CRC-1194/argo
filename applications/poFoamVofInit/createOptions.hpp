    argList::addNote
    (
        "Initialize a volume fraction field from a triangulated surface using a polynomial approximation."
    );

    argList::addOption
    (
        "fieldName", 
        "alpha.water",
        "Name of the volume fraction field." 
    ); 

    argList::addOption
    (
        "refinementLevel", 
        "label",
        "Maximum refinement level to be used. Uses 'auto' -mode if not given."
    ); 

    argList::addOption
    (
        "surfaceFile",
        "Surface mesh file."
        "STL file containing the interface description. Requires a consistent, inward normal orientation. Use OpenFOAM's 'surfaceOrient' tool for this purpose."
    );

    argList::addBoolOption
    (
        "writeFields",
        "Write randomly placed surface and fields for each run."
    );

    argList::addBoolOption
    (
        "invert",
        "Invert the computed volume fraction field: cells outside of the surface are set to 1."
    );
