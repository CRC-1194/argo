    argList::addOption
    (
        "fieldName", 
        "alpha.water",
        "Name of the volume fraction field." 
    ); 

    argList::addOption
    (
        "sqrDistFactor", 
        "scalar",
        "Expand (>1) or contract (<1) the narrow band field around the interface."
    ); 

    argList::addOption
    (
        "surfaceFile",
        "Surface mesh file."
    );

    argList::addBoolOption
    (
        "writeGeometry",
        "Write the intersected geometry used to compute the volume fractions."
    );

    argList::addBoolOption
    (
        "checkVolume",
        "Check the volume given by volume fraction compared with the volume of the surface mesh. Works only with closed surfaces."
    );

    argList::addBoolOption
    (
        "writeAllFields", 
        "Write all fields used for the initialization. Used when debugging." 
    ); 
