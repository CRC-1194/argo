    argList::addNote
    (
        "Test the accuracy of the vof initialization passed on polynomial approximation using signed distances."
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

    argList::addOption
    (
        "dataFile", 
        "polynomialVofInitResults.csv",
        "Name of the file to write evaluation results to."
    ); 

    argList::addBoolOption
    (
        "writeFields",
        "Write randomly placed surface and fields for each run."
    );

    argList::addOption
    (
        "surfaceVolume", 
        "scalar",
        "Enclosed volume of the given surface. Used to calculate volume errors."
    ); 

    argList::addBoolOption
    (
        "keepOriginalInterfacePosition",
        "Place interface as given by the surface file."
    );

