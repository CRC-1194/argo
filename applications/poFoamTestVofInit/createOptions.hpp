    argList::addNote
    (
        "Test the accuracy of the vof initialization passed on polynomial approximation using signed distances."
    );

    argList::addOption
    (
        "fieldName", 
        "vof-field",
        "Name of the volume fraction field.\nDefault: alpha.water" 
    ); 

    argList::addOption
    (
        "refinementLevel", 
        "label",
        "Maximum refinement level to be used.\nDefault: auto (-1)"
    ); 

    argList::addOption
    (
        "surfaceFile",
        "Surface mesh file."
        "STL file containing the interface description. Requires a consistent, inward normal orientation. Use OpenFOAM's 'surfaceOrient' tool for this purpose.\nDefault: surface.stl"
    );

    argList::addOption
    (
        "dataFile", 
        "results-file",
        "Name of the file to write evaluation results to.\nDefault: polynomialVofInitResults.csv"
    ); 

    argList::addBoolOption
    (
        "writeFields",
        "Write randomly placed surface and fields for each run.\nDefault: false"
    );

    argList::addBoolOption
    (
        "random",
        "Place interface randomly in the bounding box of the domain.\nDefault: false"
    );
