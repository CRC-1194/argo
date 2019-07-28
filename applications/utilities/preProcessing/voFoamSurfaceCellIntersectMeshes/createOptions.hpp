    argList::addOption
    (
        "fieldName", 
        "alpha.water",
        "Name of the volume fraction field." 
    ); 

    argList::addOption
    (
        "distanceFactor", 
        "scalar",
        "Expand (>1) or contract (<1) the distance field narrow band."
    ); 

    argList::addOption
    (
        "radiusFactor", 
        "scalar",
        "Expand (>1) or contract or contract (<1) the radius used to find triangles for mesh intersection."
    ); 

    argList::addBoolOption
    (
        "checkVolume",
        "Compute the relative difference between the volume given by the set volume fraction and the input surface mesh volume given by a surfaceMeshCentroid triangulation."  
    );

    argList::addBoolOption
    (
        "fixNormals",
        "Make the normals consistent. Assumes that the surfaceMeshCentroid of the surface mesh is its Stern point: use only with convex or weakly non-convex surfaces." 
    );

    argList::addBoolOption
    (
        "writeCutCells",
        "Write the cut cell triangulations in VTK for visualization."
    );

    argList::addOption
    (
        "surfaceFile",
        "Surface mesh file."
    );

    argList::addOption
    (
        "dataFile", 
        "voFoamSetSurfaceFraction.dat",
        "Name of the volume fraction field." 
    ); 

