    argList::addNote
    (
        "Intersect the base mesh with the tool mesh and store the volume fraction given by the intersection on the base mesh."
    );

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
        "Expand (>1) or contract (<1) the distance field narrow band."
    ); 

    argList::addBoolOption
    (
        "fixNormals",
        "Make the normals consistent. Assumes that the centroid of the surface mesh is its Stern point: use only with star-shaped surfaces." 
    );

    argList::addOption
    (
        "surfaceFile",
        "Surface mesh file."
    );

    argList::addOption
    (
        "dataFile", 
        "surfaceCellMeshIntersection.csv",
        "Name of the volume fraction field." 
    ); 
