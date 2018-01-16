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
        "fixNormals",
        "Make the normals consistent. Assumes that the surfaceMeshCentroid of the surface mesh is its Stern point: use only with convex or weakly non-convex surfaces." 
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
