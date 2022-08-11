
    reconstructionSchemes* surf = 
        mesh().getObjectPtr<reconstructionSchemes>("reconstructionScheme");

vector sphereCentre(0.005, 0.005, 0.005);
scalar sphereRadius = 0.002; // Sphere radius

    const volVectorField& interfaceNormals = surf->normal();

    // Mark interface markers
    volScalarField markers
    (
        IOobject
        (
            "cellMarkers",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("cellMarkers", dimless, -1)
    );

    forAll (markers, cellI)
    {
        if (mag(interfaceNormals[cellI]) != 0)
        {
            markers[cellI] = 0;
        }
    }
    markers.correctBoundaryConditions();

    const auto& own = mesh().owner();
    const auto& nei = mesh().neighbour();

    zoneDistribute& distribute = zoneDistribute::New(mesh());

    for (label i = 0; i < nPropagate_; ++i)
    {
        forAll (nei, fid)
        {
            if (markers[own[fid]] == i && markers[nei[fid]] == -1)
            {
                markers[nei[fid]] = i + 1;
            }
            if (markers[nei[fid]] == i && markers[own[fid]] == -1)
            {
                markers[own[fid]] = i + 1;
            }
        }
        markers.correctBoundaryConditions();

        auto& meshBoundary = markers.boundaryFieldRef();
        for (auto& mb : meshBoundary)
        {
            if (isA<processorFvPatch>(mb.patch()))
            {
                mb.initEvaluate();
                mb.evaluate();

                const auto& faceToCell = mb.patch().faceCells();
                const auto& neibrValue = mb.patchNeighbourField().cref();
                
                forAll(mb, j)
                {
                    if (neibrValue[j] == i && markers[faceToCell[j]] == -1)
                    {
                        markers[faceToCell[j]] = i + 1;
                    }
                }
            }
        }
        markers.correctBoundaryConditions();
    }

    const boolList& cellTooCoarse = surf->cellTooCoarse();

    volVectorField averagedNormals_ = interfaceNormals /
    (
        mag(interfaceNormals) + 
        dimensionedScalar(
            "SMALL", interfaceNormals.dimensions(), SMALL
        )
    );

    forAll (averagedNormals_, i)
    {
        if (cellTooCoarse[i])
        {
            averagedNormals_[i] = vector::zero;
            cellCurvature_[i] = 0.0;
        }
    }
/*
*/

    averagedNormals_.correctBoundaryConditions();
    cellCurvature_.correctBoundaryConditions();

    boolList updateZone(mesh().nCells(), false);
    forAll (updateZone, uzi)
    {
        if (markers[uzi] == -1) continue;
        if (markers[uzi] ==  0) continue;
        updateZone[uzi] = true;
    }
    distribute.updateStencil(updateZone);

{
    boolList zone(mesh().nCells(), false);
    forAll (zone, zi)
    {
        if (markers[zi] == 0 && cellTooCoarse[zi])
        {
            zone[zi] = true;
        }
    }

    distribute.setUpCommforZone(zone, false);
    Map<vector> mapMC = 
        distribute.getDatafromOtherProc(zone, mesh().C());
    Map<vector> mapNormals = 
        distribute.getDatafromOtherProc(zone, averagedNormals_);

    const labelListList& stencil = distribute.getStencil();

    volVectorField avgNormTmp = averagedNormals_;
    forAll (markers, cellI)
    {
        if (markers[cellI] != 0) continue;
        if (!cellTooCoarse[cellI]) continue;

        avgNormTmp[cellI] = Zero;

        point p = mesh().C()[cellI];

        DynamicField<vector> centres;
        DynamicField<scalar> valuesX;
        DynamicField<scalar> valuesY;
        DynamicField<scalar> valuesZ;

        for (const label gblIdx : stencil[cellI])
        {
            vector n = distribute.getValue(averagedNormals_, mapNormals, gblIdx);

            if (mag(n) != 0)
            {
                n /= mag(n);

                vector centre = distribute.getValue(mesh().C(), mapMC, gblIdx);

                vector dist = centre - p;
                vector distToSurf = dist & n / mag(n) * n;
                vector verticalDist = dist - distToSurf;

                vector cc = p - verticalDist;

                centres.append(cc);
                valuesX.append(n.x());
                valuesY.append(n.y());
                valuesZ.append(n.z());
            }
        }

        avgNormTmp[cellI][0] = interpolateSecondOrder(p, centres, valuesX, rr);
        avgNormTmp[cellI][1] = interpolateSecondOrder(p, centres, valuesY, rr);
        avgNormTmp[cellI][2] = interpolateSecondOrder(p, centres, valuesZ, rr);
    }
    avgNormTmp.correctBoundaryConditions();

    averagedNormals_ = avgNormTmp / 
    (
        mag(avgNormTmp) + 
        dimensionedScalar(
            "SMALL", avgNormTmp.dimensions(), SMALL
        )
    );
    averagedNormals_.correctBoundaryConditions();
}
/*
*/

    for (label i = 0; i < nPropagate_; ++i)
    {
        boolList zone(mesh().nCells(), false);
        forAll (zone, zi)
        {
            if (markers[zi] == i + 1)
            {
                zone[zi] = true;
            }
        }

        distribute.setUpCommforZone(zone, false);
        Map<vector> mapMC = 
            distribute.getDatafromOtherProc(zone, mesh().C());
        Map<vector> mapNormals = 
            distribute.getDatafromOtherProc(zone, averagedNormals_);
        Map<scalar> mapMarkers = 
            distribute.getDatafromOtherProc(zone, markers);

        const labelListList& stencil = distribute.getStencil();

        volVectorField avgNormTmp = averagedNormals_;
        forAll (markers, cellI)
        {
            if (markers[cellI] != i + 1) continue;

            avgNormTmp[cellI] = Zero;

            point p = mesh().C()[cellI];

            DynamicField<vector> centres;
            DynamicField<scalar> valuesX;
            DynamicField<scalar> valuesY;
            DynamicField<scalar> valuesZ;

            for (const label gblIdx : stencil[cellI])
            {
                //scalar mk = distribute.getValue(markers, mapMarkers, gblIdx);
                //if (mk == -1 || mk == i + 1) continue;

                vector n = distribute.getValue(averagedNormals_, mapNormals, gblIdx);

                if (mag(n) != 0)
                {
                    n /= mag(n);

                    vector centre = distribute.getValue(mesh().C(), mapMC, gblIdx);

                    vector dist = centre - p;
                    vector distToSurf = dist & n / mag(n) * n;
                    vector verticalDist = dist - distToSurf;

                    vector cc = p - verticalDist;

                    centres.append(cc);
                    valuesX.append(n.x());
                    valuesY.append(n.y());
                    valuesZ.append(n.z());
                }
            }

            avgNormTmp[cellI][0] = interpolateSecondOrder(p, centres, valuesX, rr);
            avgNormTmp[cellI][1] = interpolateSecondOrder(p, centres, valuesY, rr);
            avgNormTmp[cellI][2] = interpolateSecondOrder(p, centres, valuesZ, rr);
        }
        avgNormTmp.correctBoundaryConditions();

        averagedNormals_ = avgNormTmp / 
        (
            mag(avgNormTmp) + 
            dimensionedScalar(
                "SMALL", avgNormTmp.dimensions(), SMALL
            )
        );
        averagedNormals_.correctBoundaryConditions();
    }

