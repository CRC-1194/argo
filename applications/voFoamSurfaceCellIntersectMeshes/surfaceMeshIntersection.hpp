/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    voFoamSetSurfaceFraction 

Description

    surface / mesh intersection helping functions. FIXME: Move to src. TM. 

Author
    Tomislav Maric maric@mma.tu-darmstadt.de, tomislav@sourceflux.de

\*---------------------------------------------------------------------------*/

#ifndef surfaceMeshIntersection_H
#define surfaceMeshIntersection_H

void calcSearchFields(
    volScalarField& sqrSearchDist, 
    volScalarField& sqrCellRadius,
    const scalar distanceFactor, 
    const scalar radiusFactor 
)
{
    const fvMesh& mesh = sqrSearchDist.mesh();

    // Sum deltaCoeffs inversed.
    const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    // Sum the deltaCoeffs for the internal faces.
    forAll(own, I)
    {
        const scalar sqrDist = (1 / (deltaCoeffs[I] * deltaCoeffs[I]));

        sqrCellRadius[own[I]] = max(
            sqrCellRadius[own[I]],
            sqrDist
        );

        sqrCellRadius[nei[I]] = max(
            sqrCellRadius[nei[I]],
            sqrDist 
        );
    }

    // Sum the deltaCoeffs for the boundary faces.
    forAll(mesh.boundary(), patchI)
    {
        const fvsPatchField<scalar>& deltaCoeffsBoundary =
            deltaCoeffs.boundaryField()[patchI];

        const labelList& faceCells =
            mesh.boundary()[patchI].faceCells();

        forAll(mesh.boundary()[patchI], faceI)
        {
            sqrCellRadius[faceCells[faceI]] = max(
                sqrCellRadius[faceCells[faceI]], 
                (2 / (deltaCoeffsBoundary[faceI] * deltaCoeffsBoundary[faceI]))
            );
        }
    }

    // Extend the search distance field with the narrow band to ensure
    // smooth distance field solutions.
    sqrSearchDist == sqrCellRadius * sqr(distanceFactor);

    // Extend the cell search radius by 2 in case of a non-uniform mesh.
    sqrCellRadius *= sqr(radiusFactor); 
}

#endif
