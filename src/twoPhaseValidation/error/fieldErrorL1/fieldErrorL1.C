/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | Copyright 2011 Tomislav Maric 
     \\/     M anipulation  |
-------------------------------------------------------------------------------

License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    A class template for L1 VoF advection error calculation. 

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#include "meshMagnitude.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename GeomField>
Foam::fieldErrorL1<GeomField>::fieldErrorL1(scalar errorTolerance)
:
    fieldError<GeomField>(errorTolerance),
    EgTmp_()
{}

template<typename GeomField>
Foam::fieldErrorL1<GeomField>::fieldErrorL1(const fieldErrorL1<GeomField>& copy)
:
    fieldError<GeomField>(copy),
    EgTmp_(copy.EgTmp_->clone())
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<typename GeomField>
Foam::fieldErrorL1<GeomField>::~fieldErrorL1()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<typename GeomField>
void Foam::fieldErrorL1<GeomField>::computeError(
    const GeomField& initialField,
    const GeomField& currentField
)
{
    const scalarField& M = getMeshMagnitude(initialField);

    // For 2D simulations on an fvMesh, the L1 error is multiplied by a face
    // area, not a volume. 
    scalar cellLayerHeight = -VGREAT; 
    if (isA<fvMesh>(initialField.mesh()))
    {
        const fvMesh& mesh = initialField.mesh(); 
        const labelList& faceOwner = mesh.faceOwner(); 
        const cellList& cells = mesh.cells(); 
        const faceList& faces = mesh.faces(); 
        const pointField& points = mesh.points(); 

        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& fvp = mesh.boundary()[patchI]; 
            if (fvp.type() == "empty")
            {
                const face& firstFace = faces[fvp.start()]; 
                vector firstFaceNormal = -1 * firstFace.unitNormal(points); 
                const vector firstFaceCenter = firstFace.centre(points);  
                const cell& firstCell = cells[faceOwner[fvp.start()]]; 
                const pointField firstCellPoints = firstCell.points(faces, points);  

                forAll (firstCellPoints, pointI)
                {
                    const point& cellPoint = firstCellPoints[pointI]; 
                    scalar pointHeight = firstFaceNormal & (cellPoint - firstFaceCenter); 
                    if (pointHeight > cellLayerHeight)
                        cellLayerHeight = pointHeight;
                }
                break;
            }
        }
    }

    //if (EgTmp_.empty())
    if (!EgTmp_)
    {
        const auto& mesh = initialField.mesh();  
        const auto& runTime = mesh.time(); 
        EgTmp_ = tmp<volScalarField>(
            new volScalarField(
                IOobject(
                    "Eg", 
                    runTime.timeName(),
                    mesh, 
                    IOobject::NO_READ, 
                    IOobject::AUTO_WRITE
                ),
                mesh, 
                dimensionedScalar("0", initialField.dimensions(), 0) 
            )
        );
    }

    auto magField = mag(initialField - currentField); 

    scalarField& Eg = EgTmp_.ref(); 
    // Scale the error with the volume.
    forAll(Eg, cellI)
        Eg[cellI] = M[cellI] * magField()[cellI]; 

    // Sum the total error.
    auto errorSum = gSum(Eg);
    
    if (cellLayerHeight > 0)
        this->setErrorValue(errorSum / cellLayerHeight);
    else
    {
        this->setErrorValue(errorSum);
    }
}

// ************************************************************************* //
