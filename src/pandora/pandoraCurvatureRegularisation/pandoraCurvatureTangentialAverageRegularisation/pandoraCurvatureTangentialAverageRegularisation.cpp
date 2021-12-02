/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 AUTHOR,AFFILIATION
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

\*---------------------------------------------------------------------------*/

#include "pandoraCurvatureTangentialAverageRegularisation.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvc.H"
#include "processorFvPatchField.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraCurvatureTangentialAverageRegularisation, false);
addToRunTimeSelectionTable(pandoraCurvatureRegularisation, pandoraCurvatureTangentialAverageRegularisation, Dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pandoraCurvatureTangentialAverageRegularisation::pandoraCurvatureTangentialAverageRegularisation(const dictionary& dict)
    :
        pandoraCurvatureRegularisation{dict},
        nAveragingIterations_{dict.get<label>("nAveragingIterations")}
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void pandoraCurvatureTangentialAverageRegularisation::regularise
(
    volScalarField& curvature,
    const volScalarField& isInterfaceCell
)
{
    const auto& mesh = curvature.mesh();
    const auto& neighbour = mesh.neighbour();
    const auto& owner = mesh.owner();

    forAll(curvature, cid)
    {
        if (isInterfaceCell[cid] != 1.0)
        {
            curvature[cid] = 0.0;
        }
    }

    labelField count{curvature.size(), 0};
    scalarField curvatureSum{curvature.size(), 0.0};

    for (label iter = 0; iter != nAveragingIterations_; ++iter)
    {
        auto faceCurvatureTmp = fvc::interpolate(curvature);
        const auto& faceCurvature = faceCurvatureTmp.cref();

        forAll (faceCurvature, fid)
        {
            if ((isInterfaceCell[owner[fid]] == 1.0) && (isInterfaceCell[neighbour[fid]] == 1.0))
            {
                // Curvature okay: face shared by interface cells
                count[owner[fid]] += 1;
                curvatureSum[owner[fid]] += faceCurvature[fid];

                count[neighbour[fid]] += 1;
                curvatureSum[neighbour[fid]] += faceCurvature[fid];
            }
        }

        // Iterate processor boundaries
        const auto& isInterfaceCellBoundary = isInterfaceCell.boundaryField();
        const auto& faceCurvatureBoundary = faceCurvature.boundaryField();

        forAll(isInterfaceCellBoundary, patchID)
        {
            const auto& isInterfaceCellPatch = isInterfaceCellBoundary[patchID];
            const auto& faceCurvaturePatch = faceCurvatureBoundary[patchID];

            if (isType<processorFvPatch>(isInterfaceCellPatch.patch()))
            {
                // Values of isInterfaceCell do not change here. Ensure at computation of
                // isInterFaceCell that processor neighbour fields are up-to-date (TT)
                const auto& faceToCell = isInterfaceCellPatch.patch().faceCells();
                const auto& neighbourValues = isInterfaceCellPatch.patchNeighbourField().cref();

                forAll(isInterfaceCellPatch, I)
                {
                    if ((isInterfaceCell[faceToCell[I]] == 1.0) && (neighbourValues[I] == 1.0))
                    {
                        // Curvature okay: face shared by interface cells
                        count[faceToCell[I]] += 1;
                        curvatureSum[faceToCell[I]] += faceCurvaturePatch[I];
                    }
                }
            }
        }

        forAll(count, cid)
        {
            if (count[cid] > 0)
            {
                curvature[cid] = curvatureSum[cid]/count[cid];
            }
        }

        count = 0;
        curvatureSum = 0.0;

        // Update processor neighbour values of curvature field
        curvature.correctBoundaryConditions();
    }
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// ************************************************************************* //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
